"""Thin, polite RCSB API client (Search v2 + Data GraphQL).

No coordinates are ever fetched here — only metadata and sequences (see
PLAN.md §1.5). Endpoints are centralized so field paths live in one place, and
all requests retry with backoff on transient failures (429 / 5xx / network).
"""

from __future__ import annotations

import time
from collections.abc import Iterator

import httpx

from .config import Config

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
DATA_GRAPHQL_URL = "https://data.rcsb.org/graphql"

# Conservative page/batch sizes: large enough to keep request counts sane,
# small enough to stay friendly and within response-size limits. The Data API
# caps batch endpoints at 1000 ids.
SEARCH_PAGE_ROWS = 5000
DATA_BATCH_SIZE = 200

_RETRY_STATUS = {429, 500, 502, 503, 504}

# GraphQL: everything Stages 3-6 need, no coordinates. Validated live against
# 4HHB / 1A1F / 1IEP. Notable curation signals:
#   - rcsb_cluster_membership: precomputed cluster id per identity level
#     (30/50/70/90/95/100) for protein entities -> no cluster-file download.
#   - rcsb_entry_info.nonpolymer_bound_components: comp ids that actually contact
#     the protein (the cheap buffer-vs-ligand gate; e.g. 4HHB -> ["HEM"], its
#     PO4/Cl buffer is absent).
#   - rcsb_binding_affinity.comp_id: comps with a measured affinity (sparse but a
#     strong positive "this is a real ligand" signal).
#   - pdbx_vrpt_summary_{geometry,diffraction,em}: wwPDB validation-report metrics
#     (clashscore, Ramachandran/rotamer outliers, R-free, RSRZ). Geometry is
#     reported for X-ray AND EM; diffraction is X-ray-only; each comes back as a
#     1-element list. Metadata, not coordinates — keeps the no-download invariant.
#   - assemblies.interfaces.rcsb_interface_info.polymer_composition: RCSB-computed
#     assembly interfaces; a "Protein/NA" interface verifies a real protein<->DNA/RNA
#     contact (the holo gate for the nucleotide class). Present for X-ray AND EM.
_ENTRY_QUERY = """
query($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    exptl { method }
    rcsb_entry_info {
      resolution_combined
      deposited_polymer_monomer_count
      nonpolymer_bound_components
    }
    rcsb_accession_info { initial_release_date }
    rcsb_binding_affinity { comp_id }
    pdbx_vrpt_summary_geometry {
      clashscore
      percent_ramachandran_outliers
      percent_rotamer_outliers
    }
    pdbx_vrpt_summary_diffraction {
      DCC_Rfree
      percent_RSRZ_outliers
    }
    pdbx_vrpt_summary_em {
      atom_inclusion_backbone
    }
    polymer_entities {
      rcsb_id
      entity_poly {
        rcsb_entity_polymer_type
        pdbx_seq_one_letter_code_can
      }
      rcsb_cluster_membership { cluster_id identity }
    }
    nonpolymer_entities {
      nonpolymer_comp { chem_comp { id name formula type } }
    }
    assemblies {
      rcsb_id
      rcsb_assembly_info { polymer_monomer_count }
      interfaces {
        rcsb_interface_info { polymer_composition }
      }
    }
  }
}
"""


class RcsbError(RuntimeError):
    """Raised when an RCSB request fails after exhausting retries."""


class RcsbClient:
    """Minimal client for the two RCSB services IF-Split needs."""

    def __init__(
        self,
        *,
        timeout: float = 60.0,
        max_retries: int = 4,
        backoff_base: float = 1.5,
        sleep=time.sleep,
    ) -> None:
        self._client = httpx.Client(
            timeout=timeout,
            headers={"User-Agent": "IF-Split/0.1 (reproducible PDB splitter)"},
        )
        self._max_retries = max_retries
        self._backoff_base = backoff_base
        self._sleep = sleep

    def __enter__(self) -> RcsbClient:
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def close(self) -> None:
        self._client.close()

    # --- low-level POST with retry/backoff ---
    def _post(self, url: str, json_body: dict) -> httpx.Response:
        last_exc: Exception | None = None
        for attempt in range(self._max_retries + 1):
            try:
                resp = self._client.post(url, json=json_body)
            except httpx.HTTPError as exc:  # network/timeout
                last_exc = exc
            else:
                if resp.status_code not in _RETRY_STATUS:
                    return resp
                last_exc = RcsbError(f"{url} -> HTTP {resp.status_code}")
            if attempt < self._max_retries:
                self._sleep(self._backoff_base**attempt)
        raise RcsbError(f"request to {url} failed after retries: {last_exc}")

    # --- Search API: entry IDs matching the snapshot filters ---
    def _search_query_body(self, cfg: Config) -> dict:
        cutoff = f"{cfg.snapshot_date.isoformat()}T23:59:59Z"
        return {
            "query": {
                "type": "group",
                "logical_operator": "and",
                "nodes": [
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "exptl.method",
                            "operator": "in",
                            "value": list(cfg.experimental_methods),
                        },
                    },
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_entry_info.resolution_combined",
                            "operator": "less_or_equal",
                            "value": cfg.resolution_max_A,
                        },
                    },
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_accession_info.initial_release_date",
                            "operator": "less_or_equal",
                            "value": cutoff,
                        },
                    },
                ],
            },
            "return_type": "entry",
        }

    def count_entries(self, cfg: Config) -> int:
        """Total entries matching the snapshot filters (no paging)."""
        body = self._search_query_body(cfg)
        body["request_options"] = {"return_counts": True}
        resp = self._post(SEARCH_URL, body)
        resp.raise_for_status()
        return int(resp.json()["total_count"])

    def search_entry_ids(
        self, cfg: Config, limit: int | None = None, *, progress=None
    ) -> list[str]:
        """All matching entry IDs, sorted ascending for determinism.

        ``limit`` (dev convenience) takes the first N in sorted order, so a
        limited run is itself reproducible. ``progress`` (optional) is called with
        a status string after each page, so a full-PDB enumeration isn't silent.
        """
        ids: list[str] = []
        start = 0
        while True:
            rows = SEARCH_PAGE_ROWS
            if limit is not None:
                remaining = limit - len(ids)
                if remaining <= 0:
                    break
                rows = min(rows, remaining)
            body = self._search_query_body(cfg)
            body["request_options"] = {
                "paginate": {"start": start, "rows": rows},
                "sort": [
                    {
                        "sort_by": "rcsb_entry_container_identifiers.entry_id",
                        "direction": "asc",
                    }
                ],
            }
            resp = self._post(SEARCH_URL, body)
            if resp.status_code == 204:  # no (more) results
                break
            resp.raise_for_status()
            page = resp.json().get("result_set", [])
            if not page:
                break
            ids.extend(hit["identifier"] for hit in page)
            start += len(page)
            if progress:
                progress(f"search: {len(ids)} entry ids found...")
            if len(page) < rows:
                break
        return ids

    # --- Data API: batched metadata enrichment ---
    def fetch_entries(self, ids: list[str]) -> Iterator[dict]:
        """Yield raw Data-API entry objects for ``ids``, batched."""
        for batch in _chunks(ids, DATA_BATCH_SIZE):
            resp = self._post(
                DATA_GRAPHQL_URL,
                {"query": _ENTRY_QUERY, "variables": {"ids": batch}},
            )
            resp.raise_for_status()
            payload = resp.json()
            if payload.get("errors"):
                raise RcsbError(f"GraphQL errors: {payload['errors'][:1]}")
            for entry in payload["data"]["entries"]:
                if entry is not None:
                    yield entry


def _chunks(seq: list[str], size: int) -> Iterator[list[str]]:
    for i in range(0, len(seq), size):
        yield seq[i : i + size]
