"""Candidate-record schema + canonical serialization.

A ``CandidateRecord`` is one entry's snapshot metadata — everything Stages 3-6
need, with no coordinates. ``candidates.jsonl`` is the canonical, byte-stable
serialization of these records (sorted entries, sorted keys), which is what the
snapshot lock hashes.

PDB-ID compatibility: identifiers (entry_id, entity_id) are stored *verbatim* as
returned by the RCSB Data API in ``rcsb_id`` — never sliced, length-validated, or
case-folded. This makes the schema agnostic to legacy 4-character IDs (``4HHB``,
entity ``4HHB_1``) and the extended ``pdb_xxxxxxxx`` form alike.
"""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Iterable

from pydantic import BaseModel, ConfigDict

# Metal-name words as they appear in GO/InterPro/Pfam annotation strings, mapped to
# their CCD element symbols. A name that mentions "metal" generically (e.g.
# "transition metal ion binding") with no named element yields the GENERIC_METAL
# sentinel. Stage 4 uses this to tell a native-Ni/Co metalloenzyme (its own metal is
# annotated) from a metalloprotein whose native metal is something else (a Ni/Co
# substitute) from a protein with no metal annotation at all.
METAL_ANNOTATION_TERMS: dict[str, str] = {
    "nickel": "NI", "cobalt": "CO", "zinc": "ZN", "iron": "FE", "ferrous": "FE",
    "ferric": "FE", "magnesium": "MG", "manganese": "MN", "calcium": "CA",
    "copper": "CU", "cadmium": "CD", "molybdenum": "MO", "tungsten": "W",
    # Heavy / f-block metals that a bound-ion tiering would otherwise treat as a
    # phasing artifact: a native site is identified here so it is rescued instead.
    # Native mercury (MerR/MerA "mercuric" proteins) and lanthanide-dependent
    # dehydrogenases (XoxF/ExaF use La/Ce/Pr/Nd/...) are the real cases. Words are
    # element-specific to avoid false hits ("lead"/"gold" are deliberately omitted).
    "mercury": "HG", "mercuric": "HG", "platinum": "PT", "thallium": "TL",
    "osmium": "OS", "iridium": "IR", "lanthanum": "LA", "cerium": "CE",
    "praseodymium": "PR", "neodymium": "ND", "samarium": "SM", "europium": "EU",
    "gadolinium": "GD", "terbium": "TB", "dysprosium": "DY", "holmium": "HO",
    "erbium": "ER", "thulium": "TM", "ytterbium": "YB", "lutetium": "LU",
}  # fmt: skip
GENERIC_METAL = "METAL"


def metal_symbols_in_annotation(name: str | None) -> set[str]:
    """Metal element symbols named in an annotation string.

    ``"nickel cation binding"`` -> ``{"NI"}``; ``"transition metal ion binding"`` ->
    ``{"METAL"}`` (generic); a term with no metal -> ``set()``. A named element wins
    over the generic sentinel.
    """
    n = (name or "").lower()
    # Match whole words, not raw substrings, so "ytterbium" does not also fire
    # "terbium"/"erbium" (each is a substring of the next). The generic sentinel keeps
    # a plain substring test so "metallopeptidase" still counts as a (generic) metal.
    words = set(re.findall(r"[a-z]+", n))
    named = {sym for word, sym in METAL_ANNOTATION_TERMS.items() if word in words}
    if named:
        return named
    return {GENERIC_METAL} if "metal" in n else set()


# Structural-classification methods for fold-level leakage control (Stage 5).
# CATH's annotation_id IS the homologous-superfamily code (e.g. "1.10.490.10"), a
# ready-made grouping key. ECOD/SCOP2 expose a per-domain annotation_id, so we key
# on the (super)family *name* they carry ("Bcl-2", "Globin-like") instead.
STRUCTURAL_METHODS = ("cath", "ecod", "scop2")


def structural_families_from_instances(instances) -> dict[str, list[str]]:
    """Per-method structural (super)family keys across an entity's chain instances.

    Returns e.g. ``{"cath": ["1.10.490.10"], "scop2": ["Globin-like"]}``; only
    non-empty methods are included. A multi-domain chain contributes several keys.
    """
    fams: dict[str, set[str]] = {m: set() for m in STRUCTURAL_METHODS}
    for inst in instances or []:
        for ann in inst.get("rcsb_polymer_instance_annotation") or []:
            atype = ann.get("type")
            if atype == "CATH" and ann.get("annotation_id"):
                fams["cath"].add(ann["annotation_id"])
            elif atype == "ECOD" and ann.get("name"):
                fams["ecod"].add(ann["name"])
            elif atype == "SCOP2" and ann.get("name"):
                fams["scop2"].add(ann["name"])
    return {m: sorted(v) for m, v in fams.items() if v}


class PolymerEntity(BaseModel):
    model_config = ConfigDict(extra="forbid")

    entity_id: str  # RCSB rcsb_id, verbatim (e.g. "4HHB_1")
    polymer_type: str  # rcsb_entity_polymer_type: Protein / DNA / RNA / NA-hybrid / Other
    seq_len: int
    seq: str
    # RCSB precomputed cluster ids by identity level, e.g. {30: 48, 95: 1239}.
    # Empty for non-protein entities (RCSB clusters proteins only).
    cluster_ids: dict[int, int] = {}
    # Metal element symbols this protein is annotated (GO/InterPro/Pfam) to bind,
    # e.g. ["NI"] for a urease, ["ZN"] for a zinc enzyme, ["METAL"] for a generic
    # "metal ion binding" term. Empty when no metal annotation exists. Stage 4 uses
    # it to rescue native Ni/Co metalloenzymes from the purification-artifact demotion.
    metal_annotations: list[str] = []
    # Structural (super)family keys per classification method, e.g.
    # {"cath": ["1.10.490.10"], "ecod": ["Bcl-2"], "scop2": ["Globin-like"]}.
    # Captured from RCSB's per-instance structural annotations (metadata, no
    # coordinates). Stage 5 uses the config-selected method to union entities that
    # share a fold. A multi-domain chain carries several families (one per domain).
    # Only non-empty methods are stored; empty for non-protein / unclassified chains.
    structural_families: dict[str, list[str]] = {}

    @property
    def is_protein(self) -> bool:
        return "PROTEIN" in self.polymer_type.upper()

    @property
    def is_nucleic(self) -> bool:
        # rcsb_entity_polymer_type is one of Protein / DNA / RNA / NA-hybrid / Other.
        # "NA-hybrid" (a single mixed DNA/RNA strand — R-loops, RNA-primed DNA, chimeric
        # guides) contains neither "DNA" nor "RNA", so match "HYBRID" too, else a
        # protein/hybrid complex would silently get no nucleic_acid class.
        t = self.polymer_type.upper()
        return "DNA" in t or "RNA" in t or "HYBRID" in t


class NonpolymerComp(BaseModel):
    model_config = ConfigDict(extra="forbid")

    comp_id: str  # chem_comp id, e.g. "HEM", "ZN" (CCD codes are uppercase)
    name: str | None = None
    formula: str | None = None
    comp_type: str | None = None


class QualityMetrics(BaseModel):
    """wwPDB validation-report summary metrics (no coordinates).

    Geometry metrics (clashscore, Ramachandran, rotamer) are reported for both
    X-ray and cryo-EM; diffraction metrics (R-free, RSRZ) are X-ray only; EM
    map-fit (backbone atom inclusion) is cryo-EM only. Any field may be ``None``
    when the validation report does not provide it — a missing metric never
    penalizes an entry (see :func:`ifsplit.parse.quality_drop`).
    """

    model_config = ConfigDict(extra="forbid")

    clashscore: float | None = None  # all-atom clashes / 1000 atoms (lower better)
    ramachandran_outlier_pct: float | None = None  # % backbone Ramachandran outliers
    rotamer_outlier_pct: float | None = None  # % sidechain rotamer outliers
    rfree: float | None = None  # diffraction DCC_Rfree (X-ray only)
    rsrz_outlier_pct: float | None = None  # diffraction % RSRZ outliers (X-ray only)
    em_backbone_inclusion: float | None = None  # EM backbone atom inclusion (higher better)

    @property
    def has_report(self) -> bool:
        """True if the validation report supplied at least one metric."""
        return any(
            v is not None
            for v in (
                self.clashscore,
                self.ramachandran_outlier_pct,
                self.rotamer_outlier_pct,
                self.rfree,
                self.rsrz_outlier_pct,
                self.em_backbone_inclusion,
            )
        )


class CandidateRecord(BaseModel):
    """One PDB entry's snapshot metadata."""

    model_config = ConfigDict(extra="forbid")

    entry_id: str
    methods: list[str]
    resolution_A: float | None
    release_date: str  # YYYY-MM-DD
    deposited_residues: int | None
    assemblies: dict[str, int]  # assembly_id -> polymer_monomer_count
    polymer_entities: list[PolymerEntity]
    nonpolymer_comps: list[NonpolymerComp]
    # Curation signals (Stage 4). comp ids that actually contact the protein
    # (rcsb_entry_info.nonpolymer_bound_components) and comp ids with a measured
    # binding affinity (rcsb_binding_affinity). Both are the buffer-vs-ligand gate.
    bound_components: list[str] = []
    affinity_comp_ids: list[str] = []
    # comp ids RCSB curated as "subject of investigation" (a ligand of interest).
    # Catches tightly but non-covalently bound cofactors/inhibitors that
    # bound_components (bond-based) misses; another positive buffer-vs-ligand gate.
    investigated_comp_ids: list[str] = []
    # RCSB-computed count of protein<->nucleic-acid interface entities across the
    # entry's assemblies (rcsb_assembly_info.num_prot_na_interface_entities). > 0
    # verifies a real protein/NA contact (Stage 4 holo gate for the nucleic_acid
    # class) — distinguishing a true complex from a co-deposited oligo.
    protein_na_interface_count: int = 0
    # wwPDB validation-report summary (Stage 3 quality filters). All fields
    # optional; defaults to an empty report when RCSB has none.
    quality: QualityMetrics = QualityMetrics()

    @classmethod
    def from_data_api(cls, entry: dict) -> CandidateRecord:
        """Build a record from a raw Data-API entry object (deterministic)."""
        # Verbatim canonical id from RCSB (legacy or extended) — never reformat.
        entry_id = entry["rcsb_id"]

        methods = sorted(m["method"] for m in (entry.get("exptl") or []) if m.get("method"))

        info = entry.get("rcsb_entry_info") or {}
        res_list = info.get("resolution_combined") or []
        resolution = min(res_list) if res_list else None
        deposited = info.get("deposited_polymer_monomer_count")
        bound = sorted({c.upper() for c in (info.get("nonpolymer_bound_components") or [])})

        acc = entry.get("rcsb_accession_info") or {}
        rel = acc.get("initial_release_date")
        release_date = rel[:10] if rel else ""

        affinity = sorted(
            {
                a["comp_id"].upper()
                for a in (entry.get("rcsb_binding_affinity") or [])
                if a.get("comp_id")
            }
        )

        assemblies: dict[str, int] = {}
        prot_na_interfaces = 0
        for asm in entry.get("assemblies") or []:
            aid = asm.get("rcsb_id")
            info = asm.get("rcsb_assembly_info") or {}
            count = info.get("polymer_monomer_count")
            if aid is not None and count is not None:
                assemblies[aid] = count
            prot_na_interfaces += info.get("num_prot_na_interface_entities") or 0

        polymers: list[PolymerEntity] = []
        for p in entry.get("polymer_entities") or []:
            poly = p.get("entity_poly") or {}
            seq = poly.get("pdbx_seq_one_letter_code_can") or ""
            seq = "".join(seq.split())  # strip newlines/whitespace the API may insert
            cluster_ids: dict[int, int] = {}
            for cm in p.get("rcsb_cluster_membership") or []:
                ident = cm.get("identity")
                cid = cm.get("cluster_id")
                if ident is not None and cid is not None:
                    cluster_ids[int(ident)] = int(cid)
            metal_annots: set[str] = set()
            for ann in p.get("rcsb_polymer_entity_annotation") or []:
                metal_annots |= metal_symbols_in_annotation(ann.get("name"))
            structural = structural_families_from_instances(p.get("polymer_entity_instances"))
            polymers.append(
                PolymerEntity(
                    entity_id=p["rcsb_id"],  # verbatim
                    polymer_type=poly.get("rcsb_entity_polymer_type") or "Other",
                    seq_len=len(seq),
                    seq=seq,
                    cluster_ids=cluster_ids,
                    metal_annotations=sorted(metal_annots),
                    structural_families=structural,
                )
            )
        polymers.sort(key=lambda e: e.entity_id)

        comps: dict[str, NonpolymerComp] = {}
        investigated: set[str] = set()
        for n in entry.get("nonpolymer_entities") or []:
            cc = (n.get("nonpolymer_comp") or {}).get("chem_comp") or {}
            cid = cc.get("id")
            if cid and cid.upper() not in comps:
                comps[cid.upper()] = NonpolymerComp(
                    comp_id=cid.upper(),
                    name=cc.get("name"),
                    formula=cc.get("formula"),
                    comp_type=cc.get("type"),
                )
            # "subject of investigation" = Y on any instance -> this comp is a
            # curated ligand of interest (covalent or not).
            for inst in n.get("nonpolymer_entity_instances") or []:
                for vs in inst.get("rcsb_nonpolymer_instance_validation_score") or []:
                    if cid and vs.get("is_subject_of_investigation") == "Y":
                        investigated.add(cid.upper())
        nonpolymers = [comps[k] for k in sorted(comps)]

        # Validation-report summaries arrive as 1-element lists (or null).
        def _first(items) -> dict:
            return (items or [None])[0] or {}

        geo = _first(entry.get("pdbx_vrpt_summary_geometry"))
        dif = _first(entry.get("pdbx_vrpt_summary_diffraction"))
        em = _first(entry.get("pdbx_vrpt_summary_em"))
        quality = QualityMetrics(
            clashscore=geo.get("clashscore"),
            ramachandran_outlier_pct=geo.get("percent_ramachandran_outliers"),
            rotamer_outlier_pct=geo.get("percent_rotamer_outliers"),
            rfree=dif.get("DCC_Rfree"),
            rsrz_outlier_pct=dif.get("percent_RSRZ_outliers"),
            em_backbone_inclusion=em.get("atom_inclusion_backbone"),
        )

        return cls(
            entry_id=entry_id,
            methods=methods,
            resolution_A=resolution,
            release_date=release_date,
            deposited_residues=deposited,
            assemblies=dict(sorted(assemblies.items())),
            polymer_entities=polymers,
            nonpolymer_comps=nonpolymers,
            bound_components=bound,
            affinity_comp_ids=affinity,
            investigated_comp_ids=sorted(investigated),
            protein_na_interface_count=prot_na_interfaces,
            quality=quality,
        )

    def to_canonical_json(self) -> str:
        """Single-line, sorted-key JSON for byte-stable serialization."""
        return json.dumps(self.model_dump(mode="json"), sort_keys=True, separators=(",", ":"))


def canonical_jsonl_bytes(records: Iterable[CandidateRecord]) -> bytes:
    """Byte-stable candidates.jsonl content: records sorted by entry_id."""
    ordered = sorted(records, key=lambda r: r.entry_id)
    text = "".join(r.to_canonical_json() + "\n" for r in ordered)
    return text.encode("utf-8")


def sha256_hex(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def read_candidates_jsonl(path) -> list[CandidateRecord]:
    """Parse a candidates.jsonl file back into validated records."""
    records: list[CandidateRecord] = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line:
                records.append(CandidateRecord.model_validate_json(line))
    return records
