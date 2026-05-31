"""Stage 4 - Classify non-protein components and curate purification artifacts.

Classification is done from metadata alone (chem-comp formula + polymer types),
no coordinates:

- **nucleotide**: DNA/RNA polymer entities (typed by the Data API; the
  ATOM-record gotcha never arises on the metadata path).
- **metal**: non-polymer comps whose formula is purely metal element(s).
- **small-molecule**: remaining non-polymer comps after removing waters, the
  configured ``excluded_het``, and a crystallization-additive blacklist.

Purification-artifact curation (the LigandMPNN metal blemish): a poly-His tag
coordinating Ni/Co is an artifact of immobilized-metal affinity purification,
not a biological metal site. An entry whose *only* metal is a purification metal
(Ni/Co by default) AND that carries a His-tag is flagged; with
``exclude_purification_artifacts`` it is dropped from the metal class.
"""

from __future__ import annotations

import re

from .config import Config
from .schema import CandidateRecord, NonpolymerComp

# Ligand classes IF-Split tracks.
CLASS_METAL = "metal"
CLASS_NUCLEOTIDE = "nucleotide"
CLASS_SMALL_MOLECULE = "small_molecule"

# Elements that count as a "metal ion" when a component's formula is made up
# only of these (transition metals, alkali/alkaline-earth, and common metallic
# ions seen as ligands). A comp like HEM (C34 H32 Fe N4 O4) is NOT a metal ion —
# its formula contains C/H/N/O — so it falls through to small-molecule.
METAL_ELEMENTS: frozenset[str] = frozenset(
    {
        "LI",
        "BE",
        "NA",
        "MG",
        "AL",
        "K",
        "CA",
        "SC",
        "TI",
        "V",
        "CR",
        "MN",
        "FE",
        "CO",
        "NI",
        "CU",
        "ZN",
        "GA",
        "RB",
        "SR",
        "Y",
        "ZR",
        "NB",
        "MO",
        "TC",
        "RU",
        "RH",
        "PD",
        "AG",
        "CD",
        "IN",
        "SN",
        "CS",
        "BA",
        "LA",
        "CE",
        "PR",
        "ND",
        "PM",
        "SM",
        "EU",
        "GD",
        "TB",
        "DY",
        "HO",
        "ER",
        "TM",
        "YB",
        "LU",
        "HF",
        "TA",
        "W",
        "RE",
        "OS",
        "IR",
        "PT",
        "AU",
        "HG",
        "TL",
        "PB",
        "BI",
    }
)

# Common crystallization additives / buffer junk (not biologically meaningful
# ligands). Config-extensible via excluded_het; this is the always-on baseline.
DEFAULT_ADDITIVE_BLACKLIST: frozenset[str] = frozenset(
    {
        "HOH",
        "DOD",  # water
        "GOL",
        "EDO",
        "PEG",
        "PG4",
        "PGE",
        "1PE",
        "2PE",
        "P6G",
        "PE4",
        "MPD",
        "SO4",
        "PO4",
        "ACT",
        "ACY",
        "FMT",
        "EPE",
        "MES",
        "TRS",
        "BME",
        "DTT",
        "IMD",
        "DMS",
        "BOG",
        "OLC",
        "LDA",
        "SCN",
        "AZI",
        "NO3",
        "CO3",
        "FLC",
        "TLA",
        "CIT",
        "MLI",
        "IOD",
        "GLC",
        "BCN",
        "MRD",
        "BU3",
        "P33",
    }
)

_ELEMENT_RE = re.compile(r"[A-Z][a-z]?")


def elements_in_formula(formula: str | None) -> set[str]:
    """Element symbols present in a chem-comp formula, upper-cased.

    ``"C34 H32 Fe N4 O4"`` -> ``{"C", "H", "FE", "N", "O"}``; ``"Zn"`` -> ``{"ZN"}``.
    Charge tokens and counts are ignored.
    """
    if not formula:
        return set()
    # Drop charge markers like "1+", "2-", and standalone counts.
    cleaned = re.sub(r"[0-9]+[+-]?", " ", formula)
    return {m.group(0).upper() for m in _ELEMENT_RE.finditer(cleaned)}


def is_metal_ion(comp: NonpolymerComp) -> bool:
    """True if the component's formula is composed only of metal element(s)."""
    elems = elements_in_formula(comp.formula)
    if not elems:
        # Fall back to single-/double-letter comp ids that are bare elements.
        elems = {comp.comp_id} if comp.comp_id in METAL_ELEMENTS else set()
    return bool(elems) and elems <= METAL_ELEMENTS


def longest_residue_run(seq: str, residue: str = "H") -> int:
    """Length of the longest consecutive run of ``residue`` in ``seq``."""
    if not seq:
        return 0
    best = run = 0
    up = seq.upper()
    res = residue.upper()
    for ch in up:
        if ch == res:
            run += 1
            best = max(best, run)
        else:
            run = 0
    return best


def has_histag(seq: str, min_run: int) -> bool:
    """True if ``seq`` contains a poly-histidine run of at least ``min_run``."""
    return longest_residue_run(seq, "H") >= min_run


def metal_comps(record: CandidateRecord) -> list[NonpolymerComp]:
    return [c for c in record.nonpolymer_comps if is_metal_ion(c)]


def is_purification_artifact(
    record: CandidateRecord,
    *,
    purification_metals: set[str],
    histag_min_run: int,
) -> bool:
    """Flag the His-tag-binds-Ni/Co purification-artifact pattern.

    True only when (a) the entry has at least one metal, (b) *every* metal it has
    is a purification metal, and (c) some protein chain carries a His-tag. A real
    metal site (e.g. a catalytic Zn) present alongside a tag is NOT flagged.
    """
    if not purification_metals:
        return False
    metals = {c.comp_id for c in metal_comps(record)}
    if not metals or not metals <= purification_metals:
        return False
    return any(has_histag(e.seq, histag_min_run) for e in record.polymer_entities if e.is_protein)


def classify_components(record: CandidateRecord, cfg: Config) -> dict:
    """Classify one entry's components into ligand classes (metadata only).

    Returns a dict with the per-class component ids, the set of ligand-class
    tags, and curation flags (e.g. purification artifact).
    """
    blacklist = DEFAULT_ADDITIVE_BLACKLIST | set(cfg.excluded_het)
    purification = set(cfg.purification_metals)

    artifact = is_purification_artifact(
        record,
        purification_metals=purification,
        histag_min_run=cfg.histag_min_run,
    )

    metals: list[str] = []
    small_molecules: list[str] = []
    for comp in record.nonpolymer_comps:
        if is_metal_ion(comp):
            if artifact and cfg.exclude_purification_artifacts and comp.comp_id in purification:
                continue  # drop purification-metal from the metal class
            metals.append(comp.comp_id)
        elif comp.comp_id not in blacklist:
            small_molecules.append(comp.comp_id)

    has_nucleotide = any(e.is_nucleic for e in record.polymer_entities)

    tags: set[str] = set()
    if metals:
        tags.add(CLASS_METAL)
    if small_molecules:
        tags.add(CLASS_SMALL_MOLECULE)
    if has_nucleotide:
        tags.add(CLASS_NUCLEOTIDE)

    return {
        "entry_id": record.entry_id,
        "classes": sorted(tags),
        "metals": sorted(metals),
        "small_molecules": sorted(small_molecules),
        "has_nucleotide": has_nucleotide,
        "purification_artifact": artifact,
    }
