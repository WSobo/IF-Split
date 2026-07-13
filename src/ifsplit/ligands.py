"""Stage 4 - Tier and classify non-protein components (metadata only).

Curation comes *before* classification, and quality is **annotated, never
destroyed**: no structure is dropped for ligand-quality reasons (a protein with a
junk ion is still a good training backbone). Each non-protein component is tiered

  - ``functional``  : real, biologically meaningful ligand/site
  - ``ambiguous``   : present and contacting, but unconfirmed
  - ``artifact``    : buffer / cryoprotectant / counterion / purification tag

with a machine-readable reason. Ligand-*class* labels (metal / small_molecule /
nucleic_acid) derive from the tier via a config threshold (default: only
``functional`` sets a class label; ``ambiguous`` is reported but not labelled;
``artifact`` is excluded). A downstream featurizer reads the same per-component
tier to decide what counts as real ligand context — this is the lever that
improves *training* quality, not just test reporting.

Signals used (all from the Data API, no coordinates):
  - ``bound_components``  : the comp is covalently/coordination-bonded to the
                            protein (buffer gate; bond-based, misses non-covalent)
  - ``investigated_comp_ids`` : RCSB curated it a "subject of investigation" — a
                            ligand of interest; catches non-covalent cofactors
                            (FAD/NAD/FMN/NADP, inhibitors) bound_components misses
  - ``affinity_comp_ids`` : a measured binding affinity exists (strong positive)
  - chem-comp ``formula`` : metal-only vs organic
  - His-tag + Ni/Co       : the IMAC purification-artifact pattern (existing rule)
  - protein_na_interface_count: a protein<->nucleic-acid assembly interface (>0)
                            verifies a real contact (holo gate for nucleic_acid)
"""

from __future__ import annotations

import re

from .config import Config
from .schema import CandidateRecord, NonpolymerComp

# Ligand classes IF-Split tracks.
CLASS_METAL = "metal"
CLASS_NUCLEIC_ACID = "nucleic_acid"
CLASS_SMALL_MOLECULE = "small_molecule"

# Confidence tiers.
TIER_FUNCTIONAL = "functional"
TIER_AMBIGUOUS = "ambiguous"
TIER_ARTIFACT = "artifact"

# Elements that count as a "metal ion" when a component's formula is made up
# only of these. A comp like HEM (C34 H32 Fe N4 O4) is NOT a metal ion — its
# formula contains C/H/N/O — so it falls through to small-molecule.
METAL_ELEMENTS: frozenset[str] = frozenset(
    {
        "LI", "BE", "NA", "MG", "AL", "K", "CA", "SC", "TI", "V", "CR", "MN",
        "FE", "CO", "NI", "CU", "ZN", "GA", "RB", "SR", "Y", "ZR", "NB", "MO",
        "TC", "RU", "RH", "PD", "AG", "CD", "IN", "SN", "CS", "BA", "LA", "CE",
        "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB",
        "LU", "HF", "TA", "W", "RE", "OS", "IR", "PT", "AU", "HG", "TL", "PB",
        "BI",
    }
)  # fmt: skip

# Monatomic counterions / phasing atoms: even though some are technically metals,
# as lone ions they are crystallization/phasing additives, not biological sites.
COUNTERION_COMPS: frozenset[str] = frozenset(
    {"NA", "CL", "K", "BR", "I", "IOD", "CS", "RB", "F", "LI"}
)

# Common crystallization additives / buffer junk. Config-extensible via
# excluded_het; this is the always-on baseline.
DEFAULT_ADDITIVE_BLACKLIST: frozenset[str] = frozenset(
    {
        "HOH", "DOD",  # water
        "GOL", "EDO", "PEG", "PG4", "PGE", "1PE", "2PE", "P6G", "PE4", "MPD",
        "SO4", "PO4", "ACT", "ACY", "FMT", "EPE", "MES", "TRS", "BME", "DTT",
        "IMD", "DMS", "BOG", "OLC", "LDA", "SCN", "AZI", "NO3", "CO3", "FLC",
        "TLA", "CIT", "MLI", "IOD", "GLC", "BCN", "MRD", "BU3", "P33",
    }
)  # fmt: skip

_ELEMENT_RE = re.compile(r"[A-Z][a-z]?")


def elements_in_formula(formula: str | None) -> set[str]:
    """Element symbols present in a chem-comp formula, upper-cased.

    ``"C34 H32 Fe N4 O4"`` -> ``{"C", "H", "FE", "N", "O"}``; ``"Zn"`` -> ``{"ZN"}``.
    Charge tokens and counts are ignored.
    """
    if not formula:
        return set()
    cleaned = re.sub(r"[0-9]+[+-]?", " ", formula)
    return {m.group(0).upper() for m in _ELEMENT_RE.finditer(cleaned)}


def is_metal_ion(comp: NonpolymerComp) -> bool:
    """True if the component's formula is composed only of metal element(s)."""
    elems = elements_in_formula(comp.formula)
    if not elems:
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


def terminal_residue_run(seq: str, residue: str = "H", window: int = 20) -> int:
    """Longest run of ``residue`` within ``window`` residues of either terminus.

    A purification tag is terminal by construction, so a short His run near a
    chain end is tag-like; the same run buried internally is more likely a real
    motif. Used to recover His-tags partially unmodeled or trimmed from the
    deposited sequence (where the full 6xHis no longer appears).
    """
    if not seq:
        return 0
    s = seq.upper()
    return max(longest_residue_run(s[:window], residue), longest_residue_run(s[-window:], residue))


def has_histag(seq: str, min_run: int, terminal_min_run: int = 0) -> bool:
    """True if ``seq`` carries a poly-histidine tag.

    Positive if there is a His run of at least ``min_run`` anywhere, OR (when
    ``terminal_min_run`` > 0) a shorter His run of at least ``terminal_min_run``
    at a chain terminus — catching 6xHis tags left partial by unmodeled/trimmed
    residues without flagging internal His-rich metal sites.
    """
    if longest_residue_run(seq, "H") >= min_run:
        return True
    return bool(terminal_min_run) and terminal_residue_run(seq, "H") >= terminal_min_run


def has_protein_na_interface(record: CandidateRecord) -> bool:
    """True if the assembly has a protein<->nucleic-acid interface.

    Reads RCSB's precomputed ``num_prot_na_interface_entities`` — a metadata signal
    that the protein actually *contacts* the DNA/RNA, not just that both were
    co-deposited. No coordinates.
    """
    return record.protein_na_interface_count > 0


def metal_comps(record: CandidateRecord) -> list[NonpolymerComp]:
    return [c for c in record.nonpolymer_comps if is_metal_ion(c)]


def is_purification_artifact(
    record: CandidateRecord,
    *,
    purification_metals: set[str],
    histag_min_run: int,
    histag_terminal_min_run: int = 0,
) -> bool:
    """Flag the His-tag-binds-Ni/Co purification-artifact pattern.

    True only when (a) the entry has at least one metal, (b) *every* metal it has
    is a purification metal, and (c) some protein chain carries a His-tag (full
    run, or a short *terminal* run per ``histag_terminal_min_run``). A real metal
    site (e.g. a catalytic Zn) present alongside a tag is NOT flagged.
    """
    if not purification_metals:
        return False
    metals = {c.comp_id for c in metal_comps(record)}
    if not metals or not metals <= purification_metals:
        return False
    return any(
        has_histag(e.seq, histag_min_run, histag_terminal_min_run)
        for e in record.polymer_entities
        if e.is_protein
    )


def tier_component(
    comp: NonpolymerComp,
    record: CandidateRecord,
    cfg: Config,
    *,
    blacklist: set[str],
    purification: set[str],
    is_artifact_entry: bool,
    uncorroborated_purification: bool = False,
) -> tuple[str, str, str | None]:
    """Tier a single non-protein component.

    Returns ``(tier, reason, class_label_or_None)``. The class label is set only
    for the ``functional`` tier (the default threshold); ``ambiguous`` components
    return their would-be class as None so they are reported but not labelled.

    ``uncorroborated_purification`` is set by the caller for a purification metal
    (Ni/Co) that is the entry's *only* metal and has no measured affinity and no
    detected His-tag: ~82% of lone Ni/Co carry no detectable His-tag in the
    deposited sequence (IMAC tags are often absent from SEQRES, not just
    unmodeled), so a bare contact is weak evidence of a biological site -> it is
    demoted from functional to ambiguous (reported, never dropped). This
    deliberately over-fires on some real bare-Ni/Co enzymes; the annotate-never-
    drop tier keeps them recoverable. See scripts/audit_nico_histag.py.
    """
    cid = comp.comp_id
    bound = cid in set(record.bound_components)
    has_affinity = cid in set(record.affinity_comp_ids)
    investigated = cid in set(record.investigated_comp_ids)

    if is_metal_ion(comp):
        # Lone counterion / phasing atom -> artifact regardless of binding.
        if cid in COUNTERION_COMPS:
            return TIER_ARTIFACT, "counterion", None
        # His-tag/Ni(Co) IMAC purification artifact.
        if is_artifact_entry and cfg.exclude_purification_artifacts and cid in purification:
            return TIER_ARTIFACT, "histag_metal", None
        if has_affinity:
            return TIER_FUNCTIONAL, "metal_affinity", CLASS_METAL
        # Lone, uncorroborated purification metal: a likely IMAC artifact whose tag
        # is absent from the deposited sequence -> ambiguous, not functional.
        if uncorroborated_purification:
            return TIER_AMBIGUOUS, "purification_metal_uncorroborated", None
        if bound:
            return TIER_FUNCTIONAL, "metal_bound", CLASS_METAL
        if investigated:
            return TIER_FUNCTIONAL, "metal_investigated", CLASS_METAL
        # "Trust biological metals" but require contact: an unbound metal far from
        # the protein is adventitious -> ambiguous, not functional.
        return TIER_AMBIGUOUS, "metal_unbound", None

    # Non-metal small molecule.
    if cid in blacklist:
        return TIER_ARTIFACT, "additive", None
    if has_affinity:
        return TIER_FUNCTIONAL, "ligand_affinity", CLASS_SMALL_MOLECULE
    if bound:
        return TIER_FUNCTIONAL, "ligand_bound", CLASS_SMALL_MOLECULE
    # RCSB-curated "subject of investigation": catches non-covalently bound
    # cofactors/inhibitors (FAD/NAD/FMN/NADP, ...) that the bond-based
    # bound_components field misses.
    if investigated:
        return TIER_FUNCTIONAL, "ligand_investigated", CLASS_SMALL_MOLECULE
    return TIER_AMBIGUOUS, "ligand_unbound", None


def classify_components(record: CandidateRecord, cfg: Config) -> dict:
    """Tier + classify one entry's components into ligand classes (metadata only).

    The structure is always kept; only labels and tiers are assigned. Returns
    per-component tiers, the functional class tags, the ambiguous class tags, and
    curation flags.
    """
    blacklist = DEFAULT_ADDITIVE_BLACKLIST | set(cfg.excluded_het)
    purification = set(cfg.purification_metals)
    is_artifact_entry = is_purification_artifact(
        record,
        purification_metals=purification,
        histag_min_run=cfg.histag_min_run,
        histag_terminal_min_run=cfg.histag_terminal_min_run,
    )

    # An entry whose *only* metal is a purification metal (Ni/Co) is a lone-Ni/Co
    # case. If it's also not a detected His-tag artifact, a bare Ni/Co is still weak
    # evidence: ~82% of lone Ni/Co carry no detectable His-tag (IMAC tags are often
    # absent from the deposited SEQRES) so a scan can't recover them -> without
    # affinity/SOI, demote to ambiguous. (scripts/audit_nico_histag.py)
    entry_metals = {c.comp_id for c in metal_comps(record)}
    lone_purification_metal = (
        bool(purification)
        and bool(entry_metals)
        and entry_metals <= purification
        and not is_artifact_entry
    )

    tiers: dict[str, dict[str, str]] = {}
    functional_metals: list[str] = []
    functional_sms: list[str] = []
    ambiguous_classes: set[str] = set()

    affinity_ids = set(record.affinity_comp_ids)
    investigated_ids = set(record.investigated_comp_ids)
    for comp in record.nonpolymer_comps:
        # A lone Ni/Co is "uncorroborated" only if nothing vouches for it: no
        # measured affinity AND not curated as a subject of investigation.
        uncorroborated = (
            lone_purification_metal
            and comp.comp_id in purification
            and comp.comp_id not in affinity_ids
            and comp.comp_id not in investigated_ids
        )
        tier, reason, label = tier_component(
            comp,
            record,
            cfg,
            blacklist=blacklist,
            purification=purification,
            is_artifact_entry=is_artifact_entry,
            uncorroborated_purification=uncorroborated,
        )
        tiers[comp.comp_id] = {"tier": tier, "reason": reason}
        if tier == TIER_FUNCTIONAL and label == CLASS_METAL:
            functional_metals.append(comp.comp_id)
        elif tier == TIER_FUNCTIONAL and label == CLASS_SMALL_MOLECULE:
            functional_sms.append(comp.comp_id)
        elif tier == TIER_AMBIGUOUS:
            # Record the would-be class for reporting (metal vs small molecule).
            ambiguous_classes.add(CLASS_METAL if is_metal_ion(comp) else CLASS_SMALL_MOLECULE)

    # nucleic_acid class = the entry has DNA/RNA *polymer chains*. Functional only
    # if the protein actually *interfaces* the nucleic acid (RCSB assembly-interface
    # metadata); an NA chain with no protein/NA interface is co-deposited, not holo,
    # so the class is reported as ambiguous, not labelled. NB: this is the
    # protein/nucleic-acid *complex* class, NOT bound mononucleotide ligands
    # (ATP/GTP/NAD) -- those are handled above as small molecules.
    has_nucleic_acid = any(e.is_nucleic for e in record.polymer_entities)
    nucleic_acid_functional = has_nucleic_acid and has_protein_na_interface(record)
    if has_nucleic_acid:
        if nucleic_acid_functional:
            tiers["nucleic_acid"] = {"tier": TIER_FUNCTIONAL, "reason": "protein_na_interface"}
        else:
            ambiguous_classes.add(CLASS_NUCLEIC_ACID)
            tiers["nucleic_acid"] = {"tier": TIER_AMBIGUOUS, "reason": "no_protein_na_interface"}

    classes: set[str] = set()
    if functional_metals:
        classes.add(CLASS_METAL)
    if functional_sms:
        classes.add(CLASS_SMALL_MOLECULE)
    if nucleic_acid_functional:
        classes.add(CLASS_NUCLEIC_ACID)

    return {
        "entry_id": record.entry_id,
        "classes": sorted(classes),  # functional-tier class labels
        "ambiguous_classes": sorted(ambiguous_classes - classes),
        "metals": sorted(functional_metals),
        "small_molecules": sorted(functional_sms),
        "has_nucleic_acid": has_nucleic_acid,
        "tiers": dict(sorted(tiers.items())),
        "purification_artifact": is_artifact_entry,
    }
