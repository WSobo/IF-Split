"""IF-Split: a reproducible, ligand-aware train/val/test splitter for the PDB.

Replicates the LigandMPNN split methodology (Dauparas et al., Nature Methods
2025) on demand from a current PDB snapshot, with a manifest + lock file for
byte-for-byte reproduction. See PLAN.md for the full spec.
"""

__version__ = "0.3.0"
