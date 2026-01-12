"""
Default constants for SMARTS validation.

Provides commonly-used protecting groups and salt fragments for
automatic detection during reagent validation.

The protecting groups include SMARTS patterns for detection and
optional reaction SMARTS for deprotection. Salt fragments are
provided as SMILES patterns for common counterions and salts.

These defaults can be extended by providing additional entries
when creating a SMARTSValidator instance.
"""

from typing import List, Tuple

from .config import ProtectingGroupInfo

# Default protecting groups with detection and optional removal SMARTS
DEFAULT_PROTECTING_GROUPS: List[ProtectingGroupInfo] = [
    ProtectingGroupInfo(
        name="Boc",
        smarts="[NX3][C](=O)OC(C)(C)C",
        deprotection_smarts="[N:1][C](=O)OC(C)(C)C>>[N:1]",
    ),
    ProtectingGroupInfo(
        name="Fmoc",
        smarts="[NX3]C(=O)OCC1c2ccccc2-c2ccccc12",
        deprotection_smarts="[N:1]C(=O)OCC1c2ccccc2-c2ccccc12>>[N:1]",
    ),
    ProtectingGroupInfo(
        name="Cbz",
        smarts="[NX3]C(=O)OCc1ccccc1",
        deprotection_smarts="[N:1]C(=O)OCc1ccccc1>>[N:1]",
    ),
    ProtectingGroupInfo(
        name="Acetamide",
        smarts="[NX3][C](=O)[CH3]",
        deprotection_smarts="[N:1][C](=O)[CH3]>>[N:1]",
    ),
    ProtectingGroupInfo(
        name="TBS",
        smarts="[OX2][Si](C)(C)C(C)(C)C",
        deprotection_smarts="[O:1][Si](C)(C)C(C)(C)C>>[O:1]",
    ),
    ProtectingGroupInfo(
        name="O-Benzyl",
        smarts="[OX2]Cc1ccccc1",
        deprotection_smarts="[O:1]Cc1ccccc1>>[O:1]",
    ),
    ProtectingGroupInfo(
        name="Trityl",
        smarts="[NX3,OX2]C(c1ccccc1)(c1ccccc1)c1ccccc1",
        deprotection_smarts=None,  # Complex deprotection - not easily represented
    ),
    ProtectingGroupInfo(
        name="tBu-ester",
        smarts="[CX3](=O)OC(C)(C)C",
        deprotection_smarts="[C:1](=O)OC(C)(C)C>>[C:1](=O)O",
    ),
    ProtectingGroupInfo(
        name="Me-ester",
        smarts="[CX3](=O)O[CH3]",
        deprotection_smarts="[C:1](=O)O[CH3]>>[C:1](=O)O",
    ),
    ProtectingGroupInfo(
        name="Et-ester",
        smarts="[CX3](=O)OCC",
        deprotection_smarts="[C:1](=O)OCC>>[C:1](=O)O",
    ),
]

# Common salt fragments as (SMILES, name) tuples
# These are used to detect and remove counterions from multi-fragment SMILES
DEFAULT_SALT_FRAGMENTS: List[Tuple[str, str]] = [
    # Halides
    ("[Cl-]", "Chloride"),
    ("[Br-]", "Bromide"),
    ("[I-]", "Iodide"),
    ("[F-]", "Fluoride"),
    # Metal cations
    ("[Na+]", "Sodium"),
    ("[K+]", "Potassium"),
    ("[Li+]", "Lithium"),
    ("[Ca+2]", "Calcium"),
    ("[Mg+2]", "Magnesium"),
    # Nitrogen cations
    ("[NH4+]", "Ammonium"),
    # Organic acids/anions
    ("O=C(O)C(F)(F)F", "TFA (neutral)"),
    ("O=C([O-])C(F)(F)F", "Trifluoroacetate"),
    ("O=C([O-])C", "Acetate"),
    ("O=C([O-])c1ccccc1", "Benzoate"),
    # Sulfur-based
    ("O=S(=O)(O)O", "Sulfuric acid"),
    ("O=S(=O)([O-])O", "Bisulfate"),
    ("O=S(=O)([O-])[O-]", "Sulfate"),
    ("O=S(=O)([O-])c1ccc(C)cc1", "Tosylate"),
    ("O=S(=O)([O-])C(F)(F)F", "Triflate"),
    ("[O-]S(=O)(=O)C", "Mesylate"),
    # Phosphorus-based
    ("O=P([O-])([O-])O", "Phosphate"),
    ("O=P(O)(O)O", "Phosphoric acid"),
    # Other anions
    ("[O-]Cl(=O)(=O)=O", "Perchlorate"),
    ("[O-][N+](=O)c1ccccc1", "Nitrobenzene"),
    # Dicarboxylic acids
    ("O=C([O-])CC(=O)[O-]", "Malonate"),
    ("O=C([O-])C=CC(=O)[O-]", "Fumarate"),
    ("[O-]C(=O)C(O)C(O)C(=O)[O-]", "Tartrate"),
    ("O=C([O-])CC(O)(CC(=O)[O-])C(=O)[O-]", "Citrate"),
]

# Quick lookup for salt fragment SMILES (without names)
SALT_FRAGMENT_SMILES: List[str] = [smiles for smiles, _ in DEFAULT_SALT_FRAGMENTS]

# Protecting group name to info mapping for quick lookup
PROTECTING_GROUP_MAP: dict[str, ProtectingGroupInfo] = {
    pg.name: pg for pg in DEFAULT_PROTECTING_GROUPS
}


def get_protecting_group(name: str) -> ProtectingGroupInfo:
    """
    Get a protecting group by name.

    Args:
        name: Name of the protecting group (e.g., "Boc", "Fmoc")

    Returns:
        ProtectingGroupInfo object

    Raises:
        KeyError: If protecting group not found
    """
    if name not in PROTECTING_GROUP_MAP:
        available = ", ".join(PROTECTING_GROUP_MAP.keys())
        raise KeyError(f"Unknown protecting group '{name}'. Available: {available}")
    return PROTECTING_GROUP_MAP[name]


def get_all_protecting_group_names() -> List[str]:
    """Get list of all default protecting group names."""
    return list(PROTECTING_GROUP_MAP.keys())
