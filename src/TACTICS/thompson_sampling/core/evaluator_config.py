"""Pydantic configuration models for evaluators."""

from pydantic import BaseModel, Field
from typing import Literal, Optional, Callable, Any


class LookupEvaluatorConfig(BaseModel):
    """
    Configuration for LookupEvaluator.

    Looks up pre-computed scores from a CSV file. Primarily used for testing
    and benchmarking where scores are known in advance.
    """

    evaluator_type: Literal["lookup"] = "lookup"
    ref_filename: str = Field(..., description="Path to CSV file with pre-computed scores")
    compound_col: str = Field(default="Product_Code", description="Column name for compound identifiers")
    score_col: str = Field(default="Scores", description="Column name for scores in CSV")


class DBEvaluatorConfig(BaseModel):
    """
    Configuration for DBEvaluator.

    Looks up pre-computed scores from a SQLite database. Used for benchmarking
    with large datasets where database lookups are faster than CSV.
    """

    evaluator_type: Literal["db"] = "db"
    db_filename: str = Field(..., description="Path to SQLite database file")
    db_prefix: str = Field(default="", description="Prefix for database keys")


class FPEvaluatorConfig(BaseModel):
    """
    Configuration for FPEvaluator (Fingerprint Tanimoto similarity).

    Calculates Morgan fingerprint Tanimoto similarity to a reference molecule.
    """

    evaluator_type: Literal["fp"] = "fp"
    query_smiles: str = Field(..., description="SMILES string of reference molecule")


class MWEvaluatorConfig(BaseModel):
    """
    Configuration for MWEvaluator (Molecular Weight).

    Simple evaluator that calculates molecular weight. Primarily used for testing.
    """

    evaluator_type: Literal["mw"] = "mw"


class ROCSEvaluatorConfig(BaseModel):
    """
    Configuration for ROCSEvaluator (shape similarity).

    Calculates ROCS shape + chemistry overlay score to a reference molecule.
    Requires OpenEye toolkit.
    """

    evaluator_type: Literal["rocs"] = "rocs"
    query_molfile: str = Field(..., description="Path to reference molecule file (SDF/MOL)")
    max_confs: int = Field(default=50, gt=0, description="Maximum conformers to generate with Omega")


class FredEvaluatorConfig(BaseModel):
    """
    Configuration for FredEvaluator (docking score).

    Docks molecules using OpenEye FRED and returns docking scores.
    Requires OpenEye toolkit.
    """

    evaluator_type: Literal["fred"] = "fred"
    design_unit_file: str = Field(..., description="Path to OpenEye design unit file (.oedu)")
    max_confs: int = Field(default=50, gt=0, description="Maximum conformers to generate with Omega")


class MLClassifierEvaluatorConfig(BaseModel):
    """
    Configuration for MLClassifierEvaluator.

    Uses a trained scikit-learn classifier to predict activity scores from
    Morgan fingerprints.
    """

    evaluator_type: Literal["ml_classifier"] = "ml_classifier"
    model_filename: str = Field(..., description="Path to trained model file (joblib pickle)")

class CustomEvaluatorConfig(BaseModel):
    """
    Configuration for CustomEvaluator.

    Uses a user-provided callable scoring function that accepts an RDKit Mol
    and returns a float.
    """

    evaluator_type: Literal["custom"] = "custom"
    scoring_function: Callable[[Any], float] = Field(
        ...,
        description="Callable that accepts an RDKit Mol and returns a float score"
    )
