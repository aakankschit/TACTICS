"""
Enumeration utilities and result dataclasses.

This module provides:
- EnumerationResult: Result from product enumeration
- EnumerationError: Details about failed enumeration
- AutoDetectionResult: Results from SMARTS compatibility detection
- Utility functions for reagent processing
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Literal, Optional, Set, Tuple

import polars as pl
from rdkit import Chem


@dataclass
class EnumerationError:
    """
    Details about a failed enumeration.

    Attributes:
        step_index: Which step failed
        pattern_id: Which pattern was attempted
        error_type: Type of error
        message: Human-readable description
        reagent_smiles: Input SMILES that failed
        reagent_names: Input reagent names for traceability
        reaction_smarts: The SMARTS pattern that was attempted (if available)
    """

    step_index: int
    pattern_id: Optional[str]
    error_type: Literal[
        "no_compatible_pattern",
        "reaction_failed",
        "invalid_input",
        "deprotection_failed",
    ]
    message: str
    reagent_smiles: List[str] = field(default_factory=list)
    reagent_names: List[str] = field(default_factory=list)
    reaction_smarts: Optional[str] = None

    def __str__(self) -> str:
        return f"Step {self.step_index}: {self.error_type} - {self.message}"

    def __repr__(self) -> str:
        return self.detailed_message()

    def detailed_message(self) -> str:
        """
        Generate a detailed error message for debugging.

        Returns:
            Multi-line string with full error details
        """
        lines = [
            f"EnumerationError:",
            f"  Step: {self.step_index}",
            f"  Pattern: {self.pattern_id}",
            f"  Error Type: {self.error_type}",
            f"  Message: {self.message}",
        ]

        if self.reaction_smarts:
            # Truncate long SMARTS for display
            smarts_display = self.reaction_smarts
            if len(smarts_display) > 60:
                smarts_display = smarts_display[:57] + "..."
            lines.append(f"  Reaction SMARTS: {smarts_display}")

        if self.reagent_names or self.reagent_smiles:
            lines.append("  Reagents:")
            # Pair names with SMILES for display
            num_reagents = max(len(self.reagent_names), len(self.reagent_smiles))
            for i in range(num_reagents):
                name = (
                    self.reagent_names[i] if i < len(self.reagent_names) else "unknown"
                )
                smiles = (
                    self.reagent_smiles[i] if i < len(self.reagent_smiles) else "N/A"
                )
                lines.append(f"    [{i}] {name}: {smiles}")

        return "\n".join(lines)

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert error to dictionary for DataFrame export.

        Returns:
            Dictionary with all error fields
        """
        return {
            "step_index": self.step_index,
            "pattern_id": self.pattern_id,
            "error_type": self.error_type,
            "message": self.message,
            "reagent_smiles": "; ".join(self.reagent_smiles),
            "reagent_names": "; ".join(self.reagent_names),
            "reaction_smarts": self.reaction_smarts,
        }


@dataclass
class EnumerationResult:
    """
    Complete result from pipeline enumeration.

    Attributes:
        product: Final product molecule (None if failed)
        product_smiles: SMILES of final product
        product_name: Name of product (from reagent keys)
        patterns_used: {step_index: pattern_id} for each step
        intermediates: {step_index: mol} if store_intermediates=True
        error: EnumerationError if failed
    """

    product: Optional[Chem.Mol] = None
    product_smiles: Optional[str] = None
    product_name: Optional[str] = None
    patterns_used: Dict[int, str] = field(default_factory=dict)
    intermediates: Optional[Dict[int, Chem.Mol]] = None
    error: Optional[EnumerationError] = None

    @property
    def success(self) -> bool:
        """True if enumeration succeeded."""
        return self.product is not None and self.error is None

    def __repr__(self) -> str:
        if self.success:
            return f"EnumerationResult(success=True, smiles='{self.product_smiles}')"
        else:
            # Show detailed error information
            error = self.error
            return (
                f"EnumerationResult(success=False, product='{self.product_name}')\n"
                f"{error.detailed_message()}"
            )


@dataclass
class AutoDetectionResult:
    """
    Results from automatic SMARTS compatibility detection.

    Attributes:
        pattern_results: {step_index: {pattern_id: ValidationResult}}
        compatibility_map: {(step_index, position, reagent_key): {pattern_ids}}
        coverage_by_pattern: {step_index: {pattern_id: {position: coverage%}}}
        unmatched_reagents: {position: [reagent_names]} - Reagents matching no patterns
        warnings: List of warning messages
    """

    pattern_results: Dict[int, Dict[str, Any]] = field(default_factory=dict)
    compatibility_map: Dict[Tuple[int, int, str], Set[str]] = field(
        default_factory=dict
    )
    coverage_by_pattern: Dict[int, Dict[str, Dict[int, float]]] = field(
        default_factory=dict
    )
    unmatched_reagents: Dict[int, List[str]] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)


# === Utility Functions ===


def read_reagent_file(file_path: str) -> List[Tuple[str, str]]:
    """
    Read reagents from a SMILES file.

    Args:
        file_path: Path to .smi file with "SMILES name" format

    Returns:
        List of (smiles, name) tuples
    """
    reagents = []
    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) >= 2:
                smiles = parts[0]
                name = parts[1]
                reagents.append((smiles, name))
            elif len(parts) == 1:
                smiles = parts[0]
                reagents.append((smiles, smiles))
    return reagents


def results_to_dataframe(
    results: List[EnumerationResult],
    include_failures: bool = False,
) -> pl.DataFrame:
    """
    Convert enumeration results to Polars DataFrame.

    Args:
        results: List of EnumerationResult
        include_failures: Include failed enumerations

    Returns:
        DataFrame with columns: product_name, product_smiles, success
    """
    data = []
    for result in results:
        if result.success or include_failures:
            data.append(
                {
                    "product_name": result.product_name or "",
                    "product_smiles": result.product_smiles or "",
                    "success": result.success,
                }
            )

    if not data:
        return pl.DataFrame({"product_name": [], "product_smiles": [], "success": []})

    return pl.DataFrame(data)


def failures_to_dataframe(
    results: List[EnumerationResult],
) -> pl.DataFrame:
    """
    Convert failed enumeration results to a detailed Polars DataFrame for analysis.

    Args:
        results: List of EnumerationResult (will filter to failures only)

    Returns:
        DataFrame with columns:
            - product_name: Attempted product name
            - step_index: Which step failed
            - pattern_id: Which pattern was attempted
            - error_type: Type of error
            - message: Error message
            - reagent_names: Semicolon-separated reagent names
            - reagent_smiles: Semicolon-separated reagent SMILES
            - reaction_smarts: The SMARTS pattern attempted (if available)

    Example:
        >>> results = pipeline.enumerate_library()
        >>> failures_df = failures_to_dataframe(results)
        >>> print(failures_df.group_by("error_type").count())
    """
    failure_data = []

    for result in results:
        if not result.success and result.error:
            error = result.error
            failure_data.append(
                error.to_dict() | {"product_name": result.product_name or ""}
            )

    if not failure_data:
        return pl.DataFrame(
            {
                "product_name": [],
                "step_index": [],
                "pattern_id": [],
                "error_type": [],
                "message": [],
                "reagent_names": [],
                "reagent_smiles": [],
                "reaction_smarts": [],
            }
        )

    return pl.DataFrame(failure_data)


def summarize_failures(
    results: List[EnumerationResult],
) -> Dict[str, Any]:
    """
    Generate a summary of enumeration failures.

    Args:
        results: List of EnumerationResult

    Returns:
        Dictionary with failure statistics:
            - total: Total number of results
            - successes: Number of successful enumerations
            - failures: Number of failed enumerations
            - by_error_type: Count per error type
            - by_step: Count per step index
            - failure_rate: Percentage of failures

    Example:
        >>> results = pipeline.enumerate_library()
        >>> summary = summarize_failures(results)
        >>> print(f"Failure rate: {summary['failure_rate']:.1f}%")
    """
    total = len(results)
    failures = [r for r in results if not r.success]
    successes = total - len(failures)

    # Count by error type
    by_error_type: Dict[str, int] = {}
    by_step: Dict[int, int] = {}

    for result in failures:
        if result.error:
            error_type = result.error.error_type
            step_idx = result.error.step_index

            by_error_type[error_type] = by_error_type.get(error_type, 0) + 1
            by_step[step_idx] = by_step.get(step_idx, 0) + 1

    return {
        "total": total,
        "successes": successes,
        "failures": len(failures),
        "failure_rate": (len(failures) / total * 100) if total > 0 else 0.0,
        "by_error_type": by_error_type,
        "by_step": by_step,
    }
