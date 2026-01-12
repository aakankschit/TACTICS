"""
File writing utilities for enumerated chemical libraries.

This module provides functions to write enumeration results to various formats:
- CSV: Standard comma-separated values
- SMI: SMILES format (SMILES name)
- SDF: Structure-Data File format
"""

import os
from pathlib import Path
from typing import List, Literal, Optional

import polars as pl
from rdkit import Chem
from rdkit.Chem import SDWriter

from .enumeration_utils import EnumerationResult, results_to_dataframe


def write_enumerated_library(
    results: List[EnumerationResult],
    output_path: str,
    format: Literal["csv", "smi", "sdf"] = "csv",
    include_failures: bool = False,
) -> int:
    """
    Write enumerated library to file.

    Args:
        results: List of EnumerationResult from enumeration
        output_path: Output file path
        format: Output format (csv, smi, sdf)
        include_failures: Include failed enumerations in output

    Returns:
        Number of products written

    Example:
        >>> results = pipeline.enumerate_library()
        >>> n_written = write_enumerated_library(results, "products.csv", format="csv")
        >>> print(f"Wrote {n_written} products")
    """
    # Ensure output directory exists
    output_dir = os.path.dirname(output_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    # Filter results
    filtered = [r for r in results if r.success or include_failures]

    if not filtered:
        return 0

    if format == "csv":
        return _write_csv(filtered, output_path, include_failures)
    elif format == "smi":
        return _write_smi(filtered, output_path)
    elif format == "sdf":
        return _write_sdf(filtered, output_path)
    else:
        raise ValueError(f"Unsupported format: {format}")


def _write_csv(
    results: List[EnumerationResult],
    output_path: str,
    include_failures: bool,
) -> int:
    """Write results to CSV format."""
    df = results_to_dataframe(results, include_failures)
    df.write_csv(output_path)
    return len(df)


def _write_smi(results: List[EnumerationResult], output_path: str) -> int:
    """Write results to SMILES format."""
    count = 0
    with open(output_path, "w") as f:
        for result in results:
            if result.success and result.product_smiles:
                name = result.product_name or f"product_{count}"
                f.write(f"{result.product_smiles} {name}\n")
                count += 1
    return count


def _write_sdf(results: List[EnumerationResult], output_path: str) -> int:
    """Write results to SDF format."""
    count = 0
    writer = SDWriter(output_path)

    for result in results:
        if result.success and result.product:
            mol = result.product
            if result.product_name:
                mol.SetProp("_Name", result.product_name)
            if result.product_smiles:
                mol.SetProp("SMILES", result.product_smiles)
            writer.write(mol)
            count += 1

    writer.close()
    return count


def write_products_chunked(
    results: List[EnumerationResult],
    output_dir: str,
    products_per_file: int = 5000,
    format: Literal["csv", "smi"] = "smi",
    file_prefix: str = "products",
) -> int:
    """
    Write products to multiple files with a specified number per file.

    Args:
        results: List of EnumerationResult
        output_dir: Output directory path
        products_per_file: Maximum products per file
        format: Output format (csv or smi)
        file_prefix: Prefix for output file names

    Returns:
        Total number of products written
    """
    os.makedirs(output_dir, exist_ok=True)

    # Filter successful results
    successful = [r for r in results if r.success]

    if not successful:
        return 0

    total_written = 0
    num_files = (len(successful) + products_per_file - 1) // products_per_file

    for i in range(num_files):
        start_idx = i * products_per_file
        end_idx = min((i + 1) * products_per_file, len(successful))
        chunk = successful[start_idx:end_idx]

        extension = "csv" if format == "csv" else "smi"
        output_file = os.path.join(output_dir, f"{file_prefix}_{i + 1}.{extension}")

        if format == "csv":
            _write_csv(chunk, output_file, include_failures=False)
        else:
            _write_smi(chunk, output_file)

        total_written += len(chunk)

    return total_written


def write_products_to_files(
    df: pl.DataFrame,
    output_dir: str,
    products_per_file: int = 5000,
) -> None:
    """
    Write DataFrame of products to SMILES files.

    Legacy function for compatibility. Writes products in chunks.

    Args:
        df: DataFrame with 'Product_SMILES' and 'Product_Name' columns
        output_dir: Output directory
        products_per_file: Products per output file
    """
    os.makedirs(output_dir, exist_ok=True)

    num_files = (len(df) + products_per_file - 1) // products_per_file

    for i in range(num_files):
        start_idx = i * products_per_file
        end_idx = min((i + 1) * products_per_file, len(df))
        chunk = df[start_idx:end_idx]

        # Write to temporary CSV without header
        temp_output_file = os.path.join(output_dir, f"temp_products_{i + 1}.csv")
        chunk.write_csv(temp_output_file, include_header=False)

        # Convert commas to spaces for SMILES format
        with open(temp_output_file, "r") as temp_file:
            content = temp_file.read().replace(",", " ")

        # Write final SMILES file
        output_file = os.path.join(output_dir, f"products_{i + 1}.smi")
        with open(output_file, "w") as final_file:
            final_file.write(content)

        # Clean up temp file
        os.remove(temp_output_file)
