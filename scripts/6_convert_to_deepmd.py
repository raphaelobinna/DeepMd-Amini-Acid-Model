#!/usr/bin/env python3
"""
Convert ABACUS output files to DeePMD-kit format for training.

This script:
1. Finds all ABACUS output directories (OUT.ABACUS)
2. Converts them to DeePMD .npy format using dpdata
3. Organizes data into training/validation sets

Usage:
    python 6_convert_to_deepmd.py --abacus_dir data/abacus_inputs --output_dir data/deepmd_data
"""

import os
import sys
import shutil
from pathlib import Path
import numpy as np

try:
    import dpdata
except ImportError:
    print("Error: dpdata is not installed.")
    print("Install it with: conda install -c conda-forge dpdata")
    sys.exit(1)


def find_abacus_outputs(base_dir):
    """Find all ABACUS output directories."""
    outputs = []
    base_path = Path(base_dir)
    
    # Look for OUT.ABACUS directories
    for out_dir in base_path.rglob("OUT.ABACUS"):
        # Check if it has MD_dump (required for conversion)
        if (out_dir / "MD_dump").exists():
            outputs.append(out_dir)
    
    return sorted(outputs)


def convert_abacus_to_deepmd(abacus_dir, output_dir):
    """
    Convert a single ABACUS output directory to DeePMD format.
    
    Parameters:
    -----------
    abacus_dir : Path
        Path to ABACUS OUT.ABACUS directory
    output_dir : Path
        Path to output directory for DeePMD data
    """
    try:
        # dpdata needs the parent directory (containing INPUT, STRU, etc.)
        # not just OUT.ABACUS
        parent_dir = abacus_dir.parent
        
        # Check if MD_dump exists
        if not (abacus_dir / "MD_dump").exists():
            print(f"  ⚠ Skipping {parent_dir.name}: No MD_dump file")
            return False
        
        # Load ABACUS MD data
        print(f"  Converting {parent_dir.name}...")
        try:
            system = dpdata.LabeledSystem(str(parent_dir), fmt='abacus/md')
        except Exception as e:
            print(f"    ✗ Error loading: {e}")
            return False
        
        # Check if system has data
        if len(system) == 0:
            print(f"    ⚠ Warning: No frames found")
            return False
        
        # Create output directory
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Convert to DeePMD format
        system.to_deepmd_npy(str(output_dir))
        
        print(f"    ✓ Converted {len(system)} frames")
        return True
        
    except Exception as e:
        print(f"    ✗ Error converting {abacus_dir.parent.name}: {e}")
        return False


def split_data(data_dirs, train_ratio=0.8, output_base="data/deepmd_data"):
    """
    Split converted data into training and validation sets.
    
    Parameters:
    -----------
    data_dirs : list
        List of directories containing DeePMD .npy data
    train_ratio : float
        Ratio of data to use for training (default 0.8)
    output_base : str
        Base directory for output
    """
    output_base = Path(output_base)
    train_dir = output_base / "train"
    val_dir = output_base / "val"
    
    train_dir.mkdir(parents=True, exist_ok=True)
    val_dir.mkdir(parents=True, exist_ok=True)
    
    # Shuffle and split
    np.random.seed(42)  # For reproducibility
    indices = np.arange(len(data_dirs))
    np.random.shuffle(indices)
    
    split_idx = int(len(data_dirs) * train_ratio)
    train_indices = indices[:split_idx]
    val_indices = indices[split_idx:]
    
    # Copy training data
    print(f"\nCopying {len(train_indices)} datasets to training set...")
    for idx in train_indices:
        src = Path(data_dirs[idx])
        dst = train_dir / src.name
        if dst.exists():
            shutil.rmtree(dst)
        shutil.copytree(src, dst)
    
    # Copy validation data
    print(f"Copying {len(val_indices)} datasets to validation set...")
    for idx in val_indices:
        src = Path(data_dirs[idx])
        dst = val_dir / src.name
        if dst.exists():
            shutil.rmtree(dst)
        shutil.copytree(src, dst)
    
    print(f"\n✓ Data split complete:")
    print(f"  Training: {train_dir} ({len(train_indices)} datasets)")
    print(f"  Validation: {val_dir} ({len(val_indices)} datasets)")


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Convert ABACUS outputs to DeePMD-kit format"
    )
    parser.add_argument(
        "--abacus_dir",
        type=str,
        default="data/abacus_inputs",
        help="Directory containing ABACUS calculations (default: data/abacus_inputs)"
    )
    parser.add_argument(
        "--output",
        type=str,
        default="data/deepmd_data",
        help="Output directory for DeePMD data (default: data/deepmd_data)"
    )
    parser.add_argument(
        "--train-ratio",
        type=float,
        default=0.8,
        help="Ratio of data for training (default: 0.8)"
    )
    parser.add_argument(
        "--no-split",
        action="store_true",
        help="Don't split into train/val, just convert all data"
    )
    parser.add_argument(
        "--max-fragments",
        type=int,
        default=None,
        help="Maximum number of fragments to process (default: all)"
    )
    
    args = parser.parse_args()
    
    abacus_base = Path(args.abacus_dir)
    output_base = Path(args.output)
    
    if not abacus_base.exists():
        print(f"Error: Directory {abacus_base} does not exist")
        sys.exit(1)
    
    print("=" * 60)
    print("ABACUS to DeePMD Converter")
    print("=" * 60)
    print(f"ABACUS directory: {abacus_base}")
    print(f"Output directory: {output_base}")
    print("=" * 60)
    print()
    
    # Find all ABACUS outputs
    print("Searching for ABACUS output directories...")
    abacus_outputs = find_abacus_outputs(abacus_base)
    
    if len(abacus_outputs) == 0:
        print("No ABACUS output directories with MD_dump found!")
        sys.exit(1)
    
    total_found = len(abacus_outputs)
    
    # Limit to max_fragments if specified
    if args.max_fragments is not None:
        if args.max_fragments > 0:
            abacus_outputs = abacus_outputs[:args.max_fragments]
            print(f"Found {total_found} ABACUS output directories")
            print(f"Processing first {len(abacus_outputs)} (limited by --max-fragments)")
        else:
            print("Error: --max-fragments must be greater than 0")
            sys.exit(1)
    else:
        print(f"Found {len(abacus_outputs)} ABACUS output directories")
    
    print()
    
    # Convert each output
    converted_dirs = []
    temp_output = output_base / "converted"
    temp_output.mkdir(parents=True, exist_ok=True)
    
    print("Converting ABACUS outputs to DeePMD format...")
    for i, abacus_dir in enumerate(abacus_outputs, 1):
        dataset_name = f"conf_{abacus_dir.parent.name}"
        output_dir = temp_output / dataset_name
        
        if convert_abacus_to_deepmd(abacus_dir, output_dir):
            converted_dirs.append(output_dir)
    
    print()
    print(f"✓ Converted {len(converted_dirs)} datasets")
    
    if len(converted_dirs) == 0:
        print("No data was successfully converted!")
        sys.exit(1)
    
    # Split into train/val if requested
    if not args.no_split:
        print()
        print("Splitting data into training and validation sets...")
        split_data(converted_dirs, args.train_ratio, output_base)
        
        # Remove temporary directory
        if temp_output.exists():
            shutil.rmtree(temp_output)
    else:
        # Just rename the converted directory
        final_dir = output_base / "all_data"
        if final_dir.exists():
            shutil.rmtree(final_dir)
        temp_output.rename(final_dir)
        print(f"\n✓ All data saved to: {final_dir}")
    
    print()
    print("=" * 60)
    print("Conversion complete!")
    print("=" * 60)
    print(f"\nNext steps:")
    print(f"1. Review the converted data in: {output_base}")
    if not args.no_split:
        print(f"2. Training data: {output_base}/train")
        print(f"3. Validation data: {output_base}/val")
    print(f"4. Create input.json for DeePMD training")
    print(f"5. Run: dp train input.json")


if __name__ == "__main__":
    main()

