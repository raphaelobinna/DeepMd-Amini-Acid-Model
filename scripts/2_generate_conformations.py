#!/usr/bin/env python3
"""
Generate conformations for each fragment and save each in a separate PDB file.

This script:
1. Reads fragment PDB files
2. Generates multiple conformations using RDKit
3. Validates conformations (filters atoms too close)
4. Saves each conformation as a separate PDB file

Usage:
    python 2_generate_conformations.py --fragments_dir data/fragments --output_dir data/conformations
"""

import os
import sys
import numpy as np
from pathlib import Path
from scipy.spatial.distance import cdist
from rdkit import Chem
from rdkit.Chem import rdForceFieldHelpers, rdDistGeom


def load_pdb_file(pdb_file):
    """Load PDB file and return RDKit molecule with hydrogens."""
    with open(pdb_file, "r") as f:
        pdb_block = f.read()
    
    mol = Chem.MolFromPDBBlock(pdb_block, removeHs=False)
    return Chem.AddHs(mol, addCoords=True)


def generate_separate_conformations(
    pdb_file,
    output_dir,
    num_conf=50,
    min_distance=0.5
):
    """
    Generate conformations and save each in a separate PDB file.
    Each fragment gets its own folder containing all its conformations.
    
    Parameters:
    -----------
    pdb_file : str
        Input PDB file
    output_dir : str
        Output directory for conformations
    num_conf : int
        Number of conformations to generate
    min_distance : float
        Minimum allowed interatomic distance in Angstroms
    
    Returns:
    --------
    n_valid : int
        Number of valid conformations saved
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    fragment_name = Path(pdb_file).stem
    
    # Create a folder for this fragment
    fragment_dir = output_path / fragment_name
    fragment_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"  Loading {Path(pdb_file).name}...")
    mol = load_pdb_file(pdb_file)
    if mol is None:
        print(f"    ✗ Error: Could not load {pdb_file}")
        return 0
    
    print(f"    Atoms: {mol.GetNumAtoms()}")
    
    # Generate conformations
    conf_ids = rdDistGeom.EmbedMultipleConfs(
        mol, 
        num_conf, 
        params=rdDistGeom.ETKDGv3()
    )
    
    print(f"    Generated {len(conf_ids)} conformations")
    
    # Optimize with MMFF
    rdForceFieldHelpers.MMFFOptimizeMoleculeConfs(mol, numThreads=0)
    
    # Validate and save each conformation separately
    valid_count = 0
    invalid_count = 0
    
    for conf_id in conf_ids:
        conf = mol.GetConformer(conf_id)
        coords = np.array([conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())])
        
        # Check minimum distance
        is_valid = True
        if len(coords) > 1:
            dists = cdist(coords, coords)
            np.fill_diagonal(dists, np.inf)
            min_dist = dists.min()
            
            if min_dist < min_distance:
                is_valid = False
                invalid_count += 1
                continue
        
        # Save valid conformation to separate file in fragment's folder
        output_file = fragment_dir / f"{fragment_name}_conf_{conf_id:03d}.pdb"
        writer = Chem.PDBWriter(str(output_file))
        writer.write(mol, confId=conf_id)
        writer.close()
        
        valid_count += 1
    
    print(f"    ✓ Saved {valid_count} valid conformations")
    if invalid_count > 0:
        print(f"    ⚠ Filtered {invalid_count} invalid conformations")
    
    return valid_count


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Generate conformations for fragment PDB files"
    )
    parser.add_argument(
        "--fragments_dir",
        type=str,
        default="data/fragments",
        help="Directory containing fragment PDB files (default: data/fragments)"
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="data/conformations",
        help="Output directory for conformations (default: data/conformations)"
    )
    parser.add_argument(
        "--num_conf",
        type=int,
        default=50,
        help="Number of conformations per fragment (default: 50)"
    )
    parser.add_argument(
        "--min_distance",
        type=float,
        default=0.5,
        help="Minimum allowed interatomic distance in Angstroms (default: 0.5)"
    )
    
    args = parser.parse_args()
    
    fragments_path = Path(args.fragments_dir)
    if not fragments_path.exists():
        print(f"Error: Fragments directory not found: {fragments_path}")
        return
    
    # Find all PDB files
    pdb_files = sorted(fragments_path.glob("*.pdb"))
    if not pdb_files:
        print(f"Error: No PDB files found in {fragments_path}")
        return
    
    print(f"Generating conformations for {len(pdb_files)} fragments...")
    print(f"Output directory: {args.output_dir}\n")
    
    total_valid = 0
    for i, pdb_file in enumerate(pdb_files, 1):
        print(f"[{i}/{len(pdb_files)}] {pdb_file.name}")
        n_valid = generate_separate_conformations(
            str(pdb_file),
            args.output_dir,
            args.num_conf,
            args.min_distance
        )
        total_valid += n_valid
        print()
    
    print(f"✓ Complete! Generated {total_valid} total conformations")
    print(f"  Average: {total_valid/len(pdb_files):.1f} conformations per fragment")


if __name__ == "__main__":
    main()

