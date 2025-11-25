#!/usr/bin/env python3
"""
Extract peptide fragments from a protein PDB file.

This script extracts all consecutive residue fragments of a specified length
(e.g., dipeptides, tripeptides) from a PDB file and saves each as a separate PDB file.

Usage:
    python 1_extract_fragments.py --pdb_file protein.pdb --output_dir data/fragments --fragment_length 2
"""

import os
import argparse
import warnings
from pathlib import Path
from Bio.PDB import PDBIO, Atom, Chain, Model, PDBParser, Residue, Selection
from Bio.PDB.PDBExceptions import PDBConstructionWarning

# Suppress PDBConstructionWarning about element inference
warnings.filterwarnings('ignore', category=PDBConstructionWarning)


def create_atom_copy(atom):
    """Create a copy of an atom."""
    new_atom = Atom.Atom(
        name=atom.get_name(),
        coord=atom.get_coord(),
        bfactor=atom.get_bfactor(),
        occupancy=atom.get_occupancy(),
        serial_number=atom.get_serial_number(),
        altloc=atom.get_altloc(),
        fullname=atom.get_fullname()
    )
    return new_atom


def save_fragment(residues, fragment_index, output_dir="fragments"):
    """
    Save a peptide fragment (N consecutive residues) as a PDB file.
    
    Parameters:
    -----------
    residues : list
        List of BioPython Residue objects
    fragment_index : int
        Index for naming the fragment
    output_dir : str
        Output directory for fragments
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate filename from residue names and indices
    name_parts = []
    for res in residues:
        name_parts.append(f"{res.get_resname()}{res.get_id()[1]}")
    name_str = "_".join(name_parts)
    
    filename = f"{output_dir}/fragment_{fragment_index:03d}_{name_str}.pdb"
    
    new_model = Model.Model(0)
    new_chain = Chain.Chain('A')
    
    def clone_residue(residue):
        """Clone a residue with all its atoms."""
        new_residue = Residue.Residue(residue.get_id(), residue.get_resname(), 0)
        for atom in residue.get_atoms():
            new_atom = create_atom_copy(atom)
            new_residue.add(new_atom)
        return new_residue
    
    # Add all residues to the chain
    for res in residues:
        new_chain.add(clone_residue(res))
    
    new_model.add(new_chain)
    
    io = PDBIO()
    io.set_structure(new_model)
    io.save(filename)
    
    return filename


def extract_fragments(pdb_file, output_dir="fragments", fragment_length=2):
    """
    Extract all peptide fragments of specified length from a PDB file.
    
    Parameters:
    -----------
    pdb_file : str
        Input PDB file
    output_dir : str
        Output directory for fragments
    fragment_length : int
        Number of consecutive residues per fragment (default: 2 for dipeptides)
    
    Returns:
    --------
    fragment_files : list
        List of generated fragment file paths
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    
    fragment_files = []
    fragment_index = 0
    
    fragment_type = {2: "dipeptide", 3: "tripeptide", 4: "tetrapeptide"}.get(
        fragment_length, f"{fragment_length}-residue"
    )
    
    print(f"Extracting {fragment_type} fragments (length={fragment_length}) from {pdb_file}...")
    print(f"Output directory: {output_dir}\n")
    
    for model in structure:
        for chain in model:
            residues = [res for res in chain.get_list() if res.get_id()[0] == ' ']
            
            # Extract fragments of the specified length
            for i in range(len(residues) - fragment_length + 1):
                fragment_residues = residues[i:i + fragment_length]
                
                filename = save_fragment(fragment_residues, fragment_index, output_dir)
                fragment_files.append(filename)
                
                # Create display string
                res_names = "-".join([f"{res.get_resname()}{res.get_id()[1]}" for res in fragment_residues])
                print(f"  [{fragment_index:03d}] {res_names} → {Path(filename).name}")
                
                fragment_index += 1
    
    print(f"\n✓ Extracted {len(fragment_files)} {fragment_type} fragments")
    return fragment_files


def main():
    parser = argparse.ArgumentParser(
        description="Extract peptide fragments from a protein PDB file"
    )
    parser.add_argument(
        "--pdb_file",
        type=str,
        required=True,
        help="Input PDB file"
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="data/fragments",
        help="Output directory for fragments (default: data/fragments)"
    )
    parser.add_argument(
        "--fragment_length",
        type=int,
        default=2,
        help="Number of consecutive residues per fragment (default: 2 for dipeptides, 3 for tripeptides, etc.)"
    )
    
    args = parser.parse_args()
    
    if not os.path.exists(args.pdb_file):
        print(f"Error: PDB file not found: {args.pdb_file}")
        return
    
    if args.fragment_length < 1:
        print(f"Error: fragment_length must be at least 1, got {args.fragment_length}")
        return
    
    extract_fragments(args.pdb_file, args.output_dir, args.fragment_length)


if __name__ == "__main__":
    main()

