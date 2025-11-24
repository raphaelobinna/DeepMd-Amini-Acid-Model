#!/usr/bin/env python3
"""
Extract dipeptide fragments from a protein PDB file.

This script extracts all consecutive residue pairs (dipeptides) from a PDB file
and saves each as a separate PDB file.

Usage:
    python 1_extract_fragments.py --pdb_file protein.pdb --output_dir data/fragments
"""

import os
import argparse
from pathlib import Path
from Bio.PDB import PDBIO, Atom, Chain, Model, PDBParser, Residue, Selection


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


def save_dipeptide_fragment(res1, res2, fragment_index, output_dir="fragments"):
    """
    Save a dipeptide fragment (two consecutive residues) as a PDB file.
    
    Parameters:
    -----------
    res1 : BioPython Residue
        First residue
    res2 : BioPython Residue
        Second residue
    fragment_index : int
        Index for naming the fragment
    output_dir : str
        Output directory for fragments
    """
    os.makedirs(output_dir, exist_ok=True)
    
    name_1 = res1.get_resname()
    name_2 = res2.get_resname()
    idx1 = res1.get_id()[1]
    idx2 = res2.get_id()[1]
    
    filename = f"{output_dir}/fragment_{fragment_index:03d}_{name_1}{idx1}_{name_2}{idx2}.pdb"
    
    new_model = Model.Model(0)
    new_chain = Chain.Chain('A')
    
    def clone_residue(residue):
        """Clone a residue with all its atoms."""
        new_residue = Residue.Residue(residue.get_id(), residue.get_resname(), 0)
        for atom in residue.get_atoms():
            new_atom = create_atom_copy(atom)
            new_residue.add(new_atom)
        return new_residue
    
    new_chain.add(clone_residue(res1))
    new_chain.add(clone_residue(res2))
    
    new_model.add(new_chain)
    
    io = PDBIO()
    io.set_structure(new_model)
    io.save(filename)
    
    return filename


def extract_fragments(pdb_file, output_dir="fragments"):
    """
    Extract all dipeptide fragments from a PDB file.
    
    Parameters:
    -----------
    pdb_file : str
        Input PDB file
    output_dir : str
        Output directory for fragments
    
    Returns:
    --------
    fragment_files : list
        List of generated fragment file paths
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    
    fragment_files = []
    fragment_index = 0
    
    print(f"Extracting dipeptide fragments from {pdb_file}...")
    print(f"Output directory: {output_dir}\n")
    
    for model in structure:
        for chain in model:
            residues = [res for res in chain.get_list() if res.get_id()[0] == ' ']
            
            for i in range(len(residues) - 1):
                res_1 = residues[i]
                res_2 = residues[i + 1]
                
                filename = save_dipeptide_fragment(res_1, res_2, fragment_index, output_dir)
                fragment_files.append(filename)
                
                print(f"  [{fragment_index:03d}] {res_1.get_resname()}{res_1.get_id()[1]}-{res_2.get_resname()}{res_2.get_id()[1]} → {Path(filename).name}")
                
                fragment_index += 1
    
    print(f"\n✓ Extracted {len(fragment_files)} dipeptide fragments")
    return fragment_files


def main():
    parser = argparse.ArgumentParser(
        description="Extract dipeptide fragments from a protein PDB file"
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
    
    args = parser.parse_args()
    
    if not os.path.exists(args.pdb_file):
        print(f"Error: PDB file not found: {args.pdb_file}")
        return
    
    extract_fragments(args.pdb_file, args.output_dir)


if __name__ == "__main__":
    main()

