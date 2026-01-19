#!/usr/bin/env python3
"""
Generate all 20 Ace-X-NMe dipeptides sampled across the Ramachandran space.

This script:
1. Builds acetyl-capped dipeptides (Ace-X-NMe) for all 20 amino acids
2. Systematically samples the Ramachandran space (phi/psi backbone dihedrals)
3. Optionally adds explicit water solvation
4. Outputs PDB files ready for ABACUS DFT calculations

The Ramachandran space sampling ensures the model learns:
- Unique potential energy surfaces for every side chain
- Peptide backbone interactions across all conformational states
- Alpha-helix, beta-sheet, and other secondary structure regions

Usage:
    python 0_generate_dipeptides.py --output_dir data/dipeptides --phi_step 30 --psi_step 30
    python 0_generate_dipeptides.py --output_dir data/dipeptides --solvate --n_waters 50
"""

import os
import sys
import argparse
import numpy as np
from pathlib import Path
from typing import List, Tuple, Optional, Dict
import warnings

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolTransforms, rdForceFieldHelpers
    from rdkit.Geometry import Point3D
except ImportError:
    print("Error: RDKit is not installed.")
    print("Install with: conda install -c conda-forge rdkit")
    sys.exit(1)

# Suppress RDKit warnings
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


# ============================================================================
# AMINO ACID DEFINITIONS
# ============================================================================

# All 20 standard amino acids with their properties
AMINO_ACIDS = {
    'ALA': {'name': 'Alanine', 'code': 'A', 'smiles': 'C'},
    'ARG': {'name': 'Arginine', 'code': 'R', 'smiles': 'CCCNC(=N)N'},
    'ASN': {'name': 'Asparagine', 'code': 'N', 'smiles': 'CC(=O)N'},
    'ASP': {'name': 'Aspartic acid', 'code': 'D', 'smiles': 'CC(=O)O'},
    'CYS': {'name': 'Cysteine', 'code': 'C', 'smiles': 'CS'},
    'GLN': {'name': 'Glutamine', 'code': 'Q', 'smiles': 'CCC(=O)N'},
    'GLU': {'name': 'Glutamic acid', 'code': 'E', 'smiles': 'CCC(=O)O'},
    'GLY': {'name': 'Glycine', 'code': 'G', 'smiles': '[H]'},
    'HIS': {'name': 'Histidine', 'code': 'H', 'smiles': 'CC1=CNC=N1'},
    'ILE': {'name': 'Isoleucine', 'code': 'I', 'smiles': 'C(C)CC'},
    'LEU': {'name': 'Leucine', 'code': 'L', 'smiles': 'CC(C)C'},
    'LYS': {'name': 'Lysine', 'code': 'K', 'smiles': 'CCCCN'},
    'MET': {'name': 'Methionine', 'code': 'M', 'smiles': 'CCSC'},
    'PHE': {'name': 'Phenylalanine', 'code': 'F', 'smiles': 'Cc1ccccc1'},
    'PRO': {'name': 'Proline', 'code': 'P', 'smiles': None},  # Special - cyclic
    'SER': {'name': 'Serine', 'code': 'S', 'smiles': 'CO'},
    'THR': {'name': 'Threonine', 'code': 'T', 'smiles': 'C(C)O'},
    'TRP': {'name': 'Tryptophan', 'code': 'W', 'smiles': 'Cc1c[nH]c2ccccc12'},
    'TYR': {'name': 'Tyrosine', 'code': 'Y', 'smiles': 'Cc1ccc(O)cc1'},
    'VAL': {'name': 'Valine', 'code': 'V', 'smiles': 'C(C)C'},
}

# Full SMILES for Ace-X-NMe dipeptides (acetyl-amino_acid-N-methylamide)
# Format: Ace-[AA]-NMe where Ace = CH3CO- and NMe = -NHCH3
ACE_X_NME_SMILES = {
    'ALA': 'CC(=O)NC(C)C(=O)NC',                           # Ace-Ala-NMe
    'ARG': 'CC(=O)NC(CCCNC(=N)N)C(=O)NC',                  # Ace-Arg-NMe
    'ASN': 'CC(=O)NC(CC(=O)N)C(=O)NC',                     # Ace-Asn-NMe
    'ASP': 'CC(=O)NC(CC(=O)O)C(=O)NC',                     # Ace-Asp-NMe
    'CYS': 'CC(=O)NC(CS)C(=O)NC',                          # Ace-Cys-NMe
    'GLN': 'CC(=O)NC(CCC(=O)N)C(=O)NC',                    # Ace-Gln-NMe
    'GLU': 'CC(=O)NC(CCC(=O)O)C(=O)NC',                    # Ace-Glu-NMe
    'GLY': 'CC(=O)NCC(=O)NC',                              # Ace-Gly-NMe
    'HIS': 'CC(=O)NC(Cc1c[nH]cn1)C(=O)NC',                 # Ace-His-NMe
    'ILE': 'CC(=O)NC(C(C)CC)C(=O)NC',                      # Ace-Ile-NMe
    'LEU': 'CC(=O)NC(CC(C)C)C(=O)NC',                      # Ace-Leu-NMe
    'LYS': 'CC(=O)NC(CCCCN)C(=O)NC',                       # Ace-Lys-NMe
    'MET': 'CC(=O)NC(CCSC)C(=O)NC',                        # Ace-Met-NMe
    'PHE': 'CC(=O)NC(Cc1ccccc1)C(=O)NC',                   # Ace-Phe-NMe
    'PRO': 'CC(=O)N1CCCC1C(=O)NC',                         # Ace-Pro-NMe (cyclic)
    'SER': 'CC(=O)NC(CO)C(=O)NC',                          # Ace-Ser-NMe
    'THR': 'CC(=O)NC(C(C)O)C(=O)NC',                       # Ace-Thr-NMe
    'TRP': 'CC(=O)NC(Cc1c[nH]c2ccccc12)C(=O)NC',           # Ace-Trp-NMe
    'TYR': 'CC(=O)NC(Cc1ccc(O)cc1)C(=O)NC',                # Ace-Tyr-NMe
    'VAL': 'CC(=O)NC(C(C)C)C(=O)NC',                       # Ace-Val-NMe
}


# ============================================================================
# RAMACHANDRAN SPACE DEFINITIONS
# ============================================================================

# Key regions in Ramachandran space (phi, psi in degrees)
RAMACHANDRAN_REGIONS = {
    'alpha_helix_R': (-60, -45),      # Right-handed alpha helix
    'alpha_helix_L': (60, 45),        # Left-handed alpha helix
    'beta_sheet': (-120, 120),        # Beta sheet (extended)
    'beta_sheet_alt': (-140, 135),    # Alternative beta
    'polyproline_II': (-75, 145),     # Polyproline II helix
    'collagen': (-75, 165),           # Collagen helix
    'turn_I': (-60, -30),             # Type I turn
    'turn_II': (-60, 120),            # Type II turn
}


# ============================================================================
# DIPEPTIDE GENERATION
# ============================================================================

def build_dipeptide(aa_code: str) -> Optional[Chem.Mol]:
    """
    Build an Ace-X-NMe dipeptide from SMILES.
    
    Parameters:
    -----------
    aa_code : str
        Three-letter amino acid code (e.g., 'ALA', 'GLY')
    
    Returns:
    --------
    mol : rdkit.Chem.Mol or None
        RDKit molecule with 3D coordinates, or None if failed
    """
    if aa_code not in ACE_X_NME_SMILES:
        print(f"  Warning: Unknown amino acid code: {aa_code}")
        return None
    
    smiles = ACE_X_NME_SMILES[aa_code]
    
    # Create molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"  Warning: Could not parse SMILES for {aa_code}")
        return None
    
    # Add hydrogens
    mol = Chem.AddHs(mol)
    
    # Generate 3D coordinates
    try:
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception as e:
        print(f"  Warning: Could not generate 3D coordinates for {aa_code}: {e}")
        return None
    
    return mol


def get_backbone_atoms(mol: Chem.Mol, aa_code: str) -> Dict[str, int]:
    """
    Find backbone atom indices for phi/psi dihedral manipulation.
    
    For Ace-X-NMe: C(Ace)-N(X)-CA(X)-C(X)-N(NMe)
    - Phi: C(Ace)-N(X)-CA(X)-C(X)
    - Psi: N(X)-CA(X)-C(X)-N(NMe)
    
    Returns dict with atom indices for phi and psi dihedrals.
    """
    # This is a simplified approach - for production, use SMARTS patterns
    # We identify atoms by connectivity patterns
    
    backbone = {}
    
    # For most amino acids, find the pattern:
    # Acetyl C=O -> N -> CA -> C=O -> N-methyl
    
    if aa_code == 'GLY':
        # Glycine: CC(=O)NCC(=O)NC
        # Atoms: 0-C, 1-C(=O), 2-O, 3-N, 4-CA, 5-C(=O), 6-O, 7-N, 8-C
        backbone['phi'] = (1, 3, 4, 5)  # C-N-CA-C
        backbone['psi'] = (3, 4, 5, 7)  # N-CA-C-N
    elif aa_code == 'PRO':
        # Proline is cyclic - phi is constrained
        backbone['phi'] = None  # Fixed due to ring
        backbone['psi'] = None  # Will need special handling
    else:
        # Standard amino acids - find by SMARTS
        # Pattern: C(=O)-N-C-C(=O)-N
        pattern = Chem.MolFromSmarts('[C;$(C(=O)N)]-[N;$(NC(=O))]-[C;!$(C=O)]-[C;$(C(=O)N)]-[N;$(NC(=O))]')
        if pattern:
            matches = mol.GetSubstructMatches(pattern)
            if matches:
                m = matches[0]
                backbone['phi'] = (m[0], m[1], m[2], m[3])
                backbone['psi'] = (m[1], m[2], m[3], m[4])
    
    return backbone


def set_dihedral(mol: Chem.Mol, atom_indices: Tuple[int, int, int, int], 
                 angle_deg: float) -> bool:
    """
    Set a dihedral angle in the molecule.
    
    Parameters:
    -----------
    mol : rdkit.Chem.Mol
        Molecule to modify
    atom_indices : tuple
        Four atom indices defining the dihedral
    angle_deg : float
        Target angle in degrees
    
    Returns:
    --------
    success : bool
    """
    try:
        conf = mol.GetConformer()
        rdMolTransforms.SetDihedralDeg(conf, *atom_indices, angle_deg)
        return True
    except Exception as e:
        return False


def generate_ramachandran_conformers(mol: Chem.Mol, aa_code: str,
                                     phi_range: Tuple[float, float] = (-180, 180),
                                     psi_range: Tuple[float, float] = (-180, 180),
                                     phi_step: float = 30,
                                     psi_step: float = 30,
                                     include_key_regions: bool = True) -> List[Tuple[Chem.Mol, float, float]]:
    """
    Generate conformers by systematically sampling Ramachandran space.
    
    Parameters:
    -----------
    mol : rdkit.Chem.Mol
        Base dipeptide molecule
    aa_code : str
        Amino acid code
    phi_range : tuple
        Range of phi angles to sample (degrees)
    psi_range : tuple
        Range of psi angles to sample (degrees)
    phi_step : float
        Step size for phi sampling (degrees)
    psi_step : float
        Step size for psi sampling (degrees)
    include_key_regions : bool
        Whether to include key Ramachandran regions regardless of grid
    
    Returns:
    --------
    conformers : list
        List of (molecule, phi, psi) tuples
    """
    conformers = []
    backbone = get_backbone_atoms(mol, aa_code)
    
    # Generate grid of phi/psi angles
    phi_angles = np.arange(phi_range[0], phi_range[1], phi_step)
    psi_angles = np.arange(psi_range[0], psi_range[1], psi_step)
    
    sampled_points = set()
    
    # Sample grid points
    for phi in phi_angles:
        for psi in psi_angles:
            sampled_points.add((phi, psi))
    
    # Add key Ramachandran regions if requested
    if include_key_regions:
        for region_name, (phi, psi) in RAMACHANDRAN_REGIONS.items():
            sampled_points.add((phi, psi))
    
    # Generate conformers
    for phi, psi in sorted(sampled_points):
        # Make a copy of the molecule
        conf_mol = Chem.Mol(mol)
        
        # Skip if proline (phi is constrained)
        if aa_code == 'PRO':
            # Proline has fixed phi around -60 degrees due to ring
            phi = -60
        
        # Try to set dihedrals
        success = True
        if backbone.get('phi'):
            success = set_dihedral(conf_mol, backbone['phi'], phi)
        if success and backbone.get('psi'):
            success = set_dihedral(conf_mol, backbone['psi'], psi)
        
        if success:
            # Quick energy minimization to relax any clashes
            try:
                AllChem.MMFFOptimizeMolecule(conf_mol, maxIters=50)
            except:
                pass
            
            conformers.append((conf_mol, phi, psi))
    
    return conformers


# ============================================================================
# SOLVATION
# ============================================================================

def add_water_shell(mol: Chem.Mol, n_waters: int = 50, 
                    min_dist: float = 2.5, max_dist: float = 6.0) -> Chem.Mol:
    """
    Add explicit water molecules around the dipeptide.
    
    Parameters:
    -----------
    mol : rdkit.Chem.Mol
        Dipeptide molecule
    n_waters : int
        Number of water molecules to add
    min_dist : float
        Minimum distance from solute (Angstroms)
    max_dist : float
        Maximum distance from solute (Angstroms)
    
    Returns:
    --------
    solvated_mol : rdkit.Chem.Mol
        Solvated molecule
    """
    # Get solute coordinates
    conf = mol.GetConformer()
    solute_coords = np.array([conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())])
    center = solute_coords.mean(axis=0)
    
    # Create water molecule
    water = Chem.MolFromSmiles('O')
    water = Chem.AddHs(water)
    AllChem.EmbedMolecule(water)
    water_conf = water.GetConformer()
    
    # Generate water positions in a shell around the solute
    waters_added = 0
    max_attempts = n_waters * 100
    attempts = 0
    
    water_coords_list = []
    
    while waters_added < n_waters and attempts < max_attempts:
        attempts += 1
        
        # Random position in spherical shell
        r = np.random.uniform(min_dist, max_dist)
        theta = np.random.uniform(0, 2 * np.pi)
        phi = np.random.uniform(0, np.pi)
        
        x = center[0] + r * np.sin(phi) * np.cos(theta)
        y = center[1] + r * np.sin(phi) * np.sin(theta)
        z = center[2] + r * np.cos(phi)
        
        new_pos = np.array([x, y, z])
        
        # Check distance from solute atoms
        dists_to_solute = np.linalg.norm(solute_coords - new_pos, axis=1)
        if dists_to_solute.min() < min_dist:
            continue
        
        # Check distance from other waters
        if water_coords_list:
            water_coords = np.array(water_coords_list)
            dists_to_waters = np.linalg.norm(water_coords - new_pos, axis=1)
            if dists_to_waters.min() < 2.0:  # Min water-water distance
                continue
        
        water_coords_list.append(new_pos)
        waters_added += 1
    
    # Build combined molecule
    combined = Chem.RWMol(mol)
    
    for water_pos in water_coords_list:
        # Add oxygen
        o_idx = combined.AddAtom(Chem.Atom(8))  # Oxygen
        
        # Add hydrogens
        h1_idx = combined.AddAtom(Chem.Atom(1))
        h2_idx = combined.AddAtom(Chem.Atom(1))
        
        # Add bonds
        combined.AddBond(o_idx, h1_idx, Chem.BondType.SINGLE)
        combined.AddBond(o_idx, h2_idx, Chem.BondType.SINGLE)
        
        # Set positions (TIP3P geometry)
        combined_conf = combined.GetConformer()
        combined_conf.SetAtomPosition(o_idx, Point3D(*water_pos))
        
        # H-O-H angle ~104.5 degrees, O-H distance ~0.9572 Angstroms
        h1_pos = water_pos + np.array([0.757, 0.587, 0.0])
        h2_pos = water_pos + np.array([-0.757, 0.587, 0.0])
        combined_conf.SetAtomPosition(h1_idx, Point3D(*h1_pos))
        combined_conf.SetAtomPosition(h2_idx, Point3D(*h2_pos))
    
    return combined.GetMol()


# ============================================================================
# FILE OUTPUT
# ============================================================================

def save_dipeptide_pdb(mol: Chem.Mol, output_file: str, aa_code: str = None):
    """Save molecule to PDB format."""
    try:
        writer = Chem.PDBWriter(str(output_file))
        writer.write(mol)
        writer.close()
        return True
    except Exception as e:
        print(f"  Error saving PDB: {e}")
        return False


def center_molecule(mol: Chem.Mol):
    """Center molecule coordinates at origin."""
    conf = mol.GetConformer()
    coords = np.array([conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())])
    center = coords.mean(axis=0)
    
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, Point3D(pos.x - center[0], 
                                         pos.y - center[1], 
                                         pos.z - center[2]))


def shift_to_positive(mol: Chem.Mol, padding: float = 5.0):
    """Shift molecule so all coordinates are positive (for ABACUS box)."""
    conf = mol.GetConformer()
    coords = np.array([conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())])
    
    min_coords = coords.min(axis=0)
    shift = -min_coords + padding
    
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, Point3D(pos.x + shift[0],
                                         pos.y + shift[1],
                                         pos.z + shift[2]))


# ============================================================================
# MAIN WORKFLOW
# ============================================================================

def generate_all_dipeptides(output_dir: str,
                           amino_acids: List[str] = None,
                           phi_step: float = 30,
                           psi_step: float = 30,
                           solvate: bool = False,
                           n_waters: int = 50,
                           include_key_regions: bool = True):
    """
    Generate all Ace-X-NMe dipeptides across Ramachandran space.
    
    Parameters:
    -----------
    output_dir : str
        Output directory for PDB files
    amino_acids : list
        List of amino acid codes to generate (default: all 20)
    phi_step : float
        Step size for phi angle sampling (degrees)
    psi_step : float
        Step size for psi angle sampling (degrees)
    solvate : bool
        Whether to add explicit water solvation
    n_waters : int
        Number of water molecules for solvation
    include_key_regions : bool
        Include key Ramachandran regions regardless of grid
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    if amino_acids is None:
        amino_acids = list(AMINO_ACIDS.keys())
    
    # Calculate expected conformations
    n_phi = int(360 / phi_step)
    n_psi = int(360 / psi_step)
    n_grid = n_phi * n_psi
    n_key_regions = len(RAMACHANDRAN_REGIONS) if include_key_regions else 0
    
    print("=" * 60)
    print("Ace-X-NMe Dipeptide Generator")
    print("=" * 60)
    print(f"Output directory: {output_dir}")
    print(f"Amino acids: {len(amino_acids)}")
    print(f"Phi step: {phi_step}°")
    print(f"Psi step: {psi_step}°")
    print(f"Grid points: ~{n_grid} per amino acid")
    print(f"Key regions: {n_key_regions}")
    print(f"Solvation: {'Yes (' + str(n_waters) + ' waters)' if solvate else 'No'}")
    print(f"Estimated total: ~{len(amino_acids) * (n_grid + n_key_regions)} structures")
    print("=" * 60)
    print()
    
    total_generated = 0
    
    for aa_idx, aa_code in enumerate(amino_acids, 1):
        aa_info = AMINO_ACIDS[aa_code]
        print(f"[{aa_idx}/{len(amino_acids)}] {aa_code} ({aa_info['name']})...")
        
        # Create amino acid directory
        aa_dir = output_path / aa_code
        aa_dir.mkdir(exist_ok=True)
        
        # Build base dipeptide
        mol = build_dipeptide(aa_code)
        if mol is None:
            print(f"  ✗ Failed to build {aa_code}")
            continue
        
        print(f"  Atoms: {mol.GetNumAtoms()}")
        
        # Generate Ramachandran conformers
        conformers = generate_ramachandran_conformers(
            mol, aa_code,
            phi_step=phi_step,
            psi_step=psi_step,
            include_key_regions=include_key_regions
        )
        
        print(f"  Generated {len(conformers)} conformations")
        
        # Save each conformer
        saved = 0
        for conf_idx, (conf_mol, phi, psi) in enumerate(conformers):
            # Center and shift to positive coordinates
            center_molecule(conf_mol)
            shift_to_positive(conf_mol)
            
            # Add solvation if requested
            if solvate:
                conf_mol = add_water_shell(conf_mol, n_waters)
            
            # Generate filename
            phi_str = f"p{int(phi):+04d}".replace('-', 'm')
            psi_str = f"s{int(psi):+04d}".replace('-', 'm')
            filename = f"Ace_{aa_code}_NMe_{phi_str}_{psi_str}.pdb"
            
            output_file = aa_dir / filename
            if save_dipeptide_pdb(conf_mol, output_file, aa_code):
                saved += 1
        
        print(f"  ✓ Saved {saved} structures to {aa_dir.name}/")
        total_generated += saved
        print()
    
    print("=" * 60)
    print(f"✓ Complete! Generated {total_generated} total structures")
    print(f"  Output: {output_dir}")
    print("=" * 60)
    print()
    print("Next steps:")
    print("  1. Prepare ABACUS inputs:")
    print(f"     python scripts/3_prepare_abacus_inputs.py \\")
    print(f"         --conformations_dir {output_dir} \\")
    print(f"         --output_dir data/abacus_inputs")
    print()
    print("  2. Run ABACUS calculations:")
    print("     ./scripts/5_batch_run_abacus.sh data/abacus_inputs")
    print()
    print("  3. Convert to DeePMD format:")
    print("     python scripts/6_convert_to_deepmd.py")


def main():
    parser = argparse.ArgumentParser(
        description="Generate Ace-X-NMe dipeptides for all 20 amino acids across Ramachandran space",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate all 20 amino acids with 30° grid
  python 0_generate_dipeptides.py --output_dir data/dipeptides

  # Finer sampling (15° grid)
  python 0_generate_dipeptides.py --output_dir data/dipeptides --phi_step 15 --psi_step 15

  # With explicit water solvation
  python 0_generate_dipeptides.py --output_dir data/dipeptides --solvate --n_waters 50

  # Only specific amino acids
  python 0_generate_dipeptides.py --output_dir data/dipeptides --amino_acids ALA GLY VAL

  # Quick test with coarse sampling
  python 0_generate_dipeptides.py --output_dir data/dipeptides_test --phi_step 60 --psi_step 60 --amino_acids ALA GLY
        """
    )
    
    parser.add_argument(
        "--output_dir",
        type=str,
        default="data/dipeptides",
        help="Output directory for dipeptide PDB files (default: data/dipeptides)"
    )
    parser.add_argument(
        "--amino_acids",
        type=str,
        nargs="+",
        default=None,
        help="Specific amino acid codes to generate (default: all 20)"
    )
    parser.add_argument(
        "--phi_step",
        type=float,
        default=30,
        help="Step size for phi angle in degrees (default: 30)"
    )
    parser.add_argument(
        "--psi_step",
        type=float,
        default=30,
        help="Step size for psi angle in degrees (default: 30)"
    )
    parser.add_argument(
        "--solvate",
        action="store_true",
        help="Add explicit water solvation"
    )
    parser.add_argument(
        "--n_waters",
        type=int,
        default=50,
        help="Number of water molecules for solvation (default: 50)"
    )
    parser.add_argument(
        "--no_key_regions",
        action="store_true",
        help="Skip adding key Ramachandran regions (alpha helix, beta sheet, etc.)"
    )
    
    args = parser.parse_args()
    
    # Validate amino acid codes if provided
    if args.amino_acids:
        valid_codes = set(AMINO_ACIDS.keys())
        for code in args.amino_acids:
            if code.upper() not in valid_codes:
                print(f"Error: Unknown amino acid code: {code}")
                print(f"Valid codes: {', '.join(sorted(valid_codes))}")
                sys.exit(1)
        args.amino_acids = [code.upper() for code in args.amino_acids]
    
    generate_all_dipeptides(
        output_dir=args.output_dir,
        amino_acids=args.amino_acids,
        phi_step=args.phi_step,
        psi_step=args.psi_step,
        solvate=args.solvate,
        n_waters=args.n_waters,
        include_key_regions=not args.no_key_regions
    )


if __name__ == "__main__":
    main()

