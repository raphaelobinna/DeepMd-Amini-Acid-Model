# DeepMD Potential for Amino Acid Dipeptides

A deep learning-based machine learning potential for accurate and efficient simulations of amino acid dipeptides, trained using the DeePMD-kit framework with DFT reference data from ABACUS.

## Model Overview

This repository contains a DeepMD potential trained on all 20 standard amino acid dipeptides (Ace-X-NMe) across diverse backbone conformations. The model achieves near-chemical accuracy while providing ~1000× speedup compared to DFT calculations.

### Key Features

- **Coverage**: All 20 standard amino acids
- **Conformational diversity**: Multiple φ/ψ backbone angles per residue
- **Elements**: C, H, N, O, S
- **Accuracy**: ~2.2 meV/atom MAE on validation set
- **Speed**: Orders of magnitude faster than DFT

## Model Performance

### Validation Set Results

| Metric | Value |
|--------|-------|
| Mean Absolute Error (MAE) | **2.181 meV/atom** |
| Root Mean Square Error (RMSE) | **2.295 meV/atom** |
| Validation systems | 66 |

### Training Convergence (1,000,000 steps)

| Metric | Final Value |
|--------|-------------|
| Force RMSE (validation) | 0.056 eV/Å |
| Force RMSE (training) | 0.071 eV/Å |
| Energy RMSE (validation) | 0.66 meV |
| Energy RMSE (training) | 2.02 meV |

## Dataset

### Composition

| Set | Systems | Description |
|-----|---------|-------------|
| Training | 262 | 80% of conformations |
| Validation | 66 | 20% of conformations |

### Amino Acids Included

All 20 standard amino acids in dipeptide form (Acetyl-X-N-methylamide):

| Non-polar | Polar | Charged | Aromatic |
|-----------|-------|---------|----------|
| ALA, VAL, LEU, ILE | SER, THR, ASN, GLN | ASP, GLU, LYS, ARG | PHE, TYR, TRP |
| GLY, PRO, MET, CYS | | HIS | |

### Backbone Conformations

Each dipeptide was sampled at multiple (φ, ψ) dihedral angles covering:
- α-helix region (φ ≈ -60°, ψ ≈ -45°)
- β-sheet region (φ ≈ -120°, ψ ≈ +120°)
- Polyproline II (φ ≈ -75°, ψ ≈ +145°)
- Additional conformational basins

## Model Architecture

```
Descriptor: se_e2_a (smooth edition, two-body embedding)
├── Cutoff radius: 6.0 Å
├── Smoothing cutoff: 0.5 Å
├── Embedding network: [25, 50, 100]
├── Axis neurons: 16
└── Selection: [60, 120, 40, 40, 10] (C, H, N, O, S)

Fitting Network: [240, 240, 240]
├── ResNet connections: enabled
└── Output: atomic energies
```

## Training Details

| Parameter | Value |
|-----------|-------|
| Training steps | 1,000,000 |
| Initial learning rate | 1.0 × 10⁻³ |
| Final learning rate | 3.5 × 10⁻⁸ |
| Learning rate schedule | Exponential decay |
| Batch size | Auto |
| Loss prefactors | Energy: 0.02→1, Force: 1000→1 |

### DFT Reference Calculations

| Parameter | Value |
|-----------|-------|
| Software | ABACUS |
| Exchange-correlation | PBE (GGA) |
| Basis set | Numerical atomic orbitals |
| Cutoff energy | 100 Ry |
| Pseudopotentials | ONCV (PBE) |

## Usage

### Installation

```bash
# Create conda environment
conda create -n deepmd python=3.12
conda activate deepmd

# Install DeePMD-kit
conda install -c conda-forge deepmd-kit=*=*cpu*
```

### Running Inference

```python
from deepmd.infer import DeepPot

# Load the model
dp = DeepPot("frozen_model.pth")

# Predict energy and forces
coord = ...  # atomic coordinates (N, 3) in Angstroms
cell = ...   # cell vectors (3, 3) in Angstroms  
atype = ...  # atom types [0=C, 1=H, 2=N, 3=O, 4=S]

e, f, v = dp.eval(coord, cell, atype)
```

### LAMMPS Integration

```lammps
units           metal
atom_style      atomic

pair_style      deepmd frozen_model.pth
pair_coeff      * *

# Run MD simulation
velocity        all create 300.0 12345
fix             1 all nvt temp 300.0 300.0 0.1
timestep        0.001
run             100000
```

## Repository Structure

```
deepmd-amino-group/
├── data/
│   ├── dipeptides/          # Initial PDB structures
│   │   ├── ALA/
│   │   ├── ARG/
│   │   └── ...
│   ├── abacus_inputs/       # DFT calculation inputs/outputs
│   └── deepmd_data/         # Training data & model
│       ├── train/           # Training systems
│       ├── val/             # Validation systems
│       ├── frozen_model.pth # Production model
│       ├── input.json       # Training configuration
│       ├── lcurve.out       # Learning curve
│       └── evaluation/      # Evaluation plots
├── resources/               # Pseudopotentials & basis sets
├── scripts/                 # Data generation & training scripts
└── README.md
```

## Scripts

| Script | Description |
|--------|-------------|
| `0_generate_dipeptides.py` | Generate dipeptide PDB structures |
| `1_extract_fragments.py` | Extract molecular fragments |
| `2_generate_conformations.py` | Generate backbone conformations |
| `3_prepare_abacus_inputs.py` | Prepare DFT input files |
| `4_run_abacus_docker.sh` | Run ABACUS in Docker |
| `5_batch_run_abacus.sh` | Batch DFT calculations |
| `6_convert_to_deepmd.py` | Convert to DeepMD format |
| `7_evaluate_model.py` | Evaluate model accuracy |

## Citation

If you use this potential in your research, please cite:

```bibtex
@article{amino_dipeptide_deepmd_2026,
  title={Deep Learning Potential for Amino Acid Dipeptides},
  author={[Author Names]},
  journal={[Journal Name]},
  year={2026},
  doi={[DOI]}
}
```

## References

1. Zhang, L., Han, J., Wang, H., Car, R., & E, W. (2018). Deep potential molecular dynamics: a scalable model with the accuracy of quantum mechanics. *Physical Review Letters*, 120(14), 143001.

2. Wang, H., Zhang, L., Han, J., & E, W. (2018). DeePMD-kit: A deep learning package for many-body potential energy representation and molecular dynamics. *Computer Physics Communications*, 228, 178-184.

3. Chen, M., Guo, G. C., & He, L. (2010). Systematically improvable optimized atomic basis sets for ab initio calculations. *Journal of Physics: Condensed Matter*, 22(44), 445501.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- DeePMD-kit development team
- ABACUS development team

