# Protein Fragment Workflow for DeePMD-kit Training

A complete, reusable workflow for generating training data from protein fragments using ABACUS DFT calculations.

## Overview

This workflow extracts dipeptide fragments from a protein PDB file, generates multiple conformations for each fragment, runs ABACUS DFT calculations, and converts the results to DeePMD-kit format for training.

## Directory Structure

```
protein_fragment_workflow/
├── README.md                    # This file
├── scripts/                     # All workflow scripts
│   ├── 1_extract_fragments.py   # Extract dipeptides from PDB
│   ├── 2_generate_conformations.py  # Generate conformations
│   ├── 3_prepare_abacus_inputs.py  # Prepare ABACUS inputs
│   ├── 4_run_abacus_docker.sh   # Run single ABACUS calculation
│   ├── 5_batch_run_abacus.sh    # Batch run ABACUS
│   └── 6_convert_to_deepmd.py   # Convert to DeePMD format
├── data/                        # Data directories (created during workflow)
│   ├── fragments/               # Extracted dipeptide fragments
│   ├── conformations/           # Generated conformations
│   ├── abacus_inputs/           # ABACUS input files
│   ├── abacus_outputs/          # ABACUS calculation results
│   └── deepmd_data/             # DeePMD format data
└── resources/                   # Pseudopotentials and orbitals
    ├── *.upf                    # Pseudopotential files
    └── *.orb                    # Orbital files
```

## Prerequisites

### Software Requirements

1. **Python 3.7+** with packages:
   ```bash
   conda install -c conda-forge biopython rdkit scipy numpy
   ```

2. **Docker** (for running ABACUS):
   ```bash
   # Install Docker Desktop or Docker Engine
   # Verify: docker --version
   ```

3. **dpdata** (for conversion):
   ```bash
   conda install -c conda-forge dpdata
   ```

4. **DeePMD-kit** (for training):
   ```bash
   conda install -c conda-forge deepmd-kit
   ```

### Resources

- **Pseudopotentials and orbitals**: Place `.upf` and `.orb` files in `resources/` directory
- **Protein PDB file**: Your input protein structure

## Complete Workflow

### Step 1: Extract Fragments

Extract all dipeptide fragments from your protein PDB file.

```bash
cd protein_fragment_workflow
python scripts/1_extract_fragments.py \
    --pdb_file /path/to/your/protein.pdb \
    --output_dir data/fragments
```

**What it does:**
- Parses the PDB file
- Extracts all consecutive residue pairs (dipeptides)
- Saves each as a separate PDB file in `data/fragments/`

**Output:**
```
data/fragments/
├── fragment_000_GLY1_SER2.pdb
├── fragment_001_SER2_PRO1.pdb
└── ...
```

### Step 2: Generate Conformations

Generate multiple conformations for each fragment using RDKit.

```bash
python scripts/2_generate_conformations.py \
    --fragments_dir data/fragments \
    --output_dir data/conformations \
    --num_conf 50 \
    --min_distance 0.5
```

**What it does:**
- Generates multiple conformations for each fragment
- Validates conformations (filters atoms too close)
- Saves each conformation as a separate PDB file

**Parameters:**
- `--num_conf`: Number of conformations per fragment (default: 50)
- `--min_distance`: Minimum interatomic distance in Angstroms (default: 0.5)

**Output:**
```
data/conformations/
├── fragment_000_GLY1_SER2_conf_000.pdb
├── fragment_000_GLY1_SER2_conf_001.pdb
├── fragment_001_SER2_PRO1_conf_000.pdb
└── ...
```

### Step 3: Prepare ABACUS Inputs

Prepare ABACUS input files (INPUT, STRU, KPT) for each conformation.

```bash
python scripts/3_prepare_abacus_inputs.py \
    --conformations_dir data/conformations \
    --output_dir data/abacus_inputs \
    --pp_dir resources \
    --calculation md \
    --md_nstep 10
```

**What it does:**
- Parses each conformation PDB file
- Generates ABACUS STRU, INPUT, and KPT files
- Copies pseudopotentials and orbitals to each directory

**Parameters:**
- `--calculation`: 'md' for MD or 'scf' for single-point (default: 'md')
- `--md_nstep`: Number of MD steps (default: 10 for testing, use 50-200 for production)

**Output:**
```
data/abacus_inputs/
├── fragment_000_GLY1_SER2_conf_000/
│   ├── INPUT
│   ├── STRU
│   ├── KPT
│   ├── H_ONCV_PBE-1.2.upf
│   ├── C_ONCV_PBE-1.0.upf
│   └── *.orb files
├── fragment_000_GLY1_SER2_conf_001/
│   └── ...
└── ...
```

### Step 4: Run ABACUS Calculations

Run ABACUS DFT calculations on all conformations using Docker.

#### Option A: Run Single Calculation (Testing)

```bash
./scripts/4_run_abacus_docker.sh data/abacus_inputs/fragment_000_GLY1_SER2_conf_000
```

#### Option B: Batch Run All Conformations

```bash
./scripts/5_batch_run_abacus.sh data/abacus_inputs
```

**For testing, limit the number:**
```bash
./scripts/5_batch_run_abacus.sh data/abacus_inputs 10
```

**What it does:**
- Runs ABACUS MD calculations in Docker containers
- Each calculation creates `OUT.ABACUS/` directory with results
- Checks for `MD_dump` file to verify success

**Output:**
```
data/abacus_inputs/
├── fragment_000_GLY1_SER2_conf_000/
│   ├── INPUT, STRU, KPT, *.upf, *.orb
│   └── OUT.ABACUS/
│       ├── MD_dump              # Trajectory data
│       ├── running_md.log       # Calculation log
│       └── ...
└── ...
```

**Note:** ABACUS calculations can take a long time. For testing, use `--md_nstep 10` in Step 3. For production, use 50-200 steps.

### Step 5: Convert to DeePMD Format

Convert ABACUS outputs to DeePMD-kit training format.

```bash
python scripts/6_convert_to_deepmd.py \
    --abacus_dir data/abacus_inputs \
    --output_dir data/deepmd_data \
    --train-ratio 0.8
```

**What it does:**
- Finds all `OUT.ABACUS` directories with `MD_dump` files
- Converts each to DeePMD `.npy` format using `dpdata`
- Splits into training (80%) and validation (20%) sets

**Output:**
```
data/deepmd_data/
├── train/                       # Training datasets
│   ├── conf_fragment_000_GLY1_SER2_conf_000/
│   │   ├── type.raw
│   │   ├── type_map.raw
│   │   └── set.000/
│   │       ├── box.npy
│   │       ├── coord.npy
│   │       ├── energy.npy
│   │       └── force.npy
│   └── ...
└── val/                         # Validation datasets
    └── ...
```

## Quick Start Example

```bash
# 1. Extract fragments
python scripts/1_extract_fragments.py \
    --pdb_file protein.pdb \
    --output_dir data/fragments

# 2. Generate conformations (use fewer for testing)
python scripts/2_generate_conformations.py \
    --fragments_dir data/fragments \
    --output_dir data/conformations \
    --num_conf 10

# 3. Prepare ABACUS inputs (use 10 steps for testing)
python scripts/3_prepare_abacus_inputs.py \
    --conformations_dir data/conformations \
    --output_dir data/abacus_inputs \
    --pp_dir resources \
    --md_nstep 10

# 4. Run ABACUS (test with 5 conformations first)
./scripts/5_batch_run_abacus.sh data/abacus_inputs 5

# 5. Convert to DeePMD format
python scripts/6_convert_to_deepmd.py \
    --abacus_dir data/abacus_inputs \
    --output_dir data/deepmd_data
```

## Configuration

### ABACUS Input Parameters

Edit `scripts/3_prepare_abacus_inputs.py` to customize:

- **MD steps**: Change `md_nstep` (default: 10 for testing)
- **Temperature**: Change `md_tfirst` (default: 300 K)
- **Cutoff**: Change `ecutwfc` (default: 100 Ry)
- **SCF threshold**: Change `scf_thr` (default: 1e-7)

### Conformation Generation

Edit `scripts/2_generate_conformations.py` to customize:

- **Number of conformations**: Change `--num_conf` (default: 50)
- **Minimum distance**: Change `--min_distance` (default: 0.5 Å)

## Troubleshooting

### "No MD_dump file created"

- Check `OUT.ABACUS/running_md.log` for errors
- Verify `md_dumpfreq=1` in INPUT file
- Check `dump_force=1`, `dump_vel=1`, `dump_virial=1` are set
- Ensure calculation completed (not just started)

### "Structure is unreasonable" warning

- Conformations may have atoms too close
- Reduce `--min_distance` or regenerate conformations
- Check `OUT.ABACUS/warning.log` for details

### Docker issues

- Ensure Docker is running: `docker info`
- Check Docker image: `docker images | grep abacus`
- Pull image manually: `docker pull registry.dp.tech/deepmodeling/abacus`

### Conversion errors

- Ensure `MD_dump` files exist in `OUT.ABACUS/` directories
- Check `dpdata` is installed: `conda install -c conda-forge dpdata`
- Verify ABACUS calculations completed successfully

## Next Steps: Training DeePMD Model

After conversion, create `input.json` for training:

```json
{
    "model": {
        "type_map": ["H", "C", "N", "O"],
        "descriptor": {
            "type": "se_e2_a",
            "sel": [60, 30, 30, 30],
            "rcut": 6.0,
            "rcut_smth": 0.5,
            "neuron": [25, 50, 100],
            "resnet_dt": false,
            "axis_neuron": 16
        },
        "fitting_net": {
            "neuron": [240, 240, 240],
            "resnet_dt": true
        }
    },
    "learning_rate": {
        "type": "exp",
        "start_lr": 0.001,
        "decay_steps": 5000,
        "decay_rate": 0.95
    },
    "loss": {
        "start_pref_e": 0.02,
        "limit_pref_e": 1,
        "start_pref_f": 1000,
        "limit_pref_f": 1
    },
    "training": {
        "training_data": {
            "systems": ["data/deepmd_data/train/*"],
            "batch_size": 1,
            "set_prefix": "set"
        },
        "validation_data": {
            "systems": ["data/deepmd_data/val/*"],
            "batch_size": 1,
            "set_prefix": "set"
        },
        "numb_steps": 1000000,
        "seed": 42,
        "disp_file": "lcurve.out",
        "disp_freq": 1000,
        "save_freq": 10000
    }
}
```

Then train:
```bash
dp train input.json
```

## Tips for Production

1. **Start small**: Test with 5-10 conformations and 10 MD steps
2. **Scale gradually**: Increase to 50 conformations, then 50-100 MD steps
3. **Monitor resources**: ABACUS calculations are CPU/memory intensive
4. **Use parallel execution**: Run multiple ABACUS calculations simultaneously if you have resources
5. **Check data quality**: Verify converted data before training

## File Naming Convention

- **Fragments**: `fragment_XXX_RES1YY_RES2ZZ.pdb`
- **Conformations**: `fragment_XXX_RES1YY_RES2ZZ_conf_NNN.pdb`
- **ABACUS directories**: `fragment_XXX_RES1YY_RES2ZZ_conf_NNN/`

## Support

For issues or questions:
- Check ABACUS documentation: https://abacus.deepmodeling.com/
- Check DeePMD-kit documentation: https://docs.deepmodeling.com/projects/deepmd/
- Check dpdata documentation: https://docs.deepmodeling.com/projects/dpdata/

## License

This workflow is provided as-is for research and educational purposes.

