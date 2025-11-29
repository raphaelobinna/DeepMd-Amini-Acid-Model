#!/usr/bin/env python3
"""
Evaluate trained DeepMD model by comparing DFT energies with predictions.

This script:
1. Loads individual training/validation systems
2. Uses the trained model to predict energies for each system
3. Plots DFT vs predicted energies for comparison (default: 3 systems)

Usage:
    python 7_evaluate_model.py --model graph.pb --data data/deepmd_data/train
    python 7_evaluate_model.py --model graph.pb --data data/deepmd_data/val --num-plots 5
"""

import argparse
import sys
from pathlib import Path

try:
    import dpdata
except ImportError:
    print("Error: dpdata is not installed.")
    print("Install it with: conda install -c conda-forge dpdata")
    sys.exit(1)

try:
    import matplotlib.pyplot as plt
    import numpy as np
except ImportError:
    print("Error: matplotlib or numpy is not installed.")
    print("Install with: conda install matplotlib numpy")
    sys.exit(1)


def find_system_dirs(data_dir):
    """Find all DeepMD system directories."""
    data_path = Path(data_dir)
    
    if not data_path.exists():
        print(f"Error: Data directory not found: {data_path}")
        sys.exit(1)
    
    # Find all system directories (those with type.raw)
    system_dirs = sorted([d for d in data_path.iterdir() if d.is_dir() and (d / "type.raw").exists()])
    
    if not system_dirs:
        # Maybe the data_dir itself is a system
        if (data_path / "type.raw").exists():
            system_dirs = [data_path]
        else:
            print(f"Error: No DeepMD systems found in {data_path}")
            sys.exit(1)
    
    return system_dirs


def evaluate_single_system(system_dir, model_path):
    """Evaluate model predictions for a single system."""
    try:
        system = dpdata.LabeledSystem(str(system_dir), fmt="deepmd/npy")
        if len(system) == 0:
            return None
        
        predictions = system.predict(str(model_path))
        
        dft_energies = system["energies"]
        pred_energies = predictions["energies"]
        n_atoms = len(system["atom_types"])
        
        # Calculate statistics
        energy_diff = pred_energies - dft_energies
        mae = np.mean(np.abs(energy_diff))
        rmse = np.sqrt(np.mean(energy_diff**2))
        mae_per_atom = mae / n_atoms
        rmse_per_atom = rmse / n_atoms
        
        return {
            "name": system_dir.name,
            "dft_energies": dft_energies,
            "pred_energies": pred_energies,
            "n_frames": len(dft_energies),
            "n_atoms": n_atoms,
            "mae": mae,
            "rmse": rmse,
            "mae_per_atom": mae_per_atom,
            "rmse_per_atom": rmse_per_atom,
        }
    except Exception as e:
        print(f"  Warning: Failed to evaluate {system_dir.name}: {e}")
        return None


def plot_system_comparison(result, output_dir, idx):
    """Plot DFT vs predicted energies for a single system."""
    plt.figure(figsize=(8, 8))
    
    dft_energies = result["dft_energies"]
    pred_energies = result["pred_energies"]
    
    plt.scatter(dft_energies, pred_energies, alpha=0.7, edgecolors='none', s=50, c='#2E86AB')
    
    # Add diagonal reference line
    all_energies = np.concatenate([dft_energies, pred_energies])
    plot_min = all_energies.min() - 0.5
    plot_max = all_energies.max() + 0.5
    plt.xlim(plot_min, plot_max)
    plt.ylim(plot_min, plot_max)
    
    x_range = np.linspace(plot_min, plot_max, 100)
    plt.plot(x_range, x_range, "r--", linewidth=1.5, label="Perfect prediction")
    
    plt.xlabel("DFT Energy (eV)", fontsize=12)
    plt.ylabel("Predicted Energy (eV)", fontsize=12)
    plt.title(f"{result['name']}\n"
              f"Frames: {result['n_frames']}, Atoms: {result['n_atoms']}\n"
              f"MAE: {result['mae_per_atom']*1000:.3f} meV/atom, "
              f"RMSE: {result['rmse_per_atom']*1000:.3f} meV/atom")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    output_file = output_dir / f"energy_comparison_{idx+1}_{result['name']}.png"
    plt.savefig(output_file, dpi=150)
    print(f"  Plot saved: {output_file}")
    plt.close()


def evaluate_model(model_path, data_dir, num_plots=3, output_dir=None):
    """Evaluate model predictions against DFT data for individual systems."""
    
    model_path = Path(model_path)
    if not model_path.exists():
        print(f"Error: Model file not found: {model_path}")
        sys.exit(1)
    
    # Set output directory
    if output_dir is None:
        output_dir = Path(".")
    else:
        output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Find all systems
    system_dirs = find_system_dirs(data_dir)
    print(f"\nFound {len(system_dirs)} systems in {data_dir}")
    
    # Limit number of systems to evaluate for plotting
    systems_to_plot = min(num_plots, len(system_dirs))
    print(f"Will generate plots for {systems_to_plot} systems")
    
    print(f"\nEvaluating with model: {model_path}")
    print("-" * 50)
    
    all_results = []
    
    for i, sys_dir in enumerate(system_dirs):
        result = evaluate_single_system(sys_dir, model_path)
        if result is not None:
            all_results.append(result)
            print(f"  [{i+1}/{len(system_dirs)}] {result['name']}: "
                  f"MAE={result['mae_per_atom']*1000:.3f} meV/atom, "
                  f"RMSE={result['rmse_per_atom']*1000:.3f} meV/atom")
    
    if not all_results:
        print("Error: No valid systems evaluated!")
        sys.exit(1)
    
    # Calculate overall statistics
    total_mae = np.mean([r["mae_per_atom"] for r in all_results])
    total_rmse = np.mean([r["rmse_per_atom"] for r in all_results])
    
    print(f"\n{'='*50}")
    print("Overall Statistics (averaged across systems)")
    print(f"{'='*50}")
    print(f"Systems evaluated: {len(all_results)}")
    print(f"Average MAE/atom:  {total_mae*1000:.3f} meV/atom")
    print(f"Average RMSE/atom: {total_rmse*1000:.3f} meV/atom")
    print(f"{'='*50}")
    
    # Generate individual plots for first N systems
    print(f"\nGenerating plots for {systems_to_plot} systems...")
    for i in range(systems_to_plot):
        plot_system_comparison(all_results[i], output_dir, i)
    
    # Generate summary plot with all systems
    print("\nGenerating summary plot...")
    plt.figure(figsize=(10, 8))
    
    colors = plt.cm.tab10(np.linspace(0, 1, len(all_results)))
    
    for i, result in enumerate(all_results):
        plt.scatter(result["dft_energies"], result["pred_energies"], 
                   alpha=0.7, s=30, c=[colors[i]], 
                   label=f"{result['name'][:30]}..." if len(result['name']) > 30 else result['name'])
    
    # Add diagonal reference line
    all_dft = np.concatenate([r["dft_energies"] for r in all_results])
    all_pred = np.concatenate([r["pred_energies"] for r in all_results])
    plot_min = min(all_dft.min(), all_pred.min()) - 0.5
    plot_max = max(all_dft.max(), all_pred.max()) + 0.5
    
    x_range = np.linspace(plot_min, plot_max, 100)
    plt.plot(x_range, x_range, "k--", linewidth=1.5, label="Perfect prediction")
    
    plt.xlabel("DFT Energy (eV)", fontsize=12)
    plt.ylabel("Predicted Energy (eV)", fontsize=12)
    plt.title(f"DeepMD Model Evaluation - All Systems\n"
              f"Avg MAE: {total_mae*1000:.3f} meV/atom, Avg RMSE: {total_rmse*1000:.3f} meV/atom")
    
    # Put legend outside if too many systems
    if len(all_results) > 5:
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    else:
        plt.legend()
    
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    summary_file = output_dir / "energy_comparison_summary.png"
    plt.savefig(summary_file, dpi=150, bbox_inches='tight')
    print(f"  Summary plot saved: {summary_file}")
    plt.close()
    
    print("\nDone!")
    
    return all_results


def main():
    parser = argparse.ArgumentParser(
        description="Evaluate DeepMD model predictions against DFT data"
    )
    parser.add_argument(
        "--model", "-m",
        default="graph.pb",
        help="Path to trained model (graph.pb). Default: graph.pb"
    )
    parser.add_argument(
        "--data", "-d",
        default="data/deepmd_data/train",
        help="Path to DeepMD data directory. Default: data/deepmd_data/train"
    )
    parser.add_argument(
        "--num-plots", "-n",
        type=int,
        default=3,
        help="Number of individual system plots to generate. Default: 3"
    )
    parser.add_argument(
        "--output-dir", "-o",
        default=None,
        help="Output directory for plots. Default: current directory"
    )
    
    args = parser.parse_args()
    
    print("=" * 50)
    print("DeepMD Model Evaluation")
    print("=" * 50)
    
    evaluate_model(args.model, args.data, args.num_plots, args.output_dir)


if __name__ == "__main__":
    main()
