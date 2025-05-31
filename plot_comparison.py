#!/usr/bin/env python3
"""
Minimal plotting script for chase escape simulation comparison.
"""

import matplotlib.pyplot as plt
import numpy as np
import glob
from pathlib import Path

# Set up matplotlib for clean styling
plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'legend.fontsize': 12,
    'figure.figsize': (10, 7),
    'savefig.dpi': 300,
    'savefig.bbox': 'tight'
})

def read_data(filename):
    """Read data from simulation output file."""
    p_values = []
    escape_probs = []
    std_errors = []
    
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) >= 3:
                p_values.append(float(parts[0]))
                escape_probs.append(float(parts[1]))
                std_errors.append(float(parts[2]))
    
    return np.array(p_values), np.array(escape_probs), np.array(std_errors)

def extract_L(filename):
    """Extract L from filename."""
    parts = Path(filename).stem.split('_')
    for part in parts:
        if part.startswith('L'):
            return int(part[1:])
    return None

def find_intersection_point(files_data):
    """Find approximate intersection point of all curves."""
    # Get common p range
    all_p = []
    for _, (p_vals, _, _) in files_data:
        all_p.extend(p_vals)
    
    p_min = max([p_vals.min() for _, (p_vals, _, _) in files_data])
    p_max = min([p_vals.max() for _, (p_vals, _, _) in files_data])
    
    # Interpolate all curves to common p grid
    p_common = np.linspace(p_min, p_max, 100)
    interpolated_probs = []
    
    for _, (p_vals, escape_probs, _) in files_data:
        interp_probs = np.interp(p_common, p_vals, escape_probs)
        interpolated_probs.append(interp_probs)
    
    # Find point where variance is minimum (curves are closest)
    variances = []
    for i in range(len(p_common)):
        probs_at_p = [curve[i] for curve in interpolated_probs]
        variances.append(np.var(probs_at_p))
    
    min_var_idx = np.argmin(variances)
    intersection_p = p_common[min_var_idx]
    intersection_prob = np.mean([curve[min_var_idx] for curve in interpolated_probs])
    
    return intersection_p, intersection_prob

def create_comparison_plot():
    """Create comparison plot with modifications."""
    
    # Find all result files
    files = glob.glob("escape_prob_*_combined.txt")
    if not files:
        print("No result files found.")
        return
    
    # Read data and sort by L value
    files_data = []
    for filename in files:
        L = extract_L(filename)
        if L is not None:
            data = read_data(filename)
            files_data.append((L, data))
    
    # Sort by L value
    files_data.sort(key=lambda x: x[0])
    
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # Modern color palette
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    
    # Plot each curve with straight lines
    for i, (L, (p_values, escape_probs, std_errors)) in enumerate(files_data):
        color = colors[i % len(colors)]
        ax.errorbar(p_values, escape_probs, yerr=std_errors,
                   fmt='o-', color=color, markersize=6, capsize=3,
                   linewidth=2, alpha=0.8, label=f'L = {L}',
                   markeredgecolor='white', markeredgewidth=0.5)
    
    # Find and plot intersection point
    intersection_p, intersection_prob = find_intersection_point(files_data)
    ax.axvline(x=intersection_p, color='red', linestyle='--', linewidth=2, alpha=0.7)
    ax.text(intersection_p + 0.002, 0.95, f'p = {intersection_p:.3f}', 
            rotation=0, verticalalignment='top', fontsize=12, color='red')
    
    # Customize plot
    ax.set_xlabel(r'$\lambda$', fontsize=16)
    ax.set_ylabel('Escape Probability', fontsize=16)
    
    # Make xticks bold
    ax.tick_params(axis='x', labelsize=12)
    for label in ax.get_xticklabels():
        label.set_fontweight('bold')
    
    # Make yticks bold
    ax.tick_params(axis='y', labelsize=12)
    for label in ax.get_yticklabels():
        label.set_fontweight('bold')
    
    
    # Grid and styling
    ax.grid(True, alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Legend ordered by L value (no border)
    ax.legend(frameon=False, fontsize=12)
    
    plt.tight_layout()
    plt.savefig('comparison_plot.png', format='png', bbox_inches='tight', dpi=300)
    
    print(f"Comparison plot saved as: comparison_plot.png")
    print(f"Intersection point: p = {intersection_p:.3f}, Escape Probability = {intersection_prob:.3f}")
    
    return fig, ax

if __name__ == "__main__":
    create_comparison_plot() 