#!/usr/bin/env python3
"""
Script to combine results from parallel chase escape simulations.
Averages escape probabilities across all array jobs.
"""

import glob
import numpy as np
import sys
import os

def combine_results(L, lambda_rate, total_jobs=1000):
    """
    Combine results from all parallel jobs into a single averaged result.
    
    Args:
        L: Lattice size
        lambda_rate: Spontaneous conversion rate
        total_jobs: Total number of array jobs (default 1000)
    """
    
    # Pattern to match output files
    pattern = f"escape_prob_L{L}_lambda{lambda_rate:.2f}_*_*.txt"
    files = glob.glob(pattern)
    
    if not files:
        print(f"No files found matching pattern: {pattern}")
        return
    
    print(f"Found {len(files)} result files")
    
    # Dictionary to store results: p_value -> list of escape probabilities
    results = {}
    total_realizations = 0
    
    for filename in files:
        try:
            with open(filename, 'r') as f:
                lines = f.readlines()
                
            # Extract realization count from header
            for line in lines:
                if line.startswith("# Realizations:"):
                    realizations = int(line.split()[2])
                    total_realizations += realizations
                    break
            
            # Read data lines (skip comments)
            for line in lines:
                if not line.startswith('#') and line.strip():
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        p_val = float(parts[0])
                        escape_prob = float(parts[1])
                        
                        if p_val not in results:
                            results[p_val] = []
                        results[p_val].append((escape_prob, realizations))
                        
        except Exception as e:
            print(f"Error reading {filename}: {e}")
            continue
    
    if not results:
        print("No valid data found in result files")
        return
    
    # Calculate weighted averages
    output_filename = f"escape_prob_L{L}_lambda{lambda_rate:.2f}_combined.txt"
    
    with open(output_filename, 'w') as f:
        f.write("# Combined results from parallel chase escape simulations\n")
        f.write(f"# Total realizations: {total_realizations}\n")
        f.write(f"# Number of jobs: {len(files)}\n")
        f.write("# p\tescape_probability\tstd_error\n")
        
        # Sort p values
        p_values = sorted(results.keys())
        
        for p_val in p_values:
            escape_probs = [prob for prob, _ in results[p_val]]
            realizations_list = [real for _, real in results[p_val]]
            
            # Weighted average (though all jobs should have same number of realizations)
            total_weight = sum(realizations_list)
            weighted_avg = sum(prob * real for prob, real in results[p_val]) / total_weight
            
            # Standard error calculation
            std_error = np.std(escape_probs) / np.sqrt(len(escape_probs))
            
            f.write(f"{p_val:.4f}\t{weighted_avg:.6f}\t{std_error:.6f}\n")
    
    print(f"Combined results saved to: {output_filename}")
    print(f"Total realizations processed: {total_realizations}")
    print(f"Expected realizations: {total_jobs * 10}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 combine_results.py <L> <lambda>")
        print("Example: python3 combine_results.py 50 0.10")
        sys.exit(1)
    
    L = int(sys.argv[1])
    lambda_rate = float(sys.argv[2])
    
    combine_results(L, lambda_rate) 