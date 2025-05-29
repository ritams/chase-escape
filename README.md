# Parallelized Chase Escape Simulation

This repository contains a parallelized version of the chase escape simulation designed for SLURM array jobs.

## Overview

- **Total realizations**: 10,000
- **Array jobs**: 1,000 (jobs 0-999)
- **Realizations per job**: 10
- **Parallelization**: Each job runs a subset of realizations independently

## Files

- `src/chase_escape_parallel.c` - Main simulation code (parallelized version)
- `run_chase_escape.sh` - SLURM array job script
- `Makefile` - Compilation instructions
- `combine_results.py` - Post-processing script to combine results
- `README.md` - This file

## Quick Start

### 1. Compile the program
```bash
make
```

### 2. Submit the array job
```bash
sbatch run_chase_escape.sh
```

### 3. Monitor job progress
```bash
squeue -u $USER
```

### 4. Combine results after all jobs complete
```bash
python3 combine_results.py 50 0.10
```

## Detailed Instructions

### Compilation
The program can be compiled using the provided Makefile:
```bash
make clean  # Remove old executable
make        # Compile with optimization flags
```

Or manually:
```bash
gcc -O3 -Wall -std=c99 -o chase_escape_parallel src/chase_escape_parallel.c -lm
```

### Running Individual Jobs (for testing)
```bash
# Job 0: realizations 0-9
./chase_escape_parallel 50 0.10 0 10

# Job 1: realizations 10-19
./chase_escape_parallel 50 0.10 10 20

# Job 999: realizations 9990-9999
./chase_escape_parallel 50 0.10 9990 10000
```

### SLURM Script Configuration

The `run_chase_escape.sh` script is configured for:
- **Array range**: 0-999 (1000 jobs)
- **Time limit**: 2 hours per job
- **Memory**: 4GB per job
- **Lattice size (L)**: 50
- **Lambda**: 0.10

To modify parameters, edit the variables at the top of `run_chase_escape.sh`:
```bash
L=50                    # Lattice size
LAMBDA=0.10            # Spontaneous conversion rate
REALIZATIONS_PER_JOB=10 # Number of realizations per job
```

### Output Files

Each job creates a file named:
```
escape_prob_L{L}_lambda{lambda}_{start_index}_{end_index}.txt
```

Examples:
- `escape_prob_L50_lambda0.10_0_10.txt` (Job 0)
- `escape_prob_L50_lambda0.10_10_20.txt` (Job 1)
- `escape_prob_L50_lambda0.10_9990_10000.txt` (Job 999)

### Post-processing

After all jobs complete, combine the results:
```bash
python3 combine_results.py 50 0.10
```

This creates:
- `escape_prob_L50_lambda0.10_combined.txt` - Final averaged results with standard errors

## File Formats

### Individual Job Output
```
# p	escape_probability
# Realizations: 10 (indices 0 to 9)
1.9500	0.100000
1.9525	0.200000
...
```

### Combined Output
```
# Combined results from parallel chase escape simulations
# Total realizations: 10000
# Number of jobs: 1000
# p	escape_probability	std_error
1.9500	0.156789	0.012345
1.9525	0.234567	0.015678
...
```

## Troubleshooting

### Check job status
```bash
squeue -u $USER                    # Running jobs
sacct -j JOBID --format=JobID,State,ExitCode  # Completed jobs
```

### Check for failed jobs
```bash
ls error_*.out | wc -l             # Count error files
grep -l "Error" error_*.out        # Find jobs with errors
```

### Rerun specific failed jobs
If job 123 failed:
```bash
sbatch --array=123 run_chase_escape.sh
```

### Check output file count
```bash
ls escape_prob_L50_lambda0.10_*_*.txt | wc -l  # Should be 1000
```

## Performance Notes

- Each job should complete in ~30-60 minutes
- Total wall time: ~1-2 hours (depending on cluster load)
- Memory usage: ~1-2GB per job
- Disk usage: ~1KB per output file (1MB total)

## Customization

To run with different parameters:

1. **Different lattice size**: Edit `L=50` in `run_chase_escape.sh`
2. **Different lambda**: Edit `LAMBDA=0.10` in `run_chase_escape.sh`
3. **More realizations**: Increase array size and adjust `REALIZATIONS_PER_JOB`
4. **Different p-range**: Modify the p-value loop in `chase_escape_parallel.c`