#!/bin/sh

#SBATCH --array=0-999
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=chase_escape
#SBATCH --error=error_%A_%a.out
#SBATCH --output=std_%A_%a.out
#SBATCH --partition=standard
#SBATCH --mail-type=END
#SBATCH --mail-user=ritam.pal@students.iiserpune.ac.in
#SBATCH --qos=array-job

# Parameters for the simulation
L=50                    # Lattice size
LAMBDA=0.10            # Spontaneous conversion rate
REALIZATIONS_PER_JOB=10 # Number of realizations per job

# Calculate start and end indices for this job
START_INDEX=$((SLURM_ARRAY_TASK_ID * REALIZATIONS_PER_JOB))
END_INDEX=$(((SLURM_ARRAY_TASK_ID + 1) * REALIZATIONS_PER_JOB))

echo "Job ID: $SLURM_ARRAY_TASK_ID"
echo "Running realizations $START_INDEX to $((END_INDEX - 1))"

# Change to submission directory
cd $SLURM_SUBMIT_DIR

# Compile the program if needed (uncomment if not already compiled)
# gcc -o chase_escape_parallel src/chase_escape_parallel.c -lm

# Run the simulation
./chase_escape_parallel $L $LAMBDA $START_INDEX $END_INDEX

echo "Job $SLURM_ARRAY_TASK_ID completed" 