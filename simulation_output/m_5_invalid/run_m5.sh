#!/bin/bash
#SBATCH --job-name=m5_8
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=2G               # total memory per node (4 GB per cpu-core is default)
#SBATCH --time=20:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin        # send mail when job begins
#SBATCH --mail-type=end          # send mail when job ends
#SBATCH --mail-type=fail         # send mail if job fails
#SBATCH --mail-user=az2747@princeton.edu #CHANGE THIS TO YOUR EMAIL
# args = c(n_subpops, n_sim, M, n_replicates)

Rscript script_per_k.R 8 500 5 5000 > m5_8.txt
