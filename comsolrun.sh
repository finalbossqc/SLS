#!/bin/bash
#SBATCH --job-name=comsol 
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --output=slurm1.out
#SBATCH --cpus-per-task=10
#SBATCH --mail-type=END
#SBATCH --mail-user=smukherjee@unc.edu

module load comsol/6.1

comsol batch -alloc scalable  -usebatchlic -forcegcc -blas mkl -keeplicenses on -np 10 -inputfile $1 -study $2
