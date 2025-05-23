#!/bin/bash
#SBATCH --job-name=comsol 
#SBATCH --time=12:00:00
#SBATCH --output=slurm4.out
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=END
#SBATCH --mail-user=smukherjee@unc.edu

module load comsol/6.1

comsol batch -alloc scalable -usebatchlic -forcegcc -blas mkl -keeplicenses on -np 20 -inputfile $1 -study $2
