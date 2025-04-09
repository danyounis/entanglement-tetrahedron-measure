#!/bin/bash
#SBATCH --partition=qot
#SBATCH --cpus-per-task=12
#SBATCH --mem=1G
#SBATCH --time=5-00:00:00
#SBATCH --job-name=tetra
#SBATCH --mail-user=ayounis@ur.rochester.edu
#SBATCH --mail-type=end

module unload gcc lapack openblas
module load gcc/11.2.0/b1 lapack/3.9.0/b2 openblas/0.3.10/b1

export OMP_NUM_THREADS=12
export PWD=/scratch/ayounis/tetrahedron/bin

$PWD/tetra_main input.deck
