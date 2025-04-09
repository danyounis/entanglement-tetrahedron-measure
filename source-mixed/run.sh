#!/bin/bash
#SBATCH --partition=qot
#SBATCH --cpus-per-task=12
#SBATCH --mem=8G
#SBATCH --time=5-00:00:00
#SBATCH --job-name=tetra-ea
#SBATCH --mail-user=ayounis@ur.rochester.edu
#SBATCH --mail-type=end

export OMP_NUM_THREADS=12

module unload anaconda3 gcc lapack openblas
module load anaconda3/2021.11 gcc/11.2.0/b1 lapack/3.9.0/b2 openblas/0.3.10/b1

export home=/scratch/ayounis/tetrahedron

cd $home/LEAP/
make setup

cd $SLURM_SUBMIT_DIR
python3 main.py serial input.deck output.h5

unset home
