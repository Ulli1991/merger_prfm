#!/bin/bash
#SBATCH -J dm_cubes
#SBATCH -o /ptmp/uli/dwarf_merger/logs/cubes_%a.log
#SBATCH -p general
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=72
#SBATCH --mem=110000
#SBATCH -t 04:00:00
#SBATCH --array=0-23
# 3 cube sizes x 8 snapshot subsets (every 4th snapshot): task = cfg*8 + k
export NUMBA_NUM_THREADS=72 OMP_NUM_THREADS=72
CFG=$((SLURM_ARRAY_TASK_ID / 8)); K=$((SLURM_ARRAY_TASK_ID % 8))
case $CFG in 0) export PRFM_CUBE=0.5 ;; 1) export PRFM_CUBE=0.25 ;; 2) export PRFM_CUBE=0.125 ;; esac
cd /raven/u/uli/dwarf_merger/prfm
SNAPS=$(seq 0 4 231 | awk -v k=$K 'NR%8==k+1 || (k==7 && NR%8==0)')
echo "cube $PRFM_CUBE snaps: $SNAPS"
/raven/u/uli/jupyter_env/bin/python cubes.py $SNAPS
