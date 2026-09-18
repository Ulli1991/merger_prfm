#!/bin/bash
#SBATCH -J dm_cubetest2
#SBATCH -o /ptmp/uli/dwarf_merger/logs/cubes_test2.log
#SBATCH -p general
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=72
#SBATCH --mem=150000
#SBATCH -t 01:30:00
export NUMBA_NUM_THREADS=72 OMP_NUM_THREADS=72 PRFM_CUBE=0.5
cd /raven/u/uli/dwarf_merger/prfm
/raven/u/uli/jupyter_env/bin/python cubes.py 100 2>&1 | grep -v Warning
