#!/bin/bash
#SBATCH -J dm_prfmz05t
#SBATCH -o /ptmp/uli/dwarf_merger/logs/prfmz05_test.log
#SBATCH -p general
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=72
#SBATCH --mem=110000
#SBATCH -t 00:40:00
export NUMBA_NUM_THREADS=72 OMP_NUM_THREADS=72 PRFM_ZCOL=0.5 PRFM_ZREF=1 PRFM_TAG=_z05
cd /raven/u/uli/dwarf_merger/prfm
/raven/u/uli/jupyter_env/bin/python patches.py 120 40
