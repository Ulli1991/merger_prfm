#!/bin/bash
#SBATCH -J dm_prfmz05
#SBATCH -o /ptmp/uli/dwarf_merger/logs/prfmz05_%a.log
#SBATCH -p general
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=72
#SBATCH --mem=110000
#SBATCH -t 04:00:00
#SBATCH --array=0-28
export NUMBA_NUM_THREADS=72 OMP_NUM_THREADS=72 PRFM_ZCOL=0.5 PRFM_ZREF=1 PRFM_TAG=_z05
cd /raven/u/uli/dwarf_merger/prfm
S0=$((SLURM_ARRAY_TASK_ID*8)); S1=$((S0+7)); [ $S1 -gt 231 ] && S1=231
/raven/u/uli/jupyter_env/bin/python patches.py $(seq $S0 $S1)
