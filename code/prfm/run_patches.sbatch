#!/bin/bash
#SBATCH -J dm_prfm
#SBATCH -o /ptmp/uli/dwarf_merger/logs/prfm_%a.log
#SBATCH -p general
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=72
#SBATCH --mem=110000
#SBATCH -t 04:00:00
export NUMBA_NUM_THREADS=72 OMP_NUM_THREADS=72
cd /raven/u/uli/dwarf_merger/prfm
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
  /raven/u/uli/jupyter_env/bin/python patches.py $SNAPS
else
  S0=$((SLURM_ARRAY_TASK_ID*8)); S1=$((S0+7)); [ $S1 -gt 231 ] && S1=231
  /raven/u/uli/jupyter_env/bin/python patches.py $(seq $S0 $S1)
fi
