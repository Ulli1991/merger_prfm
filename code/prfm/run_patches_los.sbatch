#!/bin/bash
#SBATCH -J dm_prfmlos
#SBATCH -o /ptmp/uli/dwarf_merger/logs/prfmlos_%a.log
#SBATCH -p general
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=72
#SBATCH --mem=110000
#SBATCH -t 03:00:00
#SBATCH --array=0-7
export NUMBA_NUM_THREADS=72 OMP_NUM_THREADS=72 PRFM_LOS=1
cd /raven/u/uli/dwarf_merger/prfm
# every 4th snapshot (58), task k takes those with index = k mod 8  -> patch_XXX_los.h5
SNAPS=$(seq 0 4 231 | awk -v k=$SLURM_ARRAY_TASK_ID 'NR%8==k+1 || (k==7 && NR%8==0)')
echo "task $SLURM_ARRAY_TASK_ID snaps: $SNAPS"
/raven/u/uli/jupyter_env/bin/python patches.py $SNAPS
