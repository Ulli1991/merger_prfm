#!/bin/bash
#SBATCH -J dm_prfmgrid
#SBATCH -o /ptmp/uli/dwarf_merger/logs/prfmgrid_%a.log
#SBATCH -p general
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=72
#SBATCH --mem=110000
#SBATCH -t 04:00:00
#SBATCH --array=0-39
# 5 configurations x 8 snapshot subsets (every 4th snapshot, 58 snapshots): task = cfg*8 + k
export NUMBA_NUM_THREADS=72 OMP_NUM_THREADS=72 PRFM_ZREF=1
CFG=$((SLURM_ARRAY_TASK_ID / 8)); K=$((SLURM_ARRAY_TASK_ID % 8))
case $CFG in
  0) export PRFM_LOS=1 PRFM_ZCOL=0.5  PRFM_PATCH=0.5   PRFM_TAG=_z05 ;;         # layer weight per line of sight -> patch_XXX_los_z05.h5
  1) export PRFM_ZCOL=0.25  PRFM_PATCH=0.5   PRFM_TAG=_p05z025 ;;
  2) export PRFM_ZCOL=0.25  PRFM_PATCH=0.25  PRFM_TAG=_p025z025 ;;
  3) export PRFM_ZCOL=0.125 PRFM_PATCH=0.25  PRFM_TAG=_p025z0125 ;;
  4) export PRFM_ZCOL=0.125 PRFM_PATCH=0.125 PRFM_TAG=_p0125z0125 ;;
esac
cd /raven/u/uli/dwarf_merger/prfm
SNAPS=$(seq 0 4 231 | awk -v k=$K 'NR%8==k+1 || (k==7 && NR%8==0)')
echo "cfg $CFG tag $PRFM_TAG patch $PRFM_PATCH zcol $PRFM_ZCOL los=${PRFM_LOS:-0} snaps: $SNAPS"
/raven/u/uli/jupyter_env/bin/python patches.py $SNAPS
