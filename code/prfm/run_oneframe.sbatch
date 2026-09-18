#!/bin/bash
#SBATCH -J dm_oneframe
#SBATCH -o /ptmp/uli/dwarf_merger/logs/oneframe_%a.log
#SBATCH -p general
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=72
#SBATCH --mem=110000
#SBATCH -t 05:00:00
#SBATCH --array=0-7
# single merged frame from the second passage on (t >= 106 Myr, snapshots 108..231): full column (tag _one), then the layer
# cut about its midplanes (tag _one_z05).  8 subsets of the 124 snapshots; each task does full then layer for its subset.
export NUMBA_NUM_THREADS=72 OMP_NUM_THREADS=72 PRFM_ONEFRAME=106
K=$SLURM_ARRAY_TASK_ID
cd /raven/u/uli/dwarf_merger/prfm
SNAPS=$(seq 108 231 | awk -v k=$K 'NR%8==k+1 || (k==7 && NR%8==0)')
echo "subset $K snaps: $SNAPS"
PRFM_ZCOL=1.5 PRFM_TAG=_one /raven/u/uli/jupyter_env/bin/python patches.py $SNAPS
PRFM_ZCOL=0.5 PRFM_ZREF=1 PRFM_REF_TAG=_one PRFM_TAG=_one_z05 /raven/u/uli/jupyter_env/bin/python patches.py $SNAPS
