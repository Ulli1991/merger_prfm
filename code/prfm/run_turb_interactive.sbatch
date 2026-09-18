#!/bin/bash
#SBATCH -J dm_turb
#SBATCH -o /ptmp/uli/dwarf_merger/logs/turb.log
#SBATCH -p interactive
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=72
#SBATCH --mem=250000
#SBATCH -t 01:59:00
export OMP_NUM_THREADS=72
cd /raven/u/uli/dwarf_merger/prfm
SNAPS="20 40 60 80 100 110 115 120 130 140 150 160 170 180 190 200 210 220"
/raven/u/uli/jupyter_env/bin/python turb_decomp.py $SNAPS
/raven/u/uli/jupyter_env/bin/python structure_function.py $SNAPS
