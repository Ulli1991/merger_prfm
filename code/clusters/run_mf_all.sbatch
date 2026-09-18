#!/bin/bash
#SBATCH -J dm_clmf
#SBATCH -o /ptmp/uli/dwarf_merger/logs/clmf_all.log
#SBATCH -p general
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=72
#SBATCH --mem=110000
#SBATCH -t 03:00:00
# all six cluster catalogues on one node, 12 threads each, every group energy-tested
export OMP_NUM_THREADS=12 NUMBA_NUM_THREADS=12
cd /raven/u/uli/dwarf_merger/clusters
PY=/raven/u/uli/jupyter_env/bin/python; LG=/ptmp/uli/dwarf_merger/logs
$PY cluster_mf.py --step 1 --age 10 --link 5 --nmin 25 > $LG/clmf_main.log 2>&1 &
$PY cluster_mf.py --step 1 --agemin 20 --age 50 --link 5 --nmin 25 > $LG/clmf_old.log 2>&1 &
$PY cluster_mf.py --step 10 --age 10 --link 3 --nmin 10 > $LG/clmf_l3n10.log 2>&1 &
$PY cluster_mf.py --step 10 --age 10 --link 5 --nmin 10 > $LG/clmf_l5n10.log 2>&1 &
$PY cluster_mf.py --step 10 --age 10 --link 5 --nmin 50 > $LG/clmf_l5n50.log 2>&1 &
$PY cluster_mf.py --step 10 --age 10 --link 8 --nmin 10 > $LG/clmf_l8n10.log 2>&1 &
wait; echo CLMFDONE
