#!/bin/bash
#SBATCH -J dm_clmf
#SBATCH -o /ptmp/uli/dwarf_merger/logs/clmf.log
#SBATCH -p general
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8
#SBATCH --mem=40000
#SBATCH -t 01:55:00
export OMP_NUM_THREADS=8
cd /raven/u/uli/dwarf_merger/clusters
/raven/u/uli/jupyter_env/bin/python cluster_mf.py ${ARGS}
