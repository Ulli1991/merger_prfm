#!/bin/bash
#SBATCH -J dm_clouds
#SBATCH -o /ptmp/uli/dwarf_merger/logs/clouds.log
#SBATCH -p interactive
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8
#SBATCH --mem=30000
#SBATCH -t 01:55:00
export OMP_NUM_THREADS=8
cd /raven/u/uli/dwarf_merger/clusters
for k in $(seq 1 221); do
  /raven/u/uli/jupyter_env/bin/python clouds.py $k 2>&1 | grep -v Warning
done
echo ALLDONE
