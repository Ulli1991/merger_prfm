#!/bin/bash
#SBATCH -J dm_normtest
#SBATCH -o /ptmp/uli/dwarf_merger/logs/normal_test.log
#SBATCH -p interactive
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8
#SBATCH --mem=60000
#SBATCH -t 00:40:00
cd /raven/u/uli/dwarf_merger/prfm
/raven/u/uli/jupyter_env/bin/python normal_test.py 100 2>&1 | grep -v Warning
/raven/u/uli/jupyter_env/bin/python normal_test.py 60 2>&1 | grep -v Warning
