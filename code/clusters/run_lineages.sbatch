#!/bin/bash
#SBATCH -J dm_lineages
#SBATCH -o /ptmp/uli/dwarf_merger/logs/lineages_chain.log
#SBATCH -p interactive
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4
#SBATCH --mem=24000
#SBATCH -t 01:00:00
cd /raven/u/uli/dwarf_merger/clusters
PY=/raven/u/uli/jupyter_env/bin/python
$PY lineages.py cl100 2>&1 | grep -v Warning > /ptmp/uli/dwarf_merger/logs/lineages_cl100.log
$PY lineages.py cl10  2>&1 | grep -v Warning > /ptmp/uli/dwarf_merger/logs/lineages_cl10.log
$PY eff_law2.py 2>&1 | grep -v Warning > /ptmp/uli/dwarf_merger/logs/eff_law2.log
$PY lineage_env.py cl100 2>&1 | grep -v Warning > /ptmp/uli/dwarf_merger/logs/lineage_env_cl100.log
$PY lineage_env.py cl10 2>&1 | grep -v Warning > /ptmp/uli/dwarf_merger/logs/lineage_env_cl10.log
echo CHAINDONE
