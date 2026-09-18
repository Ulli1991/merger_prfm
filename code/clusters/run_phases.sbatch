#!/bin/bash
#SBATCH -J dm_phases2
#SBATCH -o /ptmp/uli/dwarf_merger/logs/phases_chain.log
#SBATCH -p interactive
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4
#SBATCH --mem=24000
#SBATCH -t 01:30:00
cd /raven/u/uli/dwarf_merger
PY=/raven/u/uli/jupyter_env/bin/python; LG=/ptmp/uli/dwarf_merger/logs
run() { echo "== $1"; $PY $1 $2 2>&1 | grep -v Warning > $LG/$(basename $1 .py)$( [ -n "$2" ] && echo _$2 ).log; }
run prfm/burst_dist.py; run prfm/burst_kernel.py; run clusters/lowmass_convergence.py; run clusters/clouds_analysis.py 10
run clusters/lineages.py cl100; run clusters/lineages.py cl10; run clusters/lineage_env.py cl100; run clusters/lineage_env.py cl10
run clusters/step4_fit.py; run clusters/partC_mf.py; run prfm/cluster_env_layer.py; run prfm/two_criteria.py; run prfm/ceiling_test.py
$PY prfm/paper_figs.py all 2>&1 | grep -v Warning > $LG/paper_figs.log
echo PHASESDONE
