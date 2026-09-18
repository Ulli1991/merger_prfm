#!/bin/bash
#SBATCH -J dm_rebuild
#SBATCH -o /ptmp/uli/dwarf_merger/logs/rebuild_oneframe.log
#SBATCH -p interactive
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8
#SBATCH --mem=40000
#SBATCH -t 01:55:00
# rebuild every patch-based product on the mixed directory (two-frame < snap 108, one-frame >= 108) and the new catalogues
export PRFM_DIR=/ptmp/uli/dwarf_merger/prfm_one OMP_NUM_THREADS=8 NUMBA_NUM_THREADS=8
cd /raven/u/uli/dwarf_merger
PY=/raven/u/uli/jupyter_env/bin/python; LG=/ptmp/uli/dwarf_merger/logs
run() { echo "== $1 $2 $(date +%H:%M)"; $PY $1 $2 2>&1 | grep -v Warning > $LG/$(basename $1 .py)$( [ -n "$2" ] && echo _$2 )_one.log; }
run prfm/layer_verdict.py; run prfm/bfield_time.py; run prfm/gamma_time.py
run prfm/cluster_pressure.py; run prfm/cluster_link.py; run prfm/cluster_env_layer.py; run prfm/two_criteria.py; run prfm/ceiling_test.py
run prfm/burst_dist.py; run prfm/burst_kernel.py; run clusters/lowmass_convergence.py
run clusters/lineage_env.py cl100; run clusters/lineage_env.py cl10; run clusters/step4_fit.py; run clusters/partC_mf.py
$PY prfm/paper_figs.py all 2>&1 | grep -v Warning > $LG/paper_figs_one.log
echo REBUILDDONE $(date +%H:%M)
