#!/bin/bash
# Copy the working code from the cluster into this repo.  The scripts run from /raven/u/uli/dwarf_merger,
# which is where the paths in them point; this directory is the archived copy.
SRC=/raven/u/uli/dwarf_merger; DST=$(dirname "$(readlink -f "$0")")
cp $SRC/common.py $DST/
cp $SRC/prfm/*.py $SRC/prfm/*.sbatch $DST/prfm/
cp $SRC/clusters/*.py $SRC/clusters/*.sbatch $DST/clusters/
rm -f $DST/prfm/*.bak* $DST/clusters/*.bak*
echo "synced from $SRC"
