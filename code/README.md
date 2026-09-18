# Code

Everything used to produce the figures in the top-level README. The scripts run from
`/raven/u/uli/dwarf_merger` on the cluster, which is where their relative imports and output paths point;
this directory is the archived copy, refreshed with `./sync.sh`.

Python 3.13 with numpy, scipy, h5py, numba and matplotlib. Snapshots are read from the path in `common.py`
(`SNAPDIR`), intermediate products are written under `DATADIR` (`/ptmp/uli/dwarf_merger`), and the figures go to
`figs/` of this repository.

## Order of the pipeline

| step | script | what it makes | where it runs |
|---|---|---|---|
| column catalogue | `prfm/patches.py` | `patch_KKK*.h5`: per 0.5 kpc column, the pressures, the weight from the particle-mesh solve, the star formation rates, the velocity-gradient diagnostics | `prfm/run_patches*.sbatch`, one node, about 4.5 min per snapshot |
| one merged frame | same, with `PRFM_ONEFRAME=106` | `patch_KKK_one*.h5` for the snapshots after the second passage | `prfm/run_oneframe.sbatch` |
| mixed directory | symlinks | `/ptmp/uli/dwarf_merger/prfm_one`: two-frame files before snapshot 108, one-frame files after; every reading script honours `PRFM_DIR` | shell |
| gravity | `prfm/gravity.py` | the two-level particle-mesh solver used by `patches.py` | imported |
| clouds and clumps | `clusters/clouds.py`, `clusters/fof.py` | `clouds_KKK.npz`: friends-of-friends groups of cold and dense gas with their members and the stars formed from them | `clusters/run_clouds.sbatch` |
| clump histories | `clusters/lineages.py`, `clusters/lineage_env.py` | lineages linked by member overlap, joined to the column they sit in | `clusters/run_lineages.sbatch` |
| star groups | `clusters/cluster_mf.py` | `clusters_age*.npz`: friends-of-friends groups of young stars, energy-tested, with the Subfind-style self-bound mass | `clusters/run_mf_all.sbatch` |
| joins and fits | `prfm/layer_verdict.py`, `prfm/cluster_pressure.py`, `prfm/cluster_link.py`, `prfm/cluster_env_layer.py`, `clusters/step4_fit.py`, `clusters/partC_mf.py` and the rest | the `.npz` files the figures read | `prfm/run_rebuild_oneframe.sbatch` runs the whole chain |
| figures | `prfm/paper_figs.py` | every figure of the top-level README; `python paper_figs.py all` or one name | login node, a few minutes |

## Plot style and checks

`prfm/paperstyle.py` holds the style (serif, ticks inside on four sides, single-column width), the phase bands, and
`place_legend`, which puts a legend where nothing is plotted and extends the axis if it has to. Two scripts check the
result over every figure: `prfm/check_legends.py` counts plotted points under each legend, and `prfm/check_bounds.py`
reports any text or legend crossing its panel edge. Both should print zero.

## One-off tests kept for the record

`prfm/hot_profile.py` (hot gas against height), `prfm/stream_pressure.py` (how much of the turbulent pressure is
coherent streaming between the two progenitors), `clusters/nucleus_bound.py` (energy test of the largest groups),
`prfm/normal_test.py` and `prfm/cubes*.py` (the frame-free approach, shelved), `prfm/scale_los.py` (column size and
line of sight).
