# merger_prfm

Figures from the PRFM (pressure-regulated, feedback-modulated; Ostriker & Kim 2022) test of the dwarf-galaxy merger
simulation and its link to the star cluster mass function. Code: `/raven/u/uli/dwarf_merger` (prfm/paper_figs.py).
Figures are converted to paper style one by one and replaced in place; the ones below the line still have the analysis-style labels.

## Merger phases
Separation of the two nuclei; the four phases used throughout.

![separation](figs/story_0_sep.png)

## Vertical equilibrium
Total midplane pressure over the weight of the gas, for the full column and for the gas layer within 0.5 kpc of the midplane.

![equilibrium](figs/story_1_equilibrium.png)

## Pressure supplied by feedback
Fraction of the total midplane pressure accounted for by the feedback yields of Ostriker & Kim (2022) at the local star formation rate; the rest is supplied by the flow.

![feedback share](figs/story_2_feedback.png)

## Shear-driven turbulence
Turbulent pressure of patches without star formation in the last 10 Myr against the shear-work term (gas surface density times scale height times squared shear rate); one coefficient of 0.02 fits all phases.

![shear closure](figs/story_3_shear.png)

## Magnetic field
Field strength in the gas layer, mass-weighted mean and median patch. Seeded at 10 nG, it saturates near 1 microgauss by 50 Myr, jumps at coalescence and reaches 10 to 20 microgauss in the remnant, where the magnetic pressure is at or above equipartition with thermal plus turbulent.

![dynamo](figs/story_4_dynamo.png)

## Clouds
Mass function of cold clouds (T < 1000 K, n > 10 cm^-3, friends-of-friends at 3 pc) per phase: the slope of 1.6 does not change through the merger; only the top grows.

![cloud mass function](figs/clouds_mf.png)

Fraction of a cloud's mass turned into stars within 10 Myr against the cloud's mean density, for clouds that form stars. Lines are medians per density bin. The efficiency rises with density and drops 15x from the disc phase to the remnant; its 0.7 to 0.8 dex scatter is what sets the width of the burst-mass distribution.

![cloud efficiency](figs/clouds_eff.png)

---
## Not yet converted

Cluster mass function per merger phase.

![MF per phase](figs/story_5_mf.png)

Largest bound cluster vs young stellar mass of the patch.

![reservoir](figs/story_6_reservoir.png)

Largest bound cluster per 25 Myr vs the weight of the star-forming patches.

![ceiling](figs/story_7_ceiling.png)

Line-of-sight dependence of P/W.

![line of sight](figs/prfm_time_los.png)

Low-mass end convergence of the cluster MF.

![low-mass convergence](figs/lowmass_convergence.png)

Burst-mass distribution at fixed weight, duty cycle, index per bin.

![burst distribution](figs/burst_dist.png)

Mock (real weights x density-threshold duty cycle x lognormal kernel) vs observed burst-mass distribution.

![burst kernel mock](figs/burst_kernel.png)
