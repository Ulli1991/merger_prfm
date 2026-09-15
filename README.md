# merger_prfm

Key figures from the PRFM (pressure-regulated, feedback-modulated; Ostriker & Kim 2022) test of the dwarf-galaxy merger
simulation and its link to the star cluster mass function. Analysis code lives in `/raven/u/uli/dwarf_merger` (prfm/, clusters/);
data products in `/ptmp/uli/dwarf_merger`.

## figs/
| file | what it shows |
|---|---|
| story_0_sep.png | nuclear separation vs time, phase boundaries |
| story_1_equilibrium.png | P_tot / W vs time, full column vs 0.5 kpc gas layer |
| layer_verdict_z05.png | exact layer-weight verdict per phase (Sigma-weighted median P/W, ratio of sums) |
| story_2_feedback.png | fraction of the midplane pressure supplied by feedback (OK22 yields) vs by the flow |
| story_3_shear.png | patch-scale shear closure for quiescent patches, epsilon_S = 0.02 |
| story_4_dynamo.png | magnetic field growth and P_mag / (P_th + P_turb) |
| sources_dynamo.png | pressure sources and dynamo, multi-panel version |
| prfm_time_los.png | P/W for different line-of-sight normals (disc is a layer only along the angular-momentum axis) |
| story_5_mf.png, mf_phases.png | cluster mass function per merger phase |
| story_6_reservoir.png | largest cluster vs young stellar mass of the patch (M_max ~ 0.5 M_young) |
| story_7_ceiling.png, ceiling.png, cluster_link.png | largest bound cluster per 25 Myr vs layer weight of the star-forming patches |
| lowmass_convergence.png | local MF slope per phase vs N_min, linking length and bound flag (robust above 300 Msun) |
| burst_dist.png | burst mass at fixed layer weight (lognormal), duty cycle per phase, running index per bin |
| burst_kernel.png | mock (real weights x density-threshold duty cycle x lognormal kernel) vs observed burst-mass distribution |
