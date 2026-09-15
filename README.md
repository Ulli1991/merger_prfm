# merger_prfm

Does pressure-regulated, feedback-modulated star formation (PRFM; Ostriker & Kim 2022) survive a galaxy merger, and what does
it imply for the masses of the star clusters that form? These figures come from one high-resolution simulation of two
gas-rich dwarf galaxies merging (4 Msun gas and star particles, 0.4 pc softening, magnetic fields, resolved supernovae),
followed for 226 Myr through first passage, a second passage, coalescence and the settling of the remnant.

## The idea in one paragraph

PRFM says that in a galactic disc the midplane gas pressure adjusts to balance the weight of the gas column in the
gravitational field, and that star formation supplies that pressure through feedback: the weight sets the pressure, the
pressure sets the star formation rate. We test both halves in 0.5 kpc patches of the merger, then follow the consequences
downward in scale. The star formation rate of a patch is really a sequence of bursts; each burst is the collapse of a single
cold cloud; about half of each burst ends up in one star cluster. So the weight of the gas layer fixes how big the largest
star-forming event can be, and with it the top of the cluster mass function, while the slope of the mass function is set by
how efficiently clouds convert gas to stars, a quantity that the merger changes by an order of magnitude.

The figures are arranged along that chain: merger phases, equilibrium, pressure sources, clouds, clusters.

## Merger phases
Separation of the two nuclei with time. We split the run into four phases used in every other figure: before first passage,
between passages, coalescence (the first time the nuclei come within 0.5 kpc, solid line) and the remnant.

![separation](figs/story_0_sep.png)

## Vertical equilibrium
The first half of PRFM: does the midplane pressure balance the weight? The ratio of total midplane pressure (thermal plus
turbulent plus magnetic stress) to weight is shown for two definitions of the weight. Using the full 3 kpc column, the
ratio drops to 0.3 after coalescence and the model looks broken. Restricting the weight to the gas within 0.5 kpc of the
midplane, where the star formation happens, the ratio stays at one to within 5 % in the disc phase and the remnant, and
falls only 20 to 35 % short during coalescence. The apparent failure is the dark-matter weight of tidal debris streaming
through the disc at 25 to 30 km/s, which is not a layer in any orientation and does not form stars.

![equilibrium](figs/story_1_equilibrium.png)

## Pressure supplied by feedback
The second half of PRFM: is the pressure supplied by star formation? Using the feedback yields of Ostriker & Kim (2022) at
the local star formation rate, feedback accounts for 50 to 95 % of the pressure before coalescence but almost none of it
in the quiet intervals after. The gas is still in vertical balance then, so something else is supplying the pressure: the
large-scale flow of the merger. The regulation half of PRFM therefore holds only while star formation is on.

![feedback share](figs/story_2_feedback.png)

## Shear-driven turbulence
What supplies the pressure between bursts. For patches with no star formation in the last 10 Myr, the turbulent pressure
follows the shear-work term (gas surface density times scale height times the squared patch-scale shear rate) with one
coefficient of 0.02 in every phase and no offset. Feedback-only closures miss these patches by up to two orders of
magnitude. This is the missing term in the PRFM energy balance for a merger.

![shear closure](figs/story_3_shear.png)

## Magnetic field
The field is seeded at 10 nG and amplified by the turbulence: it e-folds every 10 Myr, saturates near 1 microgauss by
50 Myr, jumps by a factor of five at coalescence and grows to 10 to 20 microgauss in the remnant, where the magnetic
pressure reaches equipartition with thermal plus turbulent pressure. It is treated as a measured contribution to the
vertical support, not modelled.

![dynamo](figs/story_4_dynamo.png)

## Clouds
Stars keep the identity of the gas particle they formed from, so every young star can be traced to the cold cloud it came
from. Nearly all stars formed within 10 Myr of a snapshot come from gas already sitting in a cold cloud above 10 cm^-3,
and in 85 to 95 % of star-forming events a single cloud supplies essentially all the stars. A burst is the collapse of one
cloud.

The mass function of those clouds (T < 1000 K, n > 10 cm^-3, friends-of-friends at 3 pc) has a slope of 1.6 and does not
change through the merger. Only the most massive end grows, as the remnant assembles cold complexes of a million solar
masses and more.

![cloud mass function](figs/clouds_mf.png)

What does change is how efficiently clouds turn gas into stars. The fraction of a cloud's mass converted within 10 Myr
rises with the cloud's mean density and drops fifteen-fold from the disc phase to the remnant, where the clouds are more
turbulent and less bound. Its scatter of 0.7 to 0.8 dex is what produces the 0.4 dex spread of burst masses at fixed
environment, and through it the slope of the cluster mass function. The physics of the slope lives here, at the cloud
scale, not in the disc-scale weight distribution.

![cloud efficiency](figs/clouds_eff.png)

## Cluster mass function
Bound star clusters younger than 10 Myr, per phase, with the maximum-likelihood power-law index above 300 Msun. The slope
flattens from 1.9 before the merger to 1.55 in the remnant, and the most massive cluster grows from 1e4 to over 1e5 Msun.
Both trends are robust to the minimum group size, the linking length and the bound flag. Below 300 Msun in the coalescence
and remnant phases the counts drop because only a third of the small groups are bound, so the function is quoted from
300 Msun up. Since the cloud mass function does not change (above), the flattening and the growing top come from the
star formation efficiency of the clouds.

![MF per phase](figs/story_5_mf.png)

## One dominant cluster per burst
For every 0.5 kpc patch that formed more than 500 Msun of stars in the last 10 Myr and hosts a bound cluster, the mass of
its most massive cluster against the total mass of young stars in the patch. The median runs along half the young mass up
to about 1e5 Msun: each burst makes one dominant cluster that captures about half of the stars formed, with a handful of
smaller companions. This is the link between the star formation rate of a patch and the top of the cluster mass function,
and it is why the weight of the gas layer, which sets the size of the largest burst, also sets the mass of the largest
cluster. The most massive bursts, in the remnant nucleus, fall below the line because their stars are spread over several
clusters.

![reservoir](figs/story_6_reservoir.png)

## The weight sets the ceiling
The most massive bound cluster formed in each 25 Myr interval against W, the weight of the gas layer (within 0.5 kpc of
the midplane) of the heaviest star-forming patches in that interval (90th percentile over patches). Labels give the start of each interval in
Myr. The line is not a fit: it is half of the PRFM star formation rate at that weight, times the patch area and a 10 Myr
burst, i.e. the size of the largest burst a patch of that weight can produce, halved by the capture fraction of the
previous figure. Seven active intervals spanning a factor of 30 in weight follow it with 0.1 dex scatter. The one quiescent
interval, 160 to 185 Myr, has the highest weight of all and sits 100 times below: the weight sets the ceiling of the
cluster mass function, and the duty cycle of star formation decides whether it is reached.

![ceiling](figs/story_7_ceiling.png)

---
## Not yet converted

Line-of-sight test: the pre-merger disc is an equilibrium layer only along its angular-momentum axis; after coalescence no
orientation makes the debris a layer.

![line of sight](figs/prfm_time_los.png)

Convergence of the low-mass end of the cluster mass function with the minimum group size, linking length and bound flag.

![low-mass convergence](figs/lowmass_convergence.png)

Burst mass at fixed weight (lognormal, 0.4 dex wide), the duty cycle per phase, and the running index per bin.

![burst distribution](figs/burst_dist.png)

Mock built from the real patch weights, a density-threshold duty cycle and the lognormal kernel, against the observed
burst-mass distribution: the slope is reproduced in every phase, the top is not.

![burst kernel mock](figs/burst_kernel.png)
