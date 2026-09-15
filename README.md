# merger_prfm

Does pressure-regulated, feedback-modulated star formation (PRFM; Ostriker & Kim 2022) survive a galaxy merger, and what does
it imply for the masses of the star clusters that form? One high-resolution simulation of two gas-rich dwarf galaxies merging
(4 Msun gas and star particles, 0.4 pc softening, magnetic fields, resolved supernovae), followed for 226 Myr through first
passage, a second passage, coalescence and the settling of the remnant.

This page keeps two things apart. **Part A** is what the simulation shows: measurements, stated without a model behind them.
**Part B** is a semi-empirical model built to explain Part A: PRFM plus ingredients that are either published theory or
constants calibrated on this run, each tested against one measurement, with the outcome stated including where it fails.

---
# Part A. What the simulation shows

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

## Feedback pressure
The second half of PRFM: is the pressure supplied by star formation? Using the feedback yields of Ostriker & Kim (2022) at
the local star formation rate, feedback accounts for 50 to 95 % of the pressure before coalescence but almost none of it
in the quiet intervals after. The gas is still in vertical balance then, so something else is supplying the pressure: the
large-scale flow of the merger. The regulation half of PRFM therefore holds only while star formation is on.

![feedback share](figs/story_2_feedback.png)

## Magnetic field
The field is seeded at 10 nG and amplified by the turbulence: it e-folds every 10 Myr, saturates near 1 microgauss by
50 Myr, jumps by a factor of five at coalescence and grows to 10 to 20 microgauss in the remnant, where the magnetic
pressure reaches equipartition with thermal plus turbulent pressure. It is treated as a measured contribution to the
vertical support, not modelled.

![dynamo](figs/story_4_dynamo.png)

## Clouds
Stars keep the identity of the gas particle they formed from, so every young star can be traced to the cold gas it came
from. Nearly all stars formed within 10 Myr of a snapshot come from gas already sitting in a cold cloud above 10 cm^-3.
Those clouds are 30 to 50 pc complexes of 1e4 to 1e5 Msun, at a mean density near 10 cm^-3, the low-density giant molecular
clouds of a dwarf galaxy. In 85 to 95 % of star-forming events a single complex supplies essentially all the stars, and
inside it the stars come from one to three dense clumps above 100 cm^-3 (a few thousand Msun, radii of 4 to 8 pc). A burst
is the collapse of one such complex; the largest bursts, and the remnant nucleus, draw on many clumps.

The mass function of the complexes (T < 1000 K, n > 10 cm^-3, friends-of-friends at 3 pc) has a slope of 1.6 and does not
change through the merger. Only the most massive end grows, as the remnant assembles cold complexes of a million solar
masses and more.

![cloud mass function](figs/clouds_mf.png)

What does change is how efficiently the clumps turn gas into stars. Following each dense clump through the snapshots by its
particle IDs, the figure shows the fraction of its peak mass that it ever turned into stars, against its virial parameter
just before star formation began. More strongly bound clumps convert more, by a factor of 40 over the observed range, and
the efficiency drops fifty-fold from the disc phase to the remnant, where the clumps are more turbulent. But at the same
virial parameter the phases are offset by more than an order of magnitude, so the kinetic virial parameter is not the whole
story, and at fixed state the efficiency still scatters by 0.7 dex. That scatter is what produces the 0.4 dex spread of
burst masses and the slope of the cluster mass function. This is the one link in the chain that is measured rather than
explained.

![cloud efficiency](figs/clouds_eff.png)

## Cluster mass function
Bound star clusters younger than 10 Myr, per phase, with the maximum-likelihood power-law index above 300 Msun. The slope
flattens from 1.85 before the merger to 1.55 in the remnant, and the most massive cluster grows from 1e4 to over 1e5 Msun.
Both trends are robust to the minimum group size, the linking length and the bound flag. Below 300 Msun in the coalescence
and remnant phases the counts drop because only a third of the small groups are bound, so the function is quoted from
300 Msun up. Since the cloud mass function does not change (above), the flattening and the growing top come from the
star formation efficiency of the clouds.

![MF per phase](figs/story_5_mf.png)

## Largest cluster per burst
For every 0.5 kpc patch that formed more than 500 Msun of stars in the last 10 Myr and hosts a bound cluster, the mass of
its most massive cluster against the total mass of young stars in the patch. The median runs along half the young mass up
to about 1e5 Msun: each burst makes one dominant cluster that captures about half of the stars formed, with a handful of
smaller companions. This is the link between the star formation rate of a patch and the top of the cluster mass function,
and it is why the weight of the gas layer, which sets the size of the largest burst, also sets the mass of the largest
cluster. The most massive bursts, in the remnant nucleus, fall below the line because their stars are spread over several
clusters.

![reservoir](figs/story_6_reservoir.png)

## Every cluster against its patch
All 580 bound young clusters against the weight, total pressure and star formation rate of the 0.5 kpc patch they formed in. The
dashed line in each panel is half the burst that PRFM allows at that weight or pressure, or half the stars actually formed in the
patch. Individual cluster masses fill two decades below the line and correlate only weakly with the patch state; the line is an
upper envelope that 9 % of clusters exceed. The weight bounds the cluster mass, it does not set it cluster by cluster.

![cluster environment](figs/cluster_env.png)

## Line of sight
The equilibrium test repeated with the column taken along other directions: the disc's angular-momentum axis, and normals
tilted by 30, 60 and 90 degrees from it. Before coalescence only the disc normal gives a ratio near one, so the disc is a
layer and the chosen normal is the right one. After coalescence every direction gives 0.2 to 0.5: the deficit does not
depend on viewing angle, so the streaming debris is not a layer in any orientation.

![line of sight](figs/prfm_time_los.png)

## Low-mass end of the cluster mass function
Local power-law slope of the cluster mass function in three mass bins, per phase, for all friends-of-friends groups
(filled) and for bound groups only (open). Above 300 Msun the two agree and the slopes are stable against the minimum group
size and linking length. In the 100 to 300 Msun bin after coalescence the bound-only slope collapses, because only a third
of those small groups are bound: that end of the function is set by the bound flag, not by the physics, and is not used.

![low-mass convergence](figs/lowmass_convergence.png)


---
# Part B. A semi-empirical model and its tests

This is not a theory. It is PRFM extended with ingredients that are either taken from published work or calibrated on this
simulation, and the table says which is which. A closed theory would derive the calibrated numbers; until it does, the
model is a description with a small number of measured constants, and its value is in showing which measured constant
carries which effect.

| ingredient | what it does | source | calibrated on this run |
|---|---|---|---|
| vertical equilibrium, PRFM rate | pressure balances weight; star formation rate at a given weight | Ostriker & Kim 2022 | nothing (their fits used as published) |
| flow-driving term | turbulent pressure of patches without star formation | shear-work closure, form assumed | one coefficient, epsilon_S = 0.02 |
| magnetic support | third pressure component | measured field | nothing modelled: the field is taken from the simulation |
| density threshold | which patches form stars | logistic form assumed | threshold and width (2 numbers) |
| burst size | how much a patch forms when on | PRFM rate, lognormal scatter assumed | normalisation 0.4 and width 0.4 to 0.5 dex per phase |
| cloud mass function | population of star-forming complexes | gravitational fragmentation, slope 1.6 | measured; not derived here |
| clump efficiency | fraction of a clump turned into stars | Federrath & Klessen 2012 form under test | open: virial-parameter trend measured, phase offset unexplained |
| capture fraction | share of a burst in one cluster | measured 0.5 | measured; not derived |

Each ingredient is tested against one measurement of Part A.

## Pressure budget: the flow term
What supplies the pressure between bursts. For patches with no star formation in the last 10 Myr, the turbulent pressure
follows the shear-work term (gas surface density times scale height times the squared patch-scale shear rate) with one
coefficient of 0.02 in every phase and no offset. Feedback-only closures miss these patches by up to two orders of
magnitude. This is the missing term in the PRFM energy balance for a merger.

![shear closure](figs/story_3_shear.png)

## Rate when on: the PRFM burst as the ceiling of cluster masses
The most massive bound cluster formed in each 25 Myr interval against $\mathcal{W}$, the weight of the gas layer (within 0.5 kpc of
the midplane) of the heaviest star-forming patches in that interval (90th percentile over patches). Labels give the start of each
interval in Myr. The line is not a fit: it is half of the PRFM star formation rate at that weight, times the patch area and a 10 Myr
burst, i.e. the size of the largest burst a patch of that weight can produce, halved by the capture fraction of the previous figure.
Seven active intervals spanning a factor of 30 in weight follow it (slope 0.94 against the median weight, 0.09 dex scatter). The one
quiescent interval, 160 to 185 Myr, has the highest weight of all and sits 100 times below. The weight sets the cluster mass while
star formation is on; whether it is on is a separate question, answered below.

![ceiling](figs/story_7_ceiling.png)

## Switch: the density threshold
Probability that a 0.5 kpc patch formed more than 500 Msun of stars in the last 10 Myr, as a function of its weight $\mathcal{W}$
(gas layer within 0.5 kpc of the midplane), per phase. Points are measured fractions, lines logistic fits. The weight at
which half the patches are active rises twenty-fold through the merger. The same threshold expressed in midplane density
moves by less than a factor of three, so the physical threshold is in density: after coalescence the higher turbulent and
magnetic support means the same weight produces a lower midplane density.

![duty cycle](figs/duty_cycle.png)

## Burst sizes: threshold times lognormal
Distribution of the 10 Myr burst mass of the patches between the two passages, against a model with no free shape: each
patch bursts with the probability from the duty cycle above, and its burst mass is drawn from a lognormal of 0.5 dex width
around a mean that rises with the weight. The model reproduces the distribution and its curvature (the same holds in the
other phases). The slope of the burst-mass distribution, and through it of the cluster mass function, is the width of this
lognormal, which the cloud figures above trace to the scatter in cloud star formation efficiency.

![burst kernel](figs/burst_kernel.png)

## Efficiency law
Outcome of the test. Adding magnetic support to the clump virial parameter (kinetic times one plus the inverse squared Alfvén
Mach number) halves the phase offset of the efficiency at fixed virial parameter, from 0.7 to 0.5 dex, and does not reduce the
0.7 dex clump-to-clump scatter. The Federrath & Klessen (2012) multi-freefall efficiency evaluated on the measured clump states
(density-PDF width, virial parameter, Mach number, plasma beta) does not describe the data: rank correlation 0.2 to 0.5,
scatter 1.5 to 2 dex, and it drives the predicted efficiency to zero in the remnant. So within Part B the efficiency stays a
calibrated ingredient. Where it comes from is Part C.

---
# Part C. The physical link

Part C is what follows from physics without a constant fitted on this run, tested against the measurements of Part A. It is
one chain of three steps, and the first two hold with no phase dependence.

**1. The clumps inherit the layer's turbulence.** Following each dense clump to the snapshot before it forms stars, its
velocity dispersion is set by the effective dispersion of the 0.5 kpc patch around it, the square root of midplane pressure
over midplane density, scaled down to the clump size: the clump dispersion goes as the patch dispersion to the power 0.9
times the size ratio to the power 0.6, with 0.19 dex scatter and phase residuals within 0.02 dex. This is a turbulent cascade
from the patch scale to the clump scale, and it does not care whether the patch turbulence is driven by feedback or by the
merger flow.

**2. Therefore the clump boundedness is set by the layer's effective dispersion.** The clumps sit at the star-formation
density, so their virial parameter goes as the square of their inherited dispersion, i.e. as the patch pressure over the
midplane density. All four phases lie on one line with a slope of 1.07 and phase offsets of 0.2 dex. Vertical equilibrium is
irrelevant to it: the virial parameter does not correlate with pressure over weight at all.

![clump virial parameter](figs/partC_alpha.png)

**3. And the efficiency follows.** The fraction of a clump turned into stars falls as the layer's effective dispersion to the
power minus 1.7, with the magnetic share of the pressure adding a weaker factor; the phase offsets shrink to 0.2 to 0.3 dex and
the remaining 0.7 dex is clump-to-clump scatter that no patch or clump quantity predicts. This third step is measured, not
derived: the exponent of efficiency against virial parameter (about minus 0.8) is not given by an existing theory that we
tested.

![clump efficiency](figs/partC_eff.png)

**What this says about PRFM in a merger.** PRFM fixes the midplane pressure at the weight, and it does so in the merger too.
But star formation does not respond to the pressure; it responds to the effective dispersion, pressure over midplane density.
Under feedback regulation the yield fixes that dispersion near 9 km/s, clumps are marginally bound and convert 3 to 8 % of
their mass. When the merger flow and the saturated field supply the pressure instead, the dispersion is 19 to 21 km/s at
equal or larger weight, the clump virial parameter rises tenfold, and the efficiency drops to a tenth of a percent. The entire
effect of the merger on star formation is a factor of two in the effective dispersion, squared into the virial parameter and
then into the efficiency. The quiescent interval is the point where that dispersion is highest for its weight. On the cluster
side this is why the cluster mass function flattens and its top grows while the cloud population does not change: the same
factor of two, acting through the efficiency of every clump.
