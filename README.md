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
# Part C. Towards an analytic model (candidate, partly verified)

Part C is what can be written as equations with no exponent fitted on this run, then tested against Part A. It is kept apart
from B, which is the empirical link. Status of each step:

1. **Layer dispersion.** Feedback-regulated: the Ostriker & Kim yield. Flow-driven: mixing-length balance of shear injection
   against dissipation, sigma_eff = (2 epsilon_S)^{1/2} H S. This is the shear closure of Part B with its one constant, so it
   enters C as an input, not a prediction.
2. **Cascade to the clump.** Burgers scaling, sigma(l) = sigma_eff (l/H)^{1/2}. Verified: the clump dispersion before onset is
   predicted to a median ratio of 0.95 with 0.2 dex scatter and no phase dependence, no free parameter.
3. **Clump virial parameter.** alpha = (5/3) sigma(R)^2 R / (G M) with the measured clump mass. Verified: measured over
   predicted has a median of 0.9, 0.4 dex scatter, rank correlation 0.8, and phase medians between 0.8 and 1.1. Derived.
4. **Efficiency as the bound fraction.** For a clump with density falling as r^-2, the virial parameter scales with radius,
   so only the inner fraction 1/alpha is bound: epsilon = epsilon_0 / alpha, exponent minus one derived. Measured exponent
   minus 0.8 to minus 0.9; epsilon_0 = 0.03; but epsilon_0 drifts by a factor of 4 to 8 from the first disc phase to the
   remnant, and the clump-to-clump scatter is 0.75 dex. Right form, one constant, normalisation not derived.

So C is analytic from the layer's effective dispersion down to the clump virial parameter, and analytic in form but not in
normalisation for the efficiency. What it explains in Part A: the tenfold rise of the clump virial parameter through the
merger and, through it, the bulk of the fifty-fold efficiency drop; what it does not: a residual factor of 4 to 8 between
phases and the intrinsic scatter.

![clump virial parameter](figs/partC_alpha.png)

![clump efficiency](figs/partC_eff.png)

The physical content, stated without claiming more than the verified steps: PRFM fixes the pressure at the weight, in the
merger as in the disc. Star formation responds to the effective dispersion, pressure over midplane density, because the
clumps inherit it through the cascade and are unbound in proportion to its square. Feedback regulation holds that dispersion
near 9 km/s; when the merger flow and the field carry the pressure it is near 20 km/s at equal or larger weight, the clumps
are ten times less bound, and most of the efficiency drop follows.
