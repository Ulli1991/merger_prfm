# merger_prfm

Does pressure-regulated, feedback-modulated star formation (PRFM; Ostriker & Kim 2022) survive a galaxy merger, and what does
it imply for the masses of the star clusters that form? One high-resolution simulation of two gas-rich dwarf galaxies merging
(4 Msun gas and star particles, 0.4 pc softening, magnetic fields, resolved supernovae), followed for 226 Myr through first
passage, a second passage, coalescence and the settling of the remnant.

This page keeps two things apart. **Part A** is what the simulation shows: measurements, stated without a model behind them.
**Part B** is a semi-empirical model built to explain Part A: PRFM plus ingredients that are either published theory or
constants calibrated on this run, each tested against one measurement, with the outcome stated including where it fails;
it ends with a physical interpretation of the chain, still built on measured inputs. **Part C** is reserved for a closed
calculation from the layer state to the cluster mass function, which does not exist yet; it states what it would have to
contain.

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

## Pressure per unit star formation
What is plotted: the measured total midplane pressure of the gas layer divided by the measured star formation rate
surface density of the same patches (stars younger than 10 Myr), summed over the clean patches of each snapshot, in
km/s. This is the pressure the layer holds per unit of star formation. The dashed line is the feedback yield of Ostriker &
Kim (2022) evaluated at the same pressure: the pressure per unit star formation that supernovae, winds and radiation
supply in their simulations. Where the two agree, star formation supplies the pressure. Before coalescence the measured
curve sits on the yield to within a factor of 2. After coalescence it sits 10 to 100 times above it in the quiet intervals,
reaching 100 at 160 to 185 Myr: the layer holds far more pressure than its star formation can supply, so that pressure comes
from elsewhere (the merger flow and the field, next figures), and it drops back to the yield only during the two nuclear
bursts.

![pressure per unit star formation](figs/story_2_feedback.png)

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
## B, functional form of the model for the cluster mass function

Inputs per 0.5 kpc patch: layer weight $\mathcal{W}$, midplane density $\rho_{\rm mid}$, scale height $H$, effective
dispersion $\sigma_{\rm eff} = (P/\rho_{\rm mid})^{1/2}$ with $P = \mathcal{W}$. Inputs per cold complex: mass $M$,
half-mass radius $r_h$, Alfvén Mach number $\mathcal{M}_A$.

**(1) Complex population.** Measured, invariant through the merger:
$$n(M)\,{\rm d}M \propto M^{-1.6}\,{\rm d}M \quad (M > 300\ {\rm M}_\odot),$$
with the top rising from $10^5$ to $10^7$ Msun through the merger.

**(2) Virial parameter of a complex from the layer.**
$$\alpha_{\rm tot}(M, r_h, \sigma_{\rm eff}) = \frac{5}{3}\,\frac{\sigma_{\rm eff}^2\,(r_h/H)\,r_h}{G\,M}\,\left(1 + \mathcal{M}_A^{-2}\right).$$
Measured over predicted: median 0.9, 0.4 dex scatter, rank correlation 0.8.

**(3) Efficiency of a complex.** A step in $\alpha_{\rm tot}$ with lognormal scatter:
$$\epsilon(\alpha_{\rm tot}) = \epsilon_u + \frac{\epsilon_b - \epsilon_u}{1 + (\alpha_{\rm tot}/\alpha_c)^m},\qquad
\log\epsilon \sim \mathcal{N}\big(\log\epsilon(\alpha_{\rm tot}),\ s\big).$$

| parameter | dense clumps (n > 100) | complexes (n > 10) |
|---|---|---|
| $\epsilon_b$ (bound) | 0.050 | 0.022 |
| $\epsilon_u$ (unbound) | 0.0012 | 0.0017 |
| $\alpha_c$ | 4.0 | 5.9 |
| $m$ | 4.7 | 3.1 |
| $s$ (dex) | 0.69 | 0.71 |

**(4) Stars to clusters.** The dominant cluster of a complex takes $f_c = 0.5$ of its stars; secondary clusters steepen the
mass function by about 0.2.

**(5) Cluster mass function.** The complex mass function convolved with the efficiency distribution:
$$n_{\rm cl}(M_{\rm cl}) = \int {\rm d}M\ n(M)\ \frac{1}{f_c M}\ p_\epsilon\!\left(\frac{M_{\rm cl}}{f_c M}\ \Big|\ \alpha_{\rm tot}(M,\sigma_{\rm eff})\right),$$
with $p_\epsilon$ the lognormal of (3) centred on the step. Evaluated on the measured complexes of each phase (200 draws):

| phase | $\sigma_{\rm eff}$ [km/s] | bound fraction | slope pred. | slope meas. | top pred. | top meas. |
|---|---|---|---|---|---|---|
| before 1st passage | 8.7 | 0.35 | 2.03 ± 0.09 | 1.84 | 1.0e4 | 7.4e3 |
| between passages | 9.4 | 0.11 | 1.94 ± 0.07 | 1.77 | 1.4e4 | 2.0e4 |
| coalescence | 18.3 | 0.05 | 1.83 ± 0.08 | 1.77 | 1.5e5 | 1.0e5 |
| remnant | 19.3 | 0.03 | 1.70 ± 0.06 | 1.55 | 2.8e5 | 1.5e5 |

The slope is 0.15 to 0.2 too steep in every phase; the flattening (0.33 vs 0.29) and the rise of the top (30 vs 20) are
reproduced. With the mass dependence of boundedness removed (bound flags shuffled) the slopes are 2.21, 2.26, 1.98, 1.83.

**(6) Largest cluster per interval.** Two equivalent statements: $M_{\max} \simeq 0.1\,M_{\star}(25\,{\rm Myr})$, the
largest burst of the interval; and, while star formation is on, $M_{\max} \simeq f_c\,\eta\,\Sigma_{\rm SFR}^{\rm PRFM}(\mathcal{W}_{90})\,A\,\tau$
with $\eta = 0.4$, $A = 0.25$ kpc$^2$, $\tau = 10$ Myr, slope 0.94 and 0.09 dex scatter over seven intervals.

**(7) Patch-level duty cycle** (used for the burst-mass mock of Part B above):
$P({\rm burst}\,|\,\rho_{\rm mid}) = [1 + \exp(-k(\log\rho_{\rm mid} - \log\rho_{50}))]^{-1}$, $k = 5.1$,
$\log\rho_{50} = -1.65$ (Msun pc$^{-3}$); burst mass $\log M_{\rm burst} \sim \mathcal{N}(a + b\log\mathcal{W},\ 0.4$–$0.5)$.

## B, the chain from the layer to the clusters: measured relations, with conjectured physics

This section stays inside Part B. The relations in it are measured. The physical reasoning attached to them, a turbulent
cascade from the patch to the clump, a virial criterion for the efficiency step, a mixing-length balance for the
flow-driven dispersion, is conjecture that fits the numbers and has not been established independently; it should be
read as a labelling of the correlations, not as an explanation of them. Everything below is written so it can be checked. Each step states the physics, the equation, what is assumed, and the
number it was tested against. "Measured" means taken from the simulation; "derived" means it follows from the previous
step with no constant fitted here.

### B-i. Pressure and effective dispersion of the layer

Vertical equilibrium gives the midplane pressure from the weight of the gas layer, $P = \mathcal{W}$ (Part A: holds to
5 % outside coalescence). The midplane density $\rho_{\rm mid}$ is measured. Define the effective dispersion

$$\sigma_{\rm eff}^2 \equiv P/\rho_{\rm mid}.$$

This is the total support per unit mass, thermal plus turbulent plus magnetic. Its value is not predicted by equilibrium:
equilibrium fixes $P$, and $\rho_{\rm mid}$ adjusts. What fixes $\sigma_{\rm eff}$ is the energy budget:

- feedback-regulated: the turbulent driving is set by the star formation rate through the yield, and the Ostriker & Kim
  calibration gives $\sigma_{\rm eff} \simeq 9$ km/s at these pressures (measured here: 8.7 and 9.5 km/s in the two disc
  phases);
- flow-driven: turbulent energy injected by the large-scale shear at rate $\Sigma\,\sigma\,H\,S^2$ (mixing length $H$,
  shear rate $S$) balances dissipation $\Sigma\,\sigma^3/H$, so $\sigma = \sqrt{2\epsilon_S}\,H S$. This is the shear
  closure of Part B with its one calibrated coefficient $\epsilon_S = 0.02$; measured $\sigma_{\rm eff}$ = 18 to 21 km/s
  after coalescence.

The magnetic field is measured, not modelled: it saturates near equipartition with the turbulence it is driven by.

### B-ii. Cascade from the layer to the clump

Assumption: below the driving scale $H$ the turbulence is a supersonic (Burgers) cascade,

$$\sigma(\ell) = \sigma_{\rm eff}\,(\ell/H)^{1/2}.$$

Test: for every dense clump, at the snapshot before it forms stars, the predicted $\sigma(r_h)$ against its measured 3D
velocity dispersion. Ratio measured/predicted: median 0.95, scatter 0.20 dex, per-phase medians 0.90, 0.94, 1.02, 1.03.
No constant. (A steeper exponent of 0.6 overshoots by 1.37.)

### B-iii. Virial parameter of the clump

With the clump's measured mass $M$ and half-mass radius $r_h$, and $\sigma_{\rm 1d}^2 = \sigma(r_h)^2/3$,

$$\alpha_{\rm vir} = \frac{5\,\sigma_{\rm 1d}^2\, r_h}{G M} = \frac{5}{3}\,\frac{\sigma_{\rm eff}^2\,(r_h/H)\, r_h}{G M}.$$

Magnetic support adds $E_B/E_{\rm kin} = v_A^2/\sigma^2$, so the total is

$$\alpha_{\rm tot} = \alpha_{\rm vir}\,(1 + \mathcal{M}_A^{-2}),\qquad \mathcal{M}_A = \sigma/v_A.$$

Test: predicted against measured $\alpha_{\rm vir}$ (the measured one uses the clump's own dispersion): median ratio 0.90,
scatter 0.40 dex, rank correlation 0.80, per-phase medians 0.81, 0.88, 1.03, 1.07. Derived. Note what this does and does
not use: the clump mass and radius are measured, so the mass-size relation of the clumps is an input; replacing it by a
uniform sphere at the threshold density is wrong by factors of 2 to 12 and phase dependent. The Alfvén Mach number is
measured per clump; it is 2.2 in the first disc phase and 1.0 to 1.2 afterwards.

### B-iv. Bound or not: the efficiency

The virial theorem says a clump with $2E_{\rm kin} + E_B < |E_{\rm grav}|$ collapses. In the units above that is
$\alpha_{\rm tot} \lesssim 2$ with $R$ the full radius, or $\approx 4$ with $R = r_h$. The efficiency measured against
$\alpha_{\rm tot}$ is a step at exactly that place:

$$\epsilon = \epsilon_b \ \ (\alpha_{\rm tot} < \alpha_c), \qquad \epsilon = \epsilon_u \ \ (\alpha_{\rm tot} \gg \alpha_c),$$

with the fit $\alpha_c = 4.0$, $\epsilon_b = 0.050$, $\epsilon_u = 0.0012$, a transition steepness of 4.7, and 0.69 dex
of clump-to-clump scatter around it. The location $\alpha_c$ is the virial criterion. The two levels are measured; the
bound one is consistent with about 1 % per free-fall time over the six free-fall times a clump lives, the standard number,
but is not derived here. With $\alpha_{\rm tot}$ taken not from the clump but predicted from the layer through B-ii and B-iii,
the same fit gives $\alpha_c = 3.5$, $\epsilon_b = 0.054$, $\epsilon_u = 0.0012$, and per-phase offsets of +0.26, −0.02,
−0.02, −0.10 dex. The first disc phase sits above the plateau, where the field is not yet saturated.

### B-v. Star formation of the layer

The star formation of a phase is the bound fraction times the bound efficiency, plus the floor:

$$\frac{M_\star}{M_{\rm dense}} = f_b\,\epsilon_b + (1 - f_b)\,\epsilon_u,\qquad f_b = f(\alpha_{\rm tot} < \alpha_c \mid \sigma_{\rm eff}).$$

Test per 25 Myr interval over all dense clumps: $f_b$ = 0.85, 0.52, 0.48, 0.37, 0.14, 0.25, 0.06, 0.015; measured
conversion over $f_b$ = 0.11, 0.17, 0.13, 0.08 in the disc phases, and the quiescent interval is predicted exactly
(0.0026 vs 0.0026). The two pericentre bursts are under-predicted by 3 and 70: there the compression is faster than the
clump's pre-onset state implies. The chain is quasi-static and does not contain the orbit.

Stated as a measurement about PRFM: the pressure is at the weight throughout, but the star formation rate at a given
weight is $f_b(\sigma_{\rm eff})\,\epsilon_b$ times the dense gas, and $f_b$ falls from 0.9 to 0.03 when $\sigma_{\rm eff}$
goes from 9 to 20 km/s at the same or larger weight.

### B-vi. From clumps to the cluster mass function

Measured inputs: a burst is one cold complex with one to three dense clumps (Part A), and half of a burst's stars end
in one cluster (Part A). So the cluster mass function is the distribution of $0.5\,\epsilon\,M$ over complexes, with
$\epsilon$ drawn from the bound or unbound distribution according to B-iv.

Construction: every complex of a phase, bound if $\alpha_{\rm tot} < \alpha_c$, gets an efficiency drawn from the measured
lifetime-efficiency distribution of bound (median 0.022, 0.7 dex) or unbound (0.0017, 0.85 dex) complexes; cluster mass
$0.5\,\epsilon\,M$; slope above 300 Msun and largest cluster, 200 draws.

| phase | predicted slope | measured slope | predicted top | measured top |
|---|---|---|---|---|
| before 1st passage | 2.03 ± 0.09 | 1.84 | 1.0e4 | 7.4e3 |
| between passages | 1.94 ± 0.07 | 1.77 | 1.4e4 | 2.0e4 |
| coalescence | 1.83 ± 0.08 | 1.77 | 1.5e5 | 1.0e5 |
| remnant | 1.70 ± 0.06 | 1.55 | 2.8e5 | 1.5e5 |

The flattening through the merger (0.33 predicted, 0.29 measured) and the rise of the top (factor 30 predicted, 20
measured) are reproduced. The absolute slope is 0.15 to 0.2 too steep in every phase; this is the tail of the lognormal
efficiency distribution above the 300 Msun cut, and it is a genuine shortfall of the construction. Two things make the
flattening happen in the model: in the merger phases the most massive complexes are more often bound than the small
ones (15 % above 30 000 Msun against 3 % below, in the remnant), and the mixing of bound and unbound efficiencies. Shuffling
the bound flag among the complexes of a phase, i.e. removing the mass dependence, gives 2.21, 2.26, 1.98, 1.83 and loses
most of the flattening. That mass dependence of boundedness is measured ($\alpha_{\rm tot} \propto M^{-0.26}$ in the
remnant), not derived: a uniform-density scaling would give $M^{-2/3}$ and is wrong.

### What is derived, what is measured, what is fitted

- Derived: B-ii (cascade, one assumed exponent, verified), B-iii (virial parameter), the location of the threshold in B-iv, B-v.
- Measured inputs: $\rho_{\rm mid}$, $H$, the clump masses and radii and their mass-size relation, the Alfvén Mach number,
  the two efficiency levels and their scatter, the capture fraction 0.5, the mass dependence of boundedness.
- Fitted on this run: $\epsilon_S = 0.02$ in C1 (flow-driven case only).
- Not covered: the timing of the pericentre bursts, the absolute slope offset of 0.15 to 0.2 in B-vi, the clump-to-clump
  scatter.

![clump virial parameter](figs/partC_alpha.png)

![clump efficiency](figs/partC_eff.png)

---
# Part C. The closed calculation (not done)

Part C is reserved for a calculation that takes as input only the state of the layer, its weight and what drives its
turbulence, and returns the cluster mass function, with the simulation used solely to test the result. Nothing in Part B
qualifies: every chain there reads the clump population, the efficiency levels, their scatter and the capture fraction from
the run. For C to exist it has to contain, with no quantity taken from the simulation:

1. the effective dispersion of the layer from the energy budget of the driving, feedback yield or shear injection, without
   a coefficient fitted here;
2. the population of dense clumps, mass function and mass-size relation, from fragmentation of the turbulent layer at that
   dispersion; this is also what decides whether massive clumps stay bound when the dispersion rises;
3. the fraction of each clump that turns into stars from a collapse calculation, including why bound clumps convert a few
   percent and unbound ones a tenth of that;
4. how the stars of a collapsing complex divide into clusters, in place of the measured capture fraction.

The measurements of Part A and the chain of Part B say what such a calculation must reproduce: an invariant clump mass
function of slope 1.6, a clump velocity dispersion inherited from the layer as a square-root cascade, a boundedness
threshold at the virial criterion, a bound efficiency near 5 %, a hundredfold suppression of star formation when the layer
dispersion doubles at fixed weight, and a cluster mass function that flattens from 1.85 to 1.55 with its top rising
twentyfold.
