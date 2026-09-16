# merger_prfm

Does pressure-regulated, feedback-modulated star formation (PRFM; Ostriker & Kim 2022) survive a galaxy merger, and what
does it imply for the masses of the star clusters that form? One simulation of two gas-rich dwarf galaxies merging
(4 Msun gas and star particles, 0.4 pc softening, magnetic fields, resolved supernovae), 226 Myr through three pericentre
passages to the final coalescence.

Three parts, kept apart. **A**: what the simulation shows, figures with what is plotted and the numbers. **B**: a
semi-empirical model, one set of equations with constants calibrated on the run, and the tests it passes and fails.
**C**: the closed calculation from the layer state to the cluster mass function, which does not exist yet.

---
# Part A. What the simulation shows

## Merger phases

Separation of the two nuclei with time. Pericentres at 40, 110, 169 and 217 Myr (dotted lines; the last, the final
coalescence, solid). The four phases used in every other figure are the intervals between them: before the first passage
(25 to 40 Myr), first to second passage (40 to 110), second to third passage (110 to 169), and third passage to coalescence
(169 to 226). The first 25 Myr are excluded: the initial conditions are isothermal, and until about 25 Myr the layer is held
up by that initial thermal pressure rather than by feedback (two thirds thermal, turbulent dispersion below 3 km/s, no
field); by 25 Myr feedback-driven turbulence has taken over. The run ends 9 Myr after the nuclei merge, so there is no
remnant phase. Only the 25 to 40 Myr interval is an undisturbed disc; from the first passage on the discs are tidally
disturbed even while they remain separate.

![separation](figs/story_0_sep.png)

## Vertical equilibrium

The first half of PRFM: does the midplane pressure balance the weight? The ratio of total midplane pressure (thermal plus
turbulent plus magnetic stress) to weight is shown for two definitions of the weight. Using the full 3 kpc column, the
ratio drops to 0.3 after the second passage and the model looks broken. Restricting the weight to the gas within 0.5 kpc of
the midplane, where the star formation happens, the ratio stays near one before the second passage and after the third, and
falls 20 to 35 % short between the second and third. The apparent failure is the dark-matter weight of tidal debris streaming
through the disc at 25 to 30 km/s, which is not a layer in any orientation and does not form stars.

![equilibrium](figs/story_1_equilibrium.png)

## Pressure per unit star formation

What is plotted: the measured total midplane pressure of the gas layer divided by the measured star formation rate
surface density of the same patches (stars younger than 10 Myr), summed over the clean patches of each snapshot, in
km/s. This is the pressure the layer holds per unit of star formation. The dashed line is the feedback yield of Ostriker &
Kim (2022) evaluated at the same pressure: the pressure per unit star formation that supernovae, winds and radiation
supply in their simulations. Where the two agree, star formation supplies the pressure. Before the second passage the measured
curve sits on the yield to within a factor of 2. After it, it sits 10 to 100 times above it in the quiet intervals,
reaching 100 at 160 to 185 Myr: the layer holds far more pressure than its star formation can supply, so that pressure comes
from elsewhere, and it drops back to the yield only during the two nuclear bursts. "Elsewhere" is defined in the next two
figures: the ordered large-scale motions of the merger (tidal streaming, infall, rotation, measured as the velocity
gradient across each 0.5 kpc patch) doing work on the turbulence, and the magnetic field.

![pressure per unit star formation](figs/story_2_feedback.png)

## Magnetic field

The field is seeded at 10 nG and amplified by the turbulence: it e-folds every 10 Myr, saturates near 1 microgauss by
50 Myr, jumps by a factor of five at the second passage and grows to 10 to 20 microgauss after the third, where the magnetic
pressure reaches equipartition with thermal plus turbulent pressure. It is treated as a measured contribution to the
vertical support, not modelled.

![dynamo](figs/story_4_dynamo.png)

## Clouds

Stars keep the identity of the gas particle they formed from, so every young star can be traced to the cold gas it came
from. Nearly all stars formed within 10 Myr of a snapshot come from gas already sitting in a cold cloud above 10 cm^-3.
Those clouds are 30 to 50 pc complexes of 1e4 to 1e5 Msun, at a mean density near 10 cm^-3, the low-density giant molecular
clouds of a dwarf galaxy. In 85 to 95 % of star-forming events a single complex supplies essentially all the stars, and
inside it the stars come from one to three dense clumps above 100 cm^-3 (a few thousand Msun, radii of 4 to 8 pc). A burst
is the collapse of one such complex; the largest bursts, in the merging nuclei, draw on many clumps.

The mass function of the complexes (T < 1000 K, n > 10 cm^-3, friends-of-friends at 3 pc) has a slope of 1.6 and does not
change through the merger. Only the most massive end grows, as the merging nuclei assemble cold complexes of a million solar
masses and more.

![cloud mass function](figs/clouds_mf.png)

What does change is how efficiently the clumps turn gas into stars. Following each dense clump through the snapshots by its
particle IDs, the figure shows the fraction of its peak mass that it ever turned into stars, against its virial parameter
just before star formation began. More strongly bound clumps convert more, by a factor of 40 over the observed range, and
the efficiency drops sixty-fold from before the first passage to after the third, where the clumps are more turbulent. But at the same
virial parameter the phases are offset by more than an order of magnitude, so the kinetic virial parameter is not the whole
story, and at fixed state the efficiency still scatters by 0.7 dex. That scatter is what produces the 0.4 dex spread of
burst masses and the slope of the cluster mass function. This is the one link in the chain that is measured rather than
explained.

![cloud efficiency](figs/clouds_eff.png)

## Cluster mass function

Bound star clusters younger than 10 Myr, per phase, with the maximum-likelihood power-law index above 300 Msun. The slope
flattens from 1.94 before the first passage to 1.52 after the third, and the most massive cluster grows from 7e3 to 1e5 Msun.
Both trends are robust to the minimum group size, the linking length and the bound flag. Below 300 Msun after the second
passage the counts drop because only a third of the small groups are bound, so the function is quoted from
300 Msun up. Since the cloud mass function does not change (above), the flattening and the growing top come from the
star formation efficiency of the clouds.

![MF per phase](figs/story_5_mf.png)

## Largest cluster per burst

For every 0.5 kpc patch that formed more than 500 Msun of stars in the last 10 Myr and hosts a bound cluster, the mass of
its most massive cluster against the total mass of young stars in the patch. The median runs along half the young mass up
to about 1e5 Msun: each burst makes one dominant cluster that captures about half of the stars formed, with a handful of
smaller companions. This is the link between the star formation rate of a patch and the top of the cluster mass function,
and it is why the weight of the gas layer, which sets the size of the largest burst, also sets the mass of the largest
cluster. The most massive bursts, in the merging nuclei, fall below the line because their stars are spread over several
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
tilted by 30, 60 and 90 degrees from it. Before the second passage only the disc normal gives a ratio near one, so the disc is a
layer and the chosen normal is the right one. After it every direction gives 0.2 to 0.5: the deficit does not
depend on viewing angle, so the streaming debris is not a layer in any orientation.

![line of sight](figs/prfm_time_los.png)

## Low-mass end of the cluster mass function

Local power-law slope of the cluster mass function in three mass bins, per phase, for all friends-of-friends groups
(filled) and for bound groups only (open). Above 300 Msun the two agree and the slopes are stable against the minimum group
size and linking length. In the 100 to 300 Msun bin after the second passage the bound-only slope collapses, because only a third
of those small groups are bound: that end of the function is set by the bound flag, not by the physics, and is not used.

![low-mass convergence](figs/lowmass_convergence.png)


---

---
# Part B. The model

## Functional form


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
| before 1st passage | 8.7 | 0.42 | 2.06 ± 0.11 | 1.94 | 7.4e3 | 6.7e3 |
| 1st to 2nd passage | 9.4 | 0.11 | 1.93 ± 0.06 | 1.77 | 1.9e4 | 2.0e4 |
| 2nd to 3rd passage | 18.3 | 0.05 | 1.85 ± 0.07 | 1.74 | 1.3e5 | 2.8e4 |
| 3rd passage to coalescence | 19.3 | 0.02 | 1.61 ± 0.07 | 1.52 | 2.8e5 | 1.0e5 |

The slope is 0.1 to 0.2 too steep in every phase; the flattening (0.45 vs 0.42) and the rise of the top (38 vs 15) are
reproduced. With the mass dependence of boundedness removed (bound flags shuffled) the slopes are 2.21, 2.26, 1.98, 1.83.

**(6) Largest cluster per interval.** Two equivalent statements: $M_{\max} \simeq 0.1\,M_{\star}(25\,{\rm Myr})$, the
largest burst of the interval; and, while star formation is on, $M_{\max} \simeq f_c\,\eta\,\Sigma_{\rm SFR}^{\rm PRFM}(\mathcal{W}_{90})\,A\,\tau$
with $\eta = 0.4$, $A = 0.25$ kpc$^2$, $\tau = 10$ Myr, slope 0.94 and 0.09 dex scatter over seven intervals.

**(7) Patch-level duty cycle** (used for the burst-mass mock of Part B above):
$P({\rm burst}\,|\,\rho_{\rm mid}) = [1 + \exp(-k(\log\rho_{\rm mid} - \log\rho_{50}))]^{-1}$, $k = 5.1$,
$\log\rho_{50} = -1.65$ (Msun pc$^{-3}$); burst mass $\log M_{\rm burst} \sim \mathcal{N}(a + b\log\mathcal{W},\ 0.4$–$0.5)$.

## Shear term

What "supplied by the flow" means. Each 0.5 kpc patch has, besides its turbulence, an ordered velocity field: the
differential motion across the patch from the merger's tidal streaming, infall and rotation. Its size is the shear rate S,
the norm of the velocity-gradient tensor of a linear fit to the gas velocities in the patch, in km/s per kpc. Ordered
motion with a gradient does work on the gas at the rate of a viscous stress, Sigma sigma H S^2, and that work feeds the
turbulence. The plot tests this on the patches that have had no star formation in the last 10 Myr, so feedback cannot be
the source: their turbulent pressure, after removing the linear bulk flow, is plotted against Sigma H S^2. One coefficient,
0.02, describes all four phases with no offset. Feedback-based yields miss these patches by up to two orders of magnitude.
This is what the previous figure's excess pressure is: the turbulence the merger's own velocity field keeps stirring.

![shear closure](figs/story_3_shear.png)

## Largest cluster and weight

The most massive bound cluster formed in each 25 Myr interval against $\mathcal{W}$, the weight of the gas layer (within 0.5 kpc of
the midplane) of the heaviest star-forming patches in that interval (90th percentile over patches). Labels give the start of each
interval in Myr. The line is not a fit: it is half of the PRFM star formation rate at that weight, times the patch area and a 10 Myr
burst, i.e. the size of the largest burst a patch of that weight can produce, halved by the capture fraction of the previous figure.
Seven active intervals spanning a factor of 30 in weight follow it (slope 0.94 against the median weight, 0.09 dex scatter). The one
quiescent interval, 160 to 185 Myr, has the highest weight of all and sits 100 times below. The weight sets the cluster mass while
star formation is on; whether it is on is a separate question, answered below.

![ceiling](figs/story_7_ceiling.png)

## Duty cycle

Probability that a 0.5 kpc patch formed more than 500 Msun of stars in the last 10 Myr, as a function of its weight $\mathcal{W}$
(gas layer within 0.5 kpc of the midplane), per phase. Points are measured fractions, lines logistic fits. The weight at
which half the patches are active rises twenty-fold through the merger. The same threshold expressed in midplane density
moves by less than a factor of three, so the physical threshold is in density: after the second passage the higher turbulent and
magnetic support means the same weight produces a lower midplane density.

![duty cycle](figs/duty_cycle.png)

## Burst-mass distribution

Distribution of the 10 Myr burst mass of the patches between the two passages, against a model with no free shape: each
patch bursts with the probability from the duty cycle above, and its burst mass is drawn from a lognormal of 0.5 dex width
around a mean that rises with the weight. The model reproduces the distribution and its curvature (the same holds in the
other phases). The slope of the burst-mass distribution, and through it of the cluster mass function, is the width of this
lognormal, which the cloud figures above trace to the scatter in cloud star formation efficiency.

![burst kernel](figs/burst_kernel.png)

## Clump virial parameter and efficiency
The virial parameter of each dense clump at the snapshot before it forms stars, against the effective dispersion of its
patch, sigma_eff = (P/rho_mid)^{1/2}. Points are clumps, lines medians per phase, the dashed line alpha proportional to
sigma_eff squared. All four phases follow one relation with 0.2 dex phase offsets. This is equation (2) of the functional
form with the measured clump mass and radius.

![clump virial parameter](figs/partC_alpha.png)

The lifetime-integrated efficiency of each clump against its total virial parameter, kinetic plus magnetic, before onset,
with the fitted step of equation (3): 5 % below alpha_tot = 4, 0.1 % above 16, 0.7 dex scatter. Colours are phases; the
black line is the median over all phases.

![clump efficiency](figs/partC_eff.png)

## What the model gets and does not get
- Vertical balance, the pressure per unit star formation before the second passage, the shear term for the quiet patches: yes.
- Which patches form stars: to a factor of 2 from the density threshold.
- The amount of star formation per interval from the bound fraction: to a factor of 2 before the second passage and exactly in the
  quiescent interval; wrong by 3 and 70 in the two pericentre bursts, which are compression events the model does not
  contain.
- The cluster mass function: the flattening (0.45 predicted, 0.42 measured) and the rise of the top (38 vs 15); the
  absolute slope is 0.15 to 0.2 too steep in every phase.
- The largest cluster per interval: linear in the weight of the active patches over seven intervals, 0.09 dex; the
  quiescent interval 100 times below.
- Not reproduced or not modelled: the timing of the bursts, the two efficiency levels, the 0.7 dex clump-to-clump scatter,
  the mass dependence of boundedness that drives the flattening (measured, used as input).



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
dispersion doubles at fixed weight, and a cluster mass function that flattens from 1.94 to 1.52 with its top rising
fifteenfold.

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
dispersion doubles at fixed weight, and a cluster mass function that flattens from 1.94 to 1.52 with its top rising
fifteenfold.
