## Status of the PRFM comparison

Reference: the Ostriker & Kim (2022) calibration. The run's own undisturbed disc lasts only from 25 to 40 Myr, since the
first 25 Myr are the isothermal initial state relaxing and the discs are tidally disturbed from the first passage on.

| | before 1st passage | 1st to 2nd | 2nd to 3rd | 3rd to coalescence |
|---|---|---|---|---|
| P over W, layer 1 kpc deep | 1.05 | 0.94 | 0.73 | 0.96 |
| pressure per unit SFR over the yield, 10 Myr SFR | 1.9 | 1.7 | 37 | 2.5 |
| same, 40 Myr SFR as in OK22 | 0.6 | 1.6 | 2.5 | 9.4 |
| SFR at fixed weight vs OK22, 0.5 kpc columns (dex) | -0.5 | -0.6 | -1.8 | -2.2 |
| same, 0.25 kpc columns (dex) | -0.2 | -0.5 | -1.6 | -2.1 |

Durations. Pressure over weight leaves 0.7 to 1.4 for 24 % of the run (97 to 141 Myr and short episodes at the later
pericentres) and never leaves 0.5 to 2. The pressure per unit star formation exceeds 3 times the yield for 34 % of the run
with the 40 Myr rate (107 to 125 and 158 to 207 Myr) and 10 times for 20 %; its time-median over the whole run is 1.9.

Reading. Vertical balance survives the merger, with one interval of moderate shortfall around the second passage during
which the gas is not a layer in any orientation. The supply of that pressure by star formation does not survive: after the
second passage the layer holds 10 to 100 times more pressure per unit star formation than feedback yields, and that pressure
is turbulent and magnetic, sustained by the shear of the merger flow at a fixed coefficient of 0.02. The star formation rate
at fixed weight falls below the OK22 relation by 0.5 dex in the discs and 2 dex after the second passage, at every column
size. The balance test itself is meaningful only for columns at least 1 kpc deep and within about 30 degrees of the disc
normal; the yield and rate comparisons depend on neither.

---
---
# Part A. What the simulation shows

## Definitions used by every figure

*Columns.* For each galaxy a 12 by 12 grid of columns of 0.5 kpc side is laid in the plane perpendicular to the gas
angular-momentum vector measured within 2 kpc of the galaxy centre; the centre is the median position of that galaxy's
initial stellar disc particles (membership by particle ID). One grid, centred between the two centres with the angular
momentum within 3 kpc, replaces the two whenever the centres are closer than 0.5 kpc. A *full column* spans 1.5 kpc above
and below the grid plane; a *layer column* spans 0.5 kpc above and below the midplane of that column, the mass-weighted
mean height of its gas in the full-column run. Only columns with gas surface density above 1 Msun/pc² are used, and before
the second passage only columns in which less than 10 % of the gas belongs to the other galaxy by particle ID.

*Gas quantities per column.* Sigma_gas = gas mass in the column / area. H = mass-weighted rms height of the gas about the
midplane. The *slab* is the layer of half-thickness max(25 pc, 0.25 H) about the midplane. "2p" means gas below 20 000 K.
P_th = sum over slab particles of m (gamma-1) u divided by the slab volume; P_turb = sum of m (v_n - v̄_n)² over the slab
divided by its volume, v_n the velocity component along the column normal relative to the galaxy's mean velocity and v̄_n
the mass-weighted mean over the slab (2p gas); P_mag = volume-weighted B²/8pi over the slab. P_tot = P_th + P_turb + P_mag.
All pressures are P/k_B in K cm⁻³. rho_mid = 2p gas mass in the slab / slab volume.

*Weight.* W = one half of [ sum over 2p gas above the midplane of m (-g_n) + sum over 2p gas below of m g_n ] / area, with
g_n the gravitational acceleration along the normal at each particle from a particle-mesh solve of all gas, stars and dark
matter (512³ cells over ±4 kpc nested in 512³ over ±32 kpc; median error 8 % against direct summation). In K cm⁻³.

*Star formation rate.* Sigma_SFR,10 = mass of stars younger than 10 Myr inside the column / area / 10 Myr, using only
stars formed during the run (the initial stellar particles are excluded); Sigma_SFR,40 the same over 40 Myr. Both in
Msun yr⁻¹ kpc⁻².

*Ostriker & Kim (2022) relations used.* Yield Upsilon_tot(P) = 10^(-0.212 log P + 3.86) km/s (their eq. 26c);
Sigma_SFR(W) = 10^(1.17 log W - 7.32) (eq. 28b); Sigma_SFR(P) = 10^(1.18 log P - 7.43) (eq. 28a). A yield in km/s is
converted to a pressure per unit star formation rate with 1 km/s = 4.81 × 10³ K cm⁻³ per Msun yr⁻¹ kpc⁻².

*Time-series statistic.* Where a figure shows one value per snapshot from many columns, it is the Sigma_gas-weighted median
over the columns: the value at which half of the total gas surface density lies above and half below.

*Phases.* 25 to 40 Myr (before the first passage), 40 to 110, 110 to 169, 169 to 226 Myr, split at the pericentres 40,
110, 169 Myr; the final coalescence is at 217 Myr. The first 25 Myr are excluded (isothermal initial state relaxing).

## Merger phases
x: time since the start of the run. y: distance between the two galaxy centres defined above. Dotted lines mark the
pericentres at 40, 110 and 169 Myr; the solid line the final coalescence at 217 Myr.

![separation](figs/story_0_sep.png)

## Vertical equilibrium
x: time. y: P_tot / W per snapshot, the Sigma_gas-weighted median over the columns, for the full column (light) and the
layer column (dark), both with P_tot and W from the 2p gas as defined above. Result: with the layer weight the ratio is
1.05, 0.94, 0.73, 0.96 in the four phases and never leaves 0.5 to 2; with the full column it drops to 0.4 between the
second and third passage because the column then includes the dark-matter weight of gas streaming outside the layer.

![equilibrium](figs/story_1_equilibrium.png)

## Pressure per unit star formation
x: time. y: the sum of the layer P_tot over the columns of a snapshot divided by the sum of Sigma_SFR,10 over the same
columns, converted to km/s with the factor above. Dashed: Upsilon_tot(P) evaluated at the Sigma_gas,2p-weighted mean layer
pressure of the snapshot. Result: within a factor of 2 of the yield before the second passage, 10 to 100 times above it
in the quiet intervals after, back to the yield during the two nuclear bursts.

![pressure per unit star formation](figs/story_2_feedback.png)

## Magnetic field
x: time. y: B = sqrt(8 pi P_mag) from the slab magnetic pressure of each column, in microgauss; dark line the
Sigma_gas-weighted mean over the columns, light line the median column. Result: 10 nG seed, e-folding 10 Myr, 1 to 2
microgauss from 50 Myr, a jump by 5 at the second passage, 10 to 20 microgauss after the third.

![dynamo](figs/story_4_dynamo.png)

## Clouds
*Clouds* are friends-of-friends groups, linking length 3 pc, at least 25 particles (100 Msun), of gas colder than 1000 K
and denser than 10 cm⁻³, identified on every tenth snapshot over the whole box. *Clumps* are the same with a 100 cm⁻³
threshold and 1.5 pc linking. Stars carry the ID of the gas particle they formed from, so the stars formed from any cloud
are counted exactly. Per cloud: mass, half-mass radius r_h, mass-weighted 3D velocity dispersion sigma_3d about the mean,
mean and maximum density, virial parameter alpha_vir = 5 (sigma_3d²/3) r_h / (G M), rms Alfvén speed v_A of its members.

Cloud mass function. x: cloud mass. y: number of clouds per unit mass divided by the total number of clouds in the phase;
Poisson errors. Dashed: M^-1.6 for reference. Result: the same in all four phases (index 1.57 to 1.60 above 300 Msun);
only the maximum mass grows, from 3 × 10⁵ to 10⁷ Msun.

![cloud mass function](figs/clouds_mf.png)

Efficiency. Clumps are followed from snapshot to snapshot by particle-ID overlap (a clump is linked to its successor if
each is the other's largest overlap). x: alpha_vir of the clump at the last snapshot before its first star forms. y:
epsilon_int = all stars ever formed from the lineage's members (summed over the chain plus 10 Myr after its last
snapshot) divided by the largest mass the lineage reached. Points: lineages with peak mass above 300 Msun that start after
5 Myr and end before 221 Myr; lines: medians per alpha bin per phase, phase assigned by the onset time. Result: 5 % for
bound clumps, falling to 0.1 % above alpha of about 10; the fraction of clumps below the threshold drops from 0.9 before
the first passage to 0.06 after the third.

![cloud efficiency](figs/clouds_eff.png)

## Cluster mass function
*Clusters* are friends-of-friends groups, linking length 5 pc, at least 25 particles, of stars younger than 10 Myr,
kept if energy-bound (kinetic plus Plummer-softened pairwise potential energy below zero, softening 1 pc) and not the
nucleus (fewer than 40 000 particles and r_h below 20 pc), in the clean columns, on every tenth snapshot so that no
cluster is counted twice. x: cluster mass. y: number per unit log mass per phase, Poisson errors. Lines: maximum-likelihood
power law above 300 Msun, alpha = 1 + N / sum ln(M/300). Result: alpha = 1.92, 1.77, 1.74, 1.52; the most massive cluster
6.7 × 10³, 2.0 × 10⁴, 2.8 × 10⁴, 1.0 × 10⁵ Msun. Below 300 Msun after the second passage only a third of the small groups
are bound, so the function is quoted from 300 Msun up.

![MF per phase](figs/story_5_mf.png)

## Largest cluster per burst
x: M_young, the mass of stars younger than 10 Myr in a full column. y: M_max, the mass of the most massive bound
cluster whose centre lies in that column. Clean columns with M_young above 500 Msun and at least one bound cluster, every
tenth snapshot. Line: median M_max per M_young bin; red: M_max = 0.5 M_young; dotted: equality. Result: the median follows
half the young mass up to 10⁵ Msun; the nuclear bursts above that split their stars over several clusters.

![reservoir](figs/story_6_reservoir.png)

## Every cluster against its patch
Each bound cluster of the clean sample against the state of the column it sits in (the column of the frame in which its
intruder fraction is lowest). Left x: layer W of that column. Middle x: layer P_tot. Right x: Sigma_SFR,10. y: cluster
mass. Dashed: 0.5 Sigma_SFR(W) A tau, 0.5 Sigma_SFR(P) A tau and 0.5 Sigma_SFR A tau with A = 0.25 kpc² and tau = 10 Myr,
the mass of half the stars a column forms in 10 Myr at the OK22 rate for its weight, for its pressure, and at its own rate.
Result: cluster masses fill two decades below the lines; 9 % lie above the weight line.

![cluster environment](figs/cluster_env.png)

## Line of sight
As the vertical-equilibrium figure, but the column normal is rotated: the disc normal itself, and the normal tilted by
30, 60 and 90 degrees about the first in-plane axis toward the second. Each orientation has its own full-column run and
its own layer cut about its own midplane. Every fourth snapshot; the two galaxy grids are merged with column-number
weights. y: Sigma_gas-weighted median of layer P_tot / W. Result: before the second passage only the disc normal gives one,
30 degrees gives 0.4 to 0.7; between the second and third passage all orientations give 0.6 to 1.2; after the third all
give 0.9 to 1.3.

![line of sight](figs/los_layer.png)

## Column size
x: column depth, the full height of the column: 3 kpc (full column) or 1, 0.5, 0.25 kpc (0.5, 0.25, 0.125 kpc above and
below the midplane of the 0.5 kpc reference column the point falls in). Marker: footprint, columns of 0.5, 0.25 and 0.125
kpc side on 12 by 12, 24 by 24 and 48 by 48 grids. y: for each configuration and phase, the median over the phase's
snapshots (every fourth snapshot) of the per-snapshot Sigma_gas-weighted median of P_tot / W, both from the 2p gas. Shaded:
twice the scale height. Result: the ratio depends on depth only; footprint 0.5 to 0.125 kpc at fixed depth changes nothing,
while depths below 0.5 kpc give 1.5 and 2.7.

![column depth](figs/scale_depth.png)

## Low-mass end of the cluster mass function
x: cluster mass bin. y: maximum-likelihood index of a power law restricted to that bin (likelihood maximised on a grid of
alpha; error from its curvature), for all friends-of-friends groups of the 10-particle catalogue (filled) and for the bound
ones only (open), per phase, every tenth snapshot. Result: above 300 Msun bound and all agree within the errors; in the
100 to 300 Msun bin after the second passage the bound-only index collapses because only a third of those groups are bound.

![low-mass convergence](figs/lowmass_convergence.png)

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
| before 1st passage | 8.7 | 0.37 | 2.02 ± 0.13 | 1.92 | 6.9e3 | 6.7e3 |
| 1st to 2nd passage | 9.4 | 0.11 | 1.93 ± 0.06 | 1.77 | 1.9e4 | 2.0e4 |
| 2nd to 3rd passage | 18.3 | 0.05 | 1.85 ± 0.07 | 1.74 | 1.3e5 | 2.8e4 |
| 3rd passage to coalescence | 19.3 | 0.02 | 1.61 ± 0.07 | 1.52 | 2.8e5 | 1.0e5 |

The slope is 0.1 to 0.2 too steep in every phase; the flattening (0.40 vs 0.40) and the rise of the top (45 vs 15) are
reproduced. With the mass dependence of boundedness removed (bound flags shuffled) the slopes are 2.21, 2.26, 1.98, 1.83.

**(6) Largest cluster per interval.** Two equivalent statements: $M_{\max} \simeq 0.1\,M_{\star}(25\,{\rm Myr})$, the
largest burst of the interval; and, while star formation is on, $M_{\max} \simeq f_c\,\eta\,\Sigma_{\rm SFR}^{\rm PRFM}(\mathcal{W}_{90})\,A\,\tau$
with $\eta = 0.4$, $A = 0.25$ kpc$^2$, $\tau = 10$ Myr, slope 0.94 and 0.09 dex scatter over seven intervals.

**(7) Patch-level duty cycle** (used for the burst-mass mock of Part B above):
$P({\rm burst}\,|\,\rho_{\rm mid}) = [1 + \exp(-k(\log\rho_{\rm mid} - \log\rho_{50}))]^{-1}$, $k = 5.1$,
$\log\rho_{50} = -1.65$ (Msun pc$^{-3}$); burst mass $\log M_{\rm burst} \sim \mathcal{N}(a + b\log\mathcal{W},\ 0.4$–$0.5)$.

## Shear term
Columns of the 18 snapshots with the turbulence decomposition, with Sigma_SFR,10 = 0. x: Sigma H S² in K cm⁻³, with
Sigma the column's gas surface density, H its rms height, and S the shear rate: the norm of the traceless symmetric part
of the velocity-gradient tensor of a linear fit of the gas velocities in the column, in km/s per kpc. y: P_turb of the
slab multiplied by (sigma_res/sigma_tot)², i.e. with the part of the velocity variance carried by the fitted linear flow
removed. Line: P_turb = 0.02 Sigma H S². Result: one coefficient in every phase, no offset.

![shear closure](figs/story_3_shear.png)

## Largest cluster and weight
25 Myr intervals. x: the 90th percentile, over the star-forming columns of the interval (M_young above 500 Msun), of the
layer W. y: the most massive bound cluster formed in the interval, over all its snapshots. Labels: interval start in Myr.
Line: 0.5 Sigma_SFR(W) A tau with the OK22 relation, A = 0.25 kpc², tau = 10 Myr; nothing fitted. Open circle: the interval
in which fewer than 5 % of the columns burst. Result: seven active intervals follow the line, slope 0.94 against the median
weight with 0.09 dex scatter; the quiescent interval lies 100 times below.

![ceiling](figs/story_7_ceiling.png)

## Duty cycle
x: layer W of a clean column. y: fraction of columns with M_young above 500 Msun, in bins of 0.5 dex in W (points), and
the maximum-likelihood logistic P = [1 + exp(-k (log W - log W_50))]⁻¹ per phase (lines). Result: log W_50 = 3.75, 3.86,
4.73, 4.95; the same fit against the midplane density moves by less than 0.5 dex between phases.

![duty cycle](figs/duty_cycle.png)

## Burst-mass distribution
Phase 40 to 110 Myr. Black: number per unit mass of the bursting columns (M_young above 500 Msun), all snapshots,
Poisson errors. Blue: 25 realisations of the null model: every clean column of the phase bursts with the probability of the
logistic above; a bursting column gets M = 10^(a + b log W + s N(0,1)) with b from a least-squares fit of log M_young on
log W over the phase's bursts and a, s from a lognormal likelihood truncated at 500 Msun; masses below 500 are dropped and
the rest binned as the data. All ingredients are fitted to the same bursts, so the comparison tests only whether the
histogram contains structure beyond them. Result: it does not.

![burst kernel](figs/burst_kernel.png)

## Clump virial parameter and efficiency
Left figure. x: sigma_eff = (P_tot / rho_mid)^(1/2) of the full column containing the clump at its pre-onset snapshot,
in km/s (P and rho from the 2p gas as defined above). y: alpha_vir of the clump at that snapshot. Points: lineages as in
the efficiency figure; lines: medians per sigma_eff bin per phase; dashed: alpha proportional to sigma_eff². Result: one
relation for all phases with 0.2 dex offsets.

![clump virial parameter](figs/partC_alpha.png)

Right figure. x: alpha_vir,tot = alpha_vir (1 + v_A²/sigma_3d²) of the clump at its pre-onset snapshot, v_A the rms
Alfvén speed of its members (B / sqrt(4 pi rho) per particle). y: epsilon_int as defined under Clouds. Black: median over
all phases per bin. Dashed: the two-level fit epsilon = epsilon_u + (epsilon_b - epsilon_u) / [1 + (alpha/alpha_c)^m] by
least absolute deviation in log: epsilon_b = 0.050, epsilon_u = 0.0012, alpha_c = 4.0, m = 4.7; 0.69 dex scatter about it.

![clump efficiency](figs/partC_eff.png)

## What the model gets and does not get
- Vertical balance, the pressure per unit star formation before the second passage, the shear term for the quiet patches: yes.
- Which patches form stars: to a factor of 2 from the density threshold.
- The amount of star formation per interval from the bound fraction: to a factor of 2 before the second passage and exactly in the
  quiescent interval; wrong by 3 and 70 in the two pericentre bursts, which are compression events the model does not
  contain.
- The cluster mass function: the flattening (0.40 predicted, 0.40 measured) and the rise of the top (45 vs 15); the
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
dispersion doubles at fixed weight, and a cluster mass function that flattens from 1.92 to 1.52 with its top rising
fifteenfold.

---
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
dispersion doubles at fixed weight, and a cluster mass function that flattens from 1.92 to 1.52 with its top rising
fifteenfold.
