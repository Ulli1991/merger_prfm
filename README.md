# merger_prfm

Does pressure-regulated, feedback-modulated star formation (PRFM; Ostriker & Kim 2022) survive a galaxy merger, and what
does it imply for the masses of the star clusters that form? One simulation of two gas-rich dwarf galaxies merging
(4 Msun gas and star particles, 0.4 pc softening, magnetic fields, resolved supernovae), 226 Myr through three pericentre
passages to the final coalescence.

Three parts, kept apart. **A**: what the simulation shows, figures with what is plotted and the numbers. **B**: a
semi-empirical model, one set of equations with constants calibrated on the run, and the tests it passes and fails.
**C**: the closed calculation from the layer state to the cluster mass function, which does not exist yet.

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

*Gas quantities per column.* Sums run over the gas particles i inside the column (footprint A = (0.5 kpc)², depth as
stated), with mass m_i, position z_i along the column normal n̂, velocity v_i, density rho_i, specific internal energy u_i,
temperature T_i, magnetic field B_i.

- Sigma_gas = (1/A) Σ_i m_i
- z_mid = Σ_i m_i z_i / Σ_i m_i ;  H = [ Σ_i m_i (z_i − z_mid)² / Σ_i m_i ]^(1/2)
- slab: |z_i − z_mid| < h_s with h_s = max(25 pc, 0.25 H);  V_slab = 2 h_s A
- "2p": particles with T_i < 2 × 10⁴ K
- v_n,i = (v_i − v_gal) · n̂ with v_gal the mean velocity of the galaxy's stellar particles within 1.5 kpc of its centre;
  v̄_n = Σ_{slab,2p} m_i v_n,i / Σ_{slab,2p} m_i
- P_th = (1/V_slab) Σ_{slab,2p} m_i (γ − 1) u_i ,  γ = 5/3
- P_turb = (1/V_slab) Σ_{slab,2p} m_i (v_n,i − v̄_n)²
- Π_mag = (1/V_slab) Σ_{slab,2p} (m_i / rho_i) (B_i² − 2 B_n,i²) / 8π, the vertical Maxwell stress (B_n the component along
  the normal; volume weighted; can be negative), which is the magnetic term that supports the layer against its weight
- P_mag = (1/V_slab) Σ_{slab} (m_i / rho_i) B_i² / 8π, the magnetic pressure (all slab gas), used only for the field strength
- P_tot = P_th + P_turb + Π_mag
- rho_mid = (1/V_slab) Σ_{slab,2p} m_i

Code units Msun (km/s)² kpc⁻³ are converted to K cm⁻³ (P/k_B) with the factor 4.90 × 10⁻⁶; Π_mag and P_mag are computed
in erg cm⁻³ and divided by k_B directly.

*Weight.* With g_n,i = g_i · n̂ the gravitational acceleration along the normal at particle i,

- W = (1/2A) [ Σ_{z_i > z_mid, 2p} m_i (−g_n,i) + Σ_{z_i < z_mid, 2p} m_i g_n,i ]

i.e. the mean of the weight of the gas above the midplane and of that below, each the integral of rho g_n through its half
of the column. g_i is from a particle-mesh solve of all gas, stars and dark matter on a 512³ grid over ±4 kpc nested in a
512³ grid over ±32 kpc (median error 8 % against direct summation on 64 test particles per snapshot). W in K cm⁻³ with
the same factor.

*Star formation rate.* Sigma_SFR,10 = (1/A) Σ_{stars in column, formed in the run, age < 10 Myr} m_* / (10 Myr), and
Sigma_SFR,40 the same with 40 Myr; in Msun yr⁻¹ kpc⁻². A column's M_young is the same sum without dividing by A and time.

*Ostriker & Kim (2022) relations used.* Yield Upsilon_tot(P) = 10^(-0.212 log P + 3.86) km/s (their eq. 26c);
Sigma_SFR(W) = 10^(1.17 log W - 7.32) (eq. 28b); Sigma_SFR(P) = 10^(1.18 log P - 7.43) (eq. 28a). A yield in km/s is
converted to a pressure per unit star formation rate with 1 km/s = 4.81 × 10³ K cm⁻³ per Msun yr⁻¹ kpc⁻².

*Time-series statistic.* Where a figure shows one value per snapshot from many columns, it is the Sigma_gas-weighted median
over the columns: the value at which half of the total gas surface density lies above and half below.

*Phases.* 25 to 40 Myr (before the first passage), 40 to 110, 110 to 169, 169 to 226 Myr, split at the pericentres 40,
110, 169 Myr; the final coalescence is at 217 Myr. The first 25 Myr are excluded (isothermal initial state relaxing).

## Merger phases
x: time since the start of the run. y: d(t) = |c_A − c_B|, with c_X the median position of the initial stellar disc
particles of galaxy X (galaxy A: particle IDs ≤ 26 000 000). Dotted lines: pericentres at 40, 110, 169 Myr; solid: final
coalescence at 217 Myr.

![separation](figs/story_0_sep.png)

## Vertical equilibrium
x: time. y, per snapshot: the Sigma_gas-weighted median over the clean columns c of r_c = P_tot,c / W_c, i.e. the value
r for which Σ_{c: r_c ≤ r} Sigma_gas,c = ½ Σ_c Sigma_gas,c. Light: full columns (|z| < 1.5 kpc about the grid plane); dark:
layer columns (|z − z_mid| < 0.5 kpc), with P_tot and W as defined above in both cases. Result: layer 1.05, 0.94, 0.73, 0.96
in the four phases, never outside 0.5 to 2; the full column drops to 0.4 between the second and third passage because it
then includes the dark-matter weight of gas streaming outside the layer.

![equilibrium](figs/story_1_equilibrium.png)

## Pressure per unit star formation
x: time. y, per snapshot, with c running over the clean layer columns of that snapshot:

- y = [ Σ_c P_tot,c ] / [ Σ_c Sigma_SFR,10,c ] / (4.81 × 10³)   in km/s,

where P_tot,c is the layer P_tot of column c (K cm⁻³) and Sigma_SFR,10,c = (1/A) Σ_{stars in c, formed in the run,
age < 10 Myr} m_* / (10 Myr) in Msun yr⁻¹ kpc⁻²; the divisor turns K cm⁻³ per Msun yr⁻¹ kpc⁻² into km/s. This is the
pressure the layer holds per unit star formation rate, the quantity Ostriker & Kim call the yield. Dashed: their
Upsilon_tot(P̄) = 10^(−0.212 log P̄ + 3.86) km/s at P̄ = Σ_c Sigma_gas,2p,c P_tot,c / Σ_c Sigma_gas,2p,c. Result: within a factor of 2 of the yield before the second passage, 10 to 100 times above it
in the quiet intervals after, back to the yield during the two nuclear bursts.

![pressure per unit star formation](figs/story_2_feedback.png)

## Magnetic field
x: time. y: B_c = (8π k_B P_mag,c)^(1/2) per layer column, in microgauss, with P_mag,c the slab magnetic pressure defined
above. Dark: Σ_c Sigma_gas,c B_c / Σ_c Sigma_gas,c over the clean columns; light: the median of B_c over the columns.
Result: 10 nG seed, e-folding time 10 Myr, 1 to 2 microgauss from 50 Myr, a factor 5 jump at the second passage, 10 to 20
microgauss after the third.

![dynamo](figs/story_4_dynamo.png)

## Pressure shares
x: time. y, per snapshot with c over the clean layer columns: Σ_c P_th,c / Σ_c P_tot,c (thermal), Σ_c P_turb,c / Σ_c P_tot,c
(turbulent), Σ_c Π_mag,c / Σ_c P_tot,c (Maxwell stress), with P_tot = P_th + P_turb + Π_mag so the three sum to one;
dashed: Σ_c P_mag,c / Σ_c P_tot,c, the magnetic pressure B²/8π relative to the total support, which is not part of the sum.
Result: the Maxwell stress is below 5 % of the support before the second passage, 20 % between the second and third, and
25 to 45 % from the third passage to coalescence; the magnetic pressure itself equals the total support at the second
passage and exceeds it after the third, i.e. the field is at or above equipartition even though only part of it supports
the layer vertically.

![pressure shares](figs/pressure_shares.png)

## Clouds
*Clouds* are friends-of-friends groups (linking length 3 pc, at least 25 particles = 100 Msun) of gas with T < 1000 K and
n_H > 10 cm⁻³, on every tenth snapshot, over the whole box. *Clumps* are the same with n_H > 100 cm⁻³ and 1.5 pc linking.
Stars carry the ID of the gas particle they formed from, so the stars formed from a cloud's members are counted exactly.
Per cloud, with sums over its members: M = Σ m_i; centre x_c = Σ m_i x_i / M; r_h = radius about x_c containing M/2;
v̄ = Σ m_i v_i / M; sigma_3d² = Σ m_i |v_i − v̄|² / M; alpha_vir = 5 (sigma_3d²/3) r_h / (G M); v_A² = (1/N) Σ B_i² / (4π rho_i).

Cloud mass function. x: M. y: N_k / (N_phase ΔM_k), the number of clouds in mass bin k divided by the bin width and by the
total number of clouds in the phase; error √N_k / (N_phase ΔM_k). Dashed: M^-1.6. Result: identical in all four phases
(maximum-likelihood index above 300 Msun 1.57 to 1.60); only the largest cloud grows, 3 × 10⁵ to 10⁷ Msun.

![cloud mass function](figs/clouds_mf.png)

Efficiency. Clumps are linked from snapshot k to k+1 when each is the other's largest member overlap; a chain of such links
is a lineage. For a lineage with clumps at snapshots k = k_0 … k_1: M_*(k) = mass of stars formed between snapshots k and
k+1 whose parent gas particle was a member at k; M_*,10(k_1) = the same over the 10 Myr after k_1;
epsilon_int = [ Σ_{k<k_1} M_*(k) + M_*,10(k_1) ] / max_k M(k). The pre-onset snapshot is the last k before M_*(k) > 0
(or k_0 if stars form at once). x: alpha_vir at the pre-onset snapshot. y: epsilon_int. Points: lineages with
max_k M(k) ≥ 300 Msun, k_0 after 5 Myr, k_1 before 221 Myr; lines: medians per alpha bin per phase, phase by onset time.
Result: about 5 % for alpha below 2 to 4, 0.1 % above 10; the fraction of clumps below the threshold falls from 0.9 to 0.06
across the run.

![cloud efficiency](figs/clouds_eff.png)

## Cluster mass function
*Clusters*: friends-of-friends groups (linking 5 pc, at least 25 particles) of stars younger than 10 Myr on every tenth
snapshot, kept if E_kin + E_pot < 0 with E_kin = ½ Σ m_i |v_i − v̄|² and E_pot the pairwise Plummer-softened potential
energy (softening 1 pc), and not the nucleus (fewer than 40 000 particles and r_h < 20 pc), lying in a clean column.
x: cluster mass. y: N_k / Δlog M_k per phase, error √N_k / Δlog M_k. Lines: alpha = 1 + N / Σ_{M_i ≥ 300} ln(M_i / 300)
over the clusters above 300 Msun, drawn over the fitted range. Result: alpha = 1.92, 1.77, 1.74, 1.52; largest cluster
6.7 × 10³, 2.0 × 10⁴, 2.8 × 10⁴, 1.0 × 10⁵ Msun. Below 300 Msun after the second passage only a third of the groups are
bound, so the function is quoted from 300 Msun up.

![MF per phase](figs/story_5_mf.png)

## Largest cluster per burst
x: M_young,c = Σ_{stars in full column c, formed in the run, age < 10 Myr} m_*. y: M_max,c = max over the bound clusters
whose centre lies in column c of the cluster mass. Clean columns with M_young > 500 Msun and at least one bound cluster,
every tenth snapshot. Line: median of M_max per M_young bin; red: M_max = 0.5 M_young; dotted: M_max = M_young. Result:
the median follows half the young mass up to 10⁵ Msun; the nuclear bursts above that split their stars over several
clusters.

![reservoir](figs/story_6_reservoir.png)

## Every cluster against its patch
Each bound cluster of the clean sample against the column it lies in (of the two grids, the one in which that column's
intruder fraction is lower). Left x: layer W_c. Middle x: layer P_tot,c. Right x: Sigma_SFR,10,c. y: cluster mass.
Dashed lines: 0.5 · 10^(1.17 log W − 7.32) · A · τ, 0.5 · 10^(1.18 log P − 7.43) · A · τ, and 0.5 · Sigma_SFR · A · τ,
with A = 0.25 kpc² and τ = 10⁷ yr: half the stars a column forms in 10 Myr at the OK22 rate for its weight, for its
pressure, and at its own rate. Result: cluster masses fill two decades below the lines; 9 % lie above the weight line.

![cluster environment](figs/cluster_env.png)

## Line of sight
As the vertical-equilibrium figure with the column normal replaced by n̂_θ = cos θ n̂ + sin θ ê₂, θ = 0, 30, 60, 90
degrees, ê₂ the second in-plane axis of the grid. Each θ has its own full-column run and its own layer cut about the
midplane found along n̂_θ. y: the Sigma_gas-weighted median of layer P_tot / W over the columns, computed per grid and
combined over the two grids with weights equal to their column counts. Every fourth snapshot. Result: before the second
passage θ = 0 gives one and θ = 30 gives 0.4 to 0.7; between the second and third passage all θ give 0.6 to 1.2; after the
third all θ give 0.9 to 1.3.

![line of sight](figs/los_layer.png)

## Column size
x: column depth 2 z_col: 3 kpc (full column, about the grid plane) or 1, 0.5, 0.25 kpc (z_col = 0.5, 0.25, 0.125 kpc about
the midplane of the 0.5 kpc reference column the point falls in). Marker: footprint side L = 0.5, 0.25, 0.125 kpc (grids of
12², 24², 48² columns over ±3 kpc), A = L². y: for each configuration, the median over the snapshots of the phase (every
fourth snapshot) of r(t), the per-snapshot Sigma_gas-weighted median of P_tot / W. Shaded: 2H. Result: r depends on the
depth only; L from 0.5 to 0.125 kpc at fixed depth changes nothing, depths of 0.5 and 0.25 kpc give 1.5 and 2.7.

![column depth](figs/scale_depth.png)

## Low-mass end of the cluster mass function
x: cluster mass bin [M_1, M_2). y: alpha maximising L(alpha) = −alpha Σ_i ln M_i − N ln[ (M_1^(1−alpha) − M_2^(1−alpha)) /
(alpha − 1) ] over the N groups with M_i in the bin, on a grid alpha = 0.2 … 4; error from the curvature of L at the
maximum. Filled: all friends-of-friends groups of the 10-particle catalogue (5 pc linking); open: bound ones only; per
phase, every tenth snapshot. Result: above 300 Msun bound and all agree within the errors; in the 100 to 300 Msun bin
after the second passage the bound-only index collapses because only a third of those groups are bound.

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
Columns of the 18 snapshots with the turbulence decomposition, restricted to Sigma_SFR,10 = 0. A linear velocity field
v(x) = v_0 + G · (x − x_0) is fitted by least squares to the gas of the column; S = || ½ (G + Gᵀ) − ⅓ tr(G) I ||_F is the
shear rate in km/s per kpc; sigma_tot² = variance of v_n over the column's gas and sigma_res² the variance of v_n after
subtracting the fitted field. x: Sigma_gas H S² × 4.90 × 10⁻⁶ with Sigma_gas in Msun kpc⁻² and H in kpc, i.e. in K cm⁻³.
y: P_turb (sigma_res / sigma_tot)². Line: y = 0.02 x. Result: one coefficient in every phase, no offset.

![shear closure](figs/story_3_shear.png)

## Largest cluster and weight
25 Myr intervals. x: the 90th percentile of { W_c : M_young,c > 500 Msun } over the columns and snapshots of the interval,
W_c the layer weight. y: max of the bound-cluster masses over the interval. Labels: interval start in Myr. Line:
0.5 · 10^(1.17 log W − 7.32) · A · τ, A = 0.25 kpc², τ = 10⁷ yr, nothing fitted. Open circle: the interval in which fewer
than 5 % of the clean columns have M_young > 500. Result: seven active intervals follow the line, slope 0.94 against the
median weight with 0.09 dex scatter; the quiescent interval lies 100 times below.

![ceiling](figs/story_7_ceiling.png)

## Duty cycle
x: layer W_c of a clean column. y: points, per 0.5 dex bin in W: N(M_young,c > 500 Msun and W_c in bin) / N(W_c in bin);
lines: P(W) = [1 + exp(−k (log W − log W_50))]⁻¹ with k and W_50 by maximum likelihood over the columns of the phase.
Result: log W_50 = 3.75, 3.86, 4.73, 4.95; the same fit against rho_mid moves by less than 0.5 dex between phases.

![duty cycle](figs/duty_cycle.png)

## Burst-mass distribution
Phase 40 to 110 Myr. Black: N_k / ΔM_k over the bursting columns (M_young > 500 Msun) of the phase, all snapshots,
error √N_k / ΔM_k. Blue: 25 realisations of the null model: every clean column of the phase bursts with probability P(W_c)
from the logistic above; a bursting column gets M = 10^(a + b log W_c + s ξ), ξ ~ N(0,1), with b from a least-squares fit
of log M_young on log W over the phase's bursts and (a, s) from a lognormal likelihood truncated at 500 Msun; masses below
500 are dropped and the rest binned as the data. All ingredients are fitted to the same bursts, so the comparison tests only
whether the histogram contains structure beyond them. Result: it does not.

![burst kernel](figs/burst_kernel.png)

## Clump virial parameter and efficiency
Left. x: sigma_eff = (P_tot,c / rho_mid,c)^(1/2) of the full column c containing the clump at its pre-onset snapshot, with
P_tot in Msun (km/s)² kpc⁻³ and rho_mid in Msun kpc⁻³ so that sigma_eff is in km/s. y: alpha_vir of the clump at that
snapshot (definition under Clouds). Points: the lineages of the efficiency figure; lines: medians per sigma_eff bin per
phase; dashed: alpha ∝ sigma_eff². Result: one relation for all phases with 0.2 dex offsets.

![clump virial parameter](figs/partC_alpha.png)

Right. x: alpha_vir,tot = alpha_vir (1 + v_A² / sigma_3d²) of the clump at its pre-onset snapshot, v_A² = (1/N) Σ B_i² /
(4π rho_i) over its members. y: epsilon_int. Black: median over all phases per bin. Dashed: epsilon = epsilon_u +
(epsilon_b − epsilon_u) / [1 + (alpha/alpha_c)^m] fitted by least absolute deviation in log epsilon: epsilon_b = 0.050,
epsilon_u = 0.0012, alpha_c = 4.0, m = 4.7; residual scatter 0.69 dex.

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
