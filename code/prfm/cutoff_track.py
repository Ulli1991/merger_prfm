"""Does the PRFM-driven excursion-set cutoff track the observed truncation of the
bound-cluster mass function through the merger?

Self-gravity barrier model only.  One global cluster-formation efficiency EPS,
fitted NOT to the shape but to the truncation: the likelihood of the observed
maximum cluster mass in each time bin, given N_bin draws from the theory mass
function of the clean patches in that bin (pdf of the maximum of N iid draws,
P(max<m) = F(m)^N).  Control ("null"): the same statistic with a single global
mass function (all patches, all times), so that only N_bin varies between bins,
i.e. size-of-sample without any environmental dependence.  Also reports the
Schechter-truncation MLE per epoch on both sides.

usage: python cutoff_track.py [--cs 0.2] [--p 1.0] [--mmin 300] [--dt 25]
out:   figs/cl_fig7_cutoff_{tag}.png, cutoff_track_{tag}.npz
"""
import sys, os, argparse, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import common as C
import excursion as X

ap = argparse.ArgumentParser()
ap.add_argument('--cs', type=float, default=0.6)        # cold-phase sound speed [km/s], measured (cold_cs.py); same as CS in cluster_pressure.py
ap.add_argument('--p', type=float, default=1.0)
ap.add_argument('--b2', type=float, default=0.75)
ap.add_argument('--fitstep', type=int, default=10)
ap.add_argument('--mmin', type=float, default=300.0)
ap.add_argument('--fintr', type=float, default=0.1)
ap.add_argument('--mmin_young', type=float, default=500.0)
ap.add_argument('--dt', type=float, default=25.0)       # time bin [Myr]
ap.add_argument('--nmin_bin', type=int, default=8)      # min clusters per bin for the M_max statistic
ap.add_argument('--tfb', type=float, default=0.0)       # >0: feedback-limited barrier, objects must reach rho with t_ff = tfb [Myr]
ap.add_argument('--rmax_pc', type=float, default=0.0)   # >0: largest object scale [pc] (e.g. patch half-size 250)
a = ap.parse_args()
tag = f'cs{a.cs:g}_p{a.p:g}' + (f'_fb{a.tfb:g}' if a.tfb > 0 else '') + (f'_r{a.rmax_pc:g}' if a.rmax_pc > 0 else '')
FIG = f'{C.DATADIR}/prfm/figs'; os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': True, 'grid.alpha': 0.25, 'legend.frameon': False, 'figure.dpi': 130})
KK = 1.0227e-3
PATCH_AREA = 0.5 ** 2 * 1e6
Mg = np.logspace(0, 8, 600)
EPOCHS = [(0, 50), (50, 100), (100, 150)]

D = np.load(f'{C.DATADIR}/prfm/patches_clusters.npz')
E = np.load(f'{C.DATADIR}/clusters/clusters_env.npz')
good = (D['M_young'] > a.mmin_young) & ((D['f_intruder'] < a.fintr) | C.is_merged(D['sep'], D['t'])) & (D['Sigma_gas'] > 1) & (np.mod(D['snap'], a.fitstep) == 0) \
       & np.isfinite(D['kappa']) & (D['kappa'] > 0) & (D['H'] > 0) & (D['rho_mid_2p'] > 0) & (D['Pturb_2p'] > 0) & (D['Sigma_gas_2p'] > 0.5) & (D['sigma_turb'] > a.cs)
idx = np.where(good)[0]
print(f'{len(idx)} clean patches at every {a.fitstep}th snapshot; cs = {a.cs} km/s, p = {a.p}, t_fb = {a.tfb} Myr, R_max = {a.rmax_pc} pc')

def state(i):
    return dict(Sigma=D['Sigma_gas_2p'][i], h=D['H'][i] * 1e3, vt_kms=D['sigma_turb'][i], cs_kms=a.cs,
                kappa=D['kappa'][i] * KK, Omega=max(D['Omega'][i], 1e-3) * KK, rho0=D['rho_mid_2p'][i], b2=a.b2, p=a.p)

def regrid(M, d):
    ok = (d > 0) & np.isfinite(d) & np.isfinite(M) & (M > 0)
    if ok.sum() < 5: return None
    o = np.argsort(M[ok])
    return 10 ** np.interp(np.log10(Mg), np.log10(M[ok][o]), np.log10(d[ok][o]), left=-300, right=-300) * PATCH_AREA

# ---- per-patch theory mass functions (object mass), Mach numbers consistent with cluster_pressure.py
t0 = time.time()
MF = {}; REC = []
for i in idx:
    st = state(i); kw = {}
    if a.tfb > 0: kw['rho_fb'] = X.rho_feedback(a.tfb)
    if a.rmax_pc > 0: kw['Rmax_h'] = max(a.rmax_pc / st['h'], 1.5)
    r = X.mass_function(**st, **kw); dn = regrid(r['M'], r['dndM'])
    if dn is None: continue
    MF[i] = dn
    ok = r['dndM'] > 0
    N1 = np.cumsum((dn * np.gradient(Mg))[::-1])[::-1]                       # N(>M) per patch
    M_N1 = Mg[np.argmin(np.abs(N1 - 1.0))] if N1.max() > 1 else np.nan       # object mass with one expected object above
    REC.append((D['snap'][i], D['t'][i], r['Q'], r['mach'], D['Mach_cold'][i], D['W'][i], D['sigma_turb'][i], D['H'][i] * 1e3,
                D['Sigma_gas_2p'][i], D['kappa'][i], D['M_max'][i], M_N1, r['M'][ok][np.argmax(r['M'][ok] ** 2 * r['dndM'][ok])]))
REC = np.array(REC); idx = np.array(list(MF))
print(f'grav model done ({time.time()-t0:.0f}s); theory Mach = sigma_turb/cs: median {np.median(REC[:,3]):.1f}; '
      f'Mach_cold from cluster_pressure.py: median {np.median(REC[:,4]):.1f} (identical iff cs = 0.2)')

# ---- observed clusters, same selection
sel = (E['bound'] > 0) & (E['nucleus'] == 0) & ((E['f_intruder'] < a.fintr) | C.is_merged(E['sep'], E['t'])) & (np.mod(E['snap'], a.fitstep) == 0) & (E['M'] >= a.mmin)
Mobs = E['M'][sel]; Tobs = E['t'][sel]
print(f'{len(Mobs)} bound clusters >= {a.mmin:g} Msun in clean patches; M_max = {Mobs.max():.2e}')

# ---- time bins
edges = np.arange(0, 150 + a.dt, a.dt)
bins = []
for lo, hi in zip(edges[:-1], edges[1:]):
    so = (Tobs >= lo) & (Tobs < hi); pi = idx[(D['t'][idx] >= lo) & (D['t'][idx] < hi)]
    if so.sum() >= a.nmin_bin and len(pi) > 0:
        bins.append(dict(lo=lo, hi=hi, t=0.5 * (lo + hi), N=int(so.sum()), Mmax=Mobs[so].max(), M=np.sort(Mobs[so]),
                         mf=sum(MF[i] for i in pi), npatch=len(pi),
                         W=np.median(D['W'][pi]), mach=np.median(REC[np.isin(idx, pi), 3])))
mf_all = sum(MF.values())
print(f'{len(bins)} time bins of {a.dt:g} Myr with >= {a.nmin_bin} clusters')

# ---- likelihood of the observed maximum
def cdf_pdf(mf, eps, m):
    """F(m), f(m) of the cluster-mass pdf (theory dN/dM, cluster mass = eps x object mass) truncated at mmin."""
    Mcl = Mg * eps; s = (Mcl >= a.mmin) & (mf > 0)
    x = Mcl[s]; pdf = mf[s] / eps
    cum = np.concatenate([[0], np.cumsum(0.5 * (pdf[1:] + pdf[:-1]) * np.diff(x))])
    Z = cum[-1]
    if Z <= 0: return np.nan, np.nan
    F = np.interp(m, x, cum / Z); f = np.interp(m, x, pdf / Z)
    return F, f
def lnL_max(eps, use_bin_mf=True):
    ll = 0.0
    for b in bins:
        F, f = cdf_pdf(b['mf'] if use_bin_mf else mf_all, eps, b['Mmax'])
        if not np.isfinite(F) or F <= 0 or f <= 0: return -np.inf
        ll += np.log(b['N']) + (b['N'] - 1) * np.log(F) + np.log(f)
    return ll
def max_quantiles(mf, eps, N, q=(0.05, 0.5, 0.95)):
    Mcl = Mg * eps; s = (Mcl >= a.mmin) & (mf > 0); x = Mcl[s]; pdf = mf[s] / eps
    cum = np.concatenate([[0], np.cumsum(0.5 * (pdf[1:] + pdf[:-1]) * np.diff(x))]); cum /= cum[-1]
    return [np.interp(qq ** (1.0 / N), cum, x) for qq in q]           # F(m)^N = q

grid = np.logspace(-3.5, 0.5, 200)
llT = np.array([lnL_max(e, True) for e in grid]); llN = np.array([lnL_max(e, False) for e in grid])
jT, jN = np.nanargmax(llT), np.nanargmax(llN); epsT, epsN = grid[jT], grid[jN]
# 1-sigma interval on eps
okT = llT >= llT[jT] - 0.5; okN = llN >= llN[jN] - 0.5
print(f'\nglobal efficiency from the truncation likelihood (max of N draws per bin):')
print(f'  theory (patch state varies): eps = {epsT:.3g} [{grid[okT].min():.3g}, {grid[okT].max():.3g}]  lnL = {llT[jT]:.2f}')
print(f'  null   (one global MF)     : eps = {epsN:.3g} [{grid[okN].min():.3g}, {grid[okN].max():.3g}]  lnL = {llN[jN]:.2f}')
print(f'  Delta lnL (theory - null) = {llT[jT]-llN[jN]:+.2f}   (same number of free parameters)')

print(f'\nper bin at the global eps: N, observed M_max | theory 5/50/95 % | null 5/50/95 % | median W, Mach, n_patch')
rows = []
for b in bins:
    qT = max_quantiles(b['mf'], epsT, b['N']); qN = max_quantiles(mf_all, epsN, b['N'])
    rows.append((b['t'], b['N'], b['Mmax'], *qT, *qN, b['W'], b['npatch']))
    print(f'  {b["lo"]:4.0f}-{b["hi"]:4.0f} Myr  N={b["N"]:3d}  {b["Mmax"]:8.2e} | {qT[0]:8.2e} {qT[1]:8.2e} {qT[2]:8.2e} | {qN[0]:8.2e} {qN[1]:8.2e} {qN[2]:8.2e} | W {b["W"]:8.3g}  n {b["npatch"]}')
rows = np.array(rows)
# how much does the theory cutoff itself move (independent of N)?  M at N(>M)=1 per bin-summed MF and per patch
print('\ntheory cutoff per bin (eps x object mass where M^2 dN/dM peaks, and where N(>M)=1 for the bin), vs observed M_max:')
for b in bins:
    m2 = Mg[np.argmax(Mg ** 2 * b['mf'])] * epsT
    N1 = np.cumsum((b['mf'] * np.gradient(Mg))[::-1])[::-1]; mN1 = Mg[np.argmin(np.abs(N1 - 1.0))] * epsT
    print(f'  {b["lo"]:4.0f}-{b["hi"]:4.0f} Myr: peak M^2dN/dM {m2:8.2e}   N(>M)=1: {mN1:8.2e}   observed M_max {b["Mmax"]:8.2e}')

# ---- Schechter truncation per epoch, both sides
def trunc_mle(Mv, mmin):
    Mv = Mv[Mv >= mmin]; n = len(Mv)
    al_pl = 1 + n / np.sum(np.log(Mv / mmin)); ll_pl = n * np.log(al_pl - 1) + n * (al_pl - 1) * np.log(mmin) - al_pl * np.sum(np.log(Mv))
    best = (al_pl, np.inf, ll_pl); x = np.logspace(np.log10(mmin), 8, 800)
    for lMc in np.linspace(np.log10(mmin) + 0.3, 7.5, 60):
        Mc = 10 ** lMc
        for al in np.linspace(1.0, 3.0, 41):
            Z = np.trapezoid(x ** -al * np.exp(-x / Mc), x)
            ll = np.sum(-al * np.log(Mv) - Mv / Mc) - n * np.log(Z)
            if ll > best[2]: best = (al, Mc, ll)
    return al_pl, best[0], best[1], best[2] - ll_pl
def trunc_fit_mf(mf, eps, mmin):
    """Schechter fit to a tabulated theory dN/dM (least squares in log)."""
    Mcl = Mg * eps; s = (Mcl >= mmin) & (mf > 0) & (mf > 1e-6 * mf.max()); x = Mcl[s]; y = np.log(mf[s])
    best = (np.nan, np.nan, np.inf)
    for lMc in np.linspace(np.log10(mmin) + 0.3, 7.5, 80):
        for al in np.linspace(1.0, 3.0, 41):
            model = -al * np.log(x) - x / 10 ** lMc; c = np.mean(y - model)
            r = np.sum((y - model - c) ** 2)
            if r < best[2]: best = (al, 10 ** lMc, r)
    return best[0], best[1]
print('\nSchechter truncation per epoch: observed (alpha_pl, alpha_tr, Mc, Delta lnL trunc-pl)  |  theory at global eps (alpha, Mc)')
for lo, hi in EPOCHS:
    so = (Tobs >= lo) & (Tobs < hi); pi = idx[(D['t'][idx] >= lo) & (D['t'][idx] < hi)]
    if so.sum() < 10 or len(pi) == 0: continue
    ao, at, Mc, dll = trunc_mle(Mobs[so], a.mmin); mfe = sum(MF[i] for i in pi); alt, Mct = trunc_fit_mf(mfe, epsT, a.mmin)
    print(f'  {lo:3d}-{hi:3d} Myr (N={so.sum()}): obs alpha_pl={ao:.2f} alpha_tr={at:.2f} Mc={Mc:.2e} (dlnL={dll:.2f})  |  theory alpha={alt:.2f} Mc={Mct:.2e}')

np.savez(f'{C.DATADIR}/prfm/cutoff_track_{tag}.npz', rows=rows, rec=REC, grid=grid, llT=llT, llN=llN, epsT=epsT, epsN=epsN, Mg=Mg, mf_all=mf_all)

# ---- figure
fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
ax[0].plot(grid, llT - llT[jT], 'r-', label=f'theory: patch state per bin, eps={epsT:.3g}')
ax[0].plot(grid, llN - llT[jT], 'k--', label=f'null: one global MF, eps={epsN:.3g}')
ax[0].set(xscale='log', xlabel=r'$\epsilon$ (cluster mass / object mass)', ylabel=r'$\ln L(M_{\max}\,|\,N)$ - max', ylim=(-15, 1)); ax[0].legend(fontsize=8)
ax[0].set_title(f'likelihood of the per-bin maximum; c_s={a.cs}, p={a.p}', fontsize=9)
t = rows[:, 0]
ax[1].fill_between(t, rows[:, 3], rows[:, 5], color='r', alpha=0.15); ax[1].plot(t, rows[:, 4], 'r-s', label='theory median (5-95%)')
ax[1].fill_between(t, rows[:, 6], rows[:, 8], color='grey', alpha=0.2); ax[1].plot(t, rows[:, 7], 'k--^', label='null median (5-95%)')
ax[1].plot(t, rows[:, 2], 'ko', ms=6, label=r'observed $M_{\max}$')
for tt, n in zip(t, rows[:, 1]): ax[1].annotate(f'N={int(n)}', (tt, rows[0, 2] * 0.25), fontsize=7, ha='center')
ax[1].set(yscale='log', xlabel='t [Myr]', ylabel=r'$M_{\max}$ [M$_\odot$]'); ax[1].legend(fontsize=8, loc='upper right')
ax[1].set_title(f'maximum bound cluster mass per {a.dt:g} Myr bin', fontsize=9)
ax[2].plot(t, rows[:, 9], 'b-o', label=r'median $\mathcal{W}$ of clean patches')
ax[2].set(yscale='log', xlabel='t [Myr]', ylabel=r'$\mathcal{W}/k_B$ [K cm$^{-3}$]'); ax[2].legend(fontsize=8, loc='upper left')
ax2 = ax[2].twinx(); ax2.plot(t, rows[:, 4] / rows[:, 7], 'r-s', label='theory / null median'); ax2.axhline(1, color='grey', ls=':')
ax2.set(ylabel='predicted M_max ratio', yscale='log'); ax2.grid(False); ax2.legend(fontsize=8, loc='lower right')
ax[2].set_title('environment vs shift of the predicted cutoff', fontsize=9)
fig.tight_layout(); fig.savefig(f'{FIG}/cl_fig7_cutoff_{tag}.png'); print('figure:', f'{FIG}/cl_fig7_cutoff_{tag}.png')
