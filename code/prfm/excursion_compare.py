"""Excursion-set mass functions driven by the measured PRFM patch state, summed
over clean patches, versus the simulated bound-cluster mass function from the same
patches and snapshots.  Three object models on the same log-density field:
    grav    first crossing of the self-gravity barrier (Hopkins 2012)
    thresh  first crossing of a constant density threshold rho_th = rho0 e^Bth
    hier    all barrier crossings from below (hierarchical sub-objects), Monte Carlo
Cluster mass = EPS x object mass; EPS (and Bth) fitted by maximum likelihood to the
observed cluster masses (shape only). Slopes with the same MLE estimator on both sides.

usage: python excursion_compare.py [--cs 0.5] [--p 1.0] [--fitstep 10] [--mmin 300]
out:   figs/cl_fig6_excursion_cs{cs}_p{p}.png, patches_excursion_cs{cs}_p{p}.npz
"""
import sys, os, argparse, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import common as C
import excursion as X

ap = argparse.ArgumentParser()
ap.add_argument('--cs', type=float, default=0.6)
ap.add_argument('--p', type=float, default=1.0)
ap.add_argument('--b2', type=float, default=0.75)
ap.add_argument('--fitstep', type=int, default=10)
ap.add_argument('--mmin', type=float, default=300.0)
ap.add_argument('--fintr', type=float, default=0.1)
ap.add_argument('--mmin_young', type=float, default=500.0)
ap.add_argument('--nhier', type=int, default=25)
ap.add_argument('--tmax', type=float, default=1e9)     # e.g. C.T_MERGED: restrict theory and data to the disc phase (Q>1)
ap.add_argument('--tfb', type=float, nargs='+', default=[3.0, 5.0])   # Myr, feedback-limited model       # patches per epoch for the Monte-Carlo hierarchical model
a = ap.parse_args()
tag = f'cs{a.cs:g}_p{a.p:g}' + (f'_t{a.tmax:g}' if a.tmax < 1e8 else '')
FIG = f'{C.DATADIR}/prfm/figs'; os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': True, 'grid.alpha': 0.25, 'legend.frameon': False, 'figure.dpi': 130})
KK = 1.0227e-3
EPOCHS = [(0, 50), (50, 100), (100, 150), (150, 232)]
COLS = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
PATCH_AREA = 0.5 ** 2 * 1e6
Mg = np.logspace(0, 8, 400)

D = np.load(f'{C.DATADIR}/prfm/patches_clusters.npz')
E = np.load(f'{C.DATADIR}/clusters/clusters_env.npz')
good = (D['M_young'] > a.mmin_young) & ((D['f_intruder'] < a.fintr) | C.is_merged(D['sep'], D['t'])) & (D['Sigma_gas'] > 1) & (np.mod(D['snap'], a.fitstep) == 0) \
       & (D['t'] < a.tmax) & np.isfinite(D['kappa']) & (D['kappa'] > 0) & (D['H'] > 0) & (D['rho_mid_2p'] > 0) & (D['Pturb_2p'] > 0) & (D['Sigma_gas_2p'] > 0.5) & (D['sigma_turb'] > a.cs)
idx = np.where(good)[0]
print(f'{len(idx)} clean patches at every {a.fitstep}th snapshot')
def state(i):
    return dict(Sigma=D['Sigma_gas_2p'][i], h=D['H'][i] * 1e3, vt_kms=D['sigma_turb'][i], cs_kms=a.cs,
                kappa=D['kappa'][i] * KK, Omega=max(D['Omega'][i], 1e-3) * KK, rho0=D['rho_mid_2p'][i], b2=a.b2, p=a.p)
def epoch_of(t):
    for e in EPOCHS:
        if e[0] <= t < e[1]: return e

def regrid(M, d):
    ok = (d > 0) & np.isfinite(d) & np.isfinite(M) & (M > 0)
    if ok.sum() < 5: return None
    o = np.argsort(M[ok])
    return 10 ** np.interp(np.log10(Mg), np.log10(M[ok][o]), np.log10(d[ok][o]), left=-300, right=-300) * PATCH_AREA

def summed(model, **kw):
    tot = np.zeros_like(Mg); per = {e: np.zeros_like(Mg) for e in EPOCHS}; rec = []
    for i in idx:
        st = state(i)
        if model == 'grav':
            r = X.mass_function(**st)
        elif model == 'thresh':
            r = X.mass_function_threshold(**st, Bth=kw['Bth'])
        elif model == 'gravfb':
            r = X.mass_function(**st, rho_fb=X.rho_feedback(kw['tfb']))
        dn = regrid(r['M'], r['dndM'])
        if dn is None: continue
        tot += dn; per[epoch_of(D['t'][i])] += dn
        ok = r['dndM'] > 0
        rec.append((D['snap'][i], D['t'][i], r.get('Q', np.nan), r['mach'], r['M'][ok][np.argmax(r['M'][ok] ** 2 * r['dndM'][ok])], np.trapezoid(r['f'], r['S'])))
    return tot, per, np.array(rec)

def summed_hier(nper, margin=np.log(2.0), rho_fb=None):
    """Monte-Carlo hierarchical crossings for a random subset of patches per epoch (scaled to all)."""
    rng = np.random.default_rng(0)
    tot = np.zeros_like(Mg); per = {e: np.zeros_like(Mg) for e in EPOCHS}
    for e in EPOCHS:
        ie = [i for i in idx if epoch_of(D['t'][i]) == e]
        if not ie: continue
        sub = rng.choice(ie, min(nper, len(ie)), replace=False)
        for i in sub:
            r = X.mass_function(**state(i), rho_fb=rho_fb)
            S, B, M = r['S'], r['B'], r['M']
            first, allc = X.walks_mass_function(S, B, M, nwalk=20000)
            # counts per S bin -> dN/dM per area: n_obj(M) = Sigma/M * (counts/dS) * |dS/dM| * dM  => dn/dM = Sigma/M * counts_i / |dM_i|
            dM = np.abs(np.gradient(M))
            dn = regrid(M, r['dndM'] * 0 + D['Sigma_gas_2p'][i] / M * allc / np.maximum(dM, 1e-30))
            if dn is None: continue
            w = len(ie) / len(sub)
            tot += dn * w; per[e] += dn * w
    return tot, per

t0 = time.time()
G_tot, G_per, G_rec = summed('grav'); print(f'grav model done ({time.time()-t0:.0f}s)', flush=True)
TH = {}
for Bth in [2.0, 3.0, 4.0, 5.0]:
    TH[Bth] = summed('thresh', Bth=Bth); print(f'threshold model Bth={Bth} done ({time.time()-t0:.0f}s)', flush=True)
H_tot, H_per = summed_hier(a.nhier); print(f'hierarchical MC done ({time.time()-t0:.0f}s)', flush=True)
FB = {}; HFB = {}
for tfb in a.tfb:
    FB[tfb] = summed('gravfb', tfb=tfb); print(f'feedback-limited model t_fb={tfb} Myr (rho_fb={X.rho_feedback(tfb):.1f} Msun/pc3) done ({time.time()-t0:.0f}s)', flush=True)
    HFB[tfb] = summed_hier(a.nhier, rho_fb=X.rho_feedback(tfb)); print(f'hierarchical + feedback t_fb={tfb} done ({time.time()-t0:.0f}s)', flush=True)

# ---- observed
sel = (E['bound'] > 0) & (E['nucleus'] == 0) & ((E['f_intruder'] < a.fintr) | C.is_merged(E['sep'], E['t'])) & (np.mod(E['snap'], a.fitstep) == 0) & (E['t'] < a.tmax)
Mobs = E['M'][sel]; Tobs = E['t'][sel]
def mle(Mv, mmin):
    Mv = Mv[Mv >= mmin]; return (1 + len(Mv) / np.sum(np.log(Mv / mmin)), len(Mv)) if len(Mv) > 5 else (np.nan, len(Mv))
def loglike(dNdM, eps, Mo, mmin):
    Mcl = Mg * eps; pdf = dNdM / eps; s = Mcl >= mmin
    norm = np.trapezoid(pdf[s], Mcl[s])
    if norm <= 0: return -np.inf
    return np.sum((np.interp(np.log10(Mo), np.log10(Mcl), np.log10(pdf + 1e-300)) - np.log10(norm)) * np.log(10))
def fit_eps(dNdM, Mo, mmin):
    Mo = Mo[Mo >= mmin]; grid = np.logspace(-2.5, 0.3, 140)
    ll = np.array([loglike(dNdM, e, Mo, mmin) for e in grid]); j = np.nanargmax(ll)
    return grid[j], ll[j]
def theory_alpha(dNdM, eps, mmin):
    Mcl = Mg * eps; s = (Mcl >= mmin) & (dNdM > 0)
    if s.sum() < 4: return np.nan
    w = dNdM[s] / eps * np.gradient(Mcl[s]); return 1 + w.sum() / np.sum(w * np.log(Mcl[s] / mmin))
def cutoff(dNdM, eps):
    return Mg[np.argmax(Mg ** 2 * dNdM)] * eps

models = {'grav (self-gravity barrier)': G_tot, 'hier (all crossings, MC)': H_tot}
for tfb in a.tfb: models[f'grav+fb (t_fb={tfb:g} Myr)'] = FB[tfb][0]; models[f'hier+fb (t_fb={tfb:g} Myr)'] = HFB[tfb][0]
bestB = max(TH, key=lambda b: fit_eps(TH[b][0], Mobs, a.mmin)[1]); models[f'thresh (Bth={bestB:g})'] = TH[bestB][0]
per_models = {'grav (self-gravity barrier)': G_per, 'hier (all crossings, MC)': H_per, f'thresh (Bth={bestB:g})': TH[bestB][1]}
for tfb in a.tfb: per_models[f'grav+fb (t_fb={tfb:g} Myr)'] = FB[tfb][1]; per_models[f'hier+fb (t_fb={tfb:g} Myr)'] = HFB[tfb][1]
al_o, n_o = mle(Mobs, a.mmin)
print(f'\nobserved: {n_o} bound clusters >= {a.mmin:g} Msun in clean patches; alpha = {al_o:.2f}; M_max = {Mobs.max():.2e}')
print(f'theory patches: median Q = {np.nanmedian(G_rec[:,2]):.2f}, median Mach = {np.nanmedian(G_rec[:,3]):.0f}, median f_coll = {np.nanmedian(G_rec[:,5]):.3f}')
fits = {}
for name, tot in models.items():
    eps, ll = fit_eps(tot, Mobs, a.mmin); fits[name] = (eps, ll)
    print(f'{name:32s}: eps = {eps:.3f}  lnL = {ll:9.1f}  alpha(theory) = {theory_alpha(tot, eps, a.mmin):.2f}  cutoff x eps = {cutoff(tot, eps):.2e}')
for Bth, (tot, _, _) in TH.items():
    eps, ll = fit_eps(tot, Mobs, a.mmin); print(f'   thresh Bth={Bth:g}: eps={eps:.3f} lnL={ll:.1f} alpha={theory_alpha(tot, eps, a.mmin):.2f} cutoff x eps={cutoff(tot, eps):.2e}')
print('\nper epoch (grav model): eps, alpha obs / theory, cutoff x eps / observed M_max, median Mach, Q')
erows = []
for e in EPOCHS:
    so = (Tobs >= e[0]) & (Tobs < e[1]); th = G_per[e]
    if so.sum() > 8 and th.max() > 0:
        ee, _ = fit_eps(th, Mobs[so], a.mmin); ao, no = mle(Mobs[so], a.mmin); at = theory_alpha(th, ee, a.mmin)
        pe = (G_rec[:, 1] >= e[0]) & (G_rec[:, 1] < e[1])
        erows.append((e, ee, ao, no, at, cutoff(th, ee), Mobs[so].max(), np.nanmedian(G_rec[pe, 3]), np.nanmedian(G_rec[pe, 2])))
        print(f'  {e[0]:3d}-{e[1]:3d} Myr: eps={ee:.3f}  alpha {ao:.2f} (N={no}) / {at:.2f}  cutoff {erows[-1][5]:.2e} / {erows[-1][6]:.2e}  Mach {erows[-1][7]:.0f}  Q {erows[-1][8]:.2f}')
print('\nper epoch at the global eps of each model: alpha theory / cutoff x eps   (sim alpha, M_max in brackets)')
for name in models:
    ee = fits[name][0]; line = f'  {name:28s}'
    for e in EPOCHS:
        so = (Tobs >= e[0]) & (Tobs < e[1]); th = per_models[name][e]
        if so.sum() > 8 and th.max() > 0:
            line += f' | {e[0]}-{e[1]}: {theory_alpha(th, ee, a.mmin):.2f} / {cutoff(th, ee):.1e} ({mle(Mobs[so], a.mmin)[0]:.2f}, {Mobs[so].max():.1e})'
    print(line)
mmins = np.logspace(2, 4, 13)
print('\ncurvature alpha(M>Mmin): Mmin, sim, ' + ', '.join(models))
curv = {name: [theory_alpha(tot, fits[name][0], mm) for mm in mmins] for name, tot in models.items()}
al_s = [mle(Mobs, mm)[0] for mm in mmins]; ns = [mle(Mobs, mm)[1] for mm in mmins]
for j, mm in enumerate(mmins): print(f'   {mm:7.0f}  {al_s[j]:.2f}  ' + '  '.join(f'{curv[n][j]:.2f}' for n in models))
np.savez(f'{C.DATADIR}/prfm/patches_excursion_{tag}.npz', Mg=Mg, grav=G_tot, hier=H_tot, thresh=TH[bestB][0], Bth=bestB, rec=G_rec, Mobs=Mobs, Tobs=Tobs,
         **{f'grav_{e[0]}_{e[1]}': v for e, v in G_per.items()})

# ---- figure
fig, ax = plt.subplots(1, 4, figsize=(19, 4.5))
edges = np.logspace(2, 6.5, 19); c = np.sqrt(edges[1:] * edges[:-1]); w = np.diff(edges)
hh, _ = np.histogram(Mobs, edges); ok = hh > 0
ax[0].errorbar(c[ok], hh[ok] / w[ok], yerr=np.sqrt(hh[ok]) / w[ok], fmt='o', color='k', ms=4, label=f'simulation: bound clusters, clean patches (N={len(Mobs)})')
Nabove = np.sum(hh[c >= a.mmin])
STY = [('-', 'r'), ('-.', 'm'), ('--', 'b'), ('-', 'g'), ('-.', 'g'), ('-', 'c'), ('-.', 'c'), (':', 'y')]
for (name, tot), (ls, col) in zip(models.items(), STY):
    eps, ll = fits[name]; Mcl = Mg * eps; pdf = tot / eps; s = Mcl >= a.mmin
    nrm = Nabove / max(np.trapezoid(pdf[s], Mcl[s]), 1e-300)
    ax[0].plot(Mcl, pdf * nrm, ls=ls, color=col, lw=1.8, label=rf'{name}: $\epsilon$={eps:.2f}, $\Delta$lnL={ll - max(f[1] for f in fits.values()):.0f}')
xx = np.array([a.mmin, 3e5]); A0 = (hh[ok] / w[ok])[np.argmin(np.abs(c[ok] - 1e3))] * 1e3 ** 2
ax[0].plot(xx, A0 * xx ** -2, 'k:', lw=1, label=r'$M^{-2}$')
ax[0].set(xscale='log', yscale='log', xlabel=r'$M_{\rm cl}$ [M$_\odot$]', ylabel=r'$dN/dM$', xlim=(1e2, 3e6), ylim=(1e-4, 1e2))
ax[0].set_title(f'shape-normalised above {a.mmin:g} Msun;  c_s={a.cs} km/s, p={a.p}', fontsize=9); ax[0].legend(fontsize=6.5, loc='lower left')
for (e, ee, ao, no, at, Mct, Mmx, mm_, qq), col in zip(erows, COLS):
    so = (Tobs >= e[0]) & (Tobs < e[1]); ms = np.sort(Mobs[so])[::-1]
    ax[1].plot(ms, np.arange(1, len(ms) + 1), color=col, lw=1.5, label=f'{e[0]}-{e[1]} Myr: sim α={ao:.2f}, grav α={at:.2f}')
    th = G_per[e]; Mcl_e = Mg * ee; cum = np.cumsum((th / ee * np.gradient(Mcl_e))[::-1])[::-1]
    s_ = Mcl_e >= a.mmin; ax[1].plot(Mcl_e, cum * np.sum(so & (Mobs >= a.mmin)) / max(cum[s_][0], 1e-300), color=col, ls='--', lw=1)
ax[1].set(xscale='log', yscale='log', xlabel=r'$M_{\rm cl}$ [M$_\odot$]', ylabel=r'$N(>M)$', xlim=(1e2, 3e6), ylim=(0.8, 2e3)); ax[1].legend(fontsize=7)
ax[1].set_title('per epoch: solid sim, dashed grav model', fontsize=9)
tt = [r[0][0] + 0.5 * (r[0][1] - r[0][0]) for r in erows]
ax[2].errorbar(tt, [r[2] for r in erows], [(r[2] - 1) / np.sqrt(r[3]) for r in erows], fmt='ko', capsize=3, label=r'$\alpha$ simulation')
ax[2].plot(tt, [r[4] for r in erows], 'r-s', label=r'$\alpha$ grav model'); ax[2].axhline(2, color='grey', ls=':')
ax2 = ax[2].twinx(); ax2.plot(tt, [r[5] for r in erows], 'r--^', label=r'cutoff grav $\times\epsilon$'); ax2.plot(tt, [r[6] for r in erows], 'k--v', label=r'$M_{\max}$ sim')
ax2.set(yscale='log', ylabel=r'M [M$_\odot$]'); ax2.grid(False)
ax[2].set(xlabel='t [Myr]', ylabel=r'$\alpha$', ylim=(1.2, 2.8)); ax[2].legend(fontsize=7, loc='upper left'); ax2.legend(fontsize=7, loc='upper right')
ax[3].errorbar(mmins, al_s, [(x - 1) / np.sqrt(max(n, 1)) if np.isfinite(x) else 0 for x, n in zip(al_s, ns)], fmt='ko', capsize=3, label='simulation')
for (name, cv), (ls, col) in zip(curv.items(), STY):
    ax[3].plot(mmins, cv, ls=ls, color=col, lw=1.8, label=name)
ax[3].axhline(2, color='grey', ls=':'); ax[3].set(xscale='log', xlabel=r'$M_{\min}$ [M$_\odot$]', ylabel=r'$\alpha(M>M_{\min})$', ylim=(1.0, 3.0)); ax[3].legend(fontsize=7)
ax[3].set_title('curvature: MLE slope above a running lower limit', fontsize=9)
fig.tight_layout(); fig.savefig(f'{FIG}/cl_fig6_excursion_{tag}.png'); print('figure:', f'{FIG}/cl_fig6_excursion_{tag}.png')

# ================================================================ is the cutoff real? size-of-sample and truncation tests
rng = np.random.default_rng(3)
def draw(dNdM, eps, n, mmin):
    Mcl = Mg * eps; s = (Mcl >= mmin) & (dNdM > 0)
    x = Mcl[s]; pdf = dNdM[s] / eps
    cdf = np.concatenate([[0], np.cumsum(0.5 * (pdf[1:] + pdf[:-1]) * np.diff(x))]); cdf /= cdf[-1]
    return np.interp(rng.random(n), cdf, x)
def trunc_mle(Mv, mmin):
    """MLE of (alpha, Mc) for dN/dM ~ M^-alpha exp(-M/Mc) above mmin vs pure power law; returns (alpha_pl, lnL_pl, alpha_tr, Mc, lnL_tr)."""
    Mv = Mv[Mv >= mmin]; n = len(Mv)
    al_pl = 1 + n / np.sum(np.log(Mv / mmin)); ll_pl = n * np.log(al_pl - 1) + n * (al_pl - 1) * np.log(mmin) - al_pl * np.sum(np.log(Mv))
    best = (al_pl, np.inf, ll_pl)
    x = np.logspace(np.log10(mmin), 8, 800)
    for lMc in np.linspace(np.log10(mmin) + 0.3, 7.5, 60):
        Mc = 10 ** lMc
        for al in np.linspace(1.0, 3.0, 41):
            pdf = x ** -al * np.exp(-x / Mc); Z = np.trapezoid(pdf, x)
            ll = np.sum(-al * np.log(Mv) - Mv / Mc) - n * np.log(Z)
            if ll > best[2]: best = (al, Mc, ll)
    return al_pl, ll_pl, best[0], best[1], best[2]
print('\n=== cutoff reality checks (above %g Msun) ===' % a.mmin)
al_pl, ll_pl, al_tr, Mc_tr, ll_tr = trunc_mle(Mobs, a.mmin)
print(f'observed (N={np.sum(Mobs >= a.mmin)}): pure power law alpha={al_pl:.2f}; truncated alpha={al_tr:.2f}, Mc={Mc_tr:.2e}; Delta lnL (trunc - pl) = {ll_tr - ll_pl:.2f}')
for name, tot in models.items():
    eps = fits[name][0]; n = int(np.sum(Mobs >= a.mmin))
    mx = np.array([draw(tot, eps, n, a.mmin).max() for _ in range(1500)])
    dll = []
    for _ in range(60):
        mv = draw(tot, eps, n, a.mmin); r = trunc_mle(mv, a.mmin); dll.append(r[4] - r[1])
    print(f'  {name:28s}: predicted M_max for N={n}: median {np.median(mx):.2e}, 5-95% {np.percentile(mx,5):.2e}-{np.percentile(mx,95):.2e}  (observed {Mobs.max():.2e}); '
          f'mock Delta lnL(trunc-pl): median {np.median(dll):.2f}, 95th {np.percentile(dll,95):.2f}')
print('per epoch, self-gravity model at global eps: predicted M_max range vs observed')
eps = fits['grav (self-gravity barrier)'][0]
for e in EPOCHS:
    so = (Tobs >= e[0]) & (Tobs < e[1]) & (Mobs >= a.mmin); th = G_per[e]
    if so.sum() > 8 and th.max() > 0:
        mx = np.array([draw(th, eps, int(so.sum()), a.mmin).max() for _ in range(1500)])
        print(f'  {e[0]:3d}-{e[1]:3d} Myr (N={so.sum()}): predicted {np.percentile(mx,5):.1e} - {np.median(mx):.1e} - {np.percentile(mx,95):.1e}, observed {Mobs[so].max():.1e}')
