"""Letter figures.

Fig 1  (theory): the shape of the bound-cluster mass function against the excursion-set
        prediction driven by the measured PRFM state, disc phase (t < T_MERGED).
        (a) dN/dlogM per epoch, simulation vs theory (same colour = same epoch)
        (b) curvature: MLE index above a running lower limit, simulation vs three object models
        (c) index above 300 Msun vs the epoch's median weight: pressure changes, shape does not
Fig 2  (reservoir): (a) M_max vs M_young per patch, (b) truncation and largest cluster vs
        time with the weight, (c) bound fraction vs Toomre Q per patch (the mechanism).

usage: python letter_figs.py [--cs 0.6] [--p 1.0] [--tfb 5] [--mmin 300]
out:   figs/letter_fig1_theory.png/.pdf, figs/letter_fig2_reservoir.png/.pdf
"""
import sys, os, argparse
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.stats import spearmanr
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import common as C
import excursion as X

ap = argparse.ArgumentParser()
ap.add_argument('--cs', type=float, default=0.6)
ap.add_argument('--p', type=float, default=1.0)
ap.add_argument('--tfb', type=float, default=5.0)
ap.add_argument('--mmin', type=float, default=300.0)
ap.add_argument('--fintr', type=float, default=0.1)
ap.add_argument('--mmin_young', type=float, default=500.0)
ap.add_argument('--fitstep', type=int, default=10)
ap.add_argument('--nhier', type=int, default=25)
a = ap.parse_args()
FIG = f'{C.DATADIR}/prfm/figs'; os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({'font.size': 8, 'axes.labelsize': 8.5, 'legend.fontsize': 7, 'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5,
                     'axes.spines.top': True, 'axes.spines.right': True, 'legend.frameon': False, 'figure.dpi': 200,
                     'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True, 'xtick.minor.visible': True, 'ytick.minor.visible': True,
                     'xtick.major.size': 3.5, 'ytick.major.size': 3.5, 'xtick.minor.size': 2, 'ytick.minor.size': 2, 'xtick.minor.width': 0.5, 'ytick.minor.width': 0.5,
                     'axes.linewidth': 0.7, 'xtick.major.width': 0.7, 'ytick.major.width': 0.7, 'lines.linewidth': 1.3,
                     'font.family': 'serif', 'mathtext.fontset': 'dejavuserif'})
KK = 1.0227e-3; PATCH_AREA = 0.5 ** 2 * 1e6
Mg = np.logspace(0, 8, 400)
EPOCHS = [(0, 50), (50, 100), (100, 150), (150, 232)]
ECOL = ['#1f5f9e', '#c8412f', '#2e8b57', '#8c6d31']
ELAB = ['0-50 Myr', '50-100 Myr', '100-150 Myr', '150-232 Myr']
GREY = '#555555'

D = np.load(f'{C.DATADIR}/prfm/patches_clusters.npz')
E = np.load(f'{C.DATADIR}/clusters/clusters_env.npz')
cleanD = (D['f_intruder'] < a.fintr) | C.is_merged(D['sep'], D['t'])
cleanE = (E['f_intruder'] < a.fintr) | C.is_merged(E['sep'], E['t'])

# ============================================================ theory, disc phase
good = (D['M_young'] > a.mmin_young) & cleanD & (D['Sigma_gas'] > 1) & (np.mod(D['snap'], a.fitstep) == 0) & (D['t'] < C.T_MERGED) \
       & np.isfinite(D['kappa']) & (D['kappa'] > 0) & (D['H'] > 0) & (D['rho_mid_2p'] > 0) & (D['Pturb_2p'] > 0) & (D['Sigma_gas_2p'] > 0.5) & (D['sigma_turb'] > a.cs)
idx = np.where(good)[0]
def state(i):
    return dict(Sigma=D['Sigma_gas_2p'][i], h=D['H'][i] * 1e3, vt_kms=D['sigma_turb'][i], cs_kms=a.cs,
                kappa=D['kappa'][i] * KK, Omega=max(D['Omega'][i], 1e-3) * KK, rho0=D['rho_mid_2p'][i], p=a.p, Rmax_h=max(250.0 / (D['H'][i] * 1e3), 1.5))
def epoch_of(t):
    for e in EPOCHS:
        if e[0] <= t < e[1]: return e
def regrid(M, d):
    ok = (d > 0) & np.isfinite(d) & np.isfinite(M) & (M > 0)
    if ok.sum() < 5: return None
    o = np.argsort(M[ok])
    return 10 ** np.interp(np.log10(Mg), np.log10(M[ok][o]), np.log10(d[ok][o]), left=-300, right=-300) * PATCH_AREA
PP = {}   # per-patch theory MF of the main model (for environment bins)
def summed(rho_fb=None, keep=False):
    tot = np.zeros_like(Mg); per = {e: np.zeros_like(Mg) for e in EPOCHS}
    for i in idx:
        r = X.mass_function(**state(i), rho_fb=rho_fb); dn = regrid(r['M'], r['dndM'])
        if dn is None: continue
        tot += dn; per[epoch_of(D['t'][i])] += dn
        if keep: PP[i] = dn
    return tot, per
def summed_hier(nper, rho_fb=None):
    rng = np.random.default_rng(0); tot = np.zeros_like(Mg); per = {e: np.zeros_like(Mg) for e in EPOCHS}
    for e in EPOCHS:
        ie = [i for i in idx if epoch_of(D['t'][i]) == e]
        if not ie: continue
        sub = rng.choice(ie, min(nper, len(ie)), replace=False)
        for i in sub:
            r = X.mass_function(**state(i), rho_fb=rho_fb)
            first, allc = X.walks_mass_function(r['S'], r['B'], r['M'], nwalk=20000)
            dn = regrid(r['M'], D['Sigma_gas_2p'][i] / r['M'] * allc / np.maximum(np.abs(np.gradient(r['M'])), 1e-30))
            if dn is None: continue
            w = len(ie) / len(sub); tot += dn * w; per[e] += dn * w
    return tot, per
rho_fb = X.rho_feedback(a.tfb)
FB_tot, FB_per = summed(rho_fb, keep=True)
G_tot, G_per = summed()
H_tot, H_per = summed_hier(a.nhier)
print(f'{len(idx)} disc-phase clean patches; models computed')

# ---- observed
selE = (E['bound'] > 0) & (E['nucleus'] == 0) & cleanE & (np.mod(E['snap'], a.fitstep) == 0)
Mall, Tall = E['M'][selE], E['t'][selE]
disc = Tall < C.T_MERGED
Mobs, Tobs = Mall[disc], Tall[disc]
def mle(Mv, mmin):
    Mv = Mv[Mv >= mmin]; n = len(Mv)
    return (1 + n / np.sum(np.log(Mv / mmin)), n) if n > 5 else (np.nan, n)
def loglike(dNdM, eps, Mo, mmin):
    Mcl = Mg * eps; pdf = dNdM / eps; s = Mcl >= mmin
    norm = np.trapezoid(pdf[s], Mcl[s])
    if norm <= 0: return -np.inf
    return np.sum((np.interp(np.log10(Mo), np.log10(Mcl), np.log10(pdf + 1e-300)) - np.log10(norm)) * np.log(10))
def fit_eps(dNdM, Mo, mmin):
    Mo = Mo[Mo >= mmin]; grid = np.logspace(-2.5, 0.3, 140)
    ll = np.array([loglike(dNdM, e, Mo, mmin) for e in grid]); j = np.nanargmax(ll); return grid[j]
def theory_alpha(dNdM, eps, mmin):
    Mcl = Mg * eps; s = (Mcl >= mmin) & (dNdM > 0)
    if s.sum() < 4: return np.nan
    w = dNdM[s] / eps * np.gradient(Mcl[s]); return 1 + w.sum() / np.sum(w * np.log(Mcl[s] / mmin))
models = [('PRFM', FB_tot, FB_per, '#c8412f', '-'),
          ('PRFM, no feedback limit', G_tot, G_per, GREY, '--'),
          ('PRFM, all crossings', H_tot, H_per, '#7b4ea3', ':')]
EPS = {name: fit_eps(tot, Mobs, a.mmin) for name, tot, _, _, _ in models}
for name in EPS: print(f'  eps({name}) = {EPS[name]:.3f}')
main = models[0]

# ============================================================ Fig 1
fig, ax = plt.subplots(1, 3, figsize=(7.4, 2.85))
for axi in ax: axi.set_box_aspect(1)
# (a) dN/dlogM per disc epoch
edges = np.logspace(2, 5.2, 14); cen = np.sqrt(edges[1:] * edges[:-1]); dlog = np.diff(np.log10(edges))
for e, col, lab in zip(EPOCHS[:2], ECOL, ELAB):
    so = (Tobs >= e[0]) & (Tobs < e[1]); hh, _ = np.histogram(Mobs[so], edges); ok = hh > 0
    Nab = np.sum(so & (Mobs >= a.mmin))
    ax[0].errorbar(cen[ok], hh[ok] / dlog[ok], yerr=np.sqrt(hh[ok]) / dlog[ok], fmt='o', ms=3.2, color=col, mfc=col, mec='none', ecolor=col, elinewidth=0.7, capsize=0, label=f'sim. {lab}', zorder=3)
    eps = EPS[main[0]]; th = main[2][e]; Mcl = Mg * eps; pdf = th / eps * Mcl * np.log(10)   # dN/dlogM
    s = Mcl >= a.mmin; nrm = Nab / max(np.trapezoid(th[s] / eps, Mcl[s]), 1e-300)
    ax[0].plot(Mcl, pdf * nrm, color=col, lw=1.3, label=f'PRFM {lab}', zorder=2)
xx = np.array([3e2, 3e4]); ax[0].plot(xx, 30 * (xx / 1e3) ** -1, color='k', lw=0.6, ls=(0, (1, 2)), label=r'$M^{-2}$', zorder=1)
ax[0].axvline(a.mmin, color='#999999', lw=0.5, ls=':')
ax[0].set(xscale='log', yscale='log', xlim=(1e2, 2e5), ylim=(0.3, 1e4), xlabel=r'$M_{\rm cl}$ [M$_\odot$]', ylabel=r'${\rm d}N/{\rm d}\log M$')
ax[0].legend(loc='upper right', handlelength=1.5, borderaxespad=0.4, labelspacing=0.3)
# (b) curvature
# closed form of Section 'origin of the slope': dn/dM ~ M^{-2+Delta} [ln(Mh/M)]^{-3/2} exp(-B0^2(1+p)/(2 s ln(Mh/M)))
# evaluated for the median disc-phase patch state; once with the layer mass Mh (pure barrier) and once with the
# feedback-limited mass M_fb = Mh exp[-(1+p) L_fb], L_fb = (ln(rho_fb/rho0) - B0)/(2-p), where the constant barrier takes over
st_med = {}
vals = {k: [] for k in ['mach', 'Q', 'kt', 'rho0', 'h']}
for i in idx:
    stt = state(i); vt = stt['vt_kms'] * X.KMS; cs_ = a.cs * X.KMS
    vals['mach'].append(vt / cs_); vals['kt'].append(stt['kappa'] / stt['Omega']); vals['rho0'].append(stt['rho0']); vals['h'].append(stt['h'])
    vals['Q'].append(np.sqrt(cs_ ** 2 + vt ** 2) * stt['kappa'] / (np.pi * X.G * stt['Sigma']))
st_med = {k: float(np.median(v)) for k, v in vals.items()}
pp_ = a.p; b2 = 0.75
B0 = np.log(st_med['Q'] / (2 * st_med['kt'])); Mh = 4 * np.pi / 3 * st_med['rho0'] * st_med['Q'] / (2 * st_med['kt']) * st_med['h'] ** 3
Lfb = (np.log(rho_fb / st_med['rho0']) - B0) / (2 - pp_); Mfb = Mh * np.exp(-(1 + pp_) * Lfb)
def closed_form(Mcut):
    M = np.logspace(0, np.log10(Mcut) - 0.02, 600); L = np.log(Mcut / M) / (1 + pp_)
    sL = np.log(b2 * st_med['mach'] ** 2) - pp_ * L; sL = np.maximum(sL, 0.3)
    beta = (2 - pp_) / sL + 0.5; Delta = beta ** 2 * sL / (2 * (1 + pp_))
    dn = M ** (-2 + Delta) * np.log(Mcut / M) ** -1.5 * np.exp(-B0 ** 2 * (1 + pp_) / (2 * sL * np.log(Mcut / M)))
    return M, dn, Delta
def alpha_tab(M, dn, eps, mmin):
    Mcl = M * eps; sel_ = (Mcl >= mmin) & (dn > 0)
    if sel_.sum() < 4: return np.nan
    w = dn[sel_] * np.gradient(Mcl[sel_]); return 1 + w.sum() / np.sum(w * np.log(Mcl[sel_] / mmin))
print(f'closed form, median patch: Mach {st_med["mach"]:.1f}  Q {st_med["Q"]:.2f}  kappa/Omega {st_med["kt"]:.2f}  rho0 {st_med["rho0"]:.3f}  h {st_med["h"]:.0f} pc  B0 {B0:.2f}  '
      f'M_h {Mh:.2e}  L_fb {Lfb:.2f}  M_fb {Mfb:.2e}  Delta at M_min=300: {closed_form(Mfb)[2][np.argmin(np.abs(closed_form(Mfb)[0] * EPS[main[0]] - 300))]:.2f}')
mmins = np.logspace(2, 4, 13)
al_s = np.array([mle(Mobs, mm)[0] for mm in mmins]); ns = np.array([mle(Mobs, mm)[1] for mm in mmins])
ax[1].errorbar(mmins, al_s, (al_s - 1) / np.sqrt(np.maximum(ns, 1)), fmt='o', ms=3.2, color='k', ecolor='k', elinewidth=0.7, capsize=0, label='simulation', zorder=3)
for name, tot, per, col, ls in models:
    ax[1].plot(mmins, [theory_alpha(tot, EPS[name], mm) for mm in mmins], color=col, ls=ls, label=name, zorder=2)
Mc_, dn_, _ = closed_form(Mfb); ax[1].plot(mmins, [alpha_tab(Mc_, dn_, EPS[main[0]], mm) for mm in mmins], color='k', lw=0.9, ls='-', label=r'eq. (19), $M_h \to M_{\rm fb}$', zorder=2)
Mc_, dn_, _ = closed_form(Mh); ax[1].plot(mmins, [alpha_tab(Mc_, dn_, EPS['PRFM, no feedback limit'], mm) for mm in mmins], color='k', lw=0.9, ls='-.', label=r'eq. (19), $M_h$', zorder=2)
ax[1].axhline(2, color='#999999', lw=0.5, ls=':')
ax[1].set(xscale='log', xlim=(80, 1.3e4), ylim=(1.2, 3.05), xlabel=r'$M_{\min}$ [M$_\odot$]', ylabel=r'$\alpha\,(M > M_{\min})$')
ax[1].legend(loc='upper left', handlelength=1.8, borderaxespad=0.3, labelspacing=0.3)
# (c) environment: index in weight terciles of the disc-phase patches, measured vs predicted from the same patches
Wp = D['W'][idx]; qe = np.quantile(Wp, [0, 1/3, 2/3, 1]); qe[-1] *= 1.001
selc = (E['bound'] > 0) & (E['nucleus'] == 0) & cleanE & (np.mod(E['snap'], a.fitstep) == 0) & (E['t'] < C.T_MERGED) & (E['M'] >= a.mmin)
Wc, Mc_ = E['W'][selc], E['M'][selc]
print('(c) weight terciles of disc-phase patches: W range, N clusters, alpha sim, alpha PRFM')
for j in range(3):
    pi = [i for i in idx if qe[j] <= D['W'][i] < qe[j + 1] and i in PP]
    sc_ = (Wc >= qe[j]) & (Wc < qe[j + 1]); al, n = mle(Mc_[sc_], a.mmin)
    th = sum(PP[i] for i in pi); at = theory_alpha(th, EPS[main[0]], a.mmin); Wm = np.median(Wp[(Wp >= qe[j]) & (Wp < qe[j + 1])])
    ax[2].errorbar(Wm, al, (al - 1) / np.sqrt(max(n, 1)), fmt='o', ms=4.5, color='k', ecolor='k', elinewidth=0.8, capsize=0, zorder=3, label='simulation' if j == 0 else None)
    ax[2].plot(Wm, at, marker='s', ms=5.5, color='#c8412f', mfc='none', mew=1.2, ls='none', zorder=4, label='PRFM' if j == 0 else None)
    print(f'  {qe[j]:.3g}-{qe[j+1]:.3g}: N={n}  alpha sim {al:.2f}+-{(al-1)/np.sqrt(max(n,1)):.2f}  PRFM {at:.2f}  (patches {len(pi)})')
ax[2].axhline(2, color='#999999', lw=0.5, ls=':')
ax[2].legend(loc='lower left', borderaxespad=0.4, labelspacing=0.3)
ax[2].set(xscale='log', xlim=(1.5e3, 1.5e5), ylim=(1.2, 2.3), xlabel=r'$\mathcal{W}/k_B$ of the patches [K cm$^{-3}$]', ylabel=r'$\alpha\,(M > 300\,{\rm M}_\odot)$')
for axi, lab in zip(ax, 'abc'): axi.text(0.0, 1.02, f'({lab})', transform=axi.transAxes, va='bottom', ha='left', fontsize=9)
fig.tight_layout(w_pad=1.0)
for ext in ['png', 'pdf']: fig.savefig(f'{FIG}/letter_fig1_theory.{ext}')
print('fig1 written')

# ============================================================ Fig 2
fig, ax = plt.subplots(1, 3, figsize=(7.4, 2.85))
for axi in ax: axi.set_box_aspect(1)
pp = cleanD & (D['M_young'] > a.mmin_young) & (D['M_max'] > 0) & (D['Sigma_gas'] > 1) & np.isfinite(D['W'])
My, Mx, W, tp = D['M_young'][pp], D['M_max'][pp], D['W'][pp], D['t'][pp]
o = np.argsort(W)
sc = ax[0].scatter(My[o], Mx[o], c=np.log10(W[o]), s=5, cmap='cividis', vmin=3.3, vmax=5.6, alpha=0.85, linewidths=0, rasterized=True)
xx = np.logspace(2.6, 6.1, 10); ax[0].plot(xx, xx, color='k', lw=0.6, ls=(0, (1, 2)), label=r'$M_{\max}=M_{\rm young}$')
ax[0].plot(xx, 0.5 * xx, color='k', lw=0.9, ls='--', label=r'$M_{\max}=0.5\,M_{\rm young}$')
ax[0].set(xscale='log', yscale='log', xlim=(4e2, 1.5e6), ylim=(1e2, 4e5), xlabel=r'$M_{\rm young}$ per patch [M$_\odot$]', ylabel=r'$M_{\max}$ per patch [M$_\odot$]')
ax[0].legend(loc='upper left', handlelength=1.8, borderaxespad=0.3, labelspacing=0.3)
cax = ax[0].inset_axes([0.52, 0.12, 0.42, 0.035]); cb = fig.colorbar(sc, cax=cax, orientation='horizontal'); cb.ax.xaxis.set_ticks_position('top'); cb.ax.xaxis.set_label_position('top'); cb.set_label(r'$\log\,\mathcal{W}/k_B$ [K cm$^{-3}$]', fontsize=7, labelpad=2); cb.ax.tick_params(labelsize=6.5, direction='in', length=2, pad=1.5); cb.outline.set_linewidth(0.5)
# (b) truncation vs time
def trunc_mle(Mv, mmin):
    Mv = Mv[Mv >= mmin]; n = len(Mv)
    al_pl = 1 + n / np.sum(np.log(Mv / mmin)); ll_pl = n * np.log(al_pl - 1) + n * (al_pl - 1) * np.log(mmin) - al_pl * np.sum(np.log(Mv))
    best = (al_pl, np.inf, ll_pl); x = np.logspace(np.log10(mmin), 8, 600)
    for lMc in np.linspace(np.log10(mmin) + 0.3, 7.5, 50):
        Mc = 10 ** lMc
        for al in np.linspace(1.0, 3.0, 33):
            Z = np.trapezoid(x ** -al * np.exp(-x / Mc), x); ll = np.sum(-al * np.log(Mv) - Mv / Mc) - n * np.log(Z)
            if ll > best[2]: best = (al, Mc, ll)
    return best[1], best[2] - ll_pl
selB = (E['bound'] > 0) & (E['nucleus'] == 0) & cleanE & (E['M'] >= a.mmin)
Mb, Tb = E['M'][selB], E['t'][selB]
sfD2 = cleanD & (D['M_young'] > a.mmin_young) & (D['Sigma_gas'] > 1) & np.isfinite(D['W'])
tb, mx, mc, wm, nb = [], [], [], [], []
for lo in np.arange(0, 232, 25):
    hi = lo + 25; so = (Tb >= lo) & (Tb < hi); sp = sfD2 & (D['t'] >= lo) & (D['t'] < hi)
    if so.sum() < 8 or sp.sum() == 0: continue
    Mc, dll = trunc_mle(Mb[so], a.mmin) if so.sum() >= 10 else (np.nan, 0)
    tb.append(0.5 * (lo + hi)); mx.append(Mb[so].max()); mc.append(Mc if dll > 2 else np.nan); wm.append(np.median(D['W'][sp])); nb.append(so.sum())
tb, mx, mc, wm = map(np.array, (tb, mx, mc, wm))
for m, col in [(tb < C.T_MERGED, '#1f5f9e'), (tb >= C.T_MERGED, '#c8412f')]:
    ax[1].plot(wm[m], mx[m], 'o', color=col, ms=4.5)
    ax[1].plot(wm[m], mc[m], 's', color=col, mfc='none', mew=1.2, ms=5)
hb = [Line2D([], [], marker='o', color='k', ls='none', ms=4.5, label='largest cluster'),
      Line2D([], [], marker='s', color='k', mfc='none', mew=1.2, ls='none', ms=5, label=r'Schechter $M_c$')]
ok_ = np.isfinite(mc)
print(f'(b) per 25 Myr bin: Spearman(M_max, W) = {spearmanr(mx, wm)[0]:.2f}, Spearman(M_c, W) = {spearmanr(mc[ok_], wm[ok_])[0]:.2f} (n={len(mx)}, {ok_.sum()})')
ax[1].set(xscale='log', yscale='log', xlim=(2e3, 5e5), ylim=(3e2, 2e6), xlabel=r'median $\mathcal{W}/k_B$ per 25 Myr [K cm$^{-3}$]', ylabel=r'$M$ [M$_\odot$]')
ax[1].legend(handles=hb, loc='upper left', handlelength=1.2, borderaxespad=0.4, labelspacing=0.3)
# (c) bound fraction vs Q per patch: the mechanism
sg = cleanD & (D['M_young'] > a.mmin_young) & (D['Sigma_gas'] > 1) & np.isfinite(D['Q']) & (D['Q'] > 0)
Qp, Gp, tq, Myq = D['Q'][sg], np.clip(D['Gamma'][sg], 0, 1), D['t'][sg], D['M_young'][sg]
qb = np.logspace(-1, 1.7, 10)
for m, col, lab in [(tq < C.T_MERGED, '#1f5f9e', 'disc phase'), (tq >= C.T_MERGED, '#c8412f', 'after coalescence')]:
    ax[2].scatter(Qp[m], Gp[m], s=4, color=col, alpha=0.25, linewidths=0, rasterized=True)
    qc, gm, glo, ghi = [], [], [], []
    for lo, hi in zip(qb[:-1], qb[1:]):
        s = m & (Qp >= lo) & (Qp < hi)
        if s.sum() >= 15:
            qc.append(np.sqrt(lo * hi)); w = Myq[s]; g = Gp[s]
            gm.append(np.sum(w * g) / w.sum()); glo.append(np.percentile(g, 25)); ghi.append(np.percentile(g, 75))
    ax[2].plot(qc, gm, 'o-', color=col, ms=3.5, label=lab)
print(f'  Gamma vs Q spearman (all star-forming patches): {spearmanr(Gp, Qp)[0]:+.2f}')
ax[2].axvline(1, color='#999999', lw=0.6, ls=':')
ax[2].set(xscale='log', xlim=(0.1, 50), ylim=(0, 1.5), yticks=[0, 0.2, 0.4, 0.6, 0.8, 1.0], xlabel=r'Toomre $Q$ of the patch', ylabel=r'$\Gamma$ = bound mass / $M_{\rm young}$')
ax[2].legend(loc='upper left', handlelength=1.8, borderaxespad=0.4, labelspacing=0.3, title='mass-weighted mean', title_fontsize=7)
for axi, lab in zip(ax, 'abc'): axi.text(0.0, 1.02, f'({lab})', transform=axi.transAxes, va='bottom', ha='left', fontsize=9)
fig.tight_layout(w_pad=1.0)
for ext in ['png', 'pdf']: fig.savefig(f'{FIG}/letter_fig2_reservoir.{ext}')
print('fig2 written')
