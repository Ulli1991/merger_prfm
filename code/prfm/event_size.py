"""The top of the cluster mass function as a reservoir effect.

(a) per patch: most massive bound cluster vs young stellar mass formed in the same
    10 Myr in the same 0.5 kpc column  ->  one dominant cluster per event, M_max ~ f x M_young
(b) vs time: largest bound cluster and Schechter truncation mass per time bin, with the
    median weight of the star-forming patches  ->  the truncation follows the weight
(c) bound fraction Gamma (bound cluster mass / stars formed) vs time, both age windows

usage: python event_size.py [--dt 25] [--mmin 300] [--fintr 0.1] [--mmin_young 500]
out:   figs/cl_fig8_event_size.png, event_size.npz
"""
import sys, os, argparse
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.stats import spearmanr
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import common as C

ap = argparse.ArgumentParser()
ap.add_argument('--dt', type=float, default=25.0)
ap.add_argument('--mmin', type=float, default=300.0)
ap.add_argument('--fintr', type=float, default=0.1)
ap.add_argument('--mmin_young', type=float, default=500.0)
ap.add_argument('--nmin_bin', type=int, default=8)
a = ap.parse_args()
FIG = f'{C.DATADIR}/prfm/figs'; os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': True, 'grid.alpha': 0.25, 'legend.frameon': False, 'figure.dpi': 130})

D = np.load(f'{C.DATADIR}/prfm/patches_clusters.npz')
E = np.load(f'{C.DATADIR}/clusters/clusters_env.npz')
G = np.load(f'{C.DATADIR}/clusters/gamma_time.npz')['data']      # k, t, [Mf, Mg, Mb] young, [Mf, Mg, Mb] old

cleanD = (D['f_intruder'] < a.fintr) | C.is_merged(D['sep'], D['t'])
cleanE = (E['f_intruder'] < a.fintr) | C.is_merged(E['sep'], E['t'])

# ---------------- (a) per patch
pp = cleanD & (D['M_young'] > a.mmin_young) & (D['M_max'] > 0) & (D['Sigma_gas'] > 1) & np.isfinite(D['W'])
My, Mx, W, t = D['M_young'][pp], D['M_max'][pp], D['W'][pp], D['t'][pp]
ratio = Mx / My
rho_y, _ = spearmanr(np.log10(Mx), np.log10(My)); rho_w, _ = spearmanr(np.log10(Mx), np.log10(W))
A = np.vstack([np.ones(pp.sum()), np.log10(My)]).T; coef = np.linalg.lstsq(A, np.log10(Mx), rcond=None)[0]; res = np.log10(Mx) - A @ coef
rho_res_w, _ = spearmanr(res, np.log10(W))
print(f'(a) {pp.sum()} clean patches with M_young > {a.mmin_young:g} and a bound cluster (all snapshots)')
print(f'    M_max / M_young: median {np.median(ratio):.2f}, 16-84% {np.percentile(ratio,16):.2f}-{np.percentile(ratio,84):.2f}')
print(f'    Spearman(log M_max, log M_young) = {rho_y:.2f};  Spearman(log M_max, log W) = {rho_w:.2f};  residual after M_young vs log W: {rho_res_w:+.2f}')
print(f'    fit: log M_max = {coef[0]:.2f} + {coef[1]:.2f} log M_young')
# and does M_young follow W (the PRFM link)?
sf = cleanD & (D['M_young'] > a.mmin_young) & (D['Sigma_gas'] > 1) & np.isfinite(D['W'])
rho_sw, _ = spearmanr(np.log10(D['M_young'][sf]), np.log10(D['W'][sf]))
print(f'    Spearman(log M_young, log W) over {sf.sum()} star-forming clean patches = {rho_sw:.2f}')

# ---------------- (b) vs time
def trunc_mle(Mv, mmin):
    Mv = Mv[Mv >= mmin]; n = len(Mv)
    al_pl = 1 + n / np.sum(np.log(Mv / mmin)); ll_pl = n * np.log(al_pl - 1) + n * (al_pl - 1) * np.log(mmin) - al_pl * np.sum(np.log(Mv))
    best = (al_pl, np.inf, ll_pl); x = np.logspace(np.log10(mmin), 8, 800)
    for lMc in np.linspace(np.log10(mmin) + 0.3, 7.5, 60):
        Mc = 10 ** lMc
        for al in np.linspace(1.0, 3.0, 41):
            Z = np.trapezoid(x ** -al * np.exp(-x / Mc), x); ll = np.sum(-al * np.log(Mv) - Mv / Mc) - n * np.log(Z)
            if ll > best[2]: best = (al, Mc, ll)
    return al_pl, best[0], best[1], best[2] - ll_pl
sel = (E['bound'] > 0) & (E['nucleus'] == 0) & cleanE & (E['M'] >= a.mmin)
Mo, To = E['M'][sel], E['t'][sel]
edges = np.arange(0, C.NSNAP * C.SNAP_DT + a.dt, a.dt)
rows = []
print(f'(b) per {a.dt:g} Myr bin: N, M_max, Schechter (alpha_tr, Mc, dlnL vs power law), power-law alpha, median W of SF patches')
for lo, hi in zip(edges[:-1], edges[1:]):
    so = (To >= lo) & (To < hi); sp = sf & (D['t'] >= lo) & (D['t'] < hi)
    if so.sum() < a.nmin_bin or sp.sum() == 0: continue
    al_pl, al_tr, Mc, dll = trunc_mle(Mo[so], a.mmin) if so.sum() >= 10 else (np.nan, np.nan, np.nan, np.nan)
    rows.append((0.5 * (lo + hi), so.sum(), Mo[so].max(), al_pl, al_tr, Mc, dll, np.median(D['W'][sp]), np.median(D['M_young'][sp]), sp.sum()))
    print(f'  {lo:4.0f}-{hi:4.0f}: N={so.sum():4d}  M_max={Mo[so].max():8.2e}  alpha_tr={al_tr:.2f} Mc={Mc:8.2e} (dlnL={dll:5.1f})  alpha_pl={al_pl:.2f}  W={np.median(D["W"][sp]):8.3g}  M_young={np.median(D["M_young"][sp]):8.3g}  n_patch={sp.sum()}')
rows = np.array(rows)
ok = np.isfinite(rows[:, 5]) & (rows[:, 6] > 2)
if ok.sum() > 3:
    print(f'    Spearman(Mc, W) over bins with a significant truncation: {spearmanr(rows[ok,5], rows[ok,7])[0]:.2f};  Spearman(M_max, W) all bins: {spearmanr(rows[:,2], rows[:,7])[0]:.2f}')
np.savez(f'{C.DATADIR}/prfm/event_size.npz', rows=rows, My=My, Mx=Mx, W=W, t=t)

# ---------------- (c) bound fraction
tg = G[:, 1]; Gy = G[:, 4] / np.maximum(G[:, 2], 1); Go = G[:, 7] / np.maximum(G[:, 5], 1)
print('(c) bound fraction (mass in bound clusters / stars formed in the window), per 50 Myr:')
for lo, hi in [(0, 50), (50, 100), (100, 150), (150, 232)]:
    e = (tg >= lo) & (tg < hi)
    print(f'  {lo:3d}-{hi:3d} Myr: young(0-10) {G[e,4].sum()/max(G[e,2].sum(),1):.2f}   old(20-50) {G[e,7].sum()/max(G[e,5].sum(),1):.2f}')

# ---------------- figure
fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
sc = ax[0].scatter(My, Mx, c=np.log10(W), s=10, cmap='viridis', alpha=0.8)
xx = np.logspace(2.5, 6, 10); ax[0].plot(xx, xx, 'k:', lw=1, label=r'$M_{\max}=M_{\rm young}$'); ax[0].plot(xx, 0.5 * xx, 'k--', lw=1, label=r'$M_{\max}=0.5\,M_{\rm young}$')
ax[0].set(xscale='log', yscale='log', xlabel=r'$M_{\rm young}$ (age < 10 Myr, per 0.5 kpc patch) [M$_\odot$]', ylabel=r'most massive bound cluster $M_{\max}$ [M$_\odot$]', xlim=(3e2, 1e6), ylim=(1e2, 3e5))
ax[0].legend(fontsize=8, loc='upper left'); cb = fig.colorbar(sc, ax=ax[0], pad=0.01); cb.set_label(r'$\log\,\mathcal{W}/k_B$ [K cm$^{-3}$]')
ax[0].set_title(f'one dominant cluster per event: median ratio {np.median(ratio):.2f}, $\\rho_s$={rho_y:.2f}', fontsize=9)
tt = rows[:, 0]
ax[1].plot(tt, rows[:, 2], 'ko-', label=r'largest bound cluster $M_{\max}$')
ax[1].plot(tt[ok], rows[ok, 5], 'rs--', label=r'Schechter $M_c$ (bins with $\Delta\ln L>2$)')
ax[1].set(yscale='log', xlabel='t [Myr]', ylabel=r'M [M$_\odot$]'); ax[1].legend(fontsize=8, loc='upper left')
ax1b = ax[1].twinx(); ax1b.plot(tt, rows[:, 7], 'b-^', alpha=0.7, label=r'median $\mathcal{W}$ of star-forming patches'); ax1b.set(yscale='log', ylabel=r'$\mathcal{W}/k_B$ [K cm$^{-3}$]'); ax1b.grid(False); ax1b.legend(fontsize=8, loc='lower right')
for x_, n_ in zip(tt, rows[:, 1]): ax[1].annotate(f'{int(n_)}', (x_, rows[:, 2].min() * 0.5), fontsize=6.5, ha='center')
ax[1].set_title('the truncation follows the weight through the merger', fontsize=9)
ax[2].plot(tg, Gy, '-', color='tab:red', lw=1.2, label='stars younger than 10 Myr')
ax[2].plot(tg, Go, '-', color='tab:blue', lw=1.2, label='stars 20-50 Myr old')
ax[2].set(xlabel='t [Myr]', ylabel=r'$\Gamma$ = bound cluster mass / stars formed', ylim=(0, 1.05)); ax[2].legend(fontsize=8, loc='lower left')
ax[2].axvline(C.T_MERGED, color='grey', ls=':'); ax[2].annotate('first coalescence', (C.T_MERGED, 1.0), fontsize=7, ha='right', va='top', rotation=90)
ax[2].set_title('bound fraction collapses after coalescence', fontsize=9)
fig.tight_layout(); fig.savefig(f'{FIG}/cl_fig8_event_size.png'); print('figure:', f'{FIG}/cl_fig8_event_size.png')
