"""Link between the PRFM layer state and the star clusters, using the exact layer weight (patch_XXX_z05.h5) attached to
patches_clusters.npz (rows are in grid order per snapshot and frame; frame codes 65=A, 66=B, 77=M).
(1) per patch: largest bound cluster M_max and young mass M_young vs layer weight W_L (and full W) - Spearman;
(2) per 25 Myr bin: largest bound cluster, Schechter Mc (from clusters_env.npz), vs the layer weight of the star-forming patches
    (median and 90th percentile of W_L over patches with M_young > 500) - log-log slopes;
(3) the same vs P_tot (layer) and vs Sigma_gas, to see which patch quantity tracks the truncation best.
out: prfm/cluster_link.npz, figs/cluster_link.png (+ ~/out)"""
import sys, os, numpy as np, h5py
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common as C
from scipy.stats import spearmanr
from scipy.optimize import minimize
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
D = dict(np.load(f'{C.DATADIR}/prfm/patches_clusters.npz')); E = np.load(f'{C.DATADIR}/clusters/clusters_env.npz')
snap = D['snap'].astype(int); fr = D['frame'].astype(int); n = len(snap)
WL = np.full(n, np.nan); PL = np.full(n, np.nan); SL = np.full(n, np.nan)
for k in np.unique(snap):
    fn = f'{C.PRFM_DIR}/patch_{k:03d}_z05.h5'
    if not os.path.exists(fn): continue
    with h5py.File(fn, 'r') as f:
        for code, name in ((65, 'frame_A'), (66, 'frame_B'), (77, 'frame_M')):
            rows = np.flatnonzero((snap == k) & (fr == code))
            if len(rows) == 0 or name not in f: continue
            g = f[name]; WL[rows] = g['W_2p'][:].ravel(); PL[rows] = g['Ptot_2p'][:].ravel(); SL[rows] = g['Sigma_gas'][:].ravel()
clean = (D['f_intruder'] < 0.1) | C.is_merged(D['sep'], D['t'])
base = clean & (D['Sigma_gas'] > 1) & np.isfinite(WL) & (WL > 0) & (D['t'] > 10)
print(f'{base.sum()} clean patches with layer weight ({np.isfinite(WL).sum()} of {n} rows matched)')
EP = [*C.PHASES]
# ---- (1) per patch
print('\n(1) per patch (M_young > 500, bound cluster present): Spearman of log M_max and log M_young with the weight')
print('phase                 n     M_max~W_full  M_max~W_layer  M_max~P_layer  M_max~Sigma | M_young~W_full  M_young~W_layer  M_young~P_layer  M_young~Sigma')
for t0, t1, lab in EP:
    m = base & (D['t'] >= t0) & (D['t'] < t1) & (D['M_young'] > 500) & (D['M_max'] > 0)
    if m.sum() < 20: print(lab, 'too few'); continue
    r = lambda y, x: spearmanr(np.log10(y[m]), np.log10(x[m]))[0]
    print(f'{lab:20s} {m.sum():5d}    ' + '  '.join(f'{r(D["M_max"], x):+.2f}        ' for x in (D['W'], WL, PL, D['Sigma_gas'])) + ' | ' + '  '.join(f'{r(D["M_young"], x):+.2f}          ' for x in (D['W'], WL, PL, D['Sigma_gas'])))
# ---- (2) per 25 Myr bin: truncation vs layer weight of the SF patches
cleanE = (E['f_intruder'] < 0.1) | C.is_merged(E['sep'], E['t']); cb = cleanE & (E['bound'] > 0) & (E['nucleus'] == 0)
def schechter(M, mmin=300.):
    M = M[M >= mmin]; nn = len(M); x = np.logspace(np.log10(mmin), 8, 2000)
    def nll(p):
        a, lMc = p; Mc = 10 ** lMc
        if a < 0 or a > 3.5 or lMc < np.log10(mmin) or lMc > 8: return 1e30
        Z = np.trapezoid(x ** -a * np.exp(-x / Mc), x); return -(np.sum(-a * np.log(M) - M / Mc) - nn * np.log(Z))
    best = min((minimize(nll, [a0, l0], method='Nelder-Mead') for a0 in (1.3, 2.0) for l0 in (3.5, 4.5, 5.5)), key=lambda r: r.fun)
    return best.x[0], 10 ** best.x[1], nn
edges = np.arange(10, 232, 25); rows = []
print('\n(2) 25 Myr bins: largest bound cluster and Schechter Mc vs the layer weight of star-forming patches (M_young > 500)')
print('t0-t1     n_cl   M_max     alpha   Mc      | SF patches: n   W_L median   W_L 90%   P_L median   Sigma median   W_full median')
for t0, t1 in zip(edges[:-1], edges[1:]):
    mc = cb & (E['t'] >= t0) & (E['t'] < t1); mp = base & (D['t'] >= t0) & (D['t'] < t1) & (D['M_young'] > 500)
    if mc.sum() < 30 or mp.sum() < 10: continue
    a, Mc, nn = schechter(E['M'][mc]); Mmax = E['M'][mc].max()
    rows.append((t0, t1, nn, Mmax, a, Mc, mp.sum(), np.median(WL[mp]), np.percentile(WL[mp], 90), np.median(PL[mp]), np.median(D['Sigma_gas'][mp]), np.median(D['W'][mp])))
    print(f'{t0:3.0f}-{t1:3.0f}  {nn:5d}  {Mmax:8.2e}  {a:5.2f}  {Mc:8.2e} |           {mp.sum():4d}   {np.median(WL[mp]):9.2e}  {np.percentile(WL[mp],90):9.2e}  {np.median(PL[mp]):9.2e}   {np.median(D["Sigma_gas"][mp]):6.1f}       {np.median(D["W"][mp]):9.2e}')
R = np.array(rows); np.savez(f'{C.DATADIR}/prfm/cluster_link.npz', bins=R, WL=WL, PL=PL)
print('\nlog-log slopes and Spearman over the bins:')
for yname, y in (('M_max', R[:, 3]), ('Mc', R[:, 5])):
    for xname, x in (('W_L median', R[:, 7]), ('W_L 90%', R[:, 8]), ('P_L median', R[:, 9]), ('Sigma median', R[:, 10]), ('W_full median', R[:, 11])):
        sl = np.polyfit(np.log10(x), np.log10(y), 1)[0]; rs = spearmanr(x, y)[0]
        print(f'  {yname:6s} vs {xname:14s}: slope {sl:+.2f}  rho_s {rs:+.2f}')
ink = '#1f2937'; blue = '#2f6fb2'; light = '#9dbfe3'
fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
for a in ax: a.grid(alpha=0.2); a.spines[['top', 'right']].set_visible(False); a.set_xscale('log'); a.set_yscale('log')
ax[0].scatter(R[:, 7], R[:, 3], c=blue, s=40, label='largest bound cluster'); ax[0].scatter(R[:, 7], R[:, 5], c=light, s=40, marker='s', label='Schechter $M_c$')
for i in range(len(R)): ax[0].annotate(f'{R[i,0]:.0f}', (R[i, 7], R[i, 3]), fontsize=7, color=ink, xytext=(3, 3), textcoords='offset points')
ax[0].set_xlabel('median layer weight of star-forming patches [K cm$^{-3}$]'); ax[0].set_ylabel('cluster mass [M$_\\odot$]'); ax[0].legend(frameon=False, fontsize=8); ax[0].set_title('Truncation follows the layer weight (25 Myr bins)', loc='left', color=ink, fontsize=10)
m = base & (D['M_young'] > 500) & (D['M_max'] > 0)
ax[1].scatter(WL[m], D['M_max'][m], s=4, c=light, alpha=0.5, label='patches: largest cluster'); ax[1].set_xlabel('layer weight of the patch [K cm$^{-3}$]'); ax[1].set_ylabel('largest bound cluster in the patch [M$_\\odot$]')
xb = np.logspace(3, 6, 13); med = [np.median(D['M_max'][m & (WL >= lo) & (WL < hi)]) if (m & (WL >= lo) & (WL < hi)).sum() >= 10 else np.nan for lo, hi in zip(xb[:-1], xb[1:])]
ax[1].plot(np.sqrt(xb[:-1] * xb[1:]), med, color=blue, lw=2, label='median per weight bin'); ax[1].legend(frameon=False, fontsize=8); ax[1].set_title('Per patch', loc='left', color=ink, fontsize=10)
fn = f'{C.DATADIR}/prfm/figs/cluster_link.png'; fig.savefig(fn, dpi=140, bbox_inches='tight'); os.system(f'cp {fn} ~/out/'); print('fig', fn)
