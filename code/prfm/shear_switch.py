"""Is the switch the shear?  Attach the patch-scale velocity-gradient fields (shear norm, div_v, vort; particle-scale gradients averaged
over the patch, all snapshots) to patches_clusters rows, then:
(1) per phase: partial Spearman of the burst flag with log shear, log |div_v| (convergence), log vort at fixed log rho_mid (residuals);
(2) logistic P(burst | rho) vs P(burst | rho, shear): deviance gained;
(3) among active patches: the conversion ratio M_young / (eta PRFM burst) vs shear (Spearman per phase);
(4) per 25 Myr interval: median shear, convergence and vorticity of the clean patches, and of the bursting ones - is 160-185 high-shear?"""
import sys, os, numpy as np, h5py
sys.path.insert(0, '/raven/u/uli/dwarf_merger'); sys.path.insert(0, '/raven/u/uli/dwarf_merger/prfm'); import common as C, ok22
from scipy.stats import spearmanr
from scipy.optimize import minimize
D = dict(np.load(f'{C.DATADIR}/prfm/patches_clusters.npz')); WL = np.load(f'{C.DATADIR}/prfm/cluster_link.npz')['WL']
snap = D['snap'].astype(int); fr = D['frame'].astype(int); n = len(snap)
SH = np.full(n, np.nan); DV = np.full(n, np.nan); VO = np.full(n, np.nan); VB = np.full(n, np.nan)
for k in np.unique(snap):
    fn = f'{C.PRFM_DIR}/patch_{k:03d}.h5'
    if not os.path.exists(fn): continue
    with h5py.File(fn, 'r') as f:
        for code, name in ((65, 'frame_A'), (66, 'frame_B'), (77, 'frame_M')):
            rows = np.flatnonzero((snap == k) & (fr == code))
            if len(rows) == 0 or name not in f: continue
            g = f[name]; SH[rows] = g['shear'][:].ravel(); DV[rows] = g['div_v'][:].ravel(); VO[rows] = g['vort'][:].ravel(); VB[rows] = g['vbar_z'][:].ravel()
clean = (D['f_intruder'] < 0.1) | C.is_merged(D['sep'], D['t'])
base = clean & (D['Sigma_gas'] > 1) & np.isfinite(WL) & (WL > 0) & (D['t'] > 10) & (D['rho_mid_2p'] > 0) & np.isfinite(SH) & (SH > 0)
burst = base & (D['M_young'] > 500); lr = np.log10(D['rho_mid_2p']); ls = np.log10(SH); lv = np.log10(np.maximum(VO, 1e-6)); conv = -DV
EP = [*C.PHASES]
print(f'{base.sum()} patches with shear; shear units as stored (particle-scale gradient norm, patch mean); median shear {np.nanmedian(SH[base]):.3g}, div_v {np.nanmedian(DV[base]):.3g}')
def resid(y, x):
    A = np.vstack([np.ones(len(x)), x]).T; return y - A @ np.linalg.lstsq(A, y, rcond=None)[0]
print('\n(1) partial Spearman of the burst flag with each velocity-gradient quantity at fixed log rho_mid (both residualised on log rho)')
print('phase                 n     shear    |div|(conv>0 part)  vort   vbar_z(|.|)   | raw Spearman(burst, rho)')
for t0, t1, lab in EP:
    m = base & (D['t'] >= t0) & (D['t'] < t1); y = burst[m].astype(float); ry = resid(y, lr[m])
    out = []
    for x in (ls[m], np.log10(np.maximum(conv[m], 1e-3)), lv[m], np.log10(np.maximum(np.abs(VB[m]), 1e-3))):
        out.append(spearmanr(ry, resid(x, lr[m]))[0])
    print(f'{lab:20s} {m.sum():5d}   ' + '   '.join(f'{v:+.2f}  ' for v in out) + f'      | {spearmanr(y, lr[m])[0]:+.2f}')
print('\n(2) logistic deviance explained: rho alone vs rho + shear vs rho + convergence (per phase)')
def fit(X, y):
    def nll(p):
        z = p[0] + X @ p[1:]; pr = np.clip(1 / (1 + np.exp(-z)), 1e-9, 1 - 1e-9); return -np.sum(y * np.log(pr) + (1 - y) * np.log(1 - pr))
    r = minimize(nll, np.zeros(X.shape[1] + 1), method='BFGS'); return r.fun, r.x
for t0, t1, lab in EP:
    m = base & (D['t'] >= t0) & (D['t'] < t1); y = burst[m].astype(float); p0 = y.mean(); nll0 = -m.sum() * (p0 * np.log(p0) + (1 - p0) * np.log(1 - p0))
    f1, c1 = fit(np.column_stack([lr[m]]), y); f2, c2 = fit(np.column_stack([lr[m], ls[m]]), y); f3, c3 = fit(np.column_stack([lr[m], np.log10(np.maximum(conv[m], 1e-3))]), y)
    print(f'{lab:20s} rho: {1-f1/nll0:.3f}   rho+shear: {1-f2/nll0:.3f} (shear coef {c2[2]:+.2f})   rho+conv: {1-f3/nll0:.3f} (conv coef {c3[2]:+.2f})')
print('\n(3) active patches: log[M_young / (0.4 PRFM burst)] vs shear, convergence, vorticity (Spearman) and vs rho')
prfm = lambda W: 0.4 * ok22.sfr_of_W(W) * 0.25 * 1e7
for t0, t1, lab in EP:
    m = burst & (D['t'] >= t0) & (D['t'] < t1); q = np.log10(D['M_young'][m] / prfm(WL[m]))
    print(f'{lab:20s} n={m.sum():4d}  shear {spearmanr(q, ls[m])[0]:+.2f}   conv {spearmanr(q, conv[m])[0]:+.2f}   vort {spearmanr(q, lv[m])[0]:+.2f}   rho {spearmanr(q, lr[m])[0]:+.2f}   median ratio {10**np.median(q):.2f}')
print('\n(4) per 25 Myr interval: median over clean patches of shear, convergence (-div v), vorticity, |vbar_z|; and the same over bursting patches; conversion ratio')
print('bin    shear(all)  shear(burst)   conv(all)   conv(burst)   vort(all)   |vbar_z|(all)   f_burst   M_young/(PRFM) med')
for t0 in np.arange(10, 210, 25):
    m = base & (D['t'] >= t0) & (D['t'] < t0 + 25); mb = burst & m
    q = np.log10(D['M_young'][mb] / prfm(WL[mb])) if mb.sum() else np.array([np.nan])
    print(f'{t0:3.0f}   {np.median(SH[m]):8.3g}   {np.median(SH[mb]) if mb.any() else np.nan:8.3g}     {np.median(conv[m]):+8.3g}   {np.median(conv[mb]) if mb.any() else np.nan:+8.3g}    {np.median(VO[m]):8.3g}   {np.median(np.abs(VB[m])):7.2f}       {mb.sum()/m.sum():.3f}     {10**np.nanmedian(q):.2f}')
