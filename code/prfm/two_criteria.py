"""Two-criterion model per 25 Myr interval, patch by patch:
  (1) P(burst | rho_mid,2p): logistic in log rho fitted on ALL clean patches of ALL phases (one threshold, no phase information);
  (2) burst size if it bursts: M_young = Sigma_SFR(W_L) A tau with the OK22 PRFM law (tau = 10 Myr, A = 0.25 kpc^2), x eta where eta is
      the single normalisation of the observed burst masses to the PRFM law over the active patches (printed);
  cluster: M_max = 0.5 x largest burst.
Predicted vs observed per interval: number of bursting patches, total M_young, expected largest burst (Monte Carlo over the bursting
draws), M_max; the quiescent interval 160-185 is the test.  Also the same with a W-threshold (1) to show the difference."""
import sys, numpy as np
sys.path.insert(0, '/raven/u/uli/dwarf_merger'); sys.path.insert(0, '/raven/u/uli/dwarf_merger/prfm'); import common as C, ok22
from scipy.optimize import minimize
rng = np.random.default_rng(3)
D = dict(np.load(f'{C.DATADIR}/prfm/patches_clusters.npz')); WL = np.load(f'{C.DATADIR}/prfm/cluster_link.npz')['WL']; E = np.load(f'{C.DATADIR}/clusters/clusters_env.npz')
clean = (D['f_intruder'] < 0.1) | C.is_merged(D['sep'], D['t']); base = clean & (D['Sigma_gas'] > 1) & np.isfinite(WL) & (WL > 0) & (D['t'] > 10) & (D['rho_mid_2p'] > 0)
cleanE = (E['f_intruder'] < 0.1) | C.is_merged(E['sep'], E['t']); cb = cleanE & (E['bound'] > 0) & (E['nucleus'] == 0)
burst = base & (D['M_young'] > 500); lr = np.log10(D['rho_mid_2p']); lW = np.log10(np.where(WL > 0, WL, 1))
def logistic_fit(x, y):
    def nll(p):
        pr = np.clip(1 / (1 + np.exp(-p[0] * (x - p[1]))), 1e-9, 1 - 1e-9); return -np.sum(y * np.log(pr) + (1 - y) * np.log(1 - pr))
    return min((minimize(nll, [k0, np.median(x)], method='Nelder-Mead') for k0 in (1, 4)), key=lambda r: r.fun).x
kr, r50 = logistic_fit(lr[base], burst[base].astype(float)); kw, w50 = logistic_fit(lW[base], burst[base].astype(float))
print(f'global logistic: P(burst|rho) k={kr:.2f} log rho_50={r50:.2f};  P(burst|W) k={kw:.2f} log W_50={w50:.2f}')
prfm = lambda W: ok22.sfr_of_W(W) * 0.25 * 1e7      # Msun per 10 Myr per 0.5 kpc patch
eta = 10 ** np.median(np.log10(D['M_young'][burst] / prfm(WL[burst]))); print(f'normalisation eta = median M_young / PRFM burst over active patches = {eta:.2f} (scatter {np.std(np.log10(D["M_young"][burst]/prfm(WL[burst]))):.2f} dex)')
print('\nbin   n_patch | bursts obs  pred(rho)  pred(W) | M_young obs   pred(rho)   pred(W) | M_max obs   pred(rho) [16-84%]     pred(W)')
for t0 in np.arange(10, 210, 25):
    m = base & (D['t'] >= t0) & (D['t'] < t0 + 25); idx = np.flatnonzero(m); ns = len(np.unique(D['snap'][m])) / 10.   # 25 Myr bin holds ~2.5 independent 10 Myr windows; patches of ALL snapshots are used -> scale by n_snap/10
    pr = 1 / (1 + np.exp(-kr * (lr[idx] - r50))); pw = 1 / (1 + np.exp(-kw * (lW[idx] - w50))); size = eta * prfm(WL[idx])
    mc = cb & (E['t'] >= t0) & (E['t'] < t0 + 25); Mmax_obs = E['M'][mc].max() if mc.any() else np.nan
    obs_n = burst[m].sum() / ns; obs_M = D['M_young'][m].sum() / ns
    res = {}
    for name, p in (('rho', pr), ('W', pw)):
        n_pred = p.sum() / ns; M_pred = (p * size).sum() / ns; mx = []
        for it in range(300):
            hit = rng.random(len(p)) < p; mx.append(0.5 * size[hit].max() if hit.any() else 0)
        res[name] = (n_pred, M_pred, np.median(mx), np.percentile(mx, 16), np.percentile(mx, 84))
    print(f'{t0:3.0f}   {len(idx)/ns:6.0f} |  {obs_n:6.1f}     {res["rho"][0]:6.1f}     {res["W"][0]:6.1f}  |  {obs_M:9.2e}   {res["rho"][1]:9.2e}   {res["W"][1]:9.2e} |  {Mmax_obs:8.1e}   {res["rho"][2]:8.1e} [{res["rho"][3]:7.1e}-{res["rho"][4]:7.1e}]   {res["W"][2]:8.1e}')
