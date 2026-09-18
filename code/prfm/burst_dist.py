"""Empirical P(m | W): the distribution of the 10 Myr burst mass M_young of a 0.5 kpc patch at fixed layer weight W_L.
(1) duty cycle: f_burst(W_L) = P(M_young > 500 | W_L), logistic fit per phase;
(2) shape at fixed W_L (0.5 dex bins, bursting patches): quantiles, width, MLE power law vs lognormal vs Schechter (all
    above 500), per phase and all phases;
(3) what carries the residual at fixed W_L: partial Spearman of log M_young (bursts) with every patch variable after
    regressing on log W_L; width of log M_young per bin vs the median cold Mach number of the bin;
(4) feedback-limited-collapse predictor: m_ff = eps M_gas,2p * 10 Myr / t_ff(rho_mid_2p) vs M_young at fixed W_L.
out: prfm/burst_dist.npz, figs/burst_dist.png"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common as C
from scipy.stats import spearmanr, norm
from scipy.optimize import minimize
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
D = dict(np.load(f'{C.DATADIR}/prfm/patches_clusters.npz')); WL = np.load(f'{C.DATADIR}/prfm/cluster_link.npz')['WL']
clean = (D['f_intruder'] < 0.1) | C.is_merged(D['sep'], D['t'])
base = clean & (D['Sigma_gas'] > 1) & np.isfinite(WL) & (WL > 0) & (D['t'] > 10)
MY = D['M_young']; MMIN = 500.; burst = base & (MY > MMIN); lW = np.log10(WL); lM = np.log10(np.maximum(MY, 1))
EP = [*C.PHASES, (25, 232, 'all phases')]
print(f'{base.sum()} clean patches, {burst.sum()} bursting (M_young > {MMIN:g})')

# ---------- (1) duty cycle
def logistic_fit(x, y):
    def nll(p):
        z = p[0] * (x - p[1]); pr = 1 / (1 + np.exp(-z)); pr = np.clip(pr, 1e-9, 1 - 1e-9)
        return -np.sum(y * np.log(pr) + (1 - y) * np.log(1 - pr))
    r = minimize(nll, [1.0, 4.5], method='Nelder-Mead'); return r.x
print('\n(1) duty cycle: P(burst | W_L) = 1 / (1 + exp(-k (log W_L - log W_50)))')
print('phase                  n_patch  n_burst   k     log W_50   f_burst at log W = 3.5 / 4.0 / 4.5 / 5.0 / 5.5')
dc = {}
for t0, t1, lab in EP:
    m = base & (D['t'] >= t0) & (D['t'] < t1); k, w50 = logistic_fit(lW[m], burst[m].astype(float)); dc[lab] = (k, w50)
    fb = [(burst & m & (lW >= l) & (lW < l + 0.5)).sum() / max((m & (lW >= l) & (lW < l + 0.5)).sum(), 1) for l in (3.5, 4.0, 4.5, 5.0, 5.5)]
    print(f'{lab:22s} {m.sum():6d}  {(burst & m).sum():6d}   {k:5.2f}  {w50:6.2f}     ' + ' / '.join(f'{f:.2f}' for f in fb))

# ---------- (2) shape at fixed W_L
def fit_pl(M):
    n = len(M); al = 1 + n / np.sum(np.log(M / MMIN)); ll = n * np.log(al - 1) - n * (al - 1) * np.log(MMIN) - al * np.sum(np.log(M)) + n * (al - 1) * np.log(MMIN) * 0
    ll = n * np.log((al - 1) / MMIN) - al * np.sum(np.log(M / MMIN)); return al, ll
def fit_ln(M):
    x = np.log10(M); a = np.log10(MMIN)
    def nll(p):
        mu, s = p
        if s <= 0.02: return 1e30
        Z = 1 - norm.cdf((a - mu) / s); return -(np.sum(norm.logpdf(x, mu, s)) - len(x) * np.log(max(Z, 1e-12)))
    r = min((minimize(nll, [m0, s0], method='Nelder-Mead') for m0 in (2.8, 3.3, 3.8) for s0 in (0.3, 0.6)), key=lambda r: r.fun)
    return r.x[0], r.x[1], -r.fun - np.sum(np.log(np.log(10) * M))   # convert log10-density to dM density (Jacobian) for comparability
def fit_sch(M):
    x = np.logspace(np.log10(MMIN), 8, 3000)
    def nll(p):
        al, lMc = p
        if al < -1 or al > 4 or lMc < np.log10(MMIN) or lMc > 8: return 1e30
        Z = np.trapezoid(x ** -al * np.exp(-x / 10 ** lMc), x); return -(np.sum(-al * np.log(M) - M / 10 ** lMc) - len(M) * np.log(Z))
    r = min((minimize(nll, [a0, l0], method='Nelder-Mead') for a0 in (1.0, 1.8) for l0 in (3.5, 4.5, 5.5)), key=lambda r: r.fun)
    return r.x[0], 10 ** r.x[1], -r.fun
def fit_pl_ll(M):
    n = len(M); al = 1 + n / np.sum(np.log(M / MMIN)); return al, n * np.log((al - 1) / MMIN) - al * np.sum(np.log(M / MMIN))
print('\n(2) bursting patches at fixed W_L: shape of P(M_young | W_L).  dlnL = lnL(model) - lnL(power law), same lower limit 500')
print('phase                 log W_L bin   n    med M_y   16%     84%   sig(logM)  a_PL   | LN mu  sig  dlnL | Sch a   Mc     dlnL | Mach_c med  f_burst')
rows = []
for t0, t1, lab in EP:
    for l0 in np.arange(3.0, 6.0, 0.5):
        mb = burst & (D['t'] >= t0) & (D['t'] < t1) & (lW >= l0) & (lW < l0 + 0.5); ma = base & (D['t'] >= t0) & (D['t'] < t1) & (lW >= l0) & (lW < l0 + 0.5)
        if mb.sum() < 25: continue
        M = MY[mb]; al, llp = fit_pl_ll(M); mu, s, lll = fit_ln(M); asch, Mc, lls = fit_sch(M)
        q = np.percentile(M, [16, 50, 84]); sig = np.std(np.log10(M)); mach = np.nanmedian(D['Mach_cold'][mb])
        rows.append((t0, l0, mb.sum(), q[1], q[0], q[2], sig, al, mu, s, lll - llp, asch, Mc, lls - llp, mach, mb.sum() / ma.sum()))
        print(f'{lab:20s}  {l0:.1f}-{l0+0.5:.1f}   {mb.sum():4d}  {q[1]:7.0f} {q[0]:6.0f} {q[2]:7.0f}   {sig:5.2f}    {al:4.2f}   | {mu:4.2f} {s:4.2f} {lll-llp:+6.1f} | {asch:5.2f} {Mc:8.1e} {lls-llp:+6.1f} |  {mach:5.1f}      {mb.sum()/ma.sum():.2f}')
R = np.array(rows)

# ---------- (3) residual carriers at fixed W_L
print('\n(3) partial Spearman: residual of log M_young (bursts) after linear regression on log W_L, vs log of each patch variable (per phase)')
vars_ = ['Sigma_gas', 'Sigma_gas_2p', 'rho_mid_2p', 'rho_sd', 'H', 'sigma_turb', 'sigma_z_col', 'sigma_eff', 'Mach_cold', 'fcold', 'fH2', 'Q', 'kappa', 'Omega', 'Sigma_star', 'G0_mid', 'SigSFR_40', 'Ptot_2p', 'Pturb_2p', 'P_DE', 'rho_std', 'R']
print('variable        ' + '  '.join(f'{lab[:14]:>14s}' for _, _, lab in EP))
for v in vars_:
    line = f'{v:15s} '
    for t0, t1, lab in EP:
        m = burst & (D['t'] >= t0) & (D['t'] < t1) & np.isfinite(D[v]) & (D[v] > 0)
        if m.sum() < 30: line += '            --'; continue
        A = np.vstack([np.ones(m.sum()), lW[m]]).T; res = lM[m] - A @ np.linalg.lstsq(A, lM[m], rcond=None)[0]
        # partial: also remove W_L from the variable
        xv = np.log10(D[v][m]); rx = xv - A @ np.linalg.lstsq(A, xv, rcond=None)[0]
        line += f'  {spearmanr(res, rx)[0]:+13.2f}'
    print(line)
# prior star formation as a memory term: SigSFR_40 - SigSFR_10/4 ~ SF in the previous 30 Myr
prev = D['SigSFR_40'] * 40 - D['SigSFR_10'] * 10
print('previous-30-Myr SF (SigSFR_40*40 - SigSFR_10*10):')
line = f'{"prev30":15s} '
for t0, t1, lab in EP:
    m = burst & (D['t'] >= t0) & (D['t'] < t1) & (prev > 0)
    A = np.vstack([np.ones(m.sum()), lW[m]]).T; res = lM[m] - A @ np.linalg.lstsq(A, lM[m], rcond=None)[0]
    xv = np.log10(prev[m]); rx = xv - A @ np.linalg.lstsq(A, xv, rcond=None)[0]; line += f'  {spearmanr(res, rx)[0]:+13.2f}'
print(line)
# width vs Mach
ok = R[:, 0] < 200
print(f'\nwidth sig(log M_y) per (phase, W bin) vs median cold Mach: Spearman {spearmanr(R[ok, 6], R[ok, 14])[0]:+.2f} over {ok.sum()} bins; '
      f'vs log W_L: {spearmanr(R[ok, 6], R[ok, 1])[0]:+.2f}; power-law index vs log W_L: {spearmanr(R[ok, 7], R[ok, 1])[0]:+.2f}, vs Mach: {spearmanr(R[ok, 7], R[ok, 14])[0]:+.2f}')

# ---------- (4) feedback-limited collapse predictor
# t_ff = sqrt(3 pi / 32 G rho); rho_mid_2p assumed in Msun/pc^3 if median < 10 else cm^-3 (checked below)
rho = D['rho_mid_2p']; med = np.nanmedian(rho[burst]); print(f'\n(4) rho_mid_2p median over bursts = {med:.3g} (units per patches.py)')
G = 4.30091e-3  # pc km^2/s^2 / Msun
rho_msun = rho if med < 5 else rho * 1.4 * 1.6726e-24 / (1.989e33 / (3.086e18) ** 3)   # cm^-3 -> Msun/pc^3 (mu = 1.4)
tff = np.sqrt(3 * np.pi / (32 * G * rho_msun)) * 3.086e13 / 3.156e13   # Myr
Mgas = D['Sigma_gas_2p'] * (500.) ** 2
mff = Mgas * 10. / tff   # eps = 1, 10 Myr window; compare in log
for t0, t1, lab in EP:
    m = burst & (D['t'] >= t0) & (D['t'] < t1) & np.isfinite(mff) & (mff > 0)
    A = np.vstack([np.ones(m.sum()), lW[m]]).T; res = lM[m] - A @ np.linalg.lstsq(A, lM[m], rcond=None)[0]
    xv = np.log10(mff[m]); rx = xv - A @ np.linalg.lstsq(A, xv, rcond=None)[0]
    eff = np.median(MY[m] / mff[m])
    print(f'{lab:22s} Spearman(log M_y, log M_gas/t_ff) = {spearmanr(lM[m], xv)[0]:+.2f}, partial at fixed W_L = {spearmanr(res, rx)[0]:+.2f}, median M_y/(M_gas 10Myr/t_ff) = {eff:.3f}, median t_ff = {np.median(tff[m]):.1f} Myr')
np.savez(f'{C.DATADIR}/prfm/burst_dist.npz', rows=R, duty=np.array([[dc[l][0], dc[l][1]] for _, _, l in EP]))

# ---------- figure: P(M_young | W_L) per W bin (all phases) and per phase at one W bin
ink = '#1f2937'; blue = '#2f6fb2'; light = '#9dbfe3'; orange = '#c2410c'; grey = '#6b7280'
fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
edges = np.logspace(np.log10(MMIN), 6.5, 22); c = np.sqrt(edges[1:] * edges[:-1]); w = np.diff(edges)
cols = plt.cm.viridis(np.linspace(0.1, 0.9, 6))
for i, l0 in enumerate(np.arange(3.0, 6.0, 0.5)):
    mb = burst & (lW >= l0) & (lW < l0 + 0.5)
    if mb.sum() < 25: continue
    h, _ = np.histogram(MY[mb], edges); ok = h > 0
    ax[0].plot(c[ok], h[ok] / w[ok] / mb.sum(), 'o-', ms=3, color=cols[i], lw=1.2, label=f'log W_L {l0:.1f}-{l0+0.5:.1f} (n={mb.sum()})')
ax[0].set(xscale='log', yscale='log', xlabel='$M_{\\rm young}$ (10 Myr, 0.5 kpc patch) [M$_\\odot$]', ylabel='$P(M_{\\rm young}\\,|\\,W_L)$  [M$_\\odot^{-1}$]'); ax[0].legend(fontsize=7, frameon=False); ax[0].set_title('Burst mass at fixed layer weight, all phases', loc='left', fontsize=10, color=ink)
for (t0, t1, lab), col in zip(EP[:4], (light, blue, orange, ink)):
    m = base & (D['t'] >= t0) & (D['t'] < t1); xx = np.linspace(3, 6, 50); k, w50 = dc[lab]
    ax[1].plot(xx, 1 / (1 + np.exp(-k * (xx - w50))), color=col, lw=1.8, label=lab)
    lb = np.arange(3.0, 6.0, 0.5); fb = [(burst & m & (lW >= l) & (lW < l + 0.5)).sum() / max((m & (lW >= l) & (lW < l + 0.5)).sum(), 1) if (m & (lW >= l) & (lW < l + 0.5)).sum() >= 20 else np.nan for l in lb]
    ax[1].plot(lb + 0.25, fb, 'o', color=col, ms=4)
ax[1].set(xlabel='log $W_L$ [K cm$^{-3}$]', ylabel='$P$(burst | $W_L$)', ylim=(0, 1)); ax[1].legend(fontsize=8, frameon=False); ax[1].set_title('Duty cycle', loc='left', fontsize=10, color=ink)
ok = R[:, 0] < 200
sc = ax[2].scatter(R[ok, 1] + 0.25, R[ok, 7], c=R[ok, 14], cmap='magma', s=30 + R[ok, 2] / 10, edgecolor=ink, lw=0.4)
plt.colorbar(sc, ax=ax[2], label='median cold Mach number'); ax[2].set(xlabel='log $W_L$ [K cm$^{-3}$]', ylabel='power-law index of $P(M_{\\rm young}|W_L)$ above 500'); ax[2].set_title('Index per (phase, weight) bin', loc='left', fontsize=10, color=ink)
for a in ax: a.grid(alpha=0.2); a.spines[['top', 'right']].set_visible(False)
fn = f'{C.DATADIR}/prfm/figs/burst_dist.png'; fig.savefig(fn, dpi=140, bbox_inches='tight'); os.system(f'cp {fn} ~/out/'); print('fig', fn)
