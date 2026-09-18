"""Follow-up to burst_dist.py.
(A) Is the duty cycle a threshold in W_L or in midplane density?  Logistic P(burst | x) per phase for x = W_L, rho_mid_2p,
    W_L / sigma_eff^2 (PRFM midplane density proxy), Sigma_gas_2p; compare the shift of x_50 between phases (a variable whose
    x_50 is phase-invariant is the physical threshold).
(B) Lognormal kernel: among bursts fit log M_y = a + b log W_L + N(0, s) per phase; then build the mock burst-mass distribution
    from the ACTUAL W_L values of all clean patches (draw burst with the logistic f_burst(W_L), mass from the kernel) and compare
    its running power-law index a(>300 / 1000 / 3000) and M_max with the observed M_young distribution per phase.  If they agree,
    the event-size 'slope' is the weight distribution convolved with a 0.4 dex lognormal kernel.
out: prfm/burst_kernel.npz, figs/burst_kernel.png"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common as C
from scipy.optimize import minimize
from scipy.stats import spearmanr
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
rng = np.random.default_rng(1)
D = dict(np.load(f'{C.DATADIR}/prfm/patches_clusters.npz')); WL = np.load(f'{C.DATADIR}/prfm/cluster_link.npz')['WL']
clean = (D['f_intruder'] < 0.1) | C.is_merged(D['sep'], D['t'])
base = clean & (D['Sigma_gas'] > 1) & np.isfinite(WL) & (WL > 0) & (D['t'] > 10)
MY = D['M_young']; MMIN = 500.; burst = base & (MY > MMIN); lW = np.log10(WL); lM = np.log10(np.maximum(MY, 1))
EP = [*C.PHASES]

def logistic_fit(x, y):
    def nll(p):
        pr = np.clip(1 / (1 + np.exp(-p[0] * (x - p[1]))), 1e-9, 1 - 1e-9); return -np.sum(y * np.log(pr) + (1 - y) * np.log(1 - pr))
    r = min((minimize(nll, [k0, np.median(x)], method='Nelder-Mead') for k0 in (1, 4)), key=lambda r: r.fun); return r.x, r.fun
# ---------- (A)
print('(A) logistic P(burst | log x): x_50 per phase and its spread.  Also the deviance explained (1 - nll/nll0).')
X = {'W_L': WL, 'rho_mid_2p': D['rho_mid_2p'], 'W_L/sigma_eff^2': WL / D['sigma_eff'] ** 2, 'W_L/sigma_z^2': WL / D['sigma_z_col'] ** 2,
     'Sigma_gas_2p': D['Sigma_gas_2p'], 'Sigma_gas_2p/H': D['Sigma_gas_2p'] / D['H'], 'P_DE': D['P_DE'], 'fcold*Sigma_2p': D['fcold'] * D['Sigma_gas_2p']}
print('variable            ' + '  '.join(f'{lab[:12]:>16s}' for _, _, lab in EP) + '   range of log x_50   mean k')
for name, x in X.items():
    line = f'{name:19s} '; x50 = []; ks = []
    for t0, t1, lab in EP:
        m = base & (D['t'] >= t0) & (D['t'] < t1) & np.isfinite(x) & (x > 0); lx = np.log10(x[m]); y = burst[m].astype(float)
        (k, c), nll = logistic_fit(lx, y); p0 = y.mean(); nll0 = -m.sum() * (p0 * np.log(p0) + (1 - p0) * np.log(1 - p0))
        x50.append(c); ks.append(k); line += f'  {c:6.2f} (k{k:4.1f},{1-nll/nll0:4.2f})'
    print(line + f'      {max(x50)-min(x50):5.2f}          {np.mean(ks):4.1f}')
# ---------- (B)
print('\n(B) lognormal kernel among bursts: log M_y = a + b log W_L + N(0, s)')
res = {}
for t0, t1, lab in EP:
    mb = burst & (D['t'] >= t0) & (D['t'] < t1); A = np.vstack([np.ones(mb.sum()), lW[mb]]).T
    coef, r, *_ = np.linalg.lstsq(A, lM[mb], rcond=None); s = np.std(lM[mb] - A @ coef)
    # the kernel is truncated at 500 in the data: refit with the truncated-lognormal likelihood so a, b, s are not biased by the cut
    from scipy.stats import norm
    b = coef[1]
    def nll(p):
        a, s_ = p
        if s_ < 0.05: return 1e30
        mu = a + b * lW[mb]; Z = np.clip(1 - norm.cdf((np.log10(MMIN) - mu) / s_), 1e-12, 1); return -(np.sum(norm.logpdf(lM[mb], mu, s_)) - np.sum(np.log(Z)))
    rr = min((minimize(nll, [coef[0] + d, s + ds], method='Nelder-Mead', options=dict(maxiter=4000)) for d in (-0.5, 0, 0.3) for ds in (0, 0.2)), key=lambda r: r.fun)
    a, s_ = rr.x
    m = base & (D['t'] >= t0) & (D['t'] < t1); (k, w50), _ = logistic_fit(lW[m], burst[m].astype(float))
    res[lab] = dict(ols=(coef[0], coef[1], s), tln=(a, b, s_), duty=(k, w50), mask=m)
    print(f'{lab:20s} OLS: a={coef[0]:5.2f} b={coef[1]:4.2f} s={s:4.2f} | truncated-LN MLE (b fixed): a={a:5.2f} s={s_:4.2f} | duty k={k:4.1f} logW50={w50:4.2f}   n_burst={mb.sum()}')

Mcold = D['fcold'] * D['Sigma_gas_2p'] * 500. ** 2
print('\nreservoir: M_young / M_cold,patch among bursts (M_cold = fcold x Sigma_gas,2p x (0.5 kpc)^2)')
for t0, t1, lab in EP:
    mb = burst & (D['t'] >= t0) & (D['t'] < t1) & (Mcold > 0); r = MY[mb] / Mcold[mb]
    print(f'{lab:20s} median {np.median(r):.3f}  84% {np.percentile(r,84):.3f}  95% {np.percentile(r,95):.3f}  99% {np.percentile(r,99):.3f}  max {r.max():.3f}   Spearman(M_y, M_cold) {spearmanr(np.log10(MY[mb]), np.log10(Mcold[mb]))[0]:+.2f}, partial at fixed W_L {spearmanr(np.log10(MY[mb]) - np.polyval(np.polyfit(lW[mb], np.log10(MY[mb]), 1), lW[mb]), np.log10(Mcold[mb]) - np.polyval(np.polyfit(lW[mb], np.log10(Mcold[mb]), 1), lW[mb]))[0]:+.2f}')
QCAP = 0.3
def run_index(M, mmin):
    M = M[M >= mmin]; n = len(M)
    return (1 + n / np.sum(np.log(M / mmin)), n) if n >= 5 else (np.nan, n)
print('\nmock vs observed M_young distribution (bursts, M > 500) per phase; mock = actual W_L of all clean patches x f_burst(W_L) x kernel, 200 realisations')
print('phase                 sample     n_burst   a(>500)      a(>1000)     a(>3000)     a(>10000)    M_max        f(M>1e4)')
mocks = {}
for t0, t1, lab in EP:
    m = res[lab]['mask']; a, b, s_ = res[lab]['tln']; k, w50 = res[lab]['duty']
    obs = MY[burst & m]
    o = [run_index(obs, mm)[0] for mm in (500, 1000, 3000, 10000)]
    print(f'{lab:20s}  observed    {len(obs):5d}    ' + '  '.join(f'{v:5.2f}       ' for v in o) + f'  {obs.max():8.2e}   {np.mean(obs > 1e4):.3f}')
    idx = []; mmax = []; nb = []; f4 = []; keep = []
    for it in range(200):
        pb = 1 / (1 + np.exp(-k * (lW[m] - w50))); hit = rng.random(m.sum()) < pb
        mm_ = 10 ** (a + b * lW[m][hit] + s_ * rng.standard_normal(hit.sum())); mm_ = mm_[mm_ > MMIN]
        idx.append([run_index(mm_, q)[0] for q in (500, 1000, 3000, 10000)]); mmax.append(mm_.max()); nb.append(len(mm_)); f4.append(np.mean(mm_ > 1e4))
        if it < 20: keep.append(mm_)
    idx = np.array(idx); mocks[lab] = (obs, keep)
    print(f'{lab:20s}  mock        {np.mean(nb):5.0f}    ' + '  '.join(f'{np.nanmean(idx[:, j]):5.2f}+-{np.nanstd(idx[:, j]):4.2f}' for j in range(4)) + f'  {np.median(mmax):8.2e}   {np.mean(f4):.3f}')
    # variant: kernel with W-independent mean (b = 0): does the W dependence matter for the slope?
    idx0 = []
    for it in range(200):
        pb = 1 / (1 + np.exp(-k * (lW[m] - w50))); hit = rng.random(m.sum()) < pb
        mu0 = a + b * np.median(lW[burst & m]); mm_ = 10 ** (mu0 + s_ * rng.standard_normal(hit.sum())); mm_ = mm_[mm_ > MMIN]
        idx0.append([run_index(mm_, q)[0] for q in (500, 1000, 3000, 10000)])
    idx0 = np.array(idx0)
    print(f'{lab:20s}  mock b=0             ' + '  '.join(f'{np.nanmean(idx0[:, j]):5.2f}+-{np.nanstd(idx0[:, j]):4.2f}' for j in range(4)))
    # variant: kernel capped at QCAP x M_cold of the same patch (reservoir)
    idxc = []; mmaxc = []; f4c = []
    for it in range(200):
        pb = 1 / (1 + np.exp(-k * (lW[m] - w50))); hit = rng.random(m.sum()) < pb
        mm_ = 10 ** (a + b * lW[m][hit] + s_ * rng.standard_normal(hit.sum())); mm_ = np.minimum(mm_, QCAP * np.maximum(Mcold[m][hit], 1)); mm_ = mm_[mm_ > MMIN]
        idxc.append([run_index(mm_, q)[0] for q in (500, 1000, 3000, 10000)]); mmaxc.append(mm_.max()); f4c.append(np.mean(mm_ > 1e4))
    idxc = np.array(idxc)
    print(f'{lab:20s}  mock cap {QCAP} M_cold   ' + '  '.join(f'{np.nanmean(idxc[:, j]):5.2f}+-{np.nanstd(idxc[:, j]):4.2f}' for j in range(4)) + f'  {np.median(mmaxc):8.2e}   {np.mean(f4c):.3f}')
np.savez(f'{C.DATADIR}/prfm/burst_kernel.npz', **{lab.replace(' ', '_'): np.array(res[lab]['tln'] + res[lab]['duty']) for lab in res})

ink = '#1f2937'; blue = '#2f6fb2'; light = '#9dbfe3'; orange = '#c2410c'
fig, ax = plt.subplots(1, 4, figsize=(16, 3.9), sharey=True)
edges = np.logspace(np.log10(MMIN), 6.3, 20); c = np.sqrt(edges[1:] * edges[:-1]); w = np.diff(edges)
for a_, (t0, t1, lab) in zip(ax, EP):
    obs, keep = mocks[lab]
    for mm_ in keep:
        h, _ = np.histogram(mm_, edges); ok = h > 0; a_.plot(c[ok], h[ok] / w[ok], color=light, lw=0.7, alpha=0.6)
    h, _ = np.histogram(obs, edges); ok = h > 0; a_.errorbar(c[ok], h[ok] / w[ok], yerr=np.sqrt(h[ok]) / w[ok], fmt='o', ms=4, color=ink, label='observed bursts')
    a_.plot([], [], color=light, lw=1.5, label='mock: $W_L$ distribution $\\otimes$ lognormal kernel')
    a_.set(xscale='log', yscale='log', xlabel='$M_{\\rm young}$ [M$_\\odot$]'); a_.set_title(lab, loc='left', fontsize=10, color=ink); a_.grid(alpha=0.2); a_.spines[['top', 'right']].set_visible(False)
    a, b, s_ = res[lab]['tln']; a_.text(0.97, 0.95, f'$\\mu = {a:.2f} + {b:.2f}\\,\\log W_L$\n$s = {s_:.2f}$ dex', transform=a_.transAxes, ha='right', va='top', fontsize=8, color=ink)
ax[0].set_ylabel('$dN/dM$ [M$_\\odot^{-1}$]'); ax[0].legend(fontsize=7, frameon=False, loc='lower left')
fn = f'{C.DATADIR}/prfm/figs/burst_kernel.png'; fig.savefig(fn, dpi=140, bbox_inches='tight'); os.system(f'cp {fn} ~/out/'); print('fig', fn)
