"""Is there ONE efficiency law across phases?  Pool the SF lineages (lineages_<cat>.npz), test eps_int against the pre-onset
virial parameter, sigma_3d, Sigma, n_mean: pooled Spearman, binned medians per phase (do the phases fall on one curve?), and the
Padoan+2012 form log eps = a - b sqrt(alpha_vir) fitted pooled and per phase."""
import sys, numpy as np
sys.path.insert(0, '/raven/u/uli/dwarf_merger'); import common as C
from scipy.stats import spearmanr
EP = [*C.PHASES]
for cat in ('cl100', 'cl10'):
    d = np.load(f'{C.DATADIR}/clouds/lineages_{cat}.npz'); L = d['rows']; col = {n: i for i, n in enumerate(d['names'])}
    g = (L[:, col['Mstar']] > 0) & (L[:, col['M_max']] >= 300) & (L[:, col['t_start']] > 5) & (L[:, col['t_end']] < 221) & np.isfinite(L[:, col['alpha_on']]) & (L[:, col['alpha_on']] > 0)
    le = np.log10(L[g, col['eps_int']]); al = L[g, col['alpha_on']]; ton = L[g, col['t_onset']]
    print(f'\n===== {cat}: {g.sum()} SF lineages pooled =====')
    for key in ('alpha_on', 'sig_on', 'Sigma_on', 'n_mean_on', 'n_max_on', 'M_on'):
        x = L[g, col[key]]; print(f'pooled Spearman(log eps_int, log {key:10s}) = {spearmanr(le, np.log10(x))[0]:+.2f}')
    # binned medians of eps_int vs alpha per phase: one curve?
    ab = np.array([1, 2, 4, 8, 16, 32, 64, 128, 512]); print('\nmedian eps_int per alpha_vir bin (n in brackets):')
    print('phase                 ' + '  '.join(f'{a:>4d}-{b:<4d}' for a, b in zip(ab[:-1], ab[1:])))
    for t0, t1, lab in EP:
        m = (ton >= t0) & (ton < t1); line = f'{lab:20s} '
        for a, b in zip(ab[:-1], ab[1:]):
            mm = m & (al >= a) & (al < b); line += f'  {10**np.median(le[mm]):7.4f}({mm.sum():3d})' if mm.sum() >= 5 else '        --     '
        print(line)
    line = f'{"all":20s} '
    for a, b in zip(ab[:-1], ab[1:]):
        mm = (al >= a) & (al < b); line += f'  {10**np.median(le[mm]):7.4f}({mm.sum():3d})' if mm.sum() >= 5 else '        --     '
    print(line)
    # Padoan form
    A = np.vstack([np.ones(g.sum()), np.sqrt(al)]).T; c = np.linalg.lstsq(A, le, rcond=None)[0]; res = le - A @ c
    print(f'\npooled fit log10 eps_int = {c[0]:.2f} - {-c[1]:.3f} sqrt(alpha_vir): residual scatter {np.std(res):.2f} dex (raw {np.std(le):.2f}); Padoan+12 slope would be 1.38/ln10 = 0.60')
    A2 = np.vstack([np.ones(g.sum()), np.log10(al)]).T; c2 = np.linalg.lstsq(A2, le, rcond=None)[0]; res2 = le - A2 @ c2
    print(f'pooled power law  log10 eps_int = {c2[0]:.2f} {c2[1]:+.2f} log10 alpha_vir: residual {np.std(res2):.2f} dex')
    print('per phase: sqrt-fit intercept/slope/residual, and median alpha_on, median eps_int, residual offset of the phase from the pooled fit')
    for t0, t1, lab in EP:
        m = (ton >= t0) & (ton < t1); cc = np.linalg.lstsq(A[m], le[m], rcond=None)[0]
        print(f'  {lab:20s} n={m.sum():4d}  a={cc[0]:5.2f} b={-cc[1]:5.3f} res={np.std(le[m]-A[m]@cc):.2f} | alpha med {np.median(al[m]):5.1f}  eps med {10**np.median(le[m]):.4f}  offset from pooled {np.mean(res[m]):+.2f} dex')
    # two-variable fit: alpha and Sigma
    Sg = L[g, col['Sigma_on']]; A3 = np.vstack([np.ones(g.sum()), np.sqrt(al), np.log10(Sg)]).T; c3 = np.linalg.lstsq(A3, le, rcond=None)[0]
    print(f'pooled sqrt(alpha) + log Sigma: coefficients {c3[1]:+.3f}, {c3[2]:+.2f}; residual {np.std(le - A3 @ c3):.2f} dex')
