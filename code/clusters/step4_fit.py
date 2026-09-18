"""Two-level boundedness fit of eps_int vs alpha_tot (kinetic + magnetic) per catalogue; saves step4_fit.npz used by paper_figs.partC_eff."""
import sys, numpy as np
sys.path.insert(0, '/raven/u/uli/dwarf_merger'); import common as C
from scipy.optimize import minimize
EP = [*C.PHASES]; out = {}
for cat in ('cl100', 'cl10'):
    d = np.load(f'{C.DATADIR}/clouds/lineage_env_{cat}.npz'); L = d['lineage_rows']; col = {n: i for i, n in enumerate(d['lineage_names'])}
    c = lambda k: L[:, col[k]]; al = c('alpha_on'); eps = c('eps_int'); ton = c('t_onset'); sc = c('sig_on'); vA = c('vA_on')
    ok = (al > 0) & (eps > 0) & (vA > 0); al, eps, ton, sc, vA = [a[ok] for a in (al, eps, ton, sc, vA)]; at = al * (1 + (vA / sc) ** 2); le = np.log10(eps)
    def model(p, a):
        lb, lu, lac, m = p; return lu + np.log10(1 + (10 ** (lb - lu) - 1) / (1 + (a / 10 ** lac) ** m))
    nll = lambda p: np.sum(np.abs(le - model(p, at)))
    p = min((minimize(nll, [np.log10(0.05), np.log10(0.001), np.log10(ac), m0], method='Nelder-Mead', options=dict(maxiter=5000)) for ac in (2, 4, 8) for m0 in (1, 2, 4)), key=lambda r: r.fun).x
    r = le - model(p, at); offs = [np.median(r[(ton >= t0) & (ton < t1)]) for t0, t1, _ in EP]; out[cat] = p
    print(f'{cat}: eps_b={10**p[0]:.3f} eps_u={10**p[1]:.4f} alpha_c={10**p[2]:.2f} m={p[3]:.1f} scatter {np.std(r):.2f} dex; phase medians ' + ' '.join(f'{o:+.2f}' for o in offs))
    for (t0, t1, lab) in EP:
        m = (ton >= t0) & (ton < t1); print(f'   {lab:28s} n={m.sum():4d} f(alpha_tot<alpha_c)={np.mean(at[m] < 10**p[2]):.2f} eps median {np.median(eps[m]):.4f} model median {np.median(10**model(p, at[m])):.4f}')
np.savez(f'{C.DATADIR}/clouds/step4_fit.npz', **out)
