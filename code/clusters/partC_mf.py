"""Cluster MF from bound/unbound complexes x measured efficiency distributions x 0.5, per phase, vs the measured cluster MF."""
import sys, numpy as np, glob
sys.path.insert(0, '/raven/u/uli/dwarf_merger'); import common as C
rng = np.random.default_rng(5); EP = [*C.PHASES]
E = np.load(f'{C.DATADIR}/clusters/clusters_env.npz'); cleanE = (E['f_intruder'] < 0.1) | C.is_merged(E['sep'], E['t'])
cb = cleanE & (E['bound'] > 0) & (E['nucleus'] == 0) & (np.mod(E['snap'], 10) == 0)
def ri(Mv, q):
    Mv = Mv[Mv >= q]; n = len(Mv); return (1 + n / np.sum(np.log(Mv / q)), n) if n >= 5 else (np.nan, n)
for cat, AC in (('cl10', None), ('cl100', None)):
    R = np.vstack([np.load(f)[f'{cat}_rows'] for f in sorted(glob.glob(f'{C.DATADIR}/clouds/clouds_[0-9]*.npz')) if int(f[-7:-4]) % 10 == 0])
    t = R[:, 1]; M = R[:, 2]; al = R[:, 14]; sig = R[:, 10]; vA = R[:, 19]; ms = R[:, 16]
    ok = (M >= 100) & (al > 0) & (vA > 0) & (sig > 0); t, M, al, sig, vA, ms = [a[ok] for a in (t, M, al, sig, vA, ms)]
    at = al * (1 + (vA / sig) ** 2); p = np.load(f'{C.DATADIR}/clouds/step4_fit.npz')[cat]; AC = 10 ** p[2]; bound = at < AC
    d = np.load(f'{C.DATADIR}/clouds/lineage_env_{cat}.npz'); L = d['lineage_rows']; col = {n: i for i, n in enumerate(d['lineage_names'])}
    la = L[:, col['alpha_on']]; le = L[:, col['eps_int']]; lv = L[:, col['vA_on']]; lsg = L[:, col['sig_on']]; g = (la > 0) & (le > 0) & (lv > 0) & (lsg > 0)
    lat = la[g] * (1 + (lv[g] / lsg[g]) ** 2); eb = np.log10(le[g][lat < AC]); eu = np.log10(le[g][lat >= AC])
    print(f'\n===== {cat} (alpha_c = {AC:.1f}): eps bound median {10**np.median(eb):.3f} sig {np.std(eb):.2f} (n={len(eb)}); unbound {10**np.median(eu):.4f} sig {np.std(eu):.2f} (n={len(eu)}) =====')
    print('phase                         a_all  a_bound f_bound | pred a(>300)   pred M_max              | measured clusters a(>300)  M_max   n')
    for t0, t1, lab in EP:
        m = (t >= t0) & (t < t1); mb = m & bound; mu = m & ~bound; S = []; MX = []
        for it in range(200):
            mc = 0.5 * np.concatenate([M[mb] * 10 ** rng.choice(eb, mb.sum()), M[mu] * 10 ** rng.choice(eu, mu.sum())]); S.append(ri(mc, 300)[0]); MX.append(mc.max())
        mcl = cb & (E['t'] >= t0) & (E['t'] < t1); a_m, n_m = ri(E['M'][mcl], 300)
        print(f'{lab:28s} {ri(M[m],300)[0]:.2f}   {ri(M[mb],300)[0]:.2f}   {bound[m].mean():.2f}  |  {np.nanmean(S):.2f}+-{np.nanstd(S):.2f}   {np.median(MX):.1e} [{np.percentile(MX,16):.0e}-{np.percentile(MX,84):.0e}]  |   {a_m:.2f}   {E["M"][mcl].max():.1e}   {n_m}')
