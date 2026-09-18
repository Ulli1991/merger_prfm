"""Is a burst the collapse of one cloud, and where does the 0.4 dex scatter come from?
Reads /ptmp/uli/dwarf_merger/clouds/clouds_KKK.npz (clouds.py) on every 10th snapshot (independent populations).
(1) cloud mass function per phase (cl10: cold gas n>10, cl100: n>100): running MLE index a(>300/1000/3000), M_max, and the
    same for the star-forming clouds only; compare with the burst-mass distribution indices from burst_kernel (1.65/1.56/1.54/1.39 at >500).
(2) per star-forming cloud: efficiency eps_10 = M_star(10 Myr) / M_cloud; its width; P(SF | M_cloud, alpha_vir).
(3) events: FoF of the new stars at 50 pc; per event the number of parent clouds and the mass fraction from the dominant cloud;
    the same for 5 pc groups (clusters).
(4) variance decomposition among SF clouds: var(log M_star) = var(log M_cl) + var(log eps) + 2 cov.
out: /ptmp/uli/dwarf_merger/clouds/clouds_analysis.npz, figs/clouds_mf.png"""
import sys, os, glob, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C
from fof import fof
from scipy.stats import spearmanr
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
OUT = f'{C.DATADIR}/clouds'; os.makedirs(f'{OUT}/figs', exist_ok=True)
EP = [*C.PHASES]
STEP = int(sys.argv[1]) if len(sys.argv) > 1 else 10
files = sorted(glob.glob(f'{OUT}/clouds_[0-9]*.npz')); files = [f for f in files if int(f[-7:-4]) % STEP == 0]
print(f'{len(files)} snapshots (every {STEP}th)')
cl = {'cl10': [], 'cl100': []}; ev = []; cls = []; stars_tot = 0.0; stars_from = {'cl10': 0.0, 'cl100': 0.0}; nstar_par = []
for fn in files:
    d = np.load(fn); k = int(d['snap']); t = float(d['t'])
    for cat in cl:
        R = d[f'{cat}_rows']
        if len(R): cl[cat].append(R)
    sm = d['star_mass']; stars_tot += sm.sum()
    nstar_par.append(np.column_stack([np.full(len(sm), t), d['star_par_n'], d['star_par_T']]))
    for cat in cl: stars_from[cat] += sm[d[f'{cat}_star_cl'] >= 0].sum()
    # events: FoF of the new stars at 50 pc and 5 pc; parent clouds from cl10
    if len(sm) < 10: continue
    pos = d['star_pos'].astype(np.float64); sc = d['cl10_star_cl']
    for link, store in ((50e-3, ev), (5e-3, cls)):
        lab = fof(pos, link); cnt = np.bincount(lab)
        for g in np.flatnonzero(cnt >= 25):
            m = lab == g; M = sm[m].sum(); par = sc[m]; ok = par >= 0
            if ok.sum() == 0: store.append((k, t, M, 0, 0.0, 0.0)); continue
            u, inv = np.unique(par[ok], return_inverse=True); mp = np.bincount(inv, weights=sm[m][ok])
            store.append((k, t, M, len(u), mp.max() / M, sm[m][ok].sum() / M))
for cat in cl: cl[cat] = np.vstack(cl[cat])
ev = np.array(ev); cls = np.array(cls); SP = np.vstack(nstar_par)
names = ['snap', 't', 'M', 'N', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'sig3d', 'rh_pc', 'n_mean', 'n_max', 'alpha_vir', 'Mstar_3', 'Mstar_10']
col = {n: i for i, n in enumerate(names)}
print(f'stars formed within 10 Myr of a snapshot: {stars_tot:.3g} Msun; from cl10 clouds {stars_from["cl10"]/stars_tot:.2f}, from cl100 clouds {stars_from["cl100"]/stars_tot:.2f}')
for t0, t1, lab in EP:
    m = (SP[:, 0] >= t0) & (SP[:, 0] < t1)
    print(f'   {lab:20s}: parent gas of the new stars at the snapshot: n>10 & cold {np.mean((SP[m,1]>10)&(SP[m,2]<1e3)):.2f}, n>100 {np.mean(SP[m,1]>100):.2f}, n<1 {np.mean(SP[m,1]<1):.2f}')

def run_index(M, mmin):
    M = M[M >= mmin]; n = len(M)
    return (1 + n / np.sum(np.log(M / mmin)), (1 + n / np.sum(np.log(M / mmin)) - 1) / np.sqrt(n), n) if n >= 5 else (np.nan, np.nan, n)
print('\n(1) cloud mass function per phase: running index a(>M) [all clouds | star-forming clouds (M_star,10 > 0)]; burst-mass distribution for reference')
BURST = {'before 1st passage': (1.65, 1.72, 2.08), 'between passages': (1.56, 1.72, 2.11), 'coalescence': (1.54, 1.65, 1.82), 'remnant': (1.39, 1.42, 1.37)}
print('cat    phase                 n_cl   a(>300)      a(>1000)     a(>3000)     M_max   | SF clouds: n   a(>300)  a(>1000)  a(>3000)  M_max   | bursts a(>500/1e3/3e3)')
res = {}
for cat in cl:
    R = cl[cat]
    for t0, t1, lab in EP:
        m = (R[:, 1] >= t0) & (R[:, 1] < t1); sf = m & (R[:, col['Mstar_10']] > 0); M = R[m, 2]; Ms = R[sf, 2]
        a = [run_index(M, q) for q in (300, 1000, 3000)]; b = [run_index(Ms, q) for q in (300, 1000, 3000)]
        res[(cat, lab)] = (a, b)
        print(f'{cat:6s} {lab:20s} {m.sum():5d}  ' + '  '.join(f'{x[0]:4.2f}+-{x[1]:4.2f}' for x in a) + f'  {M.max():8.2e} |       {sf.sum():4d}   ' + '  '.join(f'{x[0]:4.2f}    ' for x in b) + f'  {Ms.max() if sf.sum() else 0:8.2e} | ' + '/'.join(f'{v:.2f}' for v in BURST.get(lab, (float('nan'),)*3)))

print('\n(2) star-forming clouds (cl10): efficiency eps_10 = M_star(10 Myr)/M_cloud and what sets whether a cloud forms stars')
print('phase                 n_SF/n_cl   f_SF(M>1e3)  f_SF(M>1e4)  eps med   16-84%          sig(log eps)  Spearman(eps, M) (eps, alpha_vir) (eps, n_mean) | P(SF): rho(M) rho(alpha) rho(n_max)')
R = cl['cl10']
for t0, t1, lab in EP:
    m = (R[:, 1] >= t0) & (R[:, 1] < t1) & (R[:, 2] >= 300); sf = m & (R[:, col['Mstar_10']] > 0)
    eps = R[sf, col['Mstar_10']] / R[sf, 2]; le = np.log10(eps)
    f3 = np.mean(R[m & (R[:, 2] > 1e3), col['Mstar_10']] > 0); f4 = np.mean(R[m & (R[:, 2] > 1e4), col['Mstar_10']] > 0) if (m & (R[:, 2] > 1e4)).sum() else np.nan
    r = lambda x: spearmanr(le, np.log10(x[sf]))[0]
    pr = lambda x: spearmanr((R[m, col['Mstar_10']] > 0).astype(float), np.log10(x[m]))[0]
    print(f'{lab:20s} {sf.sum():4d}/{m.sum():5d}    {f3:.2f}         {f4:.2f}       {np.median(eps):.3f}   {np.percentile(eps,16):.3f}-{np.percentile(eps,84):.3f}     {np.std(le):.2f}          {r(R[:,2]):+.2f}          {r(R[:,col["alpha_vir"]]):+.2f}          {r(R[:,col["n_mean"]]):+.2f}    |        {pr(R[:,2]):+.2f}   {pr(R[:,col["alpha_vir"]]):+.2f}    {pr(R[:,col["n_max"]]):+.2f}')

print('\n(3) events = 50 pc FoF groups of the stars formed in 10 Myr (>= 25 stars); clusters = 5 pc groups. Parent clouds from cl10.')
print('phase                 n_ev   f(1 cloud)  f(<=2)  dominant-cloud mass fraction med 16-84%   n_clouds med   f_in_cloud | clusters: n  f(1 cloud)  dominant med')
for t0, t1, lab in EP:
    m = (ev[:, 1] >= t0) & (ev[:, 1] < t1) & (ev[:, 3] > 0); mc = (cls[:, 1] >= t0) & (cls[:, 1] < t1) & (cls[:, 3] > 0)
    print(f'{lab:20s} {m.sum():5d}   {np.mean(ev[m,3]==1):.2f}       {np.mean(ev[m,3]<=2):.2f}     {np.median(ev[m,4]):.2f}  {np.percentile(ev[m,4],16):.2f}-{np.percentile(ev[m,4],84):.2f}                  {np.median(ev[m,3]):3.0f}         {np.median(ev[m,5]):.2f}    |        {mc.sum():4d}   {np.mean(cls[mc,3]==1):.2f}       {np.median(cls[mc,4]):.2f}')
    # events above 1e4 Msun
    mb = m & (ev[:, 2] > 1e4)
    if mb.sum() >= 5: print(f'      events > 1e4 Msun: n={mb.sum()}, f(1 cloud) {np.mean(ev[mb,3]==1):.2f}, dominant fraction med {np.median(ev[mb,4]):.2f}, n_clouds med {np.median(ev[mb,3]):.0f}')

print('\n(4) variance decomposition among star-forming cl10 clouds (M_cl >= 300): log M_star = log M_cl + log eps')
print('phase                 sig(logMstar)  sig(logMcl)  sig(logeps)  corr(logMcl, logeps)   -> share of var from M_cl / eps / 2cov')
for t0, t1, lab in EP:
    sf = (R[:, 1] >= t0) & (R[:, 1] < t1) & (R[:, 2] >= 300) & (R[:, col['Mstar_10']] > 0)
    lm = np.log10(R[sf, 2]); ls = np.log10(R[sf, col['Mstar_10']]); le = ls - lm
    v = np.var(ls); c = np.cov(lm, le)[0, 1]
    print(f'{lab:20s}   {np.std(ls):.2f}          {np.std(lm):.2f}         {np.std(le):.2f}         {np.corrcoef(lm, le)[0,1]:+.2f}               {np.var(lm)/v:.2f} / {np.var(le)/v:.2f} / {2*c/v:+.2f}')
np.savez(f'{OUT}/clouds_analysis.npz', cl10=cl['cl10'], cl100=cl['cl100'], events=ev, clusters=cls)

# figure: cloud MF (cl10, all and star-forming) per phase with the burst-mass distribution overlaid
ink = '#1f2937'; blue = '#2f6fb2'; light = '#9dbfe3'; orange = '#c2410c'
D = dict(np.load(f'{C.DATADIR}/prfm/patches_clusters.npz')); clean = (D['f_intruder'] < 0.1) | C.is_merged(D['sep'], D['t'])
burst = clean & (D['Sigma_gas'] > 1) & (D['t'] > 10) & (D['M_young'] > 500) & (np.mod(D['snap'], 10) == 0)
fig, ax = plt.subplots(1, 4, figsize=(15, 3.8), sharey=True); edges = np.logspace(2, 7, 26); c = np.sqrt(edges[1:] * edges[:-1]); w = np.diff(edges)
for a_, (t0, t1, lab) in zip(ax, EP):
    m = (R[:, 1] >= t0) & (R[:, 1] < t1); sf = m & (R[:, col['Mstar_10']] > 0)
    for sel, colr, name in ((m, light, 'all cold clouds ($n > 10$)'), (sf, blue, 'star-forming clouds')):
        h, _ = np.histogram(R[sel, 2], edges); ok = h > 0; a_.plot(c[ok], h[ok] / w[ok] / sel.sum(), 'o-', ms=3, color=colr, lw=1.2, label=name)
    h, _ = np.histogram(R[sf, col['Mstar_10']], edges); ok = h > 0; a_.plot(c[ok], h[ok] / w[ok] / sf.sum(), 's-', ms=3, color=orange, lw=1.2, label='stars formed per cloud (10 Myr)')
    mb = burst & (D['t'] >= t0) & (D['t'] < t1); h, _ = np.histogram(D['M_young'][mb], edges); ok = h > 0; a_.plot(c[ok], h[ok] / w[ok] / mb.sum(), 'd--', ms=3, color=ink, lw=1.0, label='patch bursts $M_{\\rm young}$')
    a_.set(xscale='log', yscale='log', xlabel='M [M$_\\odot$]'); a_.set_title(lab, loc='left', fontsize=10, color=ink); a_.grid(alpha=0.2); a_.spines[['top', 'right']].set_visible(False)
ax[0].set_ylabel('normalised dN/dM'); ax[0].legend(fontsize=7, frameon=False)
fn = f'{OUT}/figs/clouds_mf.png'; fig.savefig(fn, dpi=140, bbox_inches='tight'); os.system(f'cp {fn} ~/out/'); print('fig', fn)
