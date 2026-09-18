"""Low-mass end of the cluster MF per phase: convergence with N_min (10/25/50 particles), FoF linking length (3/5/8 pc,
where catalogues exist) and the bound flag.  Local slope = MLE of a power law truncated to [M1, M2] (no fixed upper end),
so bins can be compared directly.  Snapshots every 10th (independent populations), nucleus excluded.
usage: python lowmass_convergence.py      out: /ptmp/uli/dwarf_merger/clusters/lowmass_convergence.npz, figs/lowmass_convergence.png"""
import sys, os, glob, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common as C
from scipy.optimize import brentq
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
OUT = f'{C.DATADIR}/clusters'
EP = [*C.PHASES]
BINS = [(100, 300), (300, 1000), (1000, 3000), (3000, 10000), (300, 1e7)]

def local_slope(M, lo, hi):
    """MLE alpha for dN/dM ~ M^-alpha on [lo, hi)."""
    M = M[(M >= lo) & (M < hi)]; n = len(M)
    if n < 8: return np.nan, np.nan, n
    L = np.log(M / lo); r = hi / lo
    if hi > 1e6:
        al = 1 + n / L.sum(); return al, (al - 1) / np.sqrt(n), n
    # d lnL / d alpha = 0 :  -sum ln(M/lo) + n * [ ln r * r^(1-a) / (1 - r^(1-a)) - 1/(a-1) ] ... solve numerically
    def dl(a):
        if abs(a - 1) < 1e-6: a = 1 + 1e-6
        q = r ** (1 - a)
        return -L.sum() + n * (np.log(r) * q / (1 - q) - 1 / (a - 1)) * (-1) - n * 0  # derivative of n*ln((1-a)/(lo^(1-a)-hi^(1-a)))... see below
    # cleaner: maximise directly on a grid then refine
    aa = np.linspace(0.2, 4.0, 381)
    def ll(a):
        if abs(a - 1) < 1e-6: a = 1 + 1e-6
        Z = (lo ** (1 - a) - hi ** (1 - a)) / (a - 1)
        return -a * np.log(M).sum() - n * np.log(Z)
    v = np.array([ll(a) for a in aa]); i = np.argmax(v); al = aa[i]
    # curvature error
    h = 0.02; c = (ll(al + h) - 2 * ll(al) + ll(al - h)) / h ** 2
    err = 1 / np.sqrt(-c) if c < 0 else np.nan
    return al, err, n

cats = {}
for fn in sorted(glob.glob(f'{OUT}/clusters_age0-10_l*_n*.npz')):
    d = np.load(fn); tag = os.path.basename(fn)[len('clusters_age0-10_'):-4]
    sel = (np.mod(d['snap'], 10) == 0) & (d['nucleus'] == 0)
    cats[tag] = dict(M=d['M'][sel], t=d['t'][sel], bound=d['bound'][sel] > 0, snap=d['snap'][sel])
    print(f'{tag}: {sel.sum()} groups on every-10th snapshots (of {len(d["M"])} total), M_min = {d["M"][sel].min():.0f}')
# clean flag for the n25 catalogue (clusters_env is row-aligned with clusters_age0-10_l5_n25)
E = np.load(f'{OUT}/clusters_env.npz'); d = np.load(f'{OUT}/clusters_age0-10_l5_n25.npz'); assert len(E['M']) == len(d['M'])
sel = (np.mod(d['snap'], 10) == 0) & (d['nucleus'] == 0)
cleanE = ((E['f_intruder'] < 0.1) | C.is_merged(E['sep'], E['t']))[sel]
cats['l5_n25_clean'] = dict(M=d['M'][sel][cleanE], t=d['t'][sel][cleanE], bound=(d['bound'][sel] > 0)[cleanE], snap=d['snap'][sel][cleanE])

order = ['l5_n10', 'l5_n25', 'l5_n25_clean', 'l5_n50'] + [k for k in cats if k not in ('l5_n10', 'l5_n25', 'l5_n25_clean', 'l5_n50')]
res = {}
for t0, t1, lab in EP:
    print(f'\n=== {lab} ({t0}-{t1} Myr) ===')
    print('catalogue        sample  ' + '  '.join(f'a[{lo:g}-{hi:g}]'.ljust(16) if hi < 1e6 else 'a(>300)'.ljust(16) for lo, hi in BINS) + '  f_bound[100-300] f_bound[300-1e3] f_bound[>1e3]')
    for tag in order:
        if tag not in cats: continue
        c = cats[tag]; m = (c['t'] >= t0) & (c['t'] < t1)
        for sname, sm in (('all', m), ('bound', m & c['bound'])):
            row = [local_slope(c['M'][sm], lo, hi) for lo, hi in BINS]
            fb = [c['bound'][m & (c['M'] >= lo) & (c['M'] < hi)].mean() if (m & (c['M'] >= lo) & (c['M'] < hi)).sum() > 0 else np.nan for lo, hi in ((100, 300), (300, 1000), (1000, 1e9))]
            res[(lab, tag, sname)] = (row, fb)
            print(f'{tag:16s} {sname:6s}  ' + '  '.join(f'{a:4.2f}+-{e:4.2f} ({n:3d})' if np.isfinite(a) else f'   --       ({n:3d})' for a, e, n in row) + '  ' + '  '.join(f'{x:5.2f}          ' for x in fb))
np.savez(f'{OUT}/lowmass_convergence.npz', keys=np.array([str(k) for k in res]), rows=np.array([np.array([r[0] for r in v[0]] + v[1]) for v in res.values()]))

# figure: local slope vs mass bin per phase, all vs bound, for the three N_min
ink = '#1f2937'; blue = '#2f6fb2'; light = '#9dbfe3'; orange = '#c2410c'
fig, ax = plt.subplots(1, 4, figsize=(15, 3.8), sharey=True)
x = np.array([np.sqrt(lo * hi) for lo, hi in BINS[:-1]])
for a_, (t0, t1, lab) in zip(ax, EP):
    for tag, col, ls in (('l5_n10', orange, '-'), ('l5_n25', blue, '-'), ('l5_n50', ink, '-')):
        if (lab, tag, 'all') not in res: continue
        r = res[(lab, tag, 'all')][0]; a_.errorbar(x * (1 + 0.08 * (['l5_n10', 'l5_n25', 'l5_n50'].index(tag) - 1)), [q[0] for q in r[:-1]], yerr=[q[1] for q in r[:-1]], color=col, ls=ls, marker='o', ms=4, lw=1.4, label=f'{tag[3:]} all')
        rb = res[(lab, tag, 'bound')][0]; a_.plot(x * (1 + 0.08 * (['l5_n10', 'l5_n25', 'l5_n50'].index(tag) - 1)), [q[0] for q in rb[:-1]], color=col, ls=':', marker='s', ms=4, mfc='none', lw=1.2, label=f'{tag[3:]} bound')
    a_.set_xscale('log'); a_.set_title(lab, loc='left', fontsize=10, color=ink); a_.set_xlabel('M [M$_\\odot$]'); a_.grid(alpha=0.2); a_.spines[['top', 'right']].set_visible(False); a_.axhline(2, color=ink, lw=0.6, ls='--')
ax[0].set_ylabel('local slope $\\alpha$'); ax[0].legend(fontsize=7, frameon=False, ncol=2); ax[0].set_ylim(0.5, 3.2)
os.makedirs(f'{OUT}/figs', exist_ok=True); fn = f'{OUT}/figs/lowmass_convergence.png'; fig.savefig(fn, dpi=140, bbox_inches='tight'); os.system(f'cp {fn} ~/out/'); print('fig', fn)
