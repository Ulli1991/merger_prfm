"""Final PRFM equilibrium table: full-column weight (patch_XXX.h5, ZCOL 1.5) vs layer weight (patch_XXX_z05.h5, column |z - zmid| < 0.5 kpc)
per phase and density class: median P/W, mass-weighted P/W, gas mass within x2 and within 25%; plus the yield ratio.
usage: python layer_verdict.py [tag]   (default _z05)   -> prints table, writes prfm/layer_verdict.npz and figs/layer_verdict.png (+ ~/out)"""
import sys, os, numpy as np, h5py
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C, ok22
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
tag = sys.argv[1] if len(sys.argv) > 1 else '_z05'
KB = 1.3807e-16; SFR_CGS = 1.989e33 / 3.156e7 / (3.0857e21)**2
def load(suffix):
    rows = []
    for k in range(C.NSNAP):
        fn = f'{C.PRFM_DIR}/patch_{k:03d}{suffix}.h5'
        if not os.path.exists(fn): continue
        with h5py.File(fn, 'r') as f:
            t = f.attrs['time_myr']
            for g in f:
                if '_' in g.replace('frame_', ''): continue
                G = f[g]; sep = G.attrs['separation_kpc']
                d = {key: G[key][:].ravel() for key in ['Sigma_gas', 'W_2p', 'Ptot_2p', 'Pturb_2p', 'Pmag_2p', 'SigSFR_10', 'f_intruder', 'H']}
                sel = (d['Sigma_gas'] > 1) & ((d['f_intruder'] < 0.1) | C.is_merged(sep, t)) & (d['W_2p'] > 0) & (d['Ptot_2p'] > 0) & (t > 10)
                for i in np.flatnonzero(sel):
                    rows.append([k, t, sep, g] + [d[key][i] for key in d])
    names = ['snap', 't', 'sep', 'frame', 'Sigma_gas', 'W_2p', 'Ptot_2p', 'Pturb_2p', 'Pmag_2p', 'SigSFR_10', 'f_intruder', 'H']
    return {nm: np.array([r[i] for r in rows]) for i, nm in enumerate(names)}
F = load(''); L = load(tag)
for R in (F, L):
    for k in ['t', 'sep', 'Sigma_gas', 'W_2p', 'Ptot_2p', 'Pturb_2p', 'Pmag_2p', 'SigSFR_10', 'H']: R[k] = R[k].astype(float)
snaps_common = sorted(set(np.unique(F['snap'].astype(int))) & set(np.unique(L['snap'].astype(int))))
print(f'full column: {len(F["t"])} patches from {len(np.unique(F["snap"]))} snapshots; layer{tag}: {len(L["t"])} patches from {len(np.unique(L["snap"]))} snapshots; common snapshots {len(snaps_common)}')
EP = [*C.PHASES]
def wmedian(x, w):
    o = np.argsort(x); cw = np.cumsum(w[o]) / w.sum(); return x[o][np.searchsorted(cw, 0.5)]
def stats(R, m):
    pw = R['Ptot_2p'][m] / R['W_2p'][m]; S = R['Sigma_gas'][m]
    # 'mw' = Sigma-weighted median of P/W (robust to the ~10% over-pressured patches); 'rs' = sum P / sum W
    return dict(med=np.median(pw), mw=wmedian(pw, S), rs=R['Ptot_2p'][m].sum() / R['W_2p'][m].sum(), f2=np.sum(S[(pw > 0.5) & (pw < 2)]) / S.sum(), f25=np.sum(S[(pw > 0.8) & (pw < 1.25)]) / S.sum(), n=m.sum())
print('\nP_tot/W by phase: full column (ZCOL 1.5 kpc) vs layer (|z - zmid| < 0.5 kpc); restricted to common snapshots')
print('phase                 | full: median  Sig-wtd-med  sumP/sumW  f(x2)  f(25%)  | layer: median  Sig-wtd-med  sumP/sumW  f(x2)  f(25%)  | layer by class (median): Sig1-10  10-30  >30')
out = {}
for t0, t1, lab in EP:
    mF = (F['t'] >= t0) & (F['t'] < t1) & np.isin(F['snap'].astype(int), snaps_common); mL = (L['t'] >= t0) & (L['t'] < t1) & np.isin(L['snap'].astype(int), snaps_common)
    if mL.sum() < 10: print(lab, 'no layer data yet'); continue
    a = stats(F, mF); b = stats(L, mL)
    cls = [np.median((L['Ptot_2p'] / L['W_2p'])[mL & (L['Sigma_gas'] >= lo) & (L['Sigma_gas'] < hi)]) if (mL & (L['Sigma_gas'] >= lo) & (L['Sigma_gas'] < hi)).sum() >= 10 else np.nan for lo, hi in [(1, 10), (10, 30), (30, 1e9)]]
    print(f'{lab:20s} |  {a["med"]:5.2f}    {a["mw"]:5.2f}      {a["rs"]:5.2f}     {a["f2"]:.2f}   {a["f25"]:.2f}   |   {b["med"]:5.2f}    {b["mw"]:5.2f}      {b["rs"]:5.2f}     {b["f2"]:.2f}   {b["f25"]:.2f}   |   ' + '  '.join(f'{c:5.2f}' for c in cls) + f'   (n={a["n"]}, {b["n"]})')
    out[lab] = dict(full=a, layer=b, cls=cls)
# per-snapshot mass-weighted P/W, both, for the figure
def per_snap(R):
    s = np.unique(R['snap'].astype(int)); T = np.array([R['t'][R['snap'].astype(int) == k][0] for k in s])
    mw = np.array([wmedian((R['Ptot_2p'] / R['W_2p'])[R['snap'].astype(int) == k], R['Sigma_gas'][R['snap'].astype(int) == k]) for k in s])
    med = np.array([np.median((R['Ptot_2p'] / R['W_2p'])[R['snap'].astype(int) == k]) for k in s]); sep = np.array([R['sep'][R['snap'].astype(int) == k][0] for k in s])
    return T, mw, med, sep
TF, mwF, medF, sepF = per_snap(F); TL, mwL, medL, sepL = per_snap(L)
np.savez(f'{C.DATADIR}/prfm/layer_verdict{tag}.npz', TF=TF, mwF=mwF, medF=medF, TL=TL, mwL=mwL, medL=medL, sepF=sepF)
ink = '#1f2937'; blue = '#2f6fb2'; light = '#9dbfe3'
fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True, gridspec_kw=dict(hspace=0.1))
for a in ax: a.grid(alpha=0.2); a.spines[['top', 'right']].set_visible(False); a.axvline(C.T_MERGED, color=ink, lw=0.8, ls=':')
ax[0].plot(TF, sepF, color=blue, lw=2); ax[0].set_ylabel('separation [kpc]'); ax[0].set_title('Vertical equilibrium: full column vs gas layer', loc='left', color=ink)
ax[1].plot(TF, mwF, color=light, lw=2, label='full column (3 kpc), $\\Sigma$-weighted median'); ax[1].plot(TL, mwL, color=blue, lw=2, label='layer (|z| < 0.5 kpc), $\\Sigma$-weighted median')
ax[1].plot(TF, medF, color=light, lw=1.2, ls='--', label='full column, median patch'); ax[1].plot(TL, medL, color=blue, lw=1.2, ls='--', label='layer, median patch')
ax[1].axhline(1, color=ink, lw=0.8); ax[1].set_yscale('log'); ax[1].set_ylim(0.1, 3); ax[1].set_ylabel('$P_{\\rm tot}/\\mathcal{W}$'); ax[1].legend(frameon=False, fontsize=8, ncol=2); ax[1].set_xlabel('t [Myr]')
fn = f'{C.DATADIR}/prfm/figs/layer_verdict{tag}.png'; fig.savefig(fn, dpi=140, bbox_inches='tight'); os.system(f'cp {fn} ~/out/'); print('fig', fn)
