"""Sensitivity of the PRFM comparison to the column size and to the projection.
Reads every available patch_KKK<tag>.h5 for the configurations below and computes per snapshot, over the clean columns with
Sigma_gas > 1 (intruder cut before the second passage): the Sigma-weighted median and the ratio of sums of P_tot,2p / W_2p, the
effective yield P / Sigma_SFR (10 and 40 Myr) against OK22, and the Sigma_SFR(W) relation (log-log slope and offset from OK22 eq 28b
over the star-forming columns).  Prints per-phase medians and saves scale_los.npz for the figures.
usage: python scale_los.py"""
import sys, os, glob, numpy as np, h5py
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C, ok22
CFG = {  # tag: (label, patch kpc, column half-depth kpc, LOS suffixes to use)
    '': ('0.5 x 3.0 kpc (full column)', 0.5, 1.5, ['']), '_z05': ('0.5 x 1.0 kpc (layer)', 0.5, 0.5, ['']),
    '_p05z025': ('0.5 x 0.5 kpc', 0.5, 0.25, ['']), '_p025z025': ('0.25 x 0.5 kpc', 0.25, 0.25, ['']),
    '_p025z0125': ('0.25 x 0.25 kpc', 0.25, 0.125, ['']), '_p0125z0125': ('0.125 x 0.25 kpc', 0.125, 0.125, ['']),
    '_los': ('0.5 x 3.0 kpc per LOS', 0.5, 1.5, ['', '_x', '_y', '_z', '_t30', '_t60', '_t90']), '_los_z05': ('0.5 x 1.0 kpc layer per LOS', 0.5, 0.5, ['', '_x', '_y', '_z', '_t30', '_t60', '_t90'])}
GAS_SPLIT = 14_500_000
def per_snapshot(fn, suffix):
    with h5py.File(fn, 'r') as f:
        t = float(f.attrs['time_myr']); rows = []
        for g in f:
            if not g.startswith('frame_'): continue
            base = g[6:]; lab = base[0]; suf = base[1:]
            if suf != suffix: continue
            G = f[g]; sep = float(G.attrs['separation_kpc']); r = lambda k: G[k][:].ravel()
            Sig = r('Sigma_gas'); W = r('W_2p'); P = r('Ptot_2p'); fi = r('f_intruder'); s10 = r('SigSFR_10'); s40 = r('SigSFR_40'); S2 = r('Sigma_gas_2p')
            ok = (Sig > 1) & (W > 0) & np.isfinite(P) & (P > 0) & ((fi < 0.1) | C.is_merged(sep, t))
            rows.append((Sig[ok], W[ok], P[ok], s10[ok], s40[ok], S2[ok]))
    if not rows: return None
    Sig, W, P, s10, s40, S2 = [np.concatenate([r[i] for r in rows]) for i in range(6)]
    if len(Sig) < 5: return None
    x = P / W; o = np.argsort(x); cw = np.cumsum(Sig[o]) / Sig.sum(); med = x[o][np.searchsorted(cw, 0.5)]
    ups10 = P.sum() / max(s10.sum(), 1e-30) / 4.81e3; ups40 = P.sum() / max(s40.sum(), 1e-30) / 4.81e3; upsok = ok22.ups_tot(np.average(P, weights=S2))
    sf = s10 > 0
    if sf.sum() >= 5:
        lw = np.log10(W[sf]); ls = np.log10(s10[sf]); sl = np.polyfit(lw, ls, 1)[0]; off = np.median(ls - np.log10(ok22.sfr_of_W(W[sf])))
    else: sl = off = np.nan
    return dict(t=t, n=len(Sig), med=med, ratio=P.sum() / W.sum(), ups10=ups10, ups40=ups40, upsok=upsok, sfr_slope=sl, sfr_off=off, fsf=sf.mean())
out = {}
for tag, (label, pk, zc, sufs) in CFG.items():
    files = sorted(glob.glob(f'{C.PRFM_DIR}/patch_[0-9][0-9][0-9]{tag}.h5'))
    if not files: continue
    for suf in sufs:
        res = [per_snapshot(fn, suf) for fn in files]; res = [r for r in res if r is not None]
        if not res: continue
        key = tag + ('|' + suf if suf else ''); T = np.array([r['t'] for r in res]); out[key] = {k: np.array([r[k] for r in res]) for k in res[0]}
        print(f'\n{label}{(" normal " + suf) if suf else ""}: {len(res)} snapshots')
        print('  phase                        n_col  P/W med   P/W sums   Ups10/OK22   Ups40/OK22   SFR(W) slope   offset from OK22 [dex]   f_SF')
        for t0, t1, lab in C.PHASES:
            m = (T >= t0) & (T < t1) & (T > 10)
            if m.sum() == 0: continue
            g = lambda k: np.nanmedian(out[key][k][m])
            print(f'  {lab:28s} {g("n"):5.0f}  {g("med"):6.2f}    {g("ratio"):6.2f}     {np.nanmedian((out[key]["ups10"] / out[key]["upsok"])[m]):7.1f}      {np.nanmedian((out[key]["ups40"] / out[key]["upsok"])[m]):7.1f}       {g("sfr_slope"):+5.2f}          {g("sfr_off"):+5.2f}               {g("fsf"):.2f}')
np.savez(f'{C.DATADIR}/prfm/scale_los.npz', keys=np.array(list(out)), **{k.replace('|', '__'): np.array([out[k][q] for q in ('t', 'n', 'med', 'ratio', 'ups10', 'ups40', 'upsok', 'sfr_slope', 'sfr_off', 'fsf')]) for k in out})
print('\nsaved scale_los.npz')
