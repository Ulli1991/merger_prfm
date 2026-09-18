"""Frame-free cubes (cubes_KKK_l<L>.h5): per snapshot, the equilibrium test conditioned on whether the local gas is a layer.
flatness q = sqrt(lambda_min/lambda_mid) of the gas within 1 kpc (disc ~0.1-0.2, blob ~0.5-1).  Per snapshot: Sigma-weighted median
P/W over all cubes, over 'layer' cubes (q < Q_LAYER), the mass fraction of gas in layer cubes, the median angle between the local
normal and the galaxy-frame normal, and the yield ratio.  usage: python cubes_analysis.py 0.5   (cube side)"""
import sys, os, glob, numpy as np, h5py
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C, ok22
L = sys.argv[1] if len(sys.argv) > 1 else '0.5'; Q_LAYER = 0.3
rows = []
for fn in sorted(glob.glob(f'{C.DATADIR}/prfm/cubes_[0-9][0-9][0-9]_l{L}.h5')):
    with h5py.File(fn, 'r') as f:
        t = float(f.attrs['time_myr']); r = lambda k: f[k][:]
        S = r('Sigma'); W = r('W_2p'); P = r('Ptot_2p'); q = r('q'); ang = r('ang_ngal'); s10 = r('SigSFR_10'); S2 = r('Sigma_2p')
    ok = (S > 1) & (W > 0) & (P > 0) & np.isfinite(q)
    if ok.sum() < 5: continue
    def wmed(x, w):
        o = np.argsort(x); cw = np.cumsum(w[o]) / w.sum(); return x[o][np.searchsorted(cw, 0.5)]
    lay = ok & (q < Q_LAYER)
    rows.append((t, ok.sum(), wmed(P[ok] / W[ok], S[ok]), P[ok].sum() / W[ok].sum(), wmed(P[lay] / W[lay], S[lay]) if lay.sum() >= 5 else np.nan, S[lay].sum() / S[ok].sum(), np.median(q[ok]), np.median(ang[ok]),
                 P[ok].sum() / max(s10[ok].sum(), 1e-30) / 4.81e3 / ok22.ups_tot(np.average(P[ok], weights=S2[ok]))))
R = np.array(rows); names = ['t', 'n', 'pw_med', 'pw_sums', 'pw_layer', 'f_layer', 'q_med', 'ang_med', 'ups_ratio']
np.savez(f'{C.DATADIR}/prfm/cubes_analysis_l{L}.npz', rows=R, names=np.array(names))
print(f'cube side {L} kpc, {len(R)} snapshots; layer = flatness q < {Q_LAYER}')
print('phase                         n_cubes  P/W med  P/W sums  P/W layer-cubes  mass frac in layers  q med  angle to galaxy normal [deg]  Ups/OK22')
for t0, t1, lab in C.PHASES:
    m = (R[:, 0] >= t0) & (R[:, 0] < t1)
    if m.sum() == 0: continue
    g = lambda i: np.nanmedian(R[m, i])
    print(f'{lab:28s} {g(1):7.0f}   {g(2):5.2f}    {g(3):5.2f}      {g(4):5.2f}            {g(5):.2f}           {g(6):.2f}       {g(7):5.0f}                    {g(8):5.1f}')
