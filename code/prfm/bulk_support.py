"""After turb_decomp: (1) how much of the 'turbulent' vertical pressure is bulk (linear) flow vs residual turbulence, by phase and density;
(2) equilibrium with the bulk part removed; (3) is the weight deficit W - P_tot located in cells with large bulk vertical motion,
large layer compression |dv_n/dz|, or large patch-scale shear; (4) shear-extended closure calibrated with the PATCH-SCALE shear:
c_s = P_turb,res / (Sigma H S_patch^2) on quiescent cells, per phase.
usage: python bulk_support.py   (reads prfm/turb_*.npz and the matching patch_XXX.h5)"""
import sys, os, glob, numpy as np, h5py
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C, ok22
from scipy.stats import spearmanr
PU = 4.903245584337325e-06
rows = []
for fn in sorted(glob.glob(f'{C.DATADIR}/prfm/turb_*.npz')):
    d = dict(np.load(fn)); k = int(d['snap'][0])
    with h5py.File(f'{C.PRFM_DIR}/patch_{k:03d}.h5', 'r') as f:
        for fr in np.unique(d['frame']):
            g = f[fr]; m = d['frame'] == fr; p = d['patch'][m].astype(int)
            H = g['H'][:].ravel()[p]; rho = g['rho_mid_2p'][:].ravel()[p] * 1e9; Pm = g['Pmag_2p'][:].ravel()[p]
            for j, i in enumerate(np.flatnonzero(m)):
                rows.append([d['t'][i], d['sep'][i], d['Sigma_gas'][i], d['W_2p'][i], d['Ptot_2p'][i], d['Pturb_2p'][i], d['f_intruder'][i], d['SigSFR_10'][i],
                             d['sig_tot'][i], d['sig_res'][i], d['sig3_tot'][i], d['sig3_res'][i], d['dvdz'][i], d['shear_patch'][i], d['div_patch'][i], H[j], rho[j], Pm[j], d['n'][i]])
R = dict(zip('t sep Sig W P Pt fi sfr st sr s3t s3r dvdz Sp divp H rho Pm n'.split(), np.array(rows).T))
clean = (R['fi'] < 0.1) | C.is_merged(R['sep'], R['t']); ok = clean & (R['W'] > 0) & (R['P'] > 0) & (R['st'] > 0) & (R['t'] > 10)
for k in R: R[k] = R[k][ok]
fres = (R['sr'] / R['st'])**2                       # residual (turbulent) fraction of the vertical variance
Pres = R['Pt'] * fres; Pbulk = R['Pt'] - Pres        # split of the measured 'turbulent' pressure
EP = [*C.PHASES]
print(f'{len(R["t"])} clean cells from {len(np.unique(R["t"]))} snapshots\n')
print('(1) residual (turbulent) fraction of the vertical velocity variance, median, by phase and density class; and of the 3D variance')
print('phase                 Sig 1-10        Sig 10-30       Sig >30        | 3D all')
for t0, t1, lab in EP:
    m = (R['t'] >= t0) & (R['t'] < t1); out = []
    for lo, hi in [(1, 10), (10, 30), (30, 1e9)]:
        q = m & (R['Sig'] >= lo) & (R['Sig'] < hi); out.append(f'{np.median(fres[q]):.2f} (n={q.sum():4d})' if q.sum() >= 10 else '    -         ')
    print(f'{lab:20s} ' + '  '.join(out) + f'  | {np.median((R["s3r"][m]/R["s3t"][m])**2):.2f}')
print('\n(2) equilibrium: median P/W with the measured pressure, and with the bulk part of the vertical variance removed (turbulence proper only)')
print('phase                 Sig 1-10: P/W  (P-Pbulk)/W    Sig 10-30: P/W  (P-Pbulk)/W    Sig>30: P/W  (P-Pbulk)/W')
for t0, t1, lab in EP:
    m = (R['t'] >= t0) & (R['t'] < t1); out = []
    for lo, hi in [(1, 10), (10, 30), (30, 1e9)]:
        q = m & (R['Sig'] >= lo) & (R['Sig'] < hi)
        out.append(f'{np.median(R["P"][q]/R["W"][q]):.2f}  {np.median((R["P"][q]-Pbulk[q])/R["W"][q]):.2f}        ' if q.sum() >= 10 else '  -      -           ')
    print(f'{lab:20s} ' + '  '.join(out))
print('\n(3) where is the deficit?  Spearman of log(W/P) with: bulk vertical dispersion sig_bulk, compression rate -dv_n/dz, patch shear S_p, Sigma  (all cells; then diffuse only)')
sb = np.sqrt(np.maximum(R['st']**2 - R['sr']**2, 0))
for t0, t1, lab in EP:
    m = (R['t'] >= t0) & (R['t'] < t1); dfs = m & (R['Sig'] < 10); y = np.log10(R['W'] / R['P'])
    c = [spearmanr(y[m], v[m])[0] for v in (sb, -R['dvdz'], R['Sp'], R['Sig'])]; cd = [spearmanr(y[dfs], v[dfs])[0] for v in (sb, -R['dvdz'], R['Sp'], R['Sig'])]
    print(f'{lab:20s} all: ' + '  '.join(f'{x:+.2f}' for x in c) + '   | diffuse: ' + '  '.join(f'{x:+.2f}' for x in cd) + f'   (n={m.sum()}, {dfs.sum()})')
print('   median -dv_n/dz [km/s/kpc] (positive = layer compressing) by phase, diffuse / dense: ' + ', '.join(
    f'{lab}: {np.median(-R["dvdz"][(R["t"]>=t0)&(R["t"]<t1)&(R["Sig"]<10)]):+.1f} / {np.median(-R["dvdz"][(R["t"]>=t0)&(R["t"]<t1)&(R["Sig"]>=10)]):+.1f}' for t0, t1, lab in EP))
print('\n(4) shear-extended closure with the patch-scale shear: quiescent cells (SFR=0): c_s = P_turb,res / (Sigma H S_p^2), log slope, sigma_res/(H S_p)')
X = R['Sig'] * 1e6 * R['H'] * R['Sp']**2 * PU
for t0, t1, lab in EP:
    q = (R['t'] >= t0) & (R['t'] < t1) & (R['sfr'] == 0) & (X > 0) & (Pres > 0)
    if q.sum() < 20: print(lab, 'too few'); continue
    b, a = np.polyfit(np.log10(X[q]), np.log10(Pres[q]), 1); cs = Pres[q] / X[q]; r = R['sr'][q] / (R['H'][q] * R['Sp'][q])
    print(f'{lab:20s} n={q.sum():4d}  c_s median {np.median(cs):6.2f} (16-84: {np.percentile(cs,16):.2f}-{np.percentile(cs,84):.2f})  slope {b:.2f}  sigma_res/(H S_p) {np.median(r):.2f}   rho_s(Pres, S_p | Sig H) = '
          f'{spearmanr(np.log10(Pres[q]) - np.log10(R["Sig"][q]*R["H"][q]), np.log10(R["Sp"][q]))[0]:+.2f}')
print('    all cells, feedback + flow: scatter of log P_turb,res about log(Ups_turb SFR + c_s Sigma H S_p^2) with c_s = 1, vs feedback-only (SF cells)')
Pfb = ok22.ups_turb(R['P']) * 1e5 * R['sfr'] * (1.989e33 / 3.156e7 / (3.0857e21)**2) / 1.3807e-16
for t0, t1, lab in EP:
    m = (R['t'] >= t0) & (R['t'] < t1) & (Pres > 0); sfc = m & (Pfb > 0)
    print(f'    {lab:20s} both: {np.std(np.log10(Pres[m]) - np.log10(Pfb[m] + X[m])):.2f} dex (offset {np.median(np.log10(Pres[m]) - np.log10(Pfb[m] + X[m])):+.2f})   feedback-only: {np.std(np.log10(Pres[sfc]) - np.log10(Pfb[sfc])):.2f} dex (offset {np.median(np.log10(Pres[sfc]) - np.log10(Pfb[sfc])):+.2f})')
