"""Every bound young cluster (0-10 Myr catalogue, every 10th snapshot, nucleus excluded) with the state of the 0.5 kpc patch it formed
in: layer weight W_L and layer pressure P_L (patch_KKK_z05.h5), full-column W, Sigma_SFR over 10 and 40 Myr, Sigma_gas, rho_mid, sigma_eff,
f_intruder.  Home frame = the frame in which the cluster's patch has the lowest intruder fraction (same rule as the fixed join).
out: /ptmp/uli/dwarf_merger/clusters/clusters_env_layer.npz; prints Spearman of log M with each quantity and the upper envelope;
figure via paper_figs.cluster_env."""
import sys, os, glob, numpy as np, h5py
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C, ok22
from scipy.stats import spearmanr
cl = np.load(f'{C.DATADIR}/clusters/clusters_age0-10_l5_n25.npz')
KEYS = ['W_2p', 'Ptot_2p', 'SigSFR_10', 'SigSFR_40', 'Sigma_gas', 'Sigma_gas_2p', 'rho_mid_2p', 'sigma_eff', 'f_intruder', 'H', 'Pturb_2p', 'Pmag_tot', 'Pth_2p']
rows = []
for k in range(0, C.NSNAP, 10):
    cs = np.flatnonzero((cl['snap'] == k) & (cl['bound'] > 0) & (cl['nucleus'] == 0))
    if len(cs) == 0 or not os.path.exists(f'{C.PRFM_DIR}/patch_{k:03d}_z05.h5'): continue
    pos = np.stack([cl['x'][cs], cl['y'][cs], cl['z'][cs]], 1)
    best = {}
    with h5py.File(f'{C.PRFM_DIR}/patch_{k:03d}_z05.h5', 'r') as fz, h5py.File(f'{C.PRFM_DIR}/patch_{k:03d}.h5', 'r') as ff:
        PATCH = float(fz.attrs['patch_kpc']); RMAX = float(fz.attrs['rmax_kpc']); ZCOL = float(fz.attrs['zcol_kpc']); nb = int(fz.attrs['nb'])
        for gname in fz:
            if not gname.startswith('frame_'): continue
            g = fz[gname]; gf = ff[gname]; cen = g.attrs['center']; e1 = g.attrs['e1']; e2 = g.attrs['e2']; nh = g.attrs['nhat']; sep = float(g.attrs['separation_kpc'])
            d = pos - cen; x = d @ e1; y = d @ e2; z = d @ nh
            ix = np.floor((x + RMAX) / PATCH).astype(int); iy = np.floor((y + RMAX) / PATCH).astype(int)
            ok = (ix >= 0) & (ix < nb) & (iy >= 0) & (iy < nb) & (np.abs(z) < 1.5)      # inside the grid and within the full column
            vz = {kk: g[kk][:].ravel() for kk in ('W_2p', 'Ptot_2p', 'Sigma_gas_2p', 'Pturb_2p', 'Pmag_tot', 'Pth_2p')}
            vf = {kk: gf[kk][:].ravel() for kk in ('W_2p', 'SigSFR_10', 'SigSFR_40', 'Sigma_gas', 'rho_mid_2p', 'sigma_eff', 'f_intruder', 'H')}
            for j in np.flatnonzero(ok):
                p = ix[j] * nb + iy[j]; fi = vf['f_intruder'][p]
                if j not in best or fi < best[j][0]:
                    best[j] = (fi, ord(gname[-1]), sep, vz['W_2p'][p], vz['Ptot_2p'][p], vf['W_2p'][p], vf['SigSFR_10'][p], vf['SigSFR_40'][p], vf['Sigma_gas'][p], vz['Sigma_gas_2p'][p], vf['rho_mid_2p'][p], vf['sigma_eff'][p], vf['H'][p], vz['Pturb_2p'][p], vz['Pmag_tot'][p], vz['Pth_2p'][p], abs(z[j]))
    for j, b in best.items():
        i = cs[j]; rows.append((k, cl['t'][i], cl['M'][i], cl['N'][i], cl['rh_pc'][i]) + b)
names = ['snap', 't', 'M', 'N', 'rh_pc', 'f_intruder', 'frame', 'sep', 'W_L', 'P_L', 'W_full', 'SigSFR_10', 'SigSFR_40', 'Sigma_gas', 'Sigma_gas_L', 'rho_mid', 'sigma_eff', 'H', 'Pturb_L', 'Pmag_L', 'Pth_L', 'zabs']
R = np.array(rows); col = {n: i for i, n in enumerate(names)}
clean = (R[:, col['f_intruder']] < 0.1) | C.is_merged(R[:, col['sep']], R[:, col['t']]); R = R[clean & (R[:, col['t']] > 10)]
np.savez(f'{C.DATADIR}/clusters/clusters_env_layer.npz', rows=R, names=np.array(names))
print(f'{len(R)} bound young clusters (every 10th snapshot, clean, t>10) with their patch state')
EP = [*C.PHASES, (25, 232, 'all')]
lM = np.log10(R[:, col['M']])
print('\nSpearman of log M_cluster with the patch quantity (per phase):')
print('phase                 n     W_L    P_L    SigSFR_10  SigSFR_40  Sigma_gas  rho_mid  sigma_eff  P_L/W_L')
for t0, t1, lab in EP:
    m = (R[:, col['t']] >= t0) & (R[:, col['t']] < t1)
    r = lambda key: spearmanr(lM[m], np.log10(np.maximum(R[m, col[key]], 1e-12)))[0]
    print(f'{lab:20s} {m.sum():4d}   {r("W_L"):+.2f}  {r("P_L"):+.2f}   {r("SigSFR_10"):+.2f}      {r("SigSFR_40"):+.2f}      {r("Sigma_gas"):+.2f}     {r("rho_mid"):+.2f}    {r("sigma_eff"):+.2f}     {spearmanr(lM[m], np.log10(R[m, col["P_L"]] / R[m, col["W_L"]]))[0]:+.2f}')
print('\nupper envelope: fraction of clusters above 0.5 x PRFM burst(W_L) [= 0.5 Sigma_SFR(W_L) A tau] and above 0.5 x M_young of the patch (= 0.5 SigSFR_10 A 10 Myr):')
A = 0.25; tau = 1e7
prfm_b = 0.5 * ok22.sfr_of_W(R[:, col['W_L']]) * A * tau; myoung = R[:, col['SigSFR_10']] * A * tau
for t0, t1, lab in EP:
    m = (R[:, col['t']] >= t0) & (R[:, col['t']] < t1)
    print(f'{lab:20s} above 0.5 PRFM burst: {np.mean(R[m, col["M"]] > prfm_b[m]):.2f}   above 0.5 M_young: {np.mean(R[m, col["M"]] > 0.5 * myoung[m]):.2f}   median M/M_young: {np.median(R[m, col["M"]] / np.maximum(myoung[m], 1)):.2f}   median M / (0.5 PRFM burst): {np.median(R[m, col["M"]] / prfm_b[m]):.3f}')
