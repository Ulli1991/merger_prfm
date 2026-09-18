"""Field strength per snapshot from the layer columns (patch_XXX_z05.h5 in C.PRFM_DIR): B_c = sqrt(8 pi k_B P_mag,c) with P_mag = B^2/8pi
over all slab particles; rows = snap, t, Sigma-weighted mean B [muG], median B, median P_mag/(P_th+P_turb), median Maxwell-stress share.
out: {C.DATADIR}/prfm/bfield_time.npz (rows)"""
import sys, os, glob, numpy as np, h5py
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); import common as C
KB = 1.380649e-16; rows = []
for fn in sorted(glob.glob(f'{C.PRFM_DIR}/patch_[0-9][0-9][0-9]_z05.h5')):
    k = int(os.path.basename(fn)[6:9]); B = []; S = []; r1 = []; r2 = []
    with h5py.File(fn, 'r') as f:
        t = float(f.attrs['time_myr'])
        for g in f:
            if not g.startswith('frame_'): continue
            G = f[g]; sep = float(G.attrs['separation_kpc']); r = lambda kk: G[kk][:].ravel()
            Sig = r('Sigma_gas'); W = r('W_2p'); fi = r('f_intruder'); ok = (Sig > 1) & (W > 0) & ((fi < 0.1) | C.is_merged(sep, t))
            Pm = r('Pmag_tot')[ok]; B.append(np.sqrt(8 * np.pi * KB * np.maximum(Pm, 0)) * 1e6); S.append(Sig[ok])
            r1.append(Pm / np.maximum(r('Pth_2p')[ok] + r('Pturb_2p')[ok], 1e-30)); r2.append(r('Pmag_2p')[ok] / np.maximum(r('Ptot_2p')[ok], 1e-30))
    if not B or sum(len(b) for b in B) < 3: continue
    B = np.concatenate(B); S = np.concatenate(S); r1 = np.concatenate(r1); r2 = np.concatenate(r2)
    rows.append((k, t, np.average(B, weights=S), np.median(B), np.median(r1), np.median(r2)))
R = np.array(rows); np.savez(f'{C.DATADIR}/prfm/bfield_time.npz', rows=R); print(f'{len(R)} snapshots -> bfield_time.npz')
