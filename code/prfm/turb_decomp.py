"""Split the vertical velocity of the two-phase gas in every 0.5 kpc patch (same slab as patches.py) into
a bulk linear part  v_n = v0 + a x + b y + c z  and a residual, mass-weighted.
sig_tot^2 = total variance (what enters P_turb), sig_res^2 = residual after the linear fit (turbulence proper),
sig_lin^2 = sig_tot^2 - sig_res^2 (tilt/warp a,b and layer compression/expansion c).  Same for the full 3D velocity.
usage: python turb_decomp.py SNAP [SNAP ...]   ->  prfm/turb_{snap:03d}.npz (one row per patch and frame)"""
import sys, os, numpy as np, h5py
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common as C
T_2P = 2.0e4; NMIN = 50
def run(k):
    g = C.read(k, 0, ['Coordinates', 'Velocities', 'Masses', 'InternalEnergy', 'ElectronAbundance'])
    pos = g['Coordinates'].astype(np.float64); vel = g['Velocities'].astype(np.float64); m = g['Masses'].astype(np.float64)
    T = C.temperature(g['InternalEnergy'], g['ElectronAbundance']); twop = T < T_2P
    rows = []
    with h5py.File(f'{C.PRFM_DIR}/patch_{k:03d}.h5', 'r') as f:
        t = f.attrs['time_myr']; PATCH = f.attrs['patch_kpc']; RMAX = f.attrs['rmax_kpc']; nb = int(f.attrs['nb']); ZCOL = f.attrs['zcol_kpc']
        for grp in f:
            G = f[grp]; a = dict(G.attrs); e1, e2, n = a['e1'], a['e2'], a['nhat']; sep = a['separation_kpc']
            zmid = G['zmid'][:].ravel(); hs = G['slab_half'][:].ravel()
            P = {key: G[key][:].ravel() for key in ['Sigma_gas', 'W_2p', 'Ptot_2p', 'Pturb_2p', 'f_intruder', 'SigSFR_10', 'div_v', 'shear']}
            d = pos - a['center']; x = d @ e1; y = d @ e2; z = d @ n; v = vel - a['vcen']; vn = v @ n
            ix = np.floor((x + RMAX) / PATCH).astype(int); iy = np.floor((y + RMAX) / PATCH).astype(int)
            ok = (ix >= 0) & (ix < nb) & (iy >= 0) & (iy < nb) & (np.abs(z) < ZCOL) & twop
            pid = np.where(ok, ix * nb + iy, -1)
            ok &= np.abs(z - zmid[np.clip(pid, 0, nb * nb - 1)]) < hs[np.clip(pid, 0, nb * nb - 1)]
            order = np.argsort(pid[ok]); pids = pid[ok][order]; idx = np.flatnonzero(ok)[order]
            starts = np.searchsorted(pids, np.arange(nb * nb)); ends = np.searchsorted(pids, np.arange(nb * nb), 'right')
            for p in range(nb * nb):
                sl = idx[starts[p]:ends[p]]
                if len(sl) < NMIN or P['Sigma_gas'][p] < 1: continue
                w = m[sl]; W = w.sum(); xc = ((p // nb) + 0.5) * PATCH - RMAX; yc = ((p % nb) + 0.5) * PATCH - RMAX
                A = np.column_stack([np.ones(len(sl)), x[sl] - xc, y[sl] - yc, z[sl] - zmid[p]])
                sw = np.sqrt(w)[:, None]
                def fit(yv):
                    coef = np.linalg.lstsq(A * sw, yv * sw[:, 0], rcond=None)[0]; r = yv - A @ coef
                    mu = np.sum(w * yv) / W
                    return np.sum(w * (yv - mu) ** 2) / W, np.sum(w * r ** 2) / W, coef
                st, sr, coef = fit(vn[sl])
                s3t = 0.0; s3r = 0.0; Gm = np.zeros((3, 3))     # Gm[i, j] = d U_i / d x_j of the patch-scale mean flow (frame axes e1, e2, n)
                vf = np.column_stack([v[sl] @ e1, v[sl] @ e2, vn[sl]])
                for c in range(3):
                    a_, b_, cf = fit(vf[:, c]); s3t += a_; s3r += b_; Gm[c] = cf[1:4]
                Sm = 0.5 * (Gm + Gm.T) - np.trace(Gm) / 3 * np.eye(3); shear_patch = np.sqrt(np.sum(Sm ** 2)); div_patch = np.trace(Gm)
                Om = 0.5 * (Gm - Gm.T); vort_patch = np.sqrt(2 * np.sum(Om ** 2))
                # two-body check: bimodality of vn via mass fraction in the minority mode of a 2-means split
                vv = vn[sl]; med = np.sum(w * vv) / W
                hi = vv > med; mh = w[hi].sum() / W; sep_v = abs(np.sum(w[hi] * vv[hi]) / max(w[hi].sum(), 1e-30) - np.sum(w[~hi] * vv[~hi]) / max(w[~hi].sum(), 1e-30))
                rows.append((k, t, sep, grp, p, len(sl), P['Sigma_gas'][p], P['W_2p'][p], P['Ptot_2p'][p], P['Pturb_2p'][p], P['f_intruder'][p], P['SigSFR_10'][p],
                             P['div_v'][p], P['shear'][p], np.sqrt(st), np.sqrt(sr), coef[1], coef[2], coef[3], np.sqrt(s3t), np.sqrt(s3r), sep_v, min(mh, 1 - mh), shear_patch, div_patch, vort_patch, Gm[0, 1], Gm[1, 0], Gm[0, 2], Gm[1, 2]))
    names = 'snap t sep frame patch n Sigma_gas W_2p Ptot_2p Pturb_2p f_intruder SigSFR_10 div_v shear sig_tot sig_res dvdx dvdy dvdz sig3_tot sig3_res dv_halves f_minor shear_patch div_patch vort_patch dUxdy dUydx dUxdz dUydz'.split()
    out = {nm: np.array([r[i] for r in rows]) for i, nm in enumerate(names)}
    np.savez(f'{C.DATADIR}/prfm/turb_{k:03d}.npz', **out)
    st = out['sig_tot']; sr = out['sig_res']
    print(f'snap {k} t={t:.1f} sep={sep:.2f}: {len(rows)} patches; median sig_z tot {np.median(st):.1f} res {np.median(sr):.1f} km/s; '
          f'median residual fraction of variance {np.median(sr**2/np.maximum(st**2,1e-30)):.2f}', flush=True)
if __name__ == '__main__':
    for a in sys.argv[1:]: run(int(a))
