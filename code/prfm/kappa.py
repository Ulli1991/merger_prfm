"""Gas rotation curve, Omega(R) and epicyclic frequency kappa(R) per galaxy frame.

Uses the L1 grid of the movie extraction (half 4 kpc, rho 512^3 / vel 256^3,
mass-weighted mean velocity) and the frame definitions (centre, vcen, nhat)
stored in the PRFM patch files.  For every snapshot and frame the azimuthal
velocity is mass-averaged in annuli of DR kpc within |z| < ZMAX kpc, then
  Omega = v_phi / R,   kappa^2 = R dOmega^2/dR + 4 Omega^2
(evaluated on a smoothed v_phi(R)).  Also Q-related: sigma_R (mass-weighted
in-plane velocity dispersion about the mean rotation) per annulus.

usage: python kappa.py [SNAP ...]     (default: all patch files present)
output: /ptmp/uli/dwarf_merger/prfm/kappa.h5  (group snap_SSS/frame_X: R, vphi, Omega, kappa, sigR, Sigma)
"""
import sys, os, glob
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, h5py
from scipy.ndimage import gaussian_filter1d
import common as C

DR = 0.2; RMAX = 3.2; ZMAX = 0.5
edges = np.arange(0, RMAX + DR / 2, DR); Rc = 0.5 * (edges[1:] + edges[:-1])

def one(k):
    pf = f'{C.DATADIR}/prfm/patch_{k:03d}.h5'
    if not os.path.exists(pf):
        return None
    with h5py.File(pf, 'r') as f:
        frames = {g: dict(f[g].attrs) for g in f if g.startswith('frame_')}
    with h5py.File(f'{C.DATADIR}/grids/grid_{k:03d}.h5', 'r') as f:
        g = f['L1']; half = float(g.attrs['half']); n = int(g.attrs['n'])
        rho = g['rho'][:].astype(np.float32); vel = g['vel'][:].astype(np.float32)
    dx = 2 * half / n
    # mass per velocity voxel (2x coarser): sum 2x2x2 rho voxels
    lin = 10.0 ** rho; lin[rho < -11] = 0
    m = lin.reshape(n // 2, 2, n // 2, 2, n // 2, 2).sum((1, 3, 5)) * (dx * 1e3) ** 3   # Msun per vel voxel
    n2 = n // 2; dx2 = 2 * half / n2
    ax = (np.arange(n2) + 0.5) * dx2 - half
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing='ij')
    P = np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1) + C.CENTER
    V = vel.reshape(-1, 3); M = m.ravel()
    ok = M > 0
    P, V, M = P[ok], V[ok], M[ok]
    out = {}
    for name, fr in frames.items():
        nh = np.asarray(fr['nhat']); e1 = np.asarray(fr['e1']); e2 = np.asarray(fr['e2'])
        d = P - fr['center']; x = d @ e1; y = d @ e2; z = d @ nh
        v = V - fr['vcen']; vx = v @ e1; vy = v @ e2
        R = np.hypot(x, y); phi = np.arctan2(y, x)
        vphi = -vx * np.sin(phi) + vy * np.cos(phi); vR = vx * np.cos(phi) + vy * np.sin(phi)
        s = (np.abs(z) < ZMAX) & (R < RMAX)
        ib = np.clip((R[s] / DR).astype(int), 0, len(Rc) - 1)
        mw = np.bincount(ib, M[s], len(Rc))
        vp = np.bincount(ib, M[s] * vphi[s], len(Rc)) / np.maximum(mw, 1e-30)
        sign = np.sign(np.sum(mw * vp)); vp *= sign        # rotate positive
        vr2 = np.bincount(ib, M[s] * vR[s] ** 2, len(Rc)) / np.maximum(mw, 1e-30)
        vp2 = np.bincount(ib, M[s] * vphi[s] ** 2, len(Rc)) / np.maximum(mw, 1e-30) - vp ** 2
        sigR = np.sqrt(np.maximum(0.5 * (vr2 + vp2), 0))
        Sig = mw / (np.pi * (edges[1:] ** 2 - edges[:-1] ** 2) * 1e6)   # Msun/pc^2 within |z|<ZMAX
        vs = gaussian_filter1d(vp, 1.0, mode='nearest')
        Om = vs / np.maximum(Rc, 1e-3)
        dOm2 = np.gradient(Om ** 2, Rc)
        kap = np.sqrt(np.maximum(Rc * dOm2 + 4 * Om ** 2, 0))
        out[name] = dict(R=Rc, vphi=vp, Omega=Om, kappa=kap, sigR=sigR, Sigma=Sig, mass=mw)
    return out

if __name__ == '__main__':
    snaps = [int(a) for a in sys.argv[1:]] or sorted(int(os.path.basename(p)[6:9]) for p in glob.glob(f'{C.DATADIR}/prfm/patch_*.h5'))
    from multiprocessing import Pool
    with h5py.File(f'{C.DATADIR}/prfm/kappa.h5', 'w') as fo, Pool(int(os.environ.get('NPROC', 8))) as pool:
        for k, res in zip(snaps, pool.imap(one, snaps)):
            if res is None: continue
            for name, d in res.items():
                g = fo.create_group(f'snap_{k:03d}/{name}')
                for kk, v in d.items(): g.create_dataset(kk, data=v)
            fa = res[list(res)[0]]
            print(f'snap {k:3d}: frames {list(res)}; v_phi(1 kpc)={np.interp(1.0, fa["R"], fa["vphi"]):.1f} km/s, '
                  f'kappa(1 kpc)={np.interp(1.0, fa["R"], fa["kappa"]):.1f} km/s/kpc', flush=True)
