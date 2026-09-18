"""Velocity structure function of the gas in the disk patches: sigma_v(l)^2 = <|v(x+l)-v(x)|^2>
(mass-weighted, in-plane separations only, |z|<ZMAX in each galaxy frame), from the L2 (5.9 pc,
vel 11.7 pc) and L1 (15.6 pc, vel 31 pc) velocity grids.  Fits sigma_v ~ l^(p/2), i.e. v_t^2 ~ l^p.

usage: python structure_function.py SNAP [SNAP ...]
out:   /ptmp/uli/dwarf_merger/prfm/structure_function.npz (per snapshot: l, sigma_v, p_fit)
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, h5py
import common as C

ZMAX = 0.3   # kpc
LAGS = [1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64]

def sf_level(k, level):
    with h5py.File(f'{C.DATADIR}/grids/grid_{k:03d}.h5', 'r') as f:
        g = f[level]; half = float(g.attrs['half']); n = int(g.attrs['n'])
        rho = g['rho'][:].astype(np.float32); vel = g['vel'][:].astype(np.float32)
    lin = 10.0 ** rho; lin[rho < -11] = 0
    m = lin.reshape(n // 2, 2, n // 2, 2, n // 2, 2).sum((1, 3, 5)); n2 = n // 2; dx2 = 2 * half / n2
    with h5py.File(f'{C.DATADIR}/prfm/patch_{k:03d}.h5', 'r') as f:
        frames = [dict(f[g].attrs) for g in f if g.startswith('frame_')]
    ax = (np.arange(n2) + 0.5) * dx2 - half
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing='ij')
    out = {}
    for fr in frames:
        d = np.stack([X, Y, Z], -1) + C.CENTER - fr['center']
        z = d @ np.asarray(fr['nhat']); mask = (np.abs(z) < ZMAX) & (m > 0)
        res = []
        for L in LAGS:
            if L >= n2 // 2: continue
            num = 0.0; den = 0.0
            for axis in range(3):
                a = np.roll(vel, -L, axis=axis); ma = np.roll(mask, -L, axis=axis); mm = np.roll(m, -L, axis=axis)
                w = (mask & ma) * np.sqrt(m * mm)
                num += np.sum(w * np.sum((vel - a) ** 2, -1)); den += np.sum(w)
            res.append((L * dx2 * 1e3, np.sqrt(num / max(den, 1e-30))))
        out[fr['label'] if 'label' in fr else str(len(out))] = np.array(res)
    return out

if __name__ == '__main__':
    snaps = [int(a) for a in sys.argv[1:]]
    allres = {}
    for k in snaps:
        for level in ['L2', 'L1']:
            try:
                r = sf_level(k, level)
            except Exception as e:
                print('skip', k, level, e); continue
            for fl, arr in r.items():
                key = f'{k:03d}_{level}_{fl}'; allres[key] = arr
                l, s = arr[:, 0], arr[:, 1]; ok = (l > 10) & (l < 400)
                pfit = 2 * np.polyfit(np.log(l[ok]), np.log(s[ok]), 1)[0] if ok.sum() > 2 else np.nan
                print(f'snap {k} {level} frame {fl}: ' + ' '.join(f'{a:.0f}pc:{b:.1f}' for a, b in arr) + f'  -> p = {pfit:.2f}', flush=True)
    np.savez(f'{C.DATADIR}/prfm/structure_function.npz', **allres)
