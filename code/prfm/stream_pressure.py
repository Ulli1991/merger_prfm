"""How much of the layer's turbulent pressure is coherent streaming between the gas of the two progenitors (ram pressure of the
collision) rather than random turbulence?  Per column the vertical velocity variance of the 2p slab gas splits exactly as
    sigma^2 = f_A sigma_A^2 + f_B sigma_B^2 + f_A f_B (dv)^2 ,
with f the mass fractions of gas born in each galaxy (particle ID) and dv the difference of their mean v_n.  The last term is
the two-stream (ram) part.  Prints the Sigma_gas-weighted median of its share of P_turb per snapshot.
usage: python stream_pressure.py 60 115 130 150 180 200"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C
from patches import galaxy_frames, ID_SPLIT, GAS_SPLIT, PATCH, RMAX, SLAB_MIN, SLAB_FRAC, T_2P
print(f'gas ID split {GAS_SPLIT}, stellar ID split {ID_SPLIT}, patch {PATCH} kpc, rmax {RMAX} kpc')
print('snap    t     frame  n_col   median f_A   median dv [km/s]   two-stream share of P_turb (Sigma-wtd median, and of all columns)')
for k in [int(x) for x in sys.argv[1:]]:
    g = C.read(k, 0, ['Coordinates', 'Velocities', 'Masses', 'InternalEnergy', 'ElectronAbundance', 'ParticleIDs'])
    s2 = C.read(k, 2, ['Coordinates', 'ParticleIDs', 'Velocities'])
    t = g['Time'] * C.MYR
    pos = g['Coordinates'].astype(np.float64); vel = g['Velocities'].astype(np.float64); m = g['Masses'].astype(np.float64)
    T = C.temperature(g['InternalEnergy'], g['ElectronAbundance']); isA = g['ParticleIDs'] <= GAS_SPLIT
    for fr in galaxy_frames(s2['Coordinates'], s2['ParticleIDs'], s2['Velocities'].astype(np.float64), pos, vel, m, t):
        e1, e2, nh = None, None, fr['nhat']
        a = np.array([0.0, 1.0, 0.0]) if abs(nh[0]) > 0.9 else np.array([1.0, 0.0, 0.0])
        e1 = np.cross(nh, a); e1 /= np.linalg.norm(e1); e2 = np.cross(nh, e1)
        d = pos - fr['cen']; x = d @ e1; y = d @ e2; z = d @ nh
        vn = (vel - fr['vcen']) @ nh
        nb = int(2 * RMAX / PATCH)
        ix = np.floor((x + RMAX) / PATCH).astype(int); iy = np.floor((y + RMAX) / PATCH).astype(int)
        ing = (ix >= 0) & (ix < nb) & (iy >= 0) & (iy < nb) & (np.abs(z) < 0.5) & (T < T_2P)
        pid = ix[ing] * nb + iy[ing]; zz = z[ing]; mm = m[ing]; vv = vn[ing]; aa = isA[ing]
        npat = nb * nb; bs = lambda p, w: np.bincount(p, weights=w, minlength=npat)
        msum = bs(pid, mm); zmid = bs(pid, mm * zz) / np.maximum(msum, 1e-30)
        hs = np.maximum(SLAB_MIN, SLAB_FRAC * np.sqrt(np.maximum(bs(pid, mm * (zz - zmid[pid]) ** 2) / np.maximum(msum, 1e-30), 0)))
        ins = np.abs(zz - zmid[pid]) < hs[pid]
        p2, m2, v2, a2 = pid[ins], mm[ins], vv[ins], aa[ins]
        mA = bs(p2[a2], m2[a2]); mB = bs(p2[~a2], m2[~a2]); mt = mA + mB
        vA = bs(p2[a2], (m2 * v2)[a2]) / np.maximum(mA, 1e-30); vB = bs(p2[~a2], (m2 * v2)[~a2]) / np.maximum(mB, 1e-30)
        vbar = bs(p2, m2 * v2) / np.maximum(mt, 1e-30)
        var = bs(p2, m2 * (v2 - vbar[p2]) ** 2) / np.maximum(mt, 1e-30)
        fA = mA / np.maximum(mt, 1e-30); fB = 1 - fA; dv = vA - vB
        stream = fA * fB * dv ** 2
        ok = (mt > 0) & (mA > 0) & (mB > 0) & (var > 0) & (bs(pid, mm) * 1e10 / (PATCH ** 2 * 1e6) > 1)
        if ok.sum() < 3: print(f'{k:4d} {t:6.1f}  {fr["label"]}   {ok.sum():4d}  (too few mixed columns)'); continue
        w = mt[ok]; sh = (stream / var)[ok]
        o = np.argsort(sh); cw = np.cumsum(w[o]) / w.sum()
        print(f'{k:4d} {t:6.1f}  {fr["label"]}   {ok.sum():4d}     {np.median(fA[ok]):6.2f}        {np.median(np.abs(dv[ok])):8.1f}          {sh[o][np.searchsorted(cw,0.5)]:.3f}   (mean {sh.mean():.3f}, 90th pct {np.percentile(sh,90):.3f})', flush=True)
