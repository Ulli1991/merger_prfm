"""Which local tracer gives the sheet normal?  At snapshot K (default 100, two clean discs) compare candidate normals on the 0.5 kpc
node grid within 2 kpc of either galaxy centre against the galaxy-frame normal (gas angular momentum within 2 kpc of the centre):
  inertia tensor (mass-weighted covariance, thinnest axis) of: (a) T < 2e4 gas, (b) T < 1e3 gas, (c) n_H > 1 gas, (d) n_H > 0.1 gas,
  (e) T < 2e4 gas weighted by density; (f) angular momentum of T < 2e4 gas; (g) of n_H > 1 gas;  each within R = 0.5 and 1 kpc.
Prints median / 84 % angle to the galaxy normal and the median flatness."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common as C
from scipy.spatial import cKDTree
k = int(sys.argv[1]) if len(sys.argv) > 1 else 100
g = C.read(k, 0, ['Coordinates', 'Velocities', 'Masses', 'Density', 'InternalEnergy', 'ElectronAbundance']); s2 = C.read(k, 2, ['Coordinates', 'ParticleIDs'])
pos = g['Coordinates'].astype(np.float64); vel = g['Velocities'].astype(np.float64); m = g['Masses'].astype(np.float64); T = C.temperature(g['InternalEnergy'], g['ElectronAbundance']); n = C.nH_cgs(g['Density'].astype(np.float64))
A = s2['ParticleIDs'] <= 26_000_000; cA = np.median(s2['Coordinates'][A], 0); cB = np.median(s2['Coordinates'][~A], 0)
def galnorm(c):
    d = pos - c; sel = (np.linalg.norm(d, axis=1) < 2.0) & (T < 2e4); vc = np.average(vel[sel], axis=0, weights=m[sel]); L = np.sum(m[sel, None] * np.cross(d[sel], vel[sel] - vc), 0); return L / np.linalg.norm(L)
nA, nB = galnorm(cA), galnorm(cB)
ax_ = np.arange(-6, 6.01, 0.5); nodes = np.array(np.meshgrid(ax_, ax_, ax_, indexing='ij')).reshape(3, -1).T + C.CENTER
dA = np.linalg.norm(nodes - cA, axis=1); dB = np.linalg.norm(nodes - cB, axis=1); nodes = nodes[(dA < 2) | (dB < 2)]; ref = np.where((dA < dB)[(dA < 2) | (dB < 2)][:, None], nA, nB)
sel2p = T < 2e4; tree = cKDTree(pos[sel2p]); idx2p = np.flatnonzero(sel2p)
cnt = tree.query_ball_point(nodes, 0.25, return_length=True); keep = cnt > 200; nodes = nodes[keep]; ref = ref[keep]
print(f'snap {k}: {len(nodes)} nodes with gas within 2 kpc of a galaxy centre')
cands = {'inertia T<2e4': (T < 2e4, None), 'inertia T<1e3': (T < 1e3, None), 'inertia n>1': (n > 1, None), 'inertia n>0.1': (n > 0.1, None), 'inertia T<2e4 rho-wtd': (T < 2e4, 'rho'), 'L T<2e4': (T < 2e4, 'L'), 'L n>1': (n > 1, 'L')}
for R in (0.5, 1.0):
    lists = tree.query_ball_point(nodes, R)
    for name, (mask, mode) in cands.items():
        angs = []; qs = []
        for i, lst in enumerate(lists):
            idx = idx2p[lst]; idx = idx[mask[idx]]
            if len(idx) < 50: continue
            p = pos[idx]; w = m[idx] * (n[idx] if mode == 'rho' else 1.0); v = vel[idx]
            com = np.average(p, axis=0, weights=w); d = p - com
            if mode == 'L':
                vc = np.average(v, axis=0, weights=w); Lv = np.sum(w[:, None] * np.cross(d, v - vc), 0); nn = Lv / max(np.linalg.norm(Lv), 1e-30); q = np.nan
            else:
                I = (d * w[:, None]).T @ d / w.sum(); ev, V = np.linalg.eigh(I); nn = V[:, 0]; q = np.sqrt(max(ev[0], 0) / max(ev[1], 1e-30))
            angs.append(np.degrees(np.arccos(min(1, abs(nn @ ref[i]))))); qs.append(q)
        angs = np.array(angs); print(f'  R={R} kpc  {name:24s} n={len(angs):4d}  angle to galaxy normal median {np.median(angs):5.1f}  84% {np.percentile(angs, 84):5.1f} deg   flatness median {np.nanmedian(qs):.2f}')
