"""Energy test for the nucleus-flagged groups (N > 40000 or r_h > 20 pc), which cluster_mf.py leaves untested.
PE from a random subsample of 4000 members against all members (Plummer 1 pc), KE about the group's mean velocity.
usage: python nucleus_bound.py  -> prints per snapshot / age window the flagged groups with M, r_h, KE, PE, bound"""
import sys, os, numpy as np, numba as nb
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C
from fof import fof
G = 4.3009e-3   # pc (km/s)^2 / Msun
@nb.njit(parallel=True, fastmath=True)
def pot_sub(pos, m, idx, soft):
    out = np.zeros(len(idx))
    for a in nb.prange(len(idx)):
        i = idx[a]; s = 0.0
        for j in range(pos.shape[0]):
            if j == i: continue
            dx = pos[j,0]-pos[i,0]; dy = pos[j,1]-pos[i,1]; dz = pos[j,2]-pos[i,2]
            s += m[j] / np.sqrt(dx*dx+dy*dy+dz*dz+soft*soft)
        out[a] = -G * s
    return out
rng = np.random.default_rng(1)
for k in (130, 150, 160, 170, 210, 220, 230):
    st = np.load(f'{C.DATADIR}/stars/stars_{k:03d}.npz'); t = float(st['time_myr']); age = t - st['tform_myr']
    for lo, hi in ((0, 10), (20, 50)):
        s = (age >= lo) & (age < hi) & (st['tform_myr'] > 0)
        pos = st['pos'][s].astype(np.float64); vel = st['vel'][s].astype(np.float64); m = st['mass'][s] * C.MSUN
        lab = fof(pos, 5e-3); ul, cnt = np.unique(lab[lab >= 0], return_counts=True)
        for g, n in zip(ul, cnt):
            if n < 25: continue
            mem = lab == g; p = (pos[mem] - pos[mem].mean(0)) * 1e3; v = vel[mem]; mm = m[mem]
            com = np.average(p, axis=0, weights=mm); r = np.linalg.norm(p - com, axis=1); rh = np.median(r)
            if not (n > 40000 or rh > 20): continue
            vc = np.average(v, axis=0, weights=mm); ke = 0.5 * np.sum(mm * np.sum((v - vc)**2, 1))
            idx = rng.choice(n, min(4000, n), replace=False); phi = pot_sub(p, mm, idx, 1.0)
            pe = 0.5 * np.sum(mm) * phi.mean()   # 0.5 sum_i m_i phi_i estimated from the subsample
            sig = np.sqrt(2 * ke / mm.sum() / 3)
            print(f'snap {k} t={t:5.1f} age {lo}-{hi}: N={n:6d} M={mm.sum():9.3g} Msun r_h={rh:5.1f} pc sigma_1d={sig:5.1f} km/s  KE={ke:9.3g} PE={pe:9.3g}  KE/|PE|={ke/abs(pe):5.2f}  bound={ke+pe<0}', flush=True)
