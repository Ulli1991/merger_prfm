"""Frame-free PRFM columns: cubes on a fixed spatial grid, each with its own vertical from the LOCAL gas sheet.

For snapshot k and cube side L (PRFM_CUBE, kpc; column half-depth PRFM_ZCOL, default L/2 -> an actual cube):
  1. normal field: on a 0.5 kpc grid of nodes, for the COLD gas (T < 1e3 K) within R_n of each node (R_n = 0.5, 1 kpc):
     inertia (mass-weighted covariance) tensor -> thinnest direction n (eigenvector of the smallest eigenvalue), flatness
     q = sqrt(lambda_min / lambda_mid) (disc ~0.1-0.2, blob ~0.5-1), the angular-momentum direction within the sphere and its
     angle to n, and the local mean velocity.
  2. cubes: all cubes of side L in a +-6 kpc box with gas mass / L^2 > 1 Msun/pc^2; each takes n, q, vcen from the nearest node
     (R_n = 1 kpc for the analysis; the other radii are stored for the sensitivity test).
  3. per cube, in its local frame (e1, e2, n), column |x|,|y| < L/2, |z| < ZCOL: Sigma, zmid, H, the weight
     W = mean of the upper and lower integrals of rho g_n from the PM gravity solve (all gas, and T < 2e4 gas), the midplane slab
     pressures (thermal, turbulent about the slab mean v_n, magnetic B^2/8pi), rho_mid, Sigma_SFR over 10 and 40 Myr from stars
     formed in the run, f_cold (T < 500 K), and the angles of n to the galaxy-frame normal of the nearest galaxy (for reference).
usage: PRFM_CUBE=0.25 python cubes.py K [K ...]      out: /ptmp/uli/dwarf_merger/prfm/cubes_KKK_l<L>.h5"""
import sys, os, time, numpy as np, h5py
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C
from gravity import two_level_accel
from scipy.spatial import cKDTree
L_CUBE = float(os.environ.get('PRFM_CUBE', 0.5)); ZCOL = float(os.environ.get('PRFM_ZCOL', 0.5)); TAG = os.environ.get('PRFM_TAG', '')   # depth +-0.5 kpc by default: shallower columns truncate the weight (notes sec 25)
BOX = 6.0; NODE = 0.5; RN = (0.5, 1.0); RN_USE = 1.0
FINE = (4.0, 512); COARSE = (32.0, 512); T_2P = 2.0e4; T_COLD = 500.; T_SHEET = 1.0e3; SLAB_MIN = 0.025; SLAB_FRAC = 0.25
P_UNIT = 1.989e33 * 1e10 / (3.085678e21 ** 3) / 1.3807e-16
ID_SPLIT = 26_000_000

def axes(n):
    n = n / np.linalg.norm(n); a = np.array([1.0, 0, 0]) if abs(n[0]) < 0.9 else np.array([0, 1.0, 0])
    e1 = np.cross(n, a); e1 /= np.linalg.norm(e1); e2 = np.cross(n, e1); return e1, e2, n

def galaxy_normals(s2pos, s2ids, gpos, gvel, gm):
    """galaxy-frame normals (as in patches.py) for the reference angle"""
    A = s2ids <= ID_SPLIT; cA = np.median(s2pos[A], 0); cB = np.median(s2pos[~A], 0); out = []
    for c in ((cA, cB) if np.linalg.norm(cA - cB) > 0.5 else (0.5 * (cA + cB),)):
        d = gpos - c; sel = np.linalg.norm(d, axis=1) < (2.0 if np.linalg.norm(cA - cB) > 0.5 else 3.0)
        vc = np.average(gvel[sel], axis=0, weights=gm[sel]); Lv = np.sum(gm[sel, None] * np.cross(d[sel], gvel[sel] - vc), 0); out.append((c, Lv / np.linalg.norm(Lv)))
    return out

def analyse(k):
    t0 = time.time(); fout = f'{C.DATADIR}/prfm/cubes_{k:03d}_l{L_CUBE:g}{TAG}.h5'
    if os.path.exists(fout): print('exists', fout); return
    g = C.read(k, 0, ['Coordinates', 'Velocities', 'Masses', 'Density', 'InternalEnergy', 'ElectronAbundance', 'MagneticField'])
    tmyr = g['Time'] * C.MYR; gpos = g['Coordinates'].astype(np.float64); gvel = g['Velocities'].astype(np.float64)
    gm = g['Masses'].astype(np.float64) * C.MSUN; T = C.temperature(g['InternalEnergy'], g['ElectronAbundance']); u = g['InternalEnergy'].astype(np.float64)
    B2 = np.sum(g['MagneticField'].astype(np.float64) ** 2, axis=1); vol = g['Masses'].astype(np.float64) / g['Density'].astype(np.float64)
    s2 = C.read(k, 2, ['Coordinates', 'Masses', 'ParticleIDs']); s4 = C.read(k, 4, ['Coordinates', 'Masses', 'StellarFormationTime']); dm = C.read(k, 1, ['Coordinates', 'Masses'])
    s2pos = s2['Coordinates'].astype(np.float64); s2m = s2['Masses'].astype(np.float64) * C.MSUN
    s4pos = s4['Coordinates'].astype(np.float64); s4m = s4['Masses'].astype(np.float64) * C.MSUN; sft = s4['StellarFormationTime'].astype(np.float64)
    age = tmyr - sft * C.MYR; formed = sft > 0
    dmpos = dm['Coordinates'].astype(np.float64); dmm = dm['Masses'].astype(np.float64) * C.MSUN
    print(f'snap {k} t={tmyr:.1f} Myr read ({time.time()-t0:.0f}s)', flush=True)
    # gravity at gas positions
    spos = np.vstack([s2pos, s4pos]); sm = np.concatenate([s2m, s4m])
    acc = two_level_accel(gpos, {'gas': (gpos, gm), 'star': (spos, sm), 'dm': (dmpos, dmm)}, C.CENTER, FINE[0], FINE[1], COARSE[0], COARSE[1])['total'].astype(np.float64)
    print(f'  gravity done ({time.time()-t0:.0f}s)', flush=True)
    # ---- 1. normal field on nodes
    twop = T < T_2P; sub = np.flatnonzero(T < T_SHEET)[::3]; tree2p = cKDTree(gpos[sub])     # the sheet is traced by the cold gas (T < 1e3 K), not by the warm envelope
    ax_ = np.arange(-BOX, BOX + 1e-6, NODE); nodes = np.array(np.meshgrid(ax_, ax_, ax_, indexing='ij')).reshape(3, -1).T + C.CENTER
    cnt = tree2p.query_ball_point(nodes, 0.25, return_length=True); nodes = nodes[cnt > 20]
    nf = {r: dict(n=np.zeros((len(nodes), 3)), q=np.full(len(nodes), np.nan), ang=np.full(len(nodes), np.nan), v=np.zeros((len(nodes), 3))) for r in RN}
    for r in RN:
        for i, lst in enumerate(tree2p.query_ball_point(nodes, r)):
            if len(lst) < 50: continue
            idx = sub[lst]; p = gpos[idx]; m = gm[idx]; v = gvel[idx]
            com = np.average(p, axis=0, weights=m); vc = np.average(v, axis=0, weights=m); d = p - com
            I = (d * m[:, None]).T @ d / m.sum(); w, V = np.linalg.eigh(I); n = V[:, 0]
            Lv = np.sum(m[:, None] * np.cross(d, v - vc), 0); Lv /= max(np.linalg.norm(Lv), 1e-30)
            nf[r]['n'][i] = n; nf[r]['q'][i] = np.sqrt(max(w[0], 0) / max(w[1], 1e-30)); nf[r]['ang'][i] = np.degrees(np.arccos(abs(n @ Lv))); nf[r]['v'][i] = vc
    print(f'  normal field on {len(nodes)} nodes ({time.time()-t0:.0f}s)', flush=True)
    # ---- 2. cubes
    ax_c = np.arange(-BOX, BOX + 1e-6, L_CUBE); edges = [ax_c + C.CENTER[i] for i in range(3)]
    Hm, _ = np.histogramdd(gpos, bins=edges, weights=gm); cid = np.argwhere(Hm / (L_CUBE * 1e3) ** 2 > 1.0)
    cen = np.stack([edges[i][cid[:, i]] + L_CUBE / 2 for i in range(3)], 1)
    ntree = cKDTree(nodes); _, jn = ntree.query(cen)
    tree = cKDTree(gpos); stree = cKDTree(s4pos[formed]); s4f = s4pos[formed]; s4mf = s4m[formed]; agef = age[formed]
    gal = galaxy_normals(s2pos, s2['ParticleIDs'], gpos, gvel, gm)
    rq = np.sqrt(2 * (L_CUBE / 2) ** 2 + ZCOL ** 2); A = L_CUBE ** 2
    keys = ['x', 'y', 'z', 'nx', 'ny', 'nz', 'q', 'ang_nL', 'ang_ngal', 'Sigma', 'Sigma_2p', 'zmid', 'H', 'W_2p', 'W_2p_up', 'W_2p_lo', 'W_tot', 'Pth_2p', 'Pturb_2p', 'Pmag_tot', 'Ptot_2p', 'rho_mid_2p', 'SigSFR_10', 'SigSFR_40', 'fcold', 'n_gas', 'sigma_z', 'vbar_n']
    out = {kk: np.full(len(cen), np.nan) for kk in keys}
    for c in range(len(cen)):
        node = jn[c]; n = nf[RN_USE]['n'][node]
        if not np.any(n): continue
        e1, e2, n = axes(n); vc = nf[RN_USE]['v'][node]
        idx = np.array(tree.query_ball_point(cen[c], rq))
        if len(idx) < 20: continue
        d = gpos[idx] - cen[c]; x = d @ e1; y = d @ e2; z = d @ n
        inc = (np.abs(x) < L_CUBE / 2) & (np.abs(y) < L_CUBE / 2) & (np.abs(z) < ZCOL)
        if inc.sum() < 20: continue
        idx = idx[inc]; z = z[inc]; m = gm[idx]; Tc = T[idx]; M = m.sum(); tp = Tc < T_2P
        zmid = np.average(z, weights=m); Hc = np.sqrt(np.average((z - zmid) ** 2, weights=m)); dz = z - zmid
        gn = acc[idx] @ n; up = dz > 0
        Wup = -(m * gn)[up].sum() / A * P_UNIT; Wlo = (m * gn)[~up].sum() / A * P_UNIT
        Wup2 = -(m * gn)[up & tp].sum() / A * P_UNIT; Wlo2 = (m * gn)[~up & tp].sum() / A * P_UNIT
        hs = min(max(SLAB_MIN, SLAB_FRAC * Hc), ZCOL); sl = np.abs(dz) < hs; Vs = 2 * hs * A
        vn = (gvel[idx] - vc) @ n; s2p = sl & tp
        if s2p.sum() < 5: continue
        vb = np.average(vn[s2p], weights=m[s2p])
        Pth = np.sum(m[s2p] * (C.GAMMA - 1) * u[idx][s2p]) / Vs * P_UNIT
        Pturb = np.sum(m[s2p] * (vn[s2p] - vb) ** 2) / Vs * P_UNIT
        Pmag = np.sum(vol[idx][sl] * B2[idx][sl] / (8 * np.pi) / C.BOLTZMANN) / Vs
        rho2p = m[s2p].sum() / Vs / 1e9
        sidx = np.array(stree.query_ball_point(cen[c], rq), dtype=int)
        sf10 = sf40 = 0.0
        if len(sidx):
            ds = s4f[sidx] - cen[c]; sx = ds @ e1; sy = ds @ e2; sz = ds @ n; ins = (np.abs(sx) < L_CUBE / 2) & (np.abs(sy) < L_CUBE / 2) & (np.abs(sz) < ZCOL)
            a_ = agef[sidx][ins]; mm = s4mf[sidx][ins]; sf10 = mm[(a_ >= 0) & (a_ < 10)].sum() / A / 1e7; sf40 = mm[(a_ >= 0) & (a_ < 40)].sum() / A / 4e7
        gc, gnrm = min(gal, key=lambda t: np.linalg.norm(cen[c] - t[0]))
        vals = dict(x=cen[c][0], y=cen[c][1], z=cen[c][2], nx=n[0], ny=n[1], nz=n[2], q=nf[RN_USE]['q'][node], ang_nL=nf[RN_USE]['ang'][node], ang_ngal=np.degrees(np.arccos(abs(n @ gnrm))),
                    Sigma=M / A / 1e6, Sigma_2p=m[tp].sum() / A / 1e6, zmid=zmid, H=Hc, W_2p=0.5 * (Wup2 + Wlo2), W_2p_up=Wup2, W_2p_lo=Wlo2, W_tot=0.5 * (Wup + Wlo),
                    Pth_2p=Pth, Pturb_2p=Pturb, Pmag_tot=Pmag, Ptot_2p=Pth + Pturb + Pmag, rho_mid_2p=rho2p, SigSFR_10=sf10, SigSFR_40=sf40,
                    fcold=m[Tc < T_COLD].sum() / M, n_gas=len(idx), sigma_z=np.sqrt(np.average((vn - np.average(vn, weights=m)) ** 2, weights=m)), vbar_n=np.average(vn, weights=m))
        for kk, v in vals.items(): out[kk][c] = v
    ok = np.isfinite(out['W_2p'])
    with h5py.File(fout + '.tmp', 'w') as f:
        f.attrs['time_myr'] = tmyr; f.attrs['cube_kpc'] = L_CUBE; f.attrs['zcol_kpc'] = ZCOL; f.attrs['rn_use'] = RN_USE; f.attrs['P_unit'] = P_UNIT
        for kk in keys: f.create_dataset(kk, data=out[kk][ok])
        gn_ = f.create_group('nodes'); gn_.create_dataset('pos', data=nodes)
        for r in RN:
            for kk in ('n', 'q', 'ang', 'v'): gn_.create_dataset(f'r{r:g}_{kk}', data=nf[r][kk])
    os.rename(fout + '.tmp', fout)
    r = out['Ptot_2p'][ok] / out['W_2p'][ok]; w = out['Sigma'][ok]; o = np.argsort(r); cw = np.cumsum(w[o]) / w.sum()
    print(f'  {ok.sum()} cubes with gas; Sigma-weighted median P/W {r[o][np.searchsorted(cw, 0.5)]:.2f}, ratio of sums {np.sum(out["Ptot_2p"][ok]) / np.sum(out["W_2p"][ok]):.2f}; flatness median {np.nanmedian(out["q"][ok]):.2f}; angle to galaxy normal median {np.nanmedian(out["ang_ngal"][ok]):.0f} deg ({time.time()-t0:.0f}s)', flush=True)

if __name__ == '__main__':
    for a in sys.argv[1:]: analyse(int(a))
