"""Cold-cloud catalogue for one snapshot and the stars each cloud forms in the next 10 Myr.

Stars inherit the ID of the gas particle they formed from (verified: all 2670 stars formed between snaps 100 and 110 have
their ID in the gas of snap 100), so the origin of every young star is known without positional matching.

For snapshot k:
  clouds  = FoF (linking length LINK) of gas with T < TMAX and n_H > NTH, groups of >= NMIN particles (4 Msun each);
            two catalogues: (NTH=10, LINK=3 pc) 'cl10' and (NTH=100, LINK=1.5 pc) 'cl100'.
  per cloud: mass, N, com, vcom, sigma_3d (mass-weighted), r_h, n_mean, n_max, alpha_vir = 5 sigma_1d^2 r_h / (G M),
            M_star_3 / M_star_10 = mass of stars formed from the cloud's members within 3 / 10 Myr after t_k.
  per new star (formed in (t_k, t_k + 10 Myr], from stars_{k+10}.npz): id, tform, pos, mass, parent gas n_H and T at k,
            cloud label in both catalogues (-1 if the parent was not in a cloud).
usage: python clouds.py K          out: /ptmp/uli/dwarf_merger/clouds/clouds_KKK.npz"""
import sys, os, time, numpy as np, h5py
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C
from fof import fof
k = int(sys.argv[1]); OUT = f'{C.DATADIR}/clouds'; os.makedirs(OUT, exist_ok=True)
TMAX = 1e3; NMIN = 25; G_PC = 4.30091e-3
CATS = {'cl10': (10., 3e-3), 'cl100': (100., 1.5e-3)}
t0 = time.time()
with h5py.File(C.snap_path(k), 'r') as f:
    g = f['PartType0']; tk = float(f['Header'].attrs['Time']) * C.MYR
    rho = g['Density'][:].astype(np.float64); u = g['InternalEnergy'][:].astype(np.float64); ne = g['ElectronAbundance'][:].astype(np.float64)
    T = C.temperature(u, ne); n = C.nH_cgs(rho); del u, ne, rho
    cold = (T < TMAX) & (n > 10.)
    ids_all = g['ParticleIDs'][:]
    pos = g['Coordinates'][:][cold].astype(np.float64); vel = g['Velocities'][:][cold].astype(np.float64); m = g['Masses'][:][cold].astype(np.float64) * C.MSUN
    Bc = np.sqrt(np.sum(g['MagneticField'][:][cold].astype(np.float64) ** 2, axis=1))      # Gauss
ids = ids_all[cold]; nc = n[cold]; Tc = T[cold]
rho_cgs = nc * C.PROTONMASS / C.XH                                                          # g cm^-3
cs_c = np.sqrt(C.BOLTZMANN * Tc / (1.27 * C.PROTONMASS)) / 1e5                                # isothermal sound speed, km/s
vA_c = Bc / np.sqrt(4 * np.pi * rho_cgs) / 1e5                                              # Alfven speed, km/s
print(f'snap {k} t={tk:.1f} Myr: {cold.sum()} cold dense particles ({m.sum():.3g} Msun) [{time.time()-t0:.0f}s]', flush=True)
# new stars in the next 10 Myr
ks = min(k + 10, C.NSNAP - 1); s = np.load(f'{C.DATADIR}/stars/stars_{ks:03d}.npz')
new = (s['tform_myr'] > tk) & (s['tform_myr'] <= tk + 10.) & (s['tform_myr'] > 0)
sid = s['ids'][new]; stf = s['tform_myr'][new]; spos = s['pos'][new]; sm = s['mass'][new].astype(np.float64) * C.MSUN
# parent gas properties at k for every new star
o = np.argsort(ids_all); j = np.searchsorted(ids_all, sid, sorter=o); j = np.minimum(j, len(o) - 1); hit = ids_all[o[j]] == sid
par_n = np.full(len(sid), np.nan); par_T = np.full(len(sid), np.nan); par_n[hit] = n[o[j[hit]]]; par_T[hit] = T[o[j[hit]]]
print(f'  {new.sum()} stars formed in (t_k, t_k+10]: parent found for {hit.sum()}, parent n>10 & cold: {np.sum((par_n > 10) & (par_T < TMAX))}, n>100: {np.sum(par_n > 100)}', flush=True)
out = dict(snap=k, t=tk, star_id=sid, star_tform=stf, star_pos=spos, star_mass=sm, star_par_n=par_n, star_par_T=par_T)
oc = np.argsort(ids)
for cat, (nth, link) in CATS.items():
    sel = nc > nth; p = pos[sel]; v = vel[sel]; mm = m[sel]; ii = ids[sel]; nn = nc[sel]; Bs = Bc[sel]; css = cs_c[sel]; vAs = vA_c[sel]; Ts = Tc[sel]
    lab = fof(p, link); cnt = np.bincount(lab); good = np.flatnonzero(cnt >= NMIN)
    relab = np.full(cnt.size, -1, np.int64); relab[good] = np.arange(len(good)); cl = relab[lab]   # -1 = not in a cloud
    rows = []
    # star -> cloud label
    o2 = np.argsort(ii); j2 = np.searchsorted(ii, sid, sorter=o2); j2 = np.minimum(j2, len(o2) - 1); h2 = ii[o2[j2]] == sid
    star_cl = np.full(len(sid), -1, np.int64); star_cl[h2] = cl[o2[j2[h2]]]
    ms3 = np.zeros(len(good)); ms10 = np.zeros(len(good))
    inc = star_cl >= 0
    np.add.at(ms10, star_cl[inc], sm[inc]); inc3 = inc & (stf <= tk + 3.); np.add.at(ms3, star_cl[inc3], sm[inc3])
    order = np.argsort(cl, kind='stable'); starts = np.searchsorted(cl[order], np.arange(len(good)))
    ends = np.append(starts[1:], len(order)) if len(good) else starts
    for c in range(len(good)):
        idx = order[starts[c]:ends[c]]; pp = p[idx]; vv = v[idx]; w = mm[idx]; M = w.sum()
        com = np.average(pp, axis=0, weights=w); vcom = np.average(vv, axis=0, weights=w)
        sig3 = np.sqrt(np.average(np.sum((vv - vcom) ** 2, axis=1), weights=w))
        r = np.linalg.norm(pp - com, axis=1) * 1e3; rs = np.sort(r); rh = np.interp(0.5, np.cumsum(w[np.argsort(r)]) / M, rs)
        avir = 5 * (sig3 ** 2 / 3) * rh / (G_PC * M)
        Brms = np.sqrt(np.mean(Bs[idx] ** 2)) * 1e6; cs = np.sqrt(np.mean(css[idx] ** 2)); vA = np.sqrt(np.mean(vAs[idx] ** 2)); sig_s = np.std(np.log(nn[idx])); Tm = np.mean(Ts[idx])
        rows.append((k, tk, M, len(idx), com[0], com[1], com[2], vcom[0], vcom[1], vcom[2], sig3, rh, nn[idx].mean(), nn[idx].max(), avir, ms3[c], ms10[c], Brms, cs, vA, sig_s, Tm))
    R = np.array(rows, dtype=np.float64).reshape(-1, 22)
    out[f'{cat}_rows'] = R; out[f'{cat}_star_cl'] = star_cl
    nneg = int(np.sum(cl < 0)); memb = order[nneg:]      # members of kept clouds, grouped by cloud (cl sorted ascending, -1 first)
    out[f'{cat}_member_ids'] = ii[memb].astype(np.uint32)
    out[f'{cat}_member_off'] = (np.append(starts, len(order)) - nneg).astype(np.int64) if len(good) else np.zeros(1, np.int64)
    if len(R):
        print(f'  {cat}: {len(R)} clouds, M_max {R[:,2].max():.3g}, total {R[:,2].sum():.3g} Msun; stars from clouds {ms10.sum():.3g} of {sm.sum():.3g} ({ms10.sum()/max(sm.sum(),1):.2f}); '
              f'clouds with SF in 10 Myr: {np.sum(ms10 > 0)} [{time.time()-t0:.0f}s]', flush=True)
out['names'] = np.array(['snap', 't', 'M', 'N', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'sig3d', 'rh_pc', 'n_mean', 'n_max', 'alpha_vir', 'Mstar_3', 'Mstar_10', 'B_rms_muG', 'cs', 'vA', 'sig_s', 'T_mean'])
np.savez(f'{OUT}/clouds_{k:03d}.npz', **out); print('done', f'{OUT}/clouds_{k:03d}.npz')
