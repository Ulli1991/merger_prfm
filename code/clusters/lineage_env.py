"""Part C test: does the SOURCE of the patch pressure set the boundedness of the clumps, and through it the efficiency?
For every star-forming lineage (lineages_<cat>.npz with k_on, x/y/z_on), take the patch it sits in at the pre-onset snapshot
(home frame = lowest intruder fraction) and record the patch state: W, P_tot, P_turb/P_tot, P_mag/P_tot, P_th/P_tot, feedback share
P_fb/P_tot with P_fb = Upsilon_OK22(P) x Sigma_SFR,40 (the previous 40 Myr, i.e. pre-existing feedback), particle-scale shear norm,
div_v, vorticity, sigma_z, H, rho_mid.
Then: (1) clump alpha_kin / alpha_tot / M_A / beta at onset vs the patch driving quantities (Spearman pooled and per phase; partial
at fixed P_tot); (2) eps_int vs the same at fixed alpha (partial); (3) the prediction 'feedback-supported patches host bound clumps,
flow-supported patches unbound ones' as medians of alpha in bins of the feedback share and of the shear.
usage: python lineage_env.py cl100|cl10      out: /ptmp/uli/dwarf_merger/clouds/lineage_env_<cat>.npz"""
import sys, os, numpy as np, h5py
sys.path.insert(0, '/raven/u/uli/dwarf_merger'); sys.path.insert(0, '/raven/u/uli/dwarf_merger/prfm'); import common as C, ok22
from scipy.stats import spearmanr
cat = sys.argv[1] if len(sys.argv) > 1 else 'cl100'
d = np.load(f'{C.DATADIR}/clouds/lineages_{cat}.npz'); L = d['rows']; col = {n: i for i, n in enumerate(d['names'])}
g = (L[:, col['Mstar']] > 0) & (L[:, col['M_max']] >= 300) & (L[:, col['t_start']] > 5) & (L[:, col['t_end']] < 221) & np.isfinite(L[:, col['k_on']])
idx = np.flatnonzero(g); ks = L[idx, col['k_on']].astype(int); pos = L[idx][:, [col['x_on'], col['y_on'], col['z_on']]]
FK = ['W_2p', 'Ptot_2p', 'Pturb_2p', 'Pmag_tot', 'Pth_2p', 'SigSFR_40', 'SigSFR_10', 'shear', 'div_v', 'vort', 'sigma_z_col', 'H', 'rho_mid_2p', 'f_intruder', 'Sigma_gas_2p']
env = np.full((len(idx), len(FK) + 2), np.nan)     # + frame code, layer W
for k in np.unique(ks):
    sel = np.flatnonzero(ks == k); fn = f'{C.PRFM_DIR}/patch_{k:03d}.h5'; fz = f'{C.PRFM_DIR}/patch_{k:03d}_z05.h5'
    if not os.path.exists(fn): continue
    with h5py.File(fn, 'r') as f, (h5py.File(fz, 'r') if os.path.exists(fz) else h5py.File(fn, 'r')) as fl:
        PATCH = float(f.attrs['patch_kpc']); RMAX = float(f.attrs['rmax_kpc']); nb = int(f.attrs['nb']); best = {}
        for gname in f:
            if not gname.startswith('frame_'): continue
            gr = f[gname]; cen = gr.attrs['center']; e1 = gr.attrs['e1']; e2 = gr.attrs['e2']; nh = gr.attrs['nhat']
            dd = pos[sel] - cen; x = dd @ e1; y = dd @ e2; z = dd @ nh
            ix = np.floor((x + RMAX) / PATCH).astype(int); iy = np.floor((y + RMAX) / PATCH).astype(int)
            ok = (ix >= 0) & (ix < nb) & (iy >= 0) & (iy < nb) & (np.abs(z) < 1.5)
            vals = {kk: gr[kk][:].ravel() for kk in FK}; WLl = fl[gname]['W_2p'][:].ravel() if gname in fl else vals['W_2p']
            for j in np.flatnonzero(ok):
                p = ix[j] * nb + iy[j]; fi = vals['f_intruder'][p]
                if j not in best or fi < best[j][0]: best[j] = (fi, [vals[kk][p] for kk in FK] + [ord(gname[-1]), WLl[p]])
        for j, (fi, row) in best.items(): env[sel[j]] = row
names = FK + ['frame', 'W_L']; ec = {n: i for i, n in enumerate(names)}
ok = np.isfinite(env[:, ec['Ptot_2p']]) & (env[:, ec['Ptot_2p']] > 0); env = env[ok]; Lg = L[idx[ok]]
np.savez(f'{C.DATADIR}/clouds/lineage_env_{cat}.npz', env=env, names=np.array(names), lineage_rows=Lg, lineage_names=d['names'])
print(f'{cat}: {len(Lg)} SF lineages joined to their pre-onset patch (of {len(idx)})')
P = env[:, ec['Ptot_2p']]; fturb = env[:, ec['Pturb_2p']] / P; fmag = env[:, ec['Pmag_tot']] / P; fth = env[:, ec['Pth_2p']] / P
Pfb = ok22.ups_tot_sfr(np.maximum(env[:, ec['SigSFR_40']], 1e-6)) * env[:, ec['SigSFR_40']] * 1e3 / 4.903245584337325e-06 * 1e-9  # K cm^-3 via P_unit? -> use ratio form below instead
# feedback share in the OK22 sense: P_fb = Upsilon_tot(P) [km/s] x Sigma_SFR [Msun/yr/kpc2] -> K cm^-3 : 1 km/s x 1 Msun/yr/kpc2 = 4.81e3 K cm^-3
ffb = np.clip(ok22.ups_tot(P) * env[:, ec['SigSFR_40']] * 4.81e3 / P, 0, 3)
shear = env[:, ec['shear']]; conv = -env[:, ec['div_v']]; vort = env[:, ec['vort']]
lc = lambda k: Lg[:, col[k]]
al = lc('alpha_on'); vA = lc('vA_on'); sig = lc('sig_on'); cs = lc('cs_on'); eps = lc('eps_int'); ton = lc('t_onset')
MA = sig / vA; beta = 2 * cs ** 2 / vA ** 2; al_tot = al * (1 + 1 / MA ** 2)
EP = [*C.PHASES]
def resid(y, x):
    A = np.vstack([np.ones(len(x)), x]).T; return y - A @ np.linalg.lstsq(A, y, rcond=None)[0]
print('\n(1) clump boundedness at onset vs the patch pressure source (Spearman; pooled over phases, and partial at fixed log P_tot)')
print('clump quantity   vs:  f_fb(40Myr)  f_turb   f_mag   f_th    shear   conv    vort    sigma_z   log P   log W_L')
for qn, q in (('alpha_kin', al), ('alpha_tot', al_tot), ('M_A', MA), ('beta', beta)):
    lq = np.log10(q); lP = np.log10(P)
    xs = [('f_fb', np.log10(np.maximum(ffb, 1e-3))), ('f_turb', fturb), ('f_mag', fmag), ('f_th', fth), ('shear', np.log10(shear)), ('conv', conv), ('vort', np.log10(np.maximum(vort, 1e-3))), ('sigma_z', np.log10(env[:, ec['sigma_z_col']])), ('logP', lP), ('logWL', np.log10(env[:, ec['W_L']]))]
    print(f'{qn:14s} pooled:  ' + '  '.join(f'{spearmanr(lq, x)[0]:+.2f}  ' for _, x in xs))
    print(f'{"":14s} at fixed P: ' + '  '.join(f'{spearmanr(resid(lq, lP), resid(x, lP))[0]:+.2f}  ' for _, x in xs[:-2]))
print('\n(2) eps_int vs the patch pressure source at fixed clump alpha_kin (partial Spearman), pooled and per phase')
le = np.log10(eps); la = np.log10(al)
print('phase                 n    f_fb    f_turb   f_mag   shear   conv    | Spearman(eps, alpha) for ref')
for t0, t1, lab in EP + [(25, 232, 'all')]:
    m = (ton >= t0) & (ton < t1)
    if m.sum() < 15: continue
    r = lambda x: spearmanr(resid(le[m], la[m]), resid(x[m], la[m]))[0]
    print(f'{lab:20s} {m.sum():4d}   {r(np.log10(np.maximum(ffb,1e-3))):+.2f}   {r(fturb):+.2f}    {r(fmag):+.2f}   {r(np.log10(shear)):+.2f}   {r(conv):+.2f}   | {spearmanr(le[m], la[m])[0]:+.2f}')
print('\n(3) median clump alpha_kin (and eps_int) in bins of the feedback share of the patch pressure, and of the magnetic share; n in brackets')
for name, x, edges in (('f_fb', ffb, [0, 0.03, 0.1, 0.3, 1, 3]), ('f_mag', fmag, [0, 0.05, 0.1, 0.2, 0.4, 1]), ('f_turb', fturb, [0, 0.3, 0.5, 0.7, 0.85, 1])):
    line = f'{name:7s} '
    for a, b in zip(edges[:-1], edges[1:]):
        m = (x >= a) & (x < b); line += f'  [{a:g},{b:g}): alpha {np.median(al[m]) if m.sum() >= 8 else np.nan:5.1f} eps {np.median(eps[m]) if m.sum() >= 8 else np.nan:.4f} ({m.sum():3d})'
    print(line)
print('\nper phase: median patch f_fb, f_turb, f_mag, shear at onset, and clump alpha, M_A, beta')
for t0, t1, lab in EP:
    m = (ton >= t0) & (ton < t1)
    print(f'{lab:20s} n={m.sum():4d}  f_fb {np.median(ffb[m]):.3f}  f_turb {np.median(fturb[m]):.2f}  f_mag {np.median(fmag[m]):.2f}  shear {np.median(shear[m]):6.0f}  | alpha {np.median(al[m]):5.1f}  M_A {np.median(MA[m]):4.2f}  beta {np.median(beta[m]):5.2f}')
