"""Cluster formation vs. the directly measured ISM weight (PRFM x clusters).

Joins, per snapshot and galaxy frame, the PRFM patch columns (patch_SSS.h5)
with the young-cluster catalogue (clusters_age0-10_l5_n25.npz), the young
stars (stars_SSS.npz, age < 10 Myr, formed in the run) and kappa(R)
(kappa.h5).  Per patch:
    M_young   mass of stars younger than AGE Myr in the column
    M_bound   mass of energetically bound FoF groups (age < AGE) in the column
    M_max     most massive bound group
    Gamma     = M_bound / M_young   (measured bound fraction, 'CFE')
and the theory:
    Gamma_K12(rho_mid, Mach)    Kruijssen 2012 without the cruel-cradle term,
                                driven by the measured midplane density
    Gamma_std(Sigma, kappa, Q)  same model with rho_ISM = phi_P Omega^2/(pi G Q^2)
                                (gas-only, standard route)
    M_T = 4 pi^5 G^2 Sigma^3 / kappa^4              Toomre mass
    M_cmax = eps Gamma f_coll M_T,  f_coll = min[1,(t_fb/t_ff)^4]   (RC&K 2017)

usage: python cluster_pressure.py [--age 10] [--tfb 3] [--fintr 0.1]
out:   /ptmp/uli/dwarf_merger/prfm/patches_clusters.npz, figs cl_fig*.png
"""
import sys, os, glob, argparse
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, h5py
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import common as C

ap = argparse.ArgumentParser()
ap.add_argument('--age', type=float, default=10.0)
ap.add_argument('--tfb', type=float, default=3.0)       # Myr, feedback time
ap.add_argument('--fintr', type=float, default=0.1)     # max intruder fraction for 'clean' patches
ap.add_argument('--mmin_young', type=float, default=500.)  # Msun of young stars for a patch to count
a = ap.parse_args()
FIG = f'{C.DATADIR}/prfm/figs'; os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.grid': True, 'grid.alpha': 0.25, 'legend.frameon': False, 'figure.dpi': 130})

G_PC = 4.30091e-3          # pc (km/s)^2 / Msun
G_MYR = 4.498e-3           # pc^3 / (Msun Myr^2)
KMS_KPC_TO_MYR = 1.0227e-3 # (km/s/kpc) -> 1/Myr
EPS_FF, EPS_CORE, B_TURB, CS, PHI_P = 0.012, 0.5, 0.5, 0.6, 3.0   # CS: measured cold-gas (T<500 K) sound speed, cold_cs.py
EPS_SFE = 0.1              # RC&K17 star formation efficiency per GMC

# ------------------------------------------------------------- theory
def gamma_k12(rho, mach, tfb):
    """Kruijssen 2012 bound fraction (no cruel cradle). rho Msun/pc^3, mach dimensionless, tfb Myr."""
    rho = np.atleast_1d(np.asarray(rho, float)); mach = np.atleast_1d(np.asarray(mach, float))
    out = np.full(np.broadcast(rho, mach).shape, np.nan)
    lx = np.linspace(-8, 12, 600); x = np.exp(lx); dlx = lx[1] - lx[0]
    for i in np.ndindex(out.shape):
        r, mm = np.broadcast_arrays(rho, mach)[0][i], np.broadcast_arrays(rho, mach)[1][i]
        if not (r > 0 and mm > 0): continue
        s2 = np.log(1 + B_TURB ** 2 * mm ** 2); mu = -0.5 * s2
        p = np.exp(-(lx - mu) ** 2 / (2 * s2)) / np.sqrt(2 * np.pi * s2)      # dP/dlnx
        tff = np.sqrt(3 * np.pi / (32 * G_MYR * r * x))
        eps = np.minimum(EPS_FF * tfb / tff, EPS_CORE)
        gam = np.minimum(eps / EPS_CORE, 1.0)
        w = x * p * eps * dlx
        out[i] = np.sum(w * gam) / max(np.sum(w), 1e-300)
    return out

def toomre_mass(Sigma, kappa):
    """Sigma Msun/pc^2, kappa km/s/kpc -> Msun"""
    k = kappa / 1e3   # km/s/pc
    return 4 * np.pi ** 5 * G_PC ** 2 * Sigma ** 3 / np.maximum(k, 1e-9) ** 4

def tff(rho):
    return np.sqrt(3 * np.pi / (32 * G_MYR * np.maximum(rho, 1e-30)))

# ------------------------------------------------------------- join
cl = np.load(f'{C.DATADIR}/clusters/clusters_age0-{a.age:g}_l5_n25.npz')
kap = h5py.File(f'{C.DATADIR}/prfm/kappa.h5', 'r')
rows = []; cl_env = {}
KEYS = ['Sigma_gas', 'Sigma_gas_2p', 'Sigma_star', 'W', 'W_2p', 'W_layer2p', 'P_DE', 'Ptot_2p', 'Pturb_2p', 'rho_mid',
        'rho_mid_2p', 'rho_sd', 'sigma_z_col', 'sigma_eff', 'H', 'f_intruder', 'SigSFR_10', 'SigSFR_40', 'fcold', 'fH2', 'G0_mid']
for pf in sorted(glob.glob(f'{C.PRFM_DIR}/patch_[0-9][0-9][0-9].h5')):
    k = int(pf[-6:-3])
    st = np.load(f'{C.DATADIR}/stars/stars_{k:03d}.npz')
    tk = float(st['time_myr']); age = tk - st['tform_myr']
    ys = (age >= 0) & (age < a.age) & (st['tform_myr'] > 0)
    ypos = st['pos'][ys].astype(np.float64); ym = st['mass'][ys].astype(np.float64) * C.MSUN
    cs = cl['snap'] == k
    cpos = np.stack([cl['x'][cs], cl['y'][cs], cl['z'][cs]], 1); cM = cl['M'][cs]
    cb = (cl['bound'][cs] > 0) & (cl['nucleus'][cs] == 0)
    with h5py.File(pf, 'r') as f:
        P_UNIT = float(f.attrs['P_unit']); PU = P_UNIT
        PATCH = float(f.attrs['patch_kpc']); RMAX = float(f.attrs['rmax_kpc']); ZCOL = float(f.attrs['zcol_kpc']); nb = int(f.attrs['nb'])
        for gname in f:
            if not gname.startswith('frame_'): continue
            g = f[gname]; cen = g.attrs['center']; e1 = g.attrs['e1']; e2 = g.attrs['e2']; nh = g.attrs['nhat']
            sep = float(g.attrs['separation_kpc'])
            def pidx(p):
                d = p - cen; x = d @ e1; y = d @ e2; z = d @ nh
                ix = np.floor((x + RMAX) / PATCH).astype(int); iy = np.floor((y + RMAX) / PATCH).astype(int)
                ok = (ix >= 0) & (ix < nb) & (iy >= 0) & (iy < nb) & (np.abs(z) < ZCOL)
                return np.where(ok, ix * nb + iy, -1)
            npatch = nb * nb
            yi = pidx(ypos); My = np.bincount(yi[yi >= 0], ym[yi >= 0], npatch)
            ci = pidx(cpos) if len(cpos) else np.zeros(0, int)
            Mcl = np.bincount(ci[ci >= 0], cM[ci >= 0], npatch)
            selb = (ci >= 0) & cb
            Mb = np.bincount(ci[selb], cM[selb], npatch); Ncl = np.bincount(ci[selb], None, npatch).astype(float)
            Mmax = np.zeros(npatch)
            if selb.any():
                np.maximum.at(Mmax, ci[selb], cM[selb])
            ix, iy = np.divmod(np.arange(npatch), nb)
            R = np.hypot((ix + 0.5) * PATCH - RMAX, (iy + 0.5) * PATCH - RMAX)
            kg = kap[f'snap_{k:03d}/{gname}'] if f'snap_{k:03d}/{gname}' in kap else None
            kappa = np.interp(R, kg['R'][:], kg['kappa'][:]) if kg is not None else np.full(npatch, np.nan)
            Omega = np.interp(R, kg['R'][:], kg['Omega'][:]) if kg is not None else np.full(npatch, np.nan)
            sigR = np.interp(R, kg['R'][:], kg['sigR'][:]) if kg is not None else np.full(npatch, np.nan)
            vals = {kk: g[kk][:].ravel() for kk in KEYS if kk in g}
            sig_t = np.sqrt(np.maximum(vals['Pturb_2p'] / P_UNIT / np.maximum(vals['rho_mid_2p'] * 1e9, 1e-30), 0))   # km/s
            for j in np.where(ci >= 0)[0]:
                gid = int(np.where(cs)[0][j]); pj = ci[j]
                # a cluster inside both galaxy grids belongs to the frame in which its patch has the LOWER intruder fraction
                # (its own galaxy); the old first-come rule gave galaxy B's clusters galaxy A's column with f_intruder ~ 1
                if gid not in cl_env or vals['f_intruder'][pj] < cl_env[gid][6]:
                    cl_env[gid] = (vals['W'][pj], vals['Sigma_gas'][pj], vals['rho_sd'][pj], vals['H'][pj], vals['sigma_z_col'][pj],
                                   sig_t[pj], vals['f_intruder'][pj], vals['rho_mid_2p'][pj], vals['P_DE'][pj], sep, ord(gname[-1]))
            for p in range(npatch):
                rows.append([k, tk, ord(gname[-1]), sep, R[p], My[p], Mcl[p], Mb[p], Mmax[p], Ncl[p], kappa[p], Omega[p], sigR[p]]
                            + [vals[kk][p] if kk in vals else np.nan for kk in KEYS])
    print(f'snap {k:3d} t={tk:6.1f}: young stars {ym.sum():9.0f} Msun, bound cluster mass {cM[cb].sum():9.0f}', flush=True)
names = ['snap', 't', 'frame', 'sep', 'R', 'M_young', 'M_cl', 'M_bound', 'M_max', 'N_cl', 'kappa', 'Omega', 'sigR'] + KEYS
env_names = ['W', 'Sigma_gas', 'rho_sd', 'H', 'sigma_z', 'sigma_turb', 'f_intruder', 'rho_mid_2p', 'P_DE', 'sep', 'frame']
gids = np.array(sorted(cl_env)); E = {n: np.array([cl_env[g][i] for g in gids]) for i, n in enumerate(env_names)}
E.update({k: cl['' + k][gids] for k in ['snap', 't', 'M', 'N', 'rh_pc', 'age_std', 'bound', 'nucleus']})
np.savez(f'{C.DATADIR}/clusters/clusters_env.npz', **E)
print(f'{len(gids)} of {len(cl["snap"])} young groups fall in a patch column')
D = {n: np.array([r[i] for r in rows], float) for i, n in enumerate(names)}

# ------------------------------------------------------------- theory per patch
D['Gamma'] = D['M_bound'] / np.maximum(D['M_young'], 1e-30)
mach = D['sigma_eff'] / CS
D['Gamma_K12'] = gamma_k12(D['rho_mid_2p'], mach, a.tfb)
# standard route: Q from Sigma_gas, kappa, sigma_eff; rho_ISM = phi_P Omega^2 / (pi G Q^2)
Q = (D['kappa'] * KMS_KPC_TO_MYR) * (D['sigma_eff'] * 1.0227) / (np.pi * G_MYR * D['Sigma_gas'])  # kappa 1/Myr, sigma pc/Myr, G pc^3/(Msun Myr^2)
D['Q'] = Q
Om_myr = D['Omega'] * KMS_KPC_TO_MYR
rho_std = PHI_P * Om_myr ** 2 / (np.pi * G_MYR * np.maximum(Q, 1e-3) ** 2)
D['rho_std'] = rho_std
D['Gamma_std'] = gamma_k12(rho_std, mach, a.tfb)
D['M_T'] = toomre_mass(D['Sigma_gas'], D['kappa'])
D['M_H'] = D['Sigma_gas'] * (D['H'] * 1e3) ** 2                                   # Sigma H^2, measured H
D['M_prfm'] = D['Sigma_gas'] * D['sigma_z_col'] ** 2 / (2 * G_PC * np.maximum(D['rho_sd'], 1e-30))   # DM-dominated PRFM H
D['sigma_turb'] = np.sqrt(np.maximum(D['Pturb_2p'] / PU / np.maximum(D['rho_mid_2p'] * 1e9, 1e-30), 0))   # km/s
D['Mach_cold'] = D['sigma_turb'] / CS; D['Mach_warm'] = D['sigma_turb'] / 8.0
fcoll = np.minimum(1.0, (a.tfb / tff(D['rho_mid_2p'])) ** 4)
D['M_cmax'] = EPS_SFE * D['Gamma_K12'] * fcoll * D['M_T']
fcoll_s = np.minimum(1.0, (a.tfb / tff(rho_std)) ** 4)
D['M_cmax_std'] = EPS_SFE * D['Gamma_std'] * fcoll_s * D['M_T']
np.savez(f'{C.DATADIR}/prfm/patches_clusters.npz', **D)

# ------------------------------------------------------------- selections
clean = (D['f_intruder'] < a.fintr) | C.is_merged(D['sep'], D['t'])      # merged phase: no line-of-sight companion
good = (D['M_young'] > a.mmin_young) & clean & (D['Sigma_gas'] > 1) & np.isfinite(D['W'])
print(f'\n{len(D["snap"])} patch rows, {good.sum()} with M_young>{a.mmin_young:g}, f_intr<{a.fintr}')
Gm = D['Gamma'][good]
print(f'measured bound fraction: median {np.median(Gm):.2f}, mass-weighted {np.sum(D["M_bound"][good])/np.sum(D["M_young"][good]):.2f}')
print(f'K12 (direct rho_mid): median {np.nanmedian(D["Gamma_K12"][good]):.2f};  K12 (standard, gas-only): median {np.nanmedian(D["Gamma_std"][good]):.2f}')

def binned(x, y, nb=10, lo=None, hi=None, stat=np.median):
    x = np.log10(x); ok = np.isfinite(x) & np.isfinite(y)
    e = np.linspace(np.nanpercentile(x[ok], 2) if lo is None else lo, np.nanpercentile(x[ok], 98) if hi is None else hi, nb + 1)
    c = 0.5 * (e[1:] + e[:-1]); m = np.full(nb, np.nan); q1 = m.copy(); q3 = m.copy()
    for i in range(nb):
        s = ok & (x >= e[i]) & (x < e[i + 1])
        if s.sum() >= 5:
            m[i] = stat(y[s]); q1[i] = np.percentile(y[s], 25); q3[i] = np.percentile(y[s], 75)
    return 10 ** c, m, q1, q3

# ---- fig 1: bound fraction vs W, Sigma_gas, P_DE
tcol = D['t'][good]
fig, ax = plt.subplots(1, 3, figsize=(15, 4.4), sharey=True)
for axi, key, lab in zip(ax, ['W', 'Sigma_gas', 'P_DE'], [r'$W/k_B$ [K cm$^{-3}$] (direct gravity)', r'$\Sigma_{\rm gas}$ [M$_\odot$ pc$^{-2}$]', r'$P_{DE}/k_B$ [K cm$^{-3}$]']):
    x = D[key][good]
    sc = axi.scatter(x, Gm, c=tcol, s=6, cmap='cividis', alpha=0.6)
    cx, m, q1, q3 = binned(x, Gm); axi.plot(cx, m, 'k-', lw=2, label='measured, median'); axi.fill_between(cx, q1, q3, color='k', alpha=0.12)
    cx, m, _, _ = binned(x, D['Gamma_K12'][good]); axi.plot(cx, m, 'r--', lw=2, label=r'K12, $\rho_{\rm mid}$ measured')
    cx, m, _, _ = binned(x, D['Gamma_std'][good]); axi.plot(cx, m, 'b:', lw=2, label=r'K12, gas-only $\rho_{\rm ISM}(\Sigma,\kappa,Q)$')
    axi.set(xscale='log', xlabel=lab, ylim=(0, 1.05))
ax[0].set_ylabel(r'bound fraction $\Gamma = M_{\rm bound}/M_{\rm young}$'); ax[0].legend(fontsize=8, loc='upper left')
plt.colorbar(sc, ax=ax[2], label='t [Myr]')
fig.suptitle(f'clean patches (f_intruder<{a.fintr}, M_young>{a.mmin_young:g} Msun), age<{a.age:g} Myr, t_fb={a.tfb:g} Myr')
fig.tight_layout(); fig.savefig(f'{FIG}/cl_fig1_cfe_vs_W.png')

# ---- fig 2: maximum cluster mass: per frame+snapshot observed vs predicted, and vs time
fs = {}
for i in np.where(good)[0]:
    key = (D['snap'][i], D['frame'][i])
    d = fs.setdefault(key, dict(t=D['t'][i], Mmax=0., Mc=0., Mc_std=0., MT=0., Mb=0., My=0., Gk=[], Gs=[], W=[]))
    d['Mmax'] = max(d['Mmax'], D['M_max'][i]); d['Mc'] = np.nanmax([d['Mc'], D['M_cmax'][i]]); d['Mc_std'] = np.nanmax([d['Mc_std'], D['M_cmax_std'][i]])
    d['MT'] = np.nanmax([d['MT'], D['M_T'][i]]); d['Mb'] += D['M_bound'][i]; d['My'] += D['M_young'][i]
    d['Gk'].append(D['Gamma_K12'][i] * D['M_young'][i]); d['Gs'].append(D['Gamma_std'][i] * D['M_young'][i])
F = {k: v for k, v in fs.items() if v['Mmax'] > 0}
T = np.array([v['t'] for v in F.values()]); Mx = np.array([v['Mmax'] for v in F.values()])
Mc = np.array([v['Mc'] for v in F.values()]); Mcs = np.array([v['Mc_std'] for v in F.values()]); MT = np.array([v['MT'] for v in F.values()])
fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
ax[0].scatter(T, Mx, s=12, c='k', label=r'$M_{\rm max}$ bound, per frame')
ax[0].scatter(T, Mc, s=10, c='r', marker='x', label=r'RC&K17 $M_{c,\max}$ ($\rho_{\rm mid}$ measured)')
ax[0].scatter(T, Mcs, s=10, c='b', marker='+', label=r'RC&K17 $M_{c,\max}$ (gas-only)')
ax[0].plot(T[np.argsort(T)], MT[np.argsort(T)], color='grey', lw=0.8, alpha=0.6, label=r'max $M_T$ (Toomre)')
ax[0].set(yscale='log', xlabel='t [Myr]', ylabel=r'M [M$_\odot$]', ylim=(1e2, 1e7)); ax[0].legend(fontsize=7)
ax[1].scatter(Mc, Mx, s=10, c=T, cmap='cividis'); ax[1].plot([1e2, 1e7], [1e2, 1e7], 'k-', lw=0.8)
ax[1].set(xscale='log', yscale='log', xlabel=r'predicted $M_{c,\max}$ ($\rho_{\rm mid}$ measured)', ylabel=r'observed $M_{\rm max}$', xlim=(1e2, 1e7), ylim=(1e2, 1e7))
ok = np.isfinite(Mc) & (Mc > 0)
ax[1].text(0.05, 0.92, f'median log(obs/pred) = {np.median(np.log10(Mx[ok] / Mc[ok])):.2f}', transform=ax[1].transAxes)
G_obs = np.array([v['Mb'] / v['My'] for v in F.values()]); G_k = np.array([np.sum(v['Gk']) / v['My'] for v in F.values()]); G_s = np.array([np.sum(v['Gs']) / v['My'] for v in F.values()])
o = np.argsort(T)
ax[2].plot(T[o], G_obs[o], 'k-', lw=1.2, label='measured (frame total)'); ax[2].plot(T[o], G_k[o], 'r--', lw=1.2, label=r'K12, $\rho_{\rm mid}$ measured')
ax[2].plot(T[o], G_s[o], 'b:', lw=1.2, label='K12, gas-only')
ax[2].set(xlabel='t [Myr]', ylabel=r'bound fraction $\Gamma$', ylim=(0, 1.05)); ax[2].legend(fontsize=8)
fig.tight_layout(); fig.savefig(f'{FIG}/cl_fig2_mmax_cfe_time.png')

# ---- fig 3: which variable does Gamma follow? partial dependence at fixed Sigma_gas
fig, ax = plt.subplots(1, 2, figsize=(10, 4.4), sharey=True)
x = D['Sigma_gas'][good]; wg = D['W'][good] / np.maximum(D['P_DE'][good], 1e-30)
for lo, hi, cc in [(1, 3, 'tab:blue'), (3, 10, 'tab:orange'), (10, 100, 'tab:red')]:
    s = (x >= lo) & (x < hi)
    if s.sum() > 20:
        cx, m, q1, q3 = binned(D['W'][good][s], Gm[s], nb=6); ax[0].plot(cx, m, '-o', color=cc, ms=3, label=rf'$\Sigma_{{gas}}$ {lo}-{hi}')
        cx, m, q1, q3 = binned(wg[s], Gm[s], nb=6); ax[1].plot(cx, m, '-o', color=cc, ms=3)
ax[0].set(xscale='log', xlabel=r'$W/k_B$ [K cm$^{-3}$]', ylabel=r'$\Gamma$', ylim=(0, 1.05)); ax[0].legend(fontsize=8)
ax[1].set(xscale='log', xlabel=r'$W / P_{DE}$ (direct weight / gas-only estimate)')
fig.tight_layout(); fig.savefig(f'{FIG}/cl_fig3_partial.png')
print('figures:', FIG)

# ---- fig 4: which cutoff predictor? frame-level M_max vs max over patches of each predictor
from scipy.stats import spearmanr
preds = [('M_T', r'Toomre $4\pi^5G^2\Sigma^3/\kappa^4$'), ('M_H', r'$\Sigma_{\rm gas} H^2$ (measured $H$)'), ('M_prfm', r'$\Sigma_{\rm gas}\sigma_z^2/(2G\rho_{sd})$ (PRFM, DM-dom.)')]
fl = {}
for i in np.where(good)[0]:
    key = (D['snap'][i], D['frame'][i]); d = fl.setdefault(key, dict(t=D['t'][i], Mmax=0., **{pn: 0. for pn, _ in preds}))
    d['Mmax'] = max(d['Mmax'], D['M_max'][i])
    for pn, _ in preds: d[pn] = np.nanmax([d[pn], D[pn][i]])
FL = {k: v for k, v in fl.items() if v['Mmax'] > 0}
Tf = np.array([v['t'] for v in FL.values()]); Mxf = np.array([v['Mmax'] for v in FL.values()])
fig, ax = plt.subplots(1, 3, figsize=(15, 4.4), sharey=True)
print('\ncutoff predictors (frame level, patch level):')
for axi, (pn, lab) in zip(ax, preds):
    X = np.array([v[pn] for v in FL.values()]); ok = np.isfinite(X) & (X > 0)
    rs = spearmanr(np.log10(X[ok]), np.log10(Mxf[ok]))[0]
    sc = axi.scatter(X[ok], Mxf[ok], c=Tf[ok], s=12, cmap='cividis')
    # patch level (grey)
    pk = good & (D['M_max'] > 0) & np.isfinite(D[pn]) & (D[pn] > 0)
    axi.scatter(D[pn][pk], D['M_max'][pk], s=3, c='grey', alpha=0.3, zorder=0)
    rp = spearmanr(np.log10(D[pn][pk]), np.log10(D['M_max'][pk]))[0]
    axi.plot([1e2, 1e8], [1e2, 1e8], 'k-', lw=0.8); axi.plot([1e2, 1e8], [1e1, 1e7], 'k:', lw=0.8)
    axi.set(xscale='log', yscale='log', xlabel=lab, xlim=(1e2, 1e8), ylim=(1e2, 1e7))
    axi.text(0.04, 0.92, f'Spearman: frames {rs:.2f}, patches {rp:.2f}', transform=axi.transAxes, fontsize=8)
    print(f'  {pn:8s}: rho_frame={rs:.2f} (N={ok.sum()}), rho_patch={rp:.2f} (N={pk.sum()}), median log(obs/pred) frames = {np.median(np.log10(Mxf[ok]/X[ok])):.2f}')
ax[0].set_ylabel(r'observed $M_{\rm max}$ (bound, age<%g Myr)' % a.age); plt.colorbar(sc, ax=ax[2], label='t [Myr]')
fig.tight_layout(); fig.savefig(f'{FIG}/cl_fig4_cutoff_predictors.png')

# ---- fig 5: slope invariance: alpha in W / Mach terciles, and per epoch vs Mach
def mle(Mv, mmin=300.0):
    Mv = Mv[Mv >= mmin]
    if len(Mv) < 8: return np.nan, np.nan, len(Mv)
    al = 1 + len(Mv) / np.sum(np.log(Mv / mmin)); return al, (al - 1) / np.sqrt(len(Mv)), len(Mv)
sel = (E['bound'] > 0) & (E['nucleus'] == 0) & ((E['f_intruder'] < a.fintr) | C.is_merged(E['sep'], E['t'])) & (np.mod(E['snap'], 10) == 0)
print(f'\nslope vs environment ({sel.sum()} bound young clusters in clean patches, every 10th snapshot):')
fig, ax = plt.subplots(1, 3, figsize=(15, 4.2), sharey=True)
for axi, key, lab in zip(ax[:2], ['W', 'sigma_turb'], [r'$W/k_B$ [K cm$^{-3}$]', r'$\sigma_{\rm turb}$ [km/s] ($\mathcal{M}_{\rm cold}=\sigma/0.2$)']):
    x = E[key][sel]; Mv = E['M'][sel]
    qs = np.nanpercentile(x, [0, 25, 50, 75, 100])
    for i in range(4):
        s_ = (x >= qs[i]) & (x <= qs[i + 1])
        al, er, n = mle(Mv[s_]); xc = np.nanmedian(x[s_])
        axi.errorbar(xc, al, er, fmt='o', color='k', capsize=3)
        print(f'  {key} quartile {i+1} (median {xc:.3g}): alpha = {al:.2f} +- {er:.2f} (N={n})')
    axi.axhline(2, color='r', ls='--', lw=1); axi.set(xscale='log', xlabel=lab, ylim=(1.2, 2.8))
ax[0].set_ylabel(r'$\alpha$ (MLE, $M>300$ M$_\odot$)')
# per epoch vs epoch-median Mach of the clean patches
for lo, hi in [(0, 50), (50, 100), (100, 150), (150, 232)]:
    s_ = sel & (E['t'] >= lo) & (E['t'] < hi); pg = good & (D['t'] >= lo) & (D['t'] < hi)
    al, er, n = mle(E['M'][s_]); mm = np.nanmedian(D['Mach_cold'][pg])
    if np.isfinite(al):
        ax[2].errorbar(mm, al, er, fmt='s', capsize=3, label=f'{lo}-{hi} Myr (N={n})')
    print(f'  epoch {lo}-{hi}: alpha={al:.2f}+-{er:.2f} (N={n}), median Mach_cold={mm:.1f}, median W={np.nanmedian(D["W"][pg]):.3g}')
ax[2].axhline(2, color='r', ls='--', lw=1); ax[2].set(xlabel=r'epoch median $\mathcal{M}_{\rm cold}$', xscale='log'); ax[2].legend(fontsize=8)
fig.tight_layout(); fig.savefig(f'{FIG}/cl_fig5_slope_vs_pressure.png')
print('figures:', FIG)
