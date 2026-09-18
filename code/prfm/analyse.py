"""Aggregate the per-snapshot patch files and make the PRFM diagnostic figures.

usage: python analyse.py            -> /ptmp/uli/dwarf_merger/prfm/patches_all.npz + figures in prfm/figs/
"""
import sys, os, glob
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, h5py
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import common as C, ok22

FIGDIR = f'{C.DATADIR}/prfm/figs'; os.makedirs(FIGDIR, exist_ok=True)
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.grid': True, 'grid.alpha': 0.25, 'grid.linewidth': 0.5, 'legend.frameon': False,
                     'figure.dpi': 130})
TCOL = 'cividis'
SIG_MIN = 1.0        # Msun/pc^2, patches considered
NSLAB_MIN = 50

def load_all():
    files = sorted(glob.glob(f'{C.DATADIR}/prfm/patch_*.h5'))
    rows = {}; meta = []
    for fn in files:
        with h5py.File(fn, 'r') as f:
            t = float(f.attrs['time_myr']); snap = int(fn[-6:-3])
            for gname in f:
                g = f[gname]
                sep = float(g.attrs['separation_kpc'])
                meta.append((snap, t, gname[-1], sep))
                n = g['Sigma_gas'].shape[0] * g['Sigma_gas'].shape[1]
                for k in g:
                    rows.setdefault(k, []).append(g[k][:].ravel())
                rows.setdefault('t', []).append(np.full(n, t))
                rows.setdefault('snap', []).append(np.full(n, snap))
                rows.setdefault('sep', []).append(np.full(n, sep))
                rows.setdefault('frame', []).append(np.full(n, ord(gname[-1])))
    d = {k: np.concatenate(v) for k, v in rows.items()}
    return d, meta

def sfr_global(times, cat):
    """global SFR (Msun/yr) from stars younger than 10 / 40 Myr, and SN rate per Myr (10 Myr window)."""
    tf = cat['tform_myr']; tsn = cat['tsn_myr']
    out = {}
    for w in (10, 40):
        out[f'sfr{w}'] = np.array([np.sum((tf > t - w) & (tf <= t)) * 3.98 / (w * 1e6) for t in times])
    out['snrate'] = np.array([np.sum((tsn > t - 10) & (tsn <= t)) / 10.0 for t in times])
    return out

def main():
    d, meta = load_all()
    np.savez_compressed(f'{C.DATADIR}/prfm/patches_all.npz', **d)
    cat = np.load(f'{C.DATADIR}/stars/catalog.npz')
    snaps = sorted(set(m[0] for m in meta)); tsn = np.array([[m[1] for m in meta if m[0] == s][0] for s in snaps])
    sep = np.array([[m[3] for m in meta if m[0] == s][0] for s in snaps])
    glob_ = sfr_global(tsn, cat)
    good = (d['Sigma_gas'] > SIG_MIN) & (d['Npart_slab'] > NSLAB_MIN) & (d['W_2p'] > 0) & (d['Ptot_2p'] > 0)
    print(f'{len(snaps)} snapshots, {good.sum()} usable patches of {len(good)}')
    W = d['W_2p']; P = d['Ptot_2p']; PDE = d['P_DE']; t = d['t']
    S40 = d['SigSFR_40']; S10 = d['SigSFR_10']

    # ---------------- Fig 1: global time series
    fig, ax = plt.subplots(4, 1, figsize=(8, 10), sharex=True)
    ax[0].plot(tsn, glob_['sfr10'], color='#1f4e79', lw=1.5, label='stars < 10 Myr')
    ax[0].plot(tsn, glob_['sfr40'], color='#7fa7d1', lw=1.5, label='stars < 40 Myr')
    ax[0].set_ylabel('SFR [Msun/yr]'); ax[0].set_yscale('log'); ax[0].legend()
    ax[1].plot(tsn, glob_['snrate'], color='#8c4a1f', lw=1.5); ax[1].set_ylabel('SN rate [Myr$^{-1}$]'); ax[1].set_yscale('log')
    ax[2].plot(tsn, sep, color='k', lw=1.5); ax[2].set_ylabel('separation [kpc]')
    # Sigma_gas-weighted mean pressures over usable patches
    def wmean(q, m):
        return np.array([np.average(q[m & (d['snap'] == s)], weights=d['Sigma_gas'][m & (d['snap'] == s)])
                         if np.any(m & (d['snap'] == s)) else np.nan for s in snaps])
    for q, lab, col in ((W, 'W (direct gravity, 2p)', 'k'), (P, 'P_tot,2p', '#c0392b'), (PDE, 'P_DE (analytic)', '#7f8c8d'),
                        (d['Pth_2p'], 'P_th,2p', '#e67e22'), (d['Pturb_2p'], 'P_turb,2p', '#2980b9'), (d['Pmag_2p'], 'Pi_mag,2p', '#8e44ad')):
        ax[3].plot(tsn, wmean(q, good), color=col, lw=1.3 if lab.startswith('P_t') or lab.startswith('W') else 1.0, label=lab,
                   ls='-' if 'W' in lab or 'tot' in lab or 'DE' in lab else '--')
    ax[3].set_yscale('log'); ax[3].set_ylabel('P/k$_B$ [K cm$^{-3}$]'); ax[3].legend(ncol=2, fontsize=8); ax[3].set_xlabel('t [Myr]')
    fig.tight_layout(); fig.savefig(f'{FIGDIR}/fig1_timeseries.png'); plt.close(fig)

    # ---------------- Fig 2: vertical equilibrium P_tot vs W ; P_DE vs W
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.4))
    for a, x, y, xl, yl in ((ax[0], W[good], P[good], 'W$_{2p}$/k$_B$ (direct)', 'P$_{tot,2p}$/k$_B$'),
                            (ax[1], W[good], PDE[good], 'W$_{2p}$/k$_B$ (direct)', 'P$_{DE}$/k$_B$ (analytic)'),
                            (ax[2], PDE[good], P[good], 'P$_{DE}$/k$_B$', 'P$_{tot,2p}$/k$_B$')):
        sc = a.scatter(x, y, c=t[good], cmap=TCOL, s=4, alpha=0.5, linewidths=0)
        lim = [1e2, 1e7]; a.plot(lim, lim, 'k-', lw=1); a.plot(lim, np.array(lim) * 3, 'k:', lw=0.8); a.plot(lim, np.array(lim) / 3, 'k:', lw=0.8)
        a.set_xscale('log'); a.set_yscale('log'); a.set_xlim(lim); a.set_ylim(lim); a.set_xlabel(xl); a.set_ylabel(yl)
        r = np.log10(y / x); a.text(0.03, 0.93, f'median log ratio = {np.nanmedian(r):.2f}, scatter {np.nanstd(r):.2f} dex', transform=a.transAxes, fontsize=8)
    fig.colorbar(sc, ax=ax, label='t [Myr]', pad=0.01)
    fig.savefig(f'{FIGDIR}/fig2_equilibrium.png', bbox_inches='tight'); plt.close(fig)

    # ---------------- Fig 3: feedback yields vs P_DE with OK22 fits
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.4))
    x = np.logspace(2.5, 6.5, 50)
    for a, Pc, fit, lab in ((ax[0], d['Pth_2p'], ok22.ups_th, r'$\Upsilon_{th}$'), (ax[1], d['Pturb_2p'], ok22.ups_turb, r'$\Upsilon_{turb}$'),
                            (ax[2], P, ok22.ups_tot, r'$\Upsilon_{tot}$')):
        m = good & (S40 > 0)
        a.scatter(PDE[m], Pc[m] / S40[m] * P_to_ups(), c=t[m], cmap=TCOL, s=4, alpha=0.5, linewidths=0)
        a.plot(x, fit(x), 'k-', lw=1.5, label='OK22 TIGRESS fit')
        a.set_xscale('log'); a.set_yscale('log'); a.set_xlabel('P$_{DE}$/k$_B$ [K cm$^{-3}$]'); a.set_ylabel(lab + ' [km/s]  (SFR: stars < 40 Myr)')
        a.set_ylim(10, 3e4); a.legend()
    fig.savefig(f'{FIGDIR}/fig3_yields.png', bbox_inches='tight'); plt.close(fig)

    # ---------------- Fig 4: Sigma_SFR vs W with OK22
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.4))
    for a, S, lab in ((ax[0], S40, '40 Myr'), (ax[1], S10, '10 Myr')):
        m = good
        a.scatter(W[m], np.maximum(S[m], 1e-6), c=t[m], cmap=TCOL, s=4, alpha=0.5, linewidths=0)
        a.plot(x, ok22.sfr_of_W(x), 'k-', lw=1.5, label='OK22: $\\Sigma_{SFR}$ = W/$\\Upsilon_{tot}$')
        a.set_xscale('log'); a.set_yscale('log'); a.set_xlabel('W$_{2p}$/k$_B$ [K cm$^{-3}$]'); a.set_ylabel(f'$\\Sigma_{{SFR}}$ [Msun/kpc$^2$/yr] (stars < {lab})')
        a.set_ylim(1e-6, 10); a.legend()
    fig.savefig(f'{FIGDIR}/fig4_sfr_vs_W.png', bbox_inches='tight'); plt.close(fig)

    # ---------------- Fig 5: ratios vs time (median + 16-84 over patches)
    fig, ax = plt.subplots(3, 1, figsize=(8, 8), sharex=True)
    def band(a, q, m, col, lab):
        med = []; lo = []; hi = []
        for s in snaps:
            sel = m & (d['snap'] == s)
            if sel.sum() < 5:
                med.append(np.nan); lo.append(np.nan); hi.append(np.nan); continue
            p = np.nanpercentile(q[sel], [16, 50, 84]); lo.append(p[0]); med.append(p[1]); hi.append(p[2])
        a.plot(tsn, med, color=col, lw=1.5, label=lab); a.fill_between(tsn, lo, hi, color=col, alpha=0.2, lw=0)
    band(ax[0], P / W, good, '#c0392b', 'P$_{tot,2p}$ / W$_{2p}$'); band(ax[0], PDE / W, good, '#7f8c8d', 'P$_{DE}$ / W$_{2p}$')
    ax[0].axhline(1, color='k', lw=0.8); ax[0].set_yscale('log'); ax[0].set_ylim(0.1, 10); ax[0].legend(); ax[0].set_ylabel('pressure / weight')
    m = good & (S40 > 0)
    band(ax[1], (P / np.maximum(S40, 1e-30) * P_to_ups()) / ok22.ups_tot(PDE), m, '#1f4e79', '$\\Upsilon_{tot}$ / OK22 (SFR 40 Myr)')
    m10 = good & (S10 > 0)
    band(ax[1], (P / np.maximum(S10, 1e-30) * P_to_ups()) / ok22.ups_tot(PDE), m10, '#7fa7d1', '$\\Upsilon_{tot}$ / OK22 (SFR 10 Myr)')
    ax[1].axhline(1, color='k', lw=0.8); ax[1].set_yscale('log'); ax[1].set_ylim(0.03, 30); ax[1].legend(); ax[1].set_ylabel('yield / TIGRESS fit')
    band(ax[2], d['Pturb_2p'] / np.maximum(P, 1e-30), good, '#2980b9', 'P$_{turb}$ / P$_{tot}$')
    band(ax[2], d['Pth_2p'] / np.maximum(P, 1e-30), good, '#e67e22', 'P$_{th}$ / P$_{tot}$')
    band(ax[2], d['Pmag_2p'] / np.maximum(P, 1e-30), good, '#8e44ad', '$\\Pi_{mag}$ / P$_{tot}$')
    ax[2].set_ylim(0, 1.05); ax[2].legend(ncol=3); ax[2].set_ylabel('fraction'); ax[2].set_xlabel('t [Myr]')
    fig.tight_layout(); fig.savefig(f'{FIGDIR}/fig5_ratios_vs_time.png'); plt.close(fig)

    # ---------------- Fig 6: turbulence vs tidal shear vs SN rate
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.4))
    m = good
    sc = ax[0].scatter(d['SNrate'][m] + 1e-3, d['Pturb_2p'][m], c=np.log10(d['shear'][m] + 1e-3), cmap='magma', s=4, alpha=0.6, linewidths=0)
    ax[0].set_xscale('log'); ax[0].set_yscale('log'); ax[0].set_xlabel('SN rate [Myr$^{-1}$ kpc$^{-2}$] (last 10 Myr)'); ax[0].set_ylabel('P$_{turb,2p}$/k$_B$')
    fig.colorbar(sc, ax=ax[0], label='log shear [km/s/kpc]')
    sc = ax[1].scatter(d['shear'][m], d['Pturb_2p'][m], c=t[m], cmap=TCOL, s=4, alpha=0.6, linewidths=0)
    ax[1].set_xscale('log'); ax[1].set_yscale('log'); ax[1].set_xlabel('shear |S| [km/s/kpc]'); ax[1].set_ylabel('P$_{turb,2p}$/k$_B$')
    fig.colorbar(sc, ax=ax[1], label='t [Myr]')
    fig.savefig(f'{FIGDIR}/fig6_turb_drivers.png', bbox_inches='tight'); plt.close(fig)
    print('figures written to', FIGDIR)

def P_to_ups():
    """K cm^-3 per (Msun/kpc^2/yr) -> km/s :  P k_B / Sigma_SFR."""
    kB = 1.3807e-16; msun = 1.989e33; kpc = 3.085678e21; yr = 3.15576e7
    return kB / (msun / kpc**2 / yr) / 1e5

if __name__ == '__main__':
    main()
