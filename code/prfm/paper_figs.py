"""Paper-style versions of the key figures, one function each.  usage: python paper_figs.py NAME [NAME ...]  (or 'all')
Outputs: /ptmp/uli/dwarf_merger/prfm/figs/paper/NAME.{png,pdf} and ~/merger_prfm/figs/NAME.{png,pdf}"""
import os, sys, glob, numpy as np, h5py
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C
import ok22
import paperstyle as P
P.use(); import matplotlib.pyplot as plt
FIGS = {}
def figure(f): FIGS[f.__name__] = f; return f

@figure
def story_0_sep():
    L = np.load(f'{C.DATADIR}/prfm/layer_verdict_z05.npz')
    fig, ax = P.fig()
    ax.plot(L['TF'], L['sepF'], color=P.blue)
    ax.set_xlim(0, 226); ax.set_ylim(0, 4.5)
    ax.set_xlabel(r'$t$ [Myr]'); ax.set_ylabel(r'$d_{\rm nuclei}$ [kpc]')
    P.phases(ax)
    P.save(fig, 'story_0_sep'); plt.close(fig)

@figure
def story_1_equilibrium():
    L = np.load(f'{C.DATADIR}/prfm/layer_verdict_z05.npz')
    fig, ax = P.fig()
    ax.plot(L['TF'], L['mwF'], color=P.light, label=r'full column ($|z| < 1.5$ kpc)')
    ax.plot(L['TL'], L['mwL'], color=P.blue, lw=1.5, label=r'gas layer ($|z| < 0.5$ kpc)')
    import os as _os
    if _os.path.exists(f'{C.DATADIR}/prfm/layer_verdict_z05_twoframe.npz'):
        L2 = np.load(f'{C.DATADIR}/prfm/layer_verdict_z05_twoframe.npz'); o2 = L2['TL'] >= 106
        ax.plot(L2['TL'][o2], L2['mwL'][o2], color=P.grey, lw=0.8, ls=':', label='layer, two grids after 2nd passage')
    ax.axhline(1, color=P.ink, lw=0.7)
    ax.set_yscale('log'); ax.set_xlim(0, 226); ax.set_ylim(0.15, 3)
    ax.set_xlabel(r'$t$ [Myr]'); ax.set_ylabel(r'$P_{\rm tot,2p}\,/\,\mathcal{W}_{\rm 2p}$')
    ax.set_yticks([0.2, 0.5, 1, 2]); ax.set_yticklabels(['0.2', '0.5', '1', '2'])
    P.phases(ax); P.place_legend(ax, loc='lower left')
    P.save(fig, 'story_1_equilibrium'); plt.close(fig)

@figure
def story_2_feedback():
    """effective yield P_tot / Sigma_SFR (measured) vs time, against the Ostriker & Kim 2022 feedback yield at the same pressure"""
    def series(suffix=''):
        D = dict(np.load(f'{C.DATADIR}/prfm/patches_clusters{suffix}.npz')); K = np.load(f'{C.DATADIR}/prfm/cluster_link{suffix}.npz'); PL = K['PL']; WL = K['WL']
        clean = (D['f_intruder'] < 0.1) | C.is_merged(D['sep'], D['t']); base = clean & (D['Sigma_gas'] > 1) & np.isfinite(PL) & (PL > 0) & np.isfinite(WL)
        T = np.unique(D['snap'][base]); t = []; ups = []; ok22y = []
        for k in T:
            m = base & (D['snap'] == k); sfr = D['SigSFR_10'][m]; Pk = PL[m]
            if sfr.sum() <= 0: continue
            t.append(D['t'][m][0]); ups.append(Pk.sum() / sfr.sum() / 4.81e3); ok22y.append(ok22.ups_tot(np.average(Pk, weights=D['Sigma_gas_2p'][m])))
        return np.array(t), np.array(ups), np.array(ok22y)
    t, ups, ok22y = series(); o = t > 10
    fig, ax = P.fig()
    import os as _os
    if _os.path.exists(f'{C.DATADIR}/prfm/cluster_link_twoframe.npz'):
        t2, ups2, _ = series('_twoframe'); o2 = t2 >= 106
        ax.plot(t2[o2], ups2[o2], color=P.grey, lw=0.8, ls=':', label='two grids after 2nd passage')
    ax.plot(t[o], ups[o], color=P.blue, lw=1.4, label=r'measured: $P_{\rm tot,2p}\,/\,\Sigma_{\rm SFR}$')
    ax.plot(t[o], ok22y[o], color=P.ink, lw=1.0, ls='--', label=r'feedback yield $\Upsilon_{\rm tot}(P)$ (OK22)')
    ax.set_yscale('log'); ax.set_xlim(10, 226); ax.set_ylim(2e2, 3e5)
    ax.set_xlabel(r'$t$ [Myr]'); ax.set_ylabel(r'$P_{\rm tot,2p}\,/\,\Sigma_{\rm SFR}$ [km s$^{-1}$]')
    ax.set_ylim(20, 1e6); P.phases(ax, box=True); P.place_legend(ax, loc='lower left')
    P.save(fig, 'story_2_feedback'); plt.close(fig)

@figure
def story_3_shear():
    rows = []
    for fn in sorted(glob.glob(f'{C.DATADIR}/prfm/turb_*.npz')):
        d = dict(np.load(fn)); k = int(d['snap'][0])
        with h5py.File(f'{C.DATADIR}/prfm/patch_{k:03d}.h5', 'r') as f:   # the turbulence decomposition was run on the two-frame grids; H from the same files
            for fr in np.unique(d['frame']):
                g = f[fr]; m = d['frame'] == fr; pi = d['patch'][m].astype(int); H = g['H'][:].ravel()[pi]
                for j, i in enumerate(np.flatnonzero(m)): rows.append([d['t'][i], d['sep'][i], d['Sigma_gas'][i], d['Pturb_2p'][i], d['f_intruder'][i], d['SigSFR_10'][i], d['sig_tot'][i], d['sig_res'][i], d['shear_patch'][i], H[j]])
    R = np.array(rows); t, sep, Sig, Pt, fi, sfr, st, sr, Sp, H = R.T
    ok = ((fi < 0.1) | C.is_merged(sep, t)) & (Pt > 0) & (st > 0) & (t > 10) & (sfr == 0) & (Sp > 0)
    X = Sig * 1e6 * H * Sp ** 2 * 4.903245584337325e-06; Pres = Pt * (sr / st) ** 2
    fig, ax = P.fig(aspect=0.8)
    ax.scatter(X[ok], Pres[ok], s=5, c=P.blue, alpha=0.4, lw=0, label='quiescent patches', rasterized=True)
    xx = np.logspace(3, 7, 50); ax.plot(xx, 0.02 * xx, color=P.ink, lw=1.2, label=r'$P_{\rm turb} = 0.02\,\Sigma H S^2$')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(1e3, 2e6); ax.set_ylim(1e1, 3e5)
    ax.set_xlabel(r'$\Sigma\,H\,S^2$ [K cm$^{-3}$]'); ax.set_ylabel(r'$P_{\rm turb,2p}$ [K cm$^{-3}$]'); P.place_legend(ax, loc='upper left')
    P.save(fig, 'story_3_shear'); plt.close(fig)

@figure
def story_4_dynamo():
    Rb = np.load(f'{C.DATADIR}/prfm/bfield_time.npz')['rows']   # snap, t, <B>_mass [muG], median B, median Pmag/(Pth+Pturb), stress share
    fig, ax = P.fig()
    ax.plot(Rb[:, 1], Rb[:, 2], color=P.blue, lw=1.5, label='mass-weighted mean')
    ax.plot(Rb[:, 1], Rb[:, 3], color=P.light, lw=1.2, label='median patch')
    ax.set_yscale('log'); ax.set_xlim(10, 226); ax.set_ylim(3e-3, 300)
    ax.set_xlabel(r'$t$ [Myr]'); ax.set_ylabel(r'$B$ [$\mu$G]')
    P.phases(ax); P.place_legend(ax, loc='lower right')
    P.save(fig, 'story_4_dynamo'); plt.close(fig)

@figure
def clouds_mf():
    d = np.load(f'{C.DATADIR}/clouds/clouds_analysis.npz'); R = d['cl10']
    fig, ax = P.fig(aspect=0.8)
    edges = np.logspace(2, 7, 21); c = np.sqrt(edges[1:] * edges[:-1]); w = np.diff(edges)
    for (t0, t1, lab), col in zip(P.PHASES, (P.light, P.blue, P.orange, P.ink)):
        m = (R[:, 1] >= t0) & (R[:, 1] < t1); h, _ = np.histogram(R[m, 2], edges); ok = h > 0
        ax.errorbar(c[ok], h[ok] / w[ok] / m.sum(), yerr=np.sqrt(h[ok]) / w[ok] / m.sum(), fmt='o-', ms=2.5, lw=1.0, color=col, label=lab.replace('\n', ' '))
    xx = np.array([3e2, 3e5]); ax.plot(xx, 3e-2 * (xx / 3e2) ** -1.6, color=P.grey, lw=0.8, ls='--', label=r'$M^{-1.6}$')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(1e2, 1e7); ax.set_ylim(1e-10, 1e-1)
    ax.set_xlabel(r'$M_{\rm cloud}$ [M$_\odot$]'); ax.set_ylabel(r'$N^{-1}\,{\rm d}N/{\rm d}M$ [M$_\odot^{-1}$]'); P.place_legend(ax, loc='lower left')
    P.save(fig, 'clouds_mf'); plt.close(fig)

@figure
def clouds_eff():
    """lifetime-integrated efficiency of dense-clump lineages vs the virial parameter before star formation starts"""
    d = np.load(f'{C.DATADIR}/clouds/lineages_cl100.npz'); L = d['rows']; col = {n: i for i, n in enumerate(d['names'])}
    g = (L[:, col['Mstar']] > 0) & (L[:, col['M_max']] >= 300) & (L[:, col['t_start']] > 5) & (L[:, col['t_end']] < 221) & (L[:, col['alpha_on']] > 0)
    e = L[:, col['eps_int']]; al = L[:, col['alpha_on']]; ton = L[:, col['t_onset']]
    fig, ax = P.fig(aspect=0.8); ab = np.array([0.5, 1, 2, 4, 8, 16, 32, 128]); xc = np.sqrt(ab[1:] * ab[:-1])
    for (t0, t1, lab), colr in zip(P.PHASES, (P.light, P.blue, P.orange, P.ink)):
        m = g & (ton >= t0) & (ton < t1)
        ax.scatter(al[m], e[m], s=4, color=colr, alpha=0.3, lw=0, rasterized=True)
        med = np.array([np.median(e[m & (al >= a) & (al < b)]) if np.sum(m & (al >= a) & (al < b)) >= 8 else np.nan for a, b in zip(ab[:-1], ab[1:])])
        ax.plot(xc, med, 'o-', ms=3, lw=1.3, color=colr, label=lab.replace('\n', ' '))
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(0.5, 200); ax.set_ylim(1e-4, 30)
    ax.set_xlabel(r'$\alpha_{\rm vir}$ before onset'); ax.set_ylabel(r'$\epsilon_{\rm int} = M_\star\,/\,M_{\rm clump}$'); P.place_legend(ax, loc='upper right')
    P.save(fig, 'clouds_eff'); plt.close(fig)

@figure
def story_5_mf():
    E = np.load(f'{C.DATADIR}/clusters/clusters_env.npz'); cleanE = (E['f_intruder'] < 0.1) | C.is_merged(E['sep'], E['t'])
    cb = cleanE & (E['bound'] > 0) & (E['nucleus'] == 0) & (np.mod(E['snap'], 10) == 0)   # independent populations, every 10th snapshot
    fig, ax = P.fig(aspect=0.8); edges = np.logspace(2, 5.6, 13); dl = np.diff(np.log10(edges)); x = np.sqrt(edges[:-1] * edges[1:])
    for (t0, t1, lab), col in zip(P.PHASES, (P.light, P.blue, P.orange, P.ink)):
        M = E['M'][cb & (E['t'] >= t0) & (E['t'] < t1)]; n, _ = np.histogram(M, edges); k = n > 0
        Mf = M[M >= 300]; al = 1 + len(Mf) / np.sum(np.log(Mf / 300.))
        ax.errorbar(x[k], n[k] / dl[k], yerr=np.sqrt(n[k]) / dl[k], fmt='o', ms=3, color=col, lw=1, label=lab.replace('\n', ' ') + rf'  ($\alpha = {al:.2f}$)')
        xx = np.logspace(np.log10(300), np.log10(M.max()), 50); ax.plot(xx, len(Mf) * (al - 1) * 300 ** (al - 1) * xx ** (1 - al) * np.log(10), color=col, lw=1.0)
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(1e2, 5e5); ax.set_ylim(0.5, 3e4)
    ax.set_xlabel(r'$M_{\rm group}$ [M$_\odot$]'); ax.set_ylabel(r'${\rm d}N/{\rm d}\log M$'); P.place_legend(ax, loc='upper right')
    P.save(fig, 'story_5_mf'); plt.close(fig)

@figure
def story_6_reservoir():
    D = dict(np.load(f'{C.DATADIR}/prfm/patches_clusters.npz')); clean = (D['f_intruder'] < 0.1) | C.is_merged(D['sep'], D['t'])
    m = clean & (D['M_young'] > 500) & (D['M_max'] > 0) & (D['Sigma_gas'] > 1) & (np.mod(D['snap'], 10) == 0)
    x = D['M_young'][m]; y = D['M_max'][m]
    fig, ax = P.fig(aspect=0.8)
    ax.scatter(x, y, s=5, c=P.blue, alpha=0.45, lw=0, rasterized=True, label='columns with a bound group')
    xb = np.logspace(np.log10(500), 6.5, 13); med = [np.median(y[(x >= a) & (x < b)]) if np.sum((x >= a) & (x < b)) >= 8 else np.nan for a, b in zip(xb[:-1], xb[1:])]
    ax.plot(np.sqrt(xb[:-1] * xb[1:]), med, 'o-', ms=3, color=P.ink, lw=1.2, label='median')
    xx = np.logspace(2.7, 6.5, 10); ax.plot(xx, 0.5 * xx, color=P.orange, lw=1.0, label=r'$M_{\max} = 0.5\,M_{\rm young}$'); ax.plot(xx, xx, color=P.grey, lw=0.8, ls=':', label=r'$M_{\max} = M_{\rm young}$')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(5e2, 3e6); ax.set_ylim(1e2, 3e6)
    ax.set_xlabel(r'$M_{\rm young}$ [M$_\odot$]'); ax.set_ylabel(r'$M_{\max}$ [M$_\odot$]'); P.place_legend(ax, loc='upper left')
    P.save(fig, 'story_6_reservoir'); plt.close(fig)

@figure
def story_7_ceiling():
    """largest bound cluster per 25 Myr vs the 90th-percentile layer weight of the star-forming patches (cluster_link.npz, corrected join)"""
    R = np.load(f'{C.DATADIR}/prfm/cluster_link.npz')['bins']          # t0, t1, n_cl, M_max, alpha, Mc, n_sf, WL med, WL 90%, PL med, Sigma med, W_full med
    D = dict(np.load(f'{C.DATADIR}/prfm/patches_clusters.npz')); WL = np.load(f'{C.DATADIR}/prfm/cluster_link.npz')['WL']
    clean = (D['f_intruder'] < 0.1) | C.is_merged(D['sep'], D['t']); base = clean & (D['Sigma_gas'] > 1) & np.isfinite(WL) & (WL > 0) & (D['t'] > 10)
    fb = np.array([np.mean(D['M_young'][base & (D['t'] >= r[0]) & (D['t'] < r[1])] > 500) for r in R]); act = fb >= 0.05
    W = R[:, 8]; M = R[:, 3]
    c = np.polyfit(np.log10(W[act]), np.log10(M[act]), 1); res = np.log10(M[act]) - np.polyval(c, np.log10(W[act]))
    print(f'ceiling: {act.sum()} active bins, log M_max = {c[1]:.2f} + {c[0]:.2f} log W_90, scatter {np.std(res):.2f} dex; quiescent bin(s): {R[~act, 0]}')
    fig, ax = P.fig(aspect=0.8)
    Wx = np.logspace(3.5, 5.8, 50); ax.plot(Wx, 0.5 * ok22.sfr_of_W(Wx) * 0.25 * 1e7, color=P.ink, lw=1.2, label=r'$0.5\,\Sigma_{\rm SFR}(\mathcal{W})\,A\,\tau$ (PRFM, no fit)')
    ax.scatter(W[act], M[act], s=28, c=P.blue, zorder=3, label='largest group per 25 Myr')
    ax.scatter(W[~act], M[~act], s=40, facecolors='none', edgecolors=P.orange, lw=1.2, zorder=3, label='quiescent interval')
    for r in R: ax.annotate(f'{r[0]:.0f}', (r[8], r[3]), xytext=(4, -2), textcoords='offset points', fontsize=6.5, color=P.grey)
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(3e3, 8e5); ax.set_ylim(5e2, 5e6)
    ax.set_xlabel(r'$\mathcal{W}_{\rm 2p}$ [K cm$^{-3}$]'); ax.set_ylabel(r'$M_{\max}$ [M$_\odot$]'); P.place_legend(ax, loc='upper left')
    P.save(fig, 'story_7_ceiling'); plt.close(fig)

@figure
def cluster_env():
    """every bound young cluster vs the state of the patch it formed in (clusters_env_layer.npz)"""
    d = np.load(f'{C.DATADIR}/clusters/clusters_env_layer.npz'); R = d['rows']; col = {n: i for i, n in enumerate(d['names'])}
    M = R[:, col['M']]; t = R[:, col['t']]
    fig, ax = plt.subplots(1, 3, figsize=(P.COL2, P.COL2 * 0.34), sharey=True)
    for (t0, t1, lab), colr in zip(P.PHASES, (P.light, P.blue, P.orange, P.ink)):
        m = (t >= t0) & (t < t1)
        ax[0].scatter(R[m, col['W_L']], M[m], s=6, color=colr, alpha=0.6, lw=0, label=lab.replace('\n', ' '), rasterized=True)
        ax[1].scatter(R[m, col['P_L']], M[m], s=6, color=colr, alpha=0.6, lw=0, rasterized=True)
        ax[2].scatter(np.maximum(R[m, col['SigSFR_10']], 1e-5), M[m], s=6, color=colr, alpha=0.6, lw=0, rasterized=True)
    Wx = np.logspace(2.5, 6, 50); ax[0].plot(Wx, 0.5 * ok22.sfr_of_W(Wx) * 0.25 * 1e7, color=P.ink, lw=1.0, ls='--', label=r'$0.5\,\Sigma_{\rm SFR}(\mathcal{W})\,A\,\tau$')
    ax[1].plot(Wx, 0.5 * ok22.sfr_of_Ptot(Wx) * 0.25 * 1e7, color=P.ink, lw=1.0, ls='--', label=r'$0.5\,\Sigma_{\rm SFR}(P)\,A\,\tau$')
    Sx = np.logspace(-5, -1, 50); ax[2].plot(Sx, 0.5 * Sx * 0.25 * 1e7, color=P.ink, lw=1.0, ls='--', label=r'$0.5\,\Sigma_{\rm SFR}\,A\,\tau$')
    for a in ax: a.set_xscale('log'); a.set_yscale('log'); a.set_ylim(80, 3e7)
    ax[0].set_xlim(3e2, 1e6); ax[1].set_xlim(3e2, 1e6); ax[2].set_xlim(1e-5, 1e-1)
    for a in ax: P.place_legend(a, loc='upper left', fontsize=6)
    ax[0].set_xlabel(r'$\mathcal{W}_{\rm 2p}$ [K cm$^{-3}$]'); ax[1].set_xlabel(r'$P_{\rm tot,2p}$ [K cm$^{-3}$]'); ax[2].set_xlabel(r'$\Sigma_{\rm SFR}$ [M$_\odot$ yr$^{-1}$ kpc$^{-2}$]'); ax[0].set_ylabel(r'$M_{\rm group}$ [M$_\odot$]')
    fig.subplots_adjust(wspace=0.06); P.save(fig, 'cluster_env'); plt.close(fig)

@figure
def prfm_time_los():
    d = np.load(f'{C.DATADIR}/prfm/prfm_time_los.npz', allow_pickle=True); fr = d['frame'].astype(str)   # NOTE: full-column LOS; layer-cut LOS version in scale_los.npz
    suf = lambda f: f.split('_', 2)[2] if f.count('_') >= 2 else 'n'
    T = np.unique(d['t'])
    fig, ax = P.fig()
    for key, col, lab in (('n', P.blue, 'disc normal'), ('t30', P.light, r'tilted $30^\circ$'), ('t60', P.orange, r'tilted $60^\circ$'), ('t90', P.ink, r'tilted $90^\circ$')):
        y = []
        for t in T:
            m = (d['t'] == t) & np.array([suf(f) == key for f in fr]) & (d['n'] > 0)
            y.append(np.sum(d['pw'][m] * d['n'][m]) / np.sum(d['n'][m]) if m.any() else np.nan)
        ok = T > 10; ax.plot(T[ok], np.array(y)[ok], 'o-', ms=2, lw=1.0, color=col, label=lab)
    ax.axhline(1, color=P.ink, lw=0.7); ax.set_yscale('log'); ax.set_xlim(10, 226); ax.set_ylim(0.04, 3)
    ax.set_yticks([0.05, 0.1, 0.2, 0.5, 1, 2]); ax.set_yticklabels(['0.05', '0.1', '0.2', '0.5', '1', '2'])
    ax.set_xlabel(r'$t$ [Myr]'); ax.set_ylabel(r'$P_{\rm tot,2p}\,/\,\mathcal{W}_{\rm 2p}$'); P.phases(ax); P.place_legend(ax, loc='lower right', ncol=2)
    P.save(fig, 'prfm_time_los'); plt.close(fig)

def _local_slope(M, lo, hi):
    M = M[(M >= lo) & (M < hi)]; n = len(M)
    if n < 8: return np.nan, np.nan
    aa = np.linspace(0.2, 4.0, 381)
    def ll(a):
        a = a if abs(a - 1) > 1e-6 else 1 + 1e-6; Z = (lo ** (1 - a) - hi ** (1 - a)) / (a - 1); return -a * np.log(M).sum() - n * np.log(Z)
    v = np.array([ll(a) for a in aa]); al = aa[np.argmax(v)]; h = 0.02; c = (ll(al + h) - 2 * ll(al) + ll(al - h)) / h ** 2
    return al, (1 / np.sqrt(-c) if c < 0 else np.nan)

@figure
def lowmass_convergence():
    d = np.load(f'{C.DATADIR}/clusters/clusters_age0-10_l5_n10.npz'); sel = (np.mod(d['snap'], 10) == 0) & (d['nucleus'] == 0)
    M = d['M'][sel]; t = d['t'][sel]; bd = d['bound'][sel] > 0
    BINS = [(100, 300), (300, 1000), (1000, 3000)]; x = np.array([np.sqrt(a * b) for a, b in BINS])
    fig, ax = P.fig(aspect=0.8)
    for i, ((t0, t1, lab), col) in enumerate(zip(P.PHASES, (P.light, P.blue, P.orange, P.ink))):
        m = (t >= t0) & (t < t1); off = 1 + 0.06 * (i - 1.5)
        r = np.array([_local_slope(M[m], a, b) for a, b in BINS]); rb = np.array([_local_slope(M[m & bd], a, b) for a, b in BINS])
        ax.errorbar(x * off, r[:, 0], yerr=r[:, 1], fmt='o-', ms=3, lw=1.0, color=col, label=lab.replace('\n', ' '))
        ax.plot(x * off, rb[:, 0], 's', ms=3.5, mfc='none', color=col)
    ax.plot([], [], 's', ms=3.5, mfc='none', color=P.grey, label='bound only')
    ax.axhline(2, color=P.grey, lw=0.6, ls='--'); ax.set_xscale('log'); ax.set_xlim(70, 5e3); ax.set_ylim(0, 4.6)
    ax.set_xlabel(r'$M_{\rm group}$ [M$_\odot$]'); ax.set_ylabel(r'local slope $\alpha$'); P.place_legend(ax, loc='upper left', ncol=2)
    P.save(fig, 'lowmass_convergence'); plt.close(fig)

@figure
def duty_cycle():
    D = dict(np.load(f'{C.DATADIR}/prfm/patches_clusters.npz')); WL = np.load(f'{C.DATADIR}/prfm/cluster_link.npz')['WL']
    clean = (D['f_intruder'] < 0.1) | C.is_merged(D['sep'], D['t']); base = clean & (D['Sigma_gas'] > 1) & np.isfinite(WL) & (WL > 0) & (D['t'] > 10)
    burst = base & (D['M_young'] > 500); lW = np.log10(np.where(WL > 0, WL, 1)); duty = np.load(f'{C.DATADIR}/prfm/burst_dist.npz')['duty']
    fig, ax = P.fig(aspect=0.8); xx = np.linspace(3, 6.2, 100); lb = np.arange(3.0, 6.0, 0.5)
    for i, ((t0, t1, lab), col) in enumerate(zip(P.PHASES, (P.light, P.blue, P.orange, P.ink))):
        m = base & (D['t'] >= t0) & (D['t'] < t1); k, w50 = duty[i]
        ax.plot(10 ** xx, 1 / (1 + np.exp(-k * (xx - w50))), color=col, lw=1.2, label=lab.replace('\n', ' '))
        fb = [(burst & m & (lW >= l) & (lW < l + 0.5)).sum() / (m & (lW >= l) & (lW < l + 0.5)).sum() if (m & (lW >= l) & (lW < l + 0.5)).sum() >= 20 else np.nan for l in lb]
        ax.plot(10 ** (lb + 0.25), fb, 'o', ms=3.5, color=col)
    ax.set_xscale('log'); ax.set_xlim(1e3, 1e6); ax.set_ylim(0, 1.45)
    ax.set_xlabel(r'$\mathcal{W}_{\rm 2p}$ [K cm$^{-3}$]'); ax.set_ylabel(r'$P(M_{\rm young} > 500\,{\rm M}_\odot)$'); P.place_legend(ax, loc='upper left', ncol=2)
    P.save(fig, 'duty_cycle'); plt.close(fig)

@figure
def burst_kernel():
    D = dict(np.load(f'{C.DATADIR}/prfm/patches_clusters.npz')); WL = np.load(f'{C.DATADIR}/prfm/cluster_link.npz')['WL']
    clean = (D['f_intruder'] < 0.1) | C.is_merged(D['sep'], D['t']); base = clean & (D['Sigma_gas'] > 1) & np.isfinite(WL) & (WL > 0) & (D['t'] > 10)
    lW = np.log10(np.where(WL > 0, WL, 1)); K = np.load(f'{C.DATADIR}/prfm/burst_kernel.npz'); rng = np.random.default_rng(2)
    t0, t1, lab = C.PHASES[1]; a, b, s_, k, w50 = K[lab.replace(' ', '_')]
    m = base & (D['t'] >= t0) & (D['t'] < t1); obs = D['M_young'][m & (D['M_young'] > 500)]
    edges = np.logspace(np.log10(500), 5.5, 16); c = np.sqrt(edges[1:] * edges[:-1]); w = np.diff(edges)
    fig, ax = P.fig(aspect=0.8)
    for it in range(25):
        pb = 1 / (1 + np.exp(-k * (lW[m] - w50))); hit = rng.random(m.sum()) < pb
        mm = 10 ** (a + b * lW[m][hit] + s_ * rng.standard_normal(hit.sum())); mm = mm[mm > 500]
        h, _ = np.histogram(mm, edges); ok = h > 0; ax.plot(c[ok], h[ok] / w[ok], color=P.light, lw=0.6, alpha=0.7)
    h, _ = np.histogram(obs, edges); ok = h > 0
    ax.errorbar(c[ok], h[ok] / w[ok], yerr=np.sqrt(h[ok]) / w[ok], fmt='o', ms=3.5, color=P.ink, label='simulation, ' + lab)
    ax.plot([], [], color=P.light, lw=1.2, label=r'model: threshold $\times$ lognormal')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(4e2, 4e5); ax.set_ylim(1e-4, 30)
    ax.set_xlabel(r'$M_{\rm young}$ [M$_\odot$]'); ax.set_ylabel(r'${\rm d}N/{\rm d}M$ [M$_\odot^{-1}$]'); P.place_legend(ax, loc='upper right')
    P.save(fig, 'burst_kernel'); plt.close(fig)

def _partC():
    d = np.load(f'{C.DATADIR}/clouds/lineage_env_cl100.npz'); env = d['env']; ec = {n: i for i, n in enumerate(d['names'])}; L = d['lineage_rows']; col = {n: i for i, n in enumerate(d['lineage_names'])}
    P = env[:, ec['Ptot_2p']]; rho = env[:, ec['rho_mid_2p']]; al = L[:, col['alpha_on']]; eps = L[:, col['eps_int']]; ton = L[:, col['t_onset']]; sc = L[:, col['sig_on']]
    ok = (rho > 0) & (P > 0) & (al > 0) & (eps > 0)
    seff = np.sqrt(P[ok] / (rho[ok] * 1e9) / 4.903245584337325e-06); return seff, al[ok], eps[ok], ton[ok], sc[ok]

@figure
def partC_alpha():
    seff, al, eps, ton, sc = _partC()
    fig, ax = P.fig(aspect=0.8); xb = np.logspace(0.6, 1.7, 9); xc = np.sqrt(xb[1:] * xb[:-1])
    for (t0, t1, lab), colr in zip(P.PHASES, (P.light, P.blue, P.orange, P.ink)):
        m = (ton >= t0) & (ton < t1); ax.scatter(seff[m], al[m], s=4, color=colr, alpha=0.3, lw=0, rasterized=True)
        med = np.array([np.median(al[m & (seff >= a) & (seff < b)]) if np.sum(m & (seff >= a) & (seff < b)) >= 8 else np.nan for a, b in zip(xb[:-1], xb[1:])])
        ax.plot(xc, med, 'o-', ms=3, lw=1.3, color=colr, label=lab.replace('\n', ' '))
    xx = np.array([4, 60]); ax.plot(xx, 1.0 * (xx / 9.) ** 2, color=P.grey, lw=0.9, ls='--', label=r'$\alpha_{\rm vir} \propto \sigma_{\rm eff}^2$')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(4, 60); ax.set_ylim(0.1, 3000)
    ax.set_xlabel(r'$\sigma_{\rm eff} = (P_{\rm tot,2p}/\rho_{\rm mid,2p})^{1/2}$ [km s$^{-1}$]'); ax.set_ylabel(r'$\alpha_{\rm vir}$ of the clump before onset'); P.place_legend(ax, loc='upper left', ncol=2)
    P.save(fig, 'partC_alpha'); plt.close(fig)

@figure
def partC_eff():
    """clump efficiency vs total (kinetic + magnetic) virial parameter before onset, with the two-level boundedness model"""
    d = np.load(f'{C.DATADIR}/clouds/lineage_env_cl100.npz'); L = d['lineage_rows']; col = {n: i for i, n in enumerate(d['lineage_names'])}
    al = L[:, col['alpha_on']]; eps = L[:, col['eps_int']]; ton = L[:, col['t_onset']]; sc = L[:, col['sig_on']]; vA = L[:, col['vA_on']]
    ok = (al > 0) & (eps > 0) & (vA > 0); al, eps, ton, sc, vA = [a[ok] for a in (al, eps, ton, sc, vA)]; at = al * (1 + (vA / sc) ** 2)
    pb, pu, lac, m_ = np.load(f'{C.DATADIR}/clouds/step4_fit.npz')['cl100']
    fig, ax = P.fig(aspect=0.8); xb = np.logspace(-0.7, 2.3, 11); xc = np.sqrt(xb[1:] * xb[:-1])
    for (t0, t1, lab), colr in zip(P.PHASES, (P.light, P.blue, P.orange, P.ink)):
        m = (ton >= t0) & (ton < t1); ax.scatter(at[m], eps[m], s=4, color=colr, alpha=0.3, lw=0, rasterized=True, label=lab.replace('\n', ' '))
    med = np.array([np.median(eps[(at >= a) & (at < b)]) if np.sum((at >= a) & (at < b)) >= 10 else np.nan for a, b in zip(xb[:-1], xb[1:])])
    ax.plot(xc, med, 'o-', ms=3.5, lw=1.4, color=P.ink, label='median, all phases')
    xx = np.logspace(-0.7, 2.3, 200); ax.plot(xx, 10 ** (pu + np.log10(1 + (10 ** (pb - pu) - 1) / (1 + (xx / 10 ** lac) ** m_))), color=P.orange, lw=1.3, ls='--', label='two-level model')
    ax.axvline(10 ** lac, color=P.grey, lw=0.6, ls=':')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(0.2, 200); ax.set_ylim(1e-4, 100)
    ax.set_xlabel(r'$\alpha_{\rm vir,tot} = \alpha_{\rm vir}\,(1 + \mathcal{M}_A^{-2})$ before onset'); ax.set_ylabel(r'$\epsilon_{\rm int} = M_\star\,/\,M_{\rm clump}$'); P.place_legend(ax, loc='upper right', ncol=2)
    P.save(fig, 'partC_eff'); plt.close(fig)

@figure
def scale_depth():
    """P/W (Sigma-weighted median over columns, phase median over snapshots) vs column half-depth for the column sizes of the grid"""
    d = np.load(f'{C.DATADIR}/prfm/scale_los.npz'); keys = list(d['keys'])
    cfg = [('', 0.5, 1.5), ('_z05', 0.5, 0.5), ('_p05z025', 0.5, 0.25), ('_p025z025', 0.25, 0.25), ('_p025z0125', 0.25, 0.125), ('_p0125z0125', 0.125, 0.125)]
    fig, ax = P.fig(aspect=0.8)
    for (t0, t1, lab), colr in zip(P.PHASES, (P.light, P.blue, P.orange, P.ink)):
        for side, mk in ((0.5, 'o'), (0.25, 's'), (0.125, '^')):
            xs = []; ys = []
            for tag, pk, zc in cfg:
                if pk != side or tag not in keys: continue
                a = d[tag.replace('|', '__') if tag else '']; T = a[0]; med = a[2]; m = (T >= t0) & (T < t1) & (T > 10)
                if m.sum() < 2: continue
                xs.append(2 * zc); ys.append(np.nanmedian(med[m]))
            if xs: ax.plot(xs, ys, marker=mk, ms=4, lw=1.0, color=colr, label=(lab.replace('\n', ' ') if side == 0.5 else None))
    for mk, side in (('o', '0.5 kpc footprint'), ('s', '0.25 kpc'), ('^', '0.125 kpc')): ax.plot([], [], marker=mk, ms=4, color=P.grey, lw=0, label=side)
    ax.axhline(1, color=P.ink, lw=0.7); ax.axvspan(0.3, 0.5, color=P.grey, alpha=0.12, lw=0); ax.text(0.39, 0.36, r'$2H$', ha='center', fontsize=7, color=P.grey)
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(0.2, 4); ax.set_ylim(0.3, 15)
    ax.set_xticks([0.25, 0.5, 1, 3]); ax.set_xticklabels(['0.25', '0.5', '1', '3']); ax.set_yticks([0.5, 1, 2, 5, 10]); ax.set_yticklabels(['0.5', '1', '2', '5', '10'])
    ax.set_xlabel('column depth [kpc]'); ax.set_ylabel(r'$P_{\rm tot,2p}\,/\,\mathcal{W}_{\rm 2p}$'); P.place_legend(ax, loc='upper right', ncol=2)
    P.save(fig, 'scale_depth'); plt.close(fig)

@figure
def los_layer():
    """P/W (Sigma-weighted median, layer-cut columns) vs time for the disc normal and tilted normals, from scale_los.npz"""
    d = np.load(f'{C.DATADIR}/prfm/scale_los.npz'); keys = list(d['keys'])
    fig, ax = P.fig()
    for suf, col, lab in (('', P.blue, 'disc normal'), ('_t30', P.light, r'tilted $30^\circ$'), ('_t60', P.orange, r'tilted $60^\circ$'), ('_t90', P.ink, r'tilted $90^\circ$')):
        k = '_los_z05' + (('|' + suf) if suf else '')
        if k not in keys: continue
        a = d[k.replace('|', '__')]; T = a[0]; med = a[2]; o = T > 10
        ax.plot(T[o], med[o], 'o-', ms=2, lw=1.0, color=col, label=lab)
    ax.axhline(1, color=P.ink, lw=0.7); ax.set_yscale('log'); ax.set_xlim(10, 226); ax.set_ylim(0.08, 8)
    ax.set_yticks([0.1, 0.2, 0.5, 1, 2, 5]); ax.set_yticklabels(['0.1', '0.2', '0.5', '1', '2', '5'])
    ax.set_xlabel(r'$t$ [Myr]'); ax.set_ylabel(r'$P_{\rm tot,2p}\,/\,\mathcal{W}_{\rm 2p}$'); P.phases(ax); P.place_legend(ax, loc='lower right', ncol=2)
    P.save(fig, 'los_layer'); plt.close(fig)

@figure
def pressure_shares():
    """fractions of the layer P_tot supplied by thermal, turbulent and Maxwell-stress terms, summed over clean layer columns per snapshot;
    dashed: the magnetic pressure B^2/8pi over P_tot"""
    import glob as _g
    T = []; F = []
    for fn in sorted(_g.glob(f'{C.PRFM_DIR}/patch_[0-9][0-9][0-9]_z05.h5')):
        with h5py.File(fn, 'r') as f:
            t = float(f.attrs['time_myr']); acc = np.zeros(5)
            for g in f:
                if not g.startswith('frame_'): continue
                G = f[g]; sep = float(G.attrs['separation_kpc']); r = lambda k: G[k][:].ravel()
                Sig = r('Sigma_gas'); W = r('W_2p'); fi = r('f_intruder'); ok = (Sig > 1) & (W > 0) & ((fi < 0.1) | C.is_merged(sep, t))
                acc += [r('Pth_2p')[ok].sum(), r('Pturb_2p')[ok].sum(), r('Pmag_2p')[ok].sum(), r('Pmag_tot')[ok].sum(), r('Ptot_2p')[ok].sum()]
        if acc[4] > 0: T.append(t); F.append(acc / acc[4])
    T = np.array(T); F = np.array(F); o = T > 10
    fig, ax = P.fig()
    ax.plot(T[o], F[o, 0], color=P.light, lw=1.3, label='thermal')
    ax.plot(T[o], F[o, 1], color=P.blue, lw=1.5, label='turbulent')
    ax.plot(T[o], F[o, 2], color=P.orange, lw=1.5, label='Maxwell stress')
    ax.plot(T[o], F[o, 3], color=P.ink, lw=1.0, ls='--', label=r'$B^2/8\pi$')
    ax.set_xlim(10, 226); ax.set_ylim(-0.1, 2.3); ax.axhline(0, color=P.grey, lw=0.5)
    ax.set_xlabel(r'$t$ [Myr]'); ax.set_ylabel(r'fraction of $P_{\rm tot,2p}$'); P.phases(ax, box=True); P.place_legend(ax, loc='upper left', ncol=2)
    P.save(fig, 'pressure_shares'); plt.close(fig)

def _frame_series(key, stat='wmed', layer=True):
    """per snapshot and per frame (A, B while two frames exist, M after): Sigma_gas-weighted median (or clean-area mean) of key over the
    clean columns.  returns dict frame -> (t, value)"""
    import glob as _g
    out = {'A': [], 'B': [], 'M': []}
    pat = 'patch_[0-9][0-9][0-9]_z05.h5' if layer else 'patch_[0-9][0-9][0-9].h5'
    for fn in sorted(_g.glob(f'{C.PRFM_DIR}/{pat}')):
        with h5py.File(fn, 'r') as f:
            t = float(f.attrs['time_myr'])
            for g in f:
                if not g.startswith('frame_'): continue
                G = f[g]; sep = float(G.attrs['separation_kpc']); r = lambda k: G[k][:].ravel()
                Sig = r('Sigma_gas'); W = r('W_2p'); fi = r('f_intruder'); ok = (Sig > 1) & (W > 0) & ((fi < 0.1) | C.is_merged(sep, t))
                if ok.sum() < 3: continue
                x = r(key)[ok]
                if stat == 'wmed':
                    o = np.argsort(x); cw = np.cumsum(Sig[ok][o]) / Sig[ok].sum(); v = x[o][np.searchsorted(cw, 0.5)]
                else: v = x.mean()
                out[g[-1]].append((t, v))
    return {k: np.array(v) for k, v in out.items() if len(v)}

def _plot_frames(ax, D, **kw):
    lab = {'A': 'galaxy A', 'B': 'galaxy B', 'M': 'merged'}; col = {'A': P.blue, 'B': P.orange, 'M': P.ink}
    for k in ('A', 'B', 'M'):
        if k in D:
            tt, vv = D[k][:, 0], D[k][:, 1].astype(float); vv[1:][np.diff(tt) > 2.5] = np.nan   # break the line across gaps (frame absent)
            o = tt > 10; ax.plot(tt[o], vv[o], color=col[k], lw=1.3, label=lab[k], **kw)

@figure
def midplane_pressure():
    """layer P_tot per snapshot: Sigma_gas-weighted median over the clean columns of each frame"""
    D = _frame_series('Ptot_2p')
    fig, ax = P.fig(); _plot_frames(ax, D)
    ax.set_yscale('log'); ax.set_xlim(10, 226); ax.set_ylim(3e2, 3e6)
    ax.set_xlabel(r'$t$ [Myr]'); ax.set_ylabel(r'$P_{\rm tot,2p}/k_{\rm B}$ [K cm$^{-3}$]'); P.phases(ax, box=True); P.place_legend(ax, loc='lower right')
    P.save(fig, 'midplane_pressure'); plt.close(fig)

@figure
def midplane_density():
    """layer rho_mid,2p per snapshot: Sigma_gas-weighted median over the clean columns of each frame, as n_H"""
    D = _frame_series('rho_mid_2p')
    f = 1.989e33 / 3.0857e18**3 * C.XH / C.PROTONMASS   # Msun/pc^3 -> H nuclei per cm^3 (= 30.8)
    D = {k: np.c_[v[:, 0], v[:, 1] * f] for k, v in D.items()}
    fig, ax = P.fig(); _plot_frames(ax, D)
    ax.set_yscale('log'); ax.set_xlim(10, 226); ax.set_ylim(3e-3, 3e1)
    ax.set_xlabel(r'$t$ [Myr]'); ax.set_ylabel(r'$n_{\rm H,mid,2p}$ [cm$^{-3}$]'); P.phases(ax, box=True); P.place_legend(ax, loc='lower right')
    P.save(fig, 'midplane_density'); plt.close(fig)

@figure
def sfr_surface_density():
    """Sigma_SFR,10 per snapshot: mean over the clean full columns of each frame (= SFR in the clean columns / their area)"""
    D = _frame_series('SigSFR_10', stat='mean', layer=False)
    fig, ax = P.fig(); _plot_frames(ax, D)
    ax.set_yscale('log'); ax.set_xlim(10, 226); ax.set_ylim(1e-6, 1e0)
    ax.set_xlabel(r'$t$ [Myr]'); ax.set_ylabel(r'$\Sigma_{\rm SFR}$ [M$_\odot$ yr$^{-1}$ kpc$^{-2}$]'); P.phases(ax, box=True); P.place_legend(ax, loc='upper left')
    P.save(fig, 'sfr_surface_density'); plt.close(fig)

@figure
def sfr_history():
    """total SFR of the run from the formation times of all stars present in the last snapshot: mass formed in (t - dt, t] / dt"""
    s4 = C.read(231, 4, ['StellarFormationTime', 'Masses']); tf = s4['StellarFormationTime'] * C.MYR; m = s4['Masses'] * 1e10
    t = np.arange(1, 226.5, 1.0); fig, ax = P.fig()
    for dt, c, lw in ((10, P.light, 1.0), (20, P.blue, 1.2), (30, P.orange, 1.2), (40, P.ink, 1.4)):
        sfr = np.array([m[(tf > ti - dt) & (tf <= ti)].sum() / (dt * 1e6) for ti in t]); o = t >= dt
        ax.plot(t[o], sfr[o], color=c, lw=lw, label=f'{dt} Myr')
    ax.set_yscale('log'); ax.set_xlim(10, 226); ax.set_ylim(3e-4, 3e0)
    ax.set_xlabel(r'$t$ [Myr]'); ax.set_ylabel(r'SFR [M$_\odot$ yr$^{-1}$]'); P.phases(ax, box=True); P.place_legend(ax, loc='upper left', ncol=2)
    P.save(fig, 'sfr_history'); plt.close(fig)

@figure
def link4_cascade():
    """link 4 test: clump 3d velocity dispersion at the pre-onset snapshot against the layer-cascade prediction
    sigma_eff (r_h/H)^(1/2) of the column the clump sits in; p fitted as the slope of log(sigma_3d/sigma_eff) on log(r_h/H)"""
    d = np.load(f'{C.DATADIR}/clouds/lineage_env_cl100.npz'); env = d['env']; ec = {n: i for i, n in enumerate(d['names'])}; L = d['lineage_rows']; col = {n: i for i, n in enumerate(d['lineage_names'])}
    Pk = env[:, ec['Ptot_2p']]; rho = env[:, ec['rho_mid_2p']]; H = env[:, ec['H']]; rh = L[:, col['rh_on']] * 1e-3; s3 = L[:, col['sig_on']]; ton = L[:, col['t_onset']]
    ok = (rho > 0) & (Pk > 0) & (H > 0) & (rh > 0) & (s3 > 0)
    seff = np.sqrt(Pk / (rho * 1e9) / 4.903245584337325e-06)
    x = np.log10(rh / H)[ok]; y = np.log10(s3 / seff)[ok]
    p, a = np.polyfit(x, y, 1); res = np.std(y - (p * x + a))
    print(f'link 4: {ok.sum()} clumps; log(sigma_3d/sigma_eff) = {a:+.2f} + {p:.2f} log(r_h/H), scatter {res:.2f} dex; median r_h/H = {10**np.median(x):.3f}')
    for t0, t1, lab in C.PHASES:
        m = (ton[ok] >= t0) & (ton[ok] < t1)
        if m.sum() > 20: pp, aa = np.polyfit(x[m], y[m], 1); print(f'   {lab:28s} n={m.sum():4d}  p={pp:.2f}  offset={aa:+.2f}  median sigma_3d/[sigma_eff (r_h/H)^0.5] = {10**np.median(y[m] - 0.5 * x[m]):.2f}')
    fig, ax = P.fig()
    xs = seff[ok] * np.sqrt(rh[ok] / H[ok]); ys = s3[ok]
    for (t0, t1, lab), c in zip(C.PHASES, (P.light, P.blue, P.orange, P.ink)):
        m = (ton[ok] >= t0) & (ton[ok] < t1); ax.scatter(xs[m], ys[m], s=4, color=c, alpha=0.6, lw=0, label=lab.replace('\n', ' '))
    g = np.logspace(-0.7, 1.5, 10); ax.plot(g, g, color=P.grey, lw=0.8, ls='--', label=r'$\sigma_{\rm 3d} = \sigma_{\rm eff}\,(r_h/H)^{1/2}$')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(0.2, 60); ax.set_ylim(0.2, 60)
    ax.set_xlabel(r'$\sigma_{\rm eff}\,(r_h/H)^{1/2}$ [km s$^{-1}$]'); ax.set_ylabel(r'$\sigma_{\rm 3d}$ of the clump [km s$^{-1}$]'); P.place_legend(ax, loc='upper left', fontsize=6.5)
    P.save(fig, 'link4_cascade'); plt.close(fig)

@figure
def bound_fraction():
    """fraction of the stellar mass formed in an age window that sits in bound groups above a mass floor, per snapshot (every 10th):
    windows 0-10 and 20-50 Myr, floors 300 and 1000 Msun.  A property of the star formation criterion, shown for completeness."""
    fig, ax = P.fig()
    for (fn, lo, hi, lab), c in ((('clusters_age0-10_l5_n25.npz', 0, 10, 'age 0 to 10 Myr'), P.blue), (('clusters_age20-50_l5_n25.npz', 20, 50, 'age 20 to 50 Myr'), P.orange)):
        cl = np.load(f'{C.DATADIR}/clusters/{fn}')
        for floor, ls in ((300, '-'), (1000, '--')):
            T = []; G = []
            for k in range(0, C.NSNAP, 10):
                st = np.load(f'{C.DATADIR}/stars/stars_{k:03d}.npz'); t = float(st['time_myr']); age = t - st['tform_myr']
                Mf = np.sum(st['mass'][(age >= lo) & (age < hi) & (st['tform_myr'] > 0)]) * C.MSUN
                if Mf < 2000 or t < lo + 5: continue
                m = (cl['snap'] == k) & (cl['bound'] > 0) & (cl['M'] >= floor); T.append(t); G.append(cl['M'][m].sum() / Mf)
            ax.plot(T, G, color=c, ls=ls, lw=1.3, label=f'{lab}, groups $\\geq$ {floor} M$_\\odot$')
    ax.set_xlim(10, 226); ax.set_ylim(0, 1.6); ax.set_xlabel(r'$t$ [Myr]'); ax.set_ylabel('bound fraction of the stars formed')
    P.phases(ax, box=True); P.place_legend(ax, loc='upper left', fontsize=6.5)
    P.save(fig, 'bound_fraction'); plt.close(fig)

@figure
def gamma_env():
    """bound fraction of the stars formed in the last 10 Myr (groups >= 300 Msun filled, >= 1000 open) against the state of the clean
    columns of the same snapshot: Sigma_SFR over the active columns, layer P_tot, sigma_eff, layer W (Sigma-weighted medians).
    Goddard+2010 observed relation on the Sigma_SFR panel."""
    d = np.load(f'{C.DATADIR}/clusters/gamma_env.npz'); R = d['rows']; c = {n: i for i, n in enumerate(d['names'])}
    o = (R[:, c['t']] > 25) & (R[:, c['M_formed']] > 2000)
    fig, axs = plt.subplots(2, 2, figsize=(P.COL2 * 0.72, P.COL2 * 0.62)); axs = axs.ravel()
    panels = [('SigSFR_active', r'$\Sigma_{\rm SFR}$ of the active columns [M$_\odot$ yr$^{-1}$ kpc$^{-2}$]', (2e-5, 1e-1)),
              ('P_L', r'$P_{\rm tot,2p}$ [K cm$^{-3}$]', (1e3, 1e6)), ('sigma_eff', r'$\sigma_{\rm eff}$ [km s$^{-1}$]', (4, 40)), ('W_L', r'$\mathcal{W}_{\rm 2p}$ [K cm$^{-3}$]', (1e3, 1e6))]
    cols = (P.light, P.blue, P.orange, P.ink)
    for ax, (key, lab, xl) in zip(axs, panels):
        x = R[:, c[key]]
        for (t0, t1, plab), col in zip(C.PHASES, cols):
            m = o & (R[:, c['t']] >= t0) & (R[:, c['t']] < t1)
            short = {'before 1st\npassage': 'before 1st', '1st to 2nd\npassage': '1st to 2nd', '2nd to 3rd\npassage': '2nd to 3rd', '3rd passage to\ncoalescence': 'after 3rd'}.get(plab, plab.replace('\n', ' '))
            ax.scatter(x[m], R[m, c['G300']], s=22, color=col, zorder=3, label=short)
            ax.scatter(x[m], R[m, c['G1000']], s=22, facecolors='none', edgecolors=col, lw=0.8, zorder=3)
        if key == 'SigSFR_active':
            g = np.logspace(-5, 0, 50); ax.plot(g, 0.29 * g**0.24, color=P.grey, lw=1.0, ls='--', label='Goddard 2010')
            ax.fill_between(g, 0.23 * g**0.28, 0.35 * g**0.20, color=P.grey, alpha=0.15, lw=0)
        ax.set_xscale('log'); ax.set_xlim(*xl); ax.set_ylim(0, 1.45); ax.set_xlabel(lab)
        if key == 'sigma_eff': ax.set_xticks([5, 10, 20, 30]); ax.set_xticklabels(['5', '10', '20', '30']); ax.xaxis.set_minor_formatter(plt.NullFormatter())
    axs[0].set_ylabel('bound fraction, age < 10 Myr'); axs[2].set_ylabel('bound fraction, age < 10 Myr')
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc='upper center', bbox_to_anchor=(0.5, 1.004), ncol=5, fontsize=6.5, frameon=False, handletextpad=0.3, columnspacing=1.0)
    fig.tight_layout(w_pad=0.6, h_pad=0.6, rect=(0, 0, 1, 0.955)); P.save(fig, 'gamma_env'); plt.close(fig)

@figure
def sfr_lag():
    """cross-correlation of log P_layer(t) and log W_layer(t) (Sigma-weighted medians over the clean layer columns) with log SFR(t + lag)
    (whole-run SFR in 2 Myr bins from the stellar formation times), over t > 25 Myr and per phase; positive lag = star formation later"""
    d = np.load(f'{C.DATADIR}/prfm/lag_series.npz'); T = d['T']; lp = np.log10(d['Pm']); lw = np.log10(d['Wm']); ls = np.log10(np.maximum(d['sfr2'], 1e-5))
    lags = np.arange(-30, 41)
    def xc(a, b, sel, lag):
        aa = a[sel]; bb = b[sel]
        if lag > 0: x = aa[:-lag]; y = bb[lag:]
        elif lag < 0: x = aa[-lag:]; y = bb[:lag]
        else: x = aa; y = bb
        if len(x) < 10: return np.nan
        x = x - x.mean(); y = y - y.mean(); return (x * y).sum() / np.sqrt((x * x).sum() * (y * y).sum())
    fig, ax = P.fig()
    ax.plot(lags, [xc(lp, ls, T > 25, l) for l in lags], color=P.ink, lw=1.6, label='$P_{\\rm tot,2p}$, whole run')
    ax.plot(lags, [xc(lw, ls, T > 25, l) for l in lags], color=P.ink, lw=1.0, ls='--', label='$\\mathcal{W}_{\\rm 2p}$, whole run')
    for (t0, t1, lab), c in zip(C.PHASES[1:], (P.blue, P.orange, P.grey)):
        ax.plot(lags, [xc(lp, ls, (T >= t0) & (T < t1), l) for l in lags], color=c, lw=1.2, label=r'$P_{\rm tot,2p}$, ' + lab.replace('\n', ' '))
    ax.axvline(0, color=P.grey, lw=0.5); ax.axhline(0, color=P.grey, lw=0.5)
    ax.set_xlim(-30, 40); ax.set_ylim(-0.8, 1.45); ax.set_xlabel('lag of the star formation rate behind the pressure [Myr]'); ax.set_ylabel('correlation coefficient')
    P.place_legend(ax, loc='upper left', fontsize=6.5)
    P.save(fig, 'sfr_lag'); plt.close(fig)

@figure
def chain_sfr():
    """star formation rate of the run (2 Myr bins) against the rate the chain gives from the clump population of the same snapshot:
    sum over the n_H > 100 clumps of eps(alpha_tot) M / tau, with the fitted efficiency step and tau the run-median of the ratio"""
    d = np.load(f'{C.DATADIR}/clouds/chain_sfr.npz'); R = d['rows']; c = {n: i for i, n in enumerate(d['names'])}
    o = R[:, c['t']] > 25; tau = np.median(R[o, c['Mstar_pred']] / np.maximum(R[o, c['sfr_now']], 1e-30))
    fig, ax = P.fig()
    ax.plot(R[:, c['t']], R[:, c['sfr_now']], color=P.ink, lw=1.4, label='measured, 2 Myr bins')
    ax.plot(R[:, c['t']], R[:, c['Mstar_pred']] / tau, color=P.orange, lw=1.4, label=r'chain: $\sum_{\rm clumps}\epsilon(\alpha_{\rm tot})\,M\,/\,\tau$')
    ax.set_yscale('log'); ax.set_xlim(10, 226); ax.set_ylim(1e-4, 3e0)
    ax.set_xlabel(r'$t$ [Myr]'); ax.set_ylabel(r'SFR [M$_\odot$ yr$^{-1}$]'); P.phases(ax, box=True); P.place_legend(ax, loc='upper left')
    print(f'chain_sfr: tau = {tau/1e6:.1f} Myr')
    P.save(fig, 'chain_sfr'); plt.close(fig)

@figure
def yield_lag():
    """pressure per unit star formation with the star formation rate taken LAG Myr later: sum_c P_c(t) / sum_c Sigma_SFR,10,c(t + lag), clean
    layer columns, against the OK22 yield at Pbar(t); lag = 12 Myr from the cross-correlation"""
    d = np.load(f'{C.DATADIR}/prfm/lag_series.npz'); T = d['T']; Psum = d['Psum']; SF = d['SFsum']; Pbar = d['Pbar']
    fig, ax = P.fig(); o = T > 10
    ups0 = Psum / np.maximum(SF, 1e-30) / 4.81e3
    ax.plot(T[o], ups0[o], color=P.light, lw=1.2, label='no lag')
    for lag, c in ((12, P.blue),):
        n = len(T) - lag; ups = Psum[:n] / np.maximum(SF[lag:], 1e-30) / 4.81e3; oo = T[:n] > 10
        ax.plot(T[:n][oo], ups[oo], color=c, lw=1.5, label=f'star formation {lag} Myr later')
    ax.plot(T[o], ok22.ups_tot(Pbar[o]), color=P.ink, lw=1.0, ls='--', label=r'feedback yield $\Upsilon_{\rm tot}(P)$ (OK22)')
    ax.set_yscale('log'); ax.set_xlim(10, 226); ax.set_ylim(20, 1e6)
    ax.set_xlabel(r'$t$ [Myr]'); ax.set_ylabel(r'$P_{\rm tot,2p}(t)\,/\,\Sigma_{\rm SFR}(t+\ell)$ [km s$^{-1}$]'); P.phases(ax, box=True); P.place_legend(ax, loc='lower left')
    for lag in (0, 12):
        n = len(T) - lag; r = Psum[:n] / np.maximum(SF[lag:], 1e-30) / 4.81e3 / ok22.ups_tot(Pbar[:n]); tt = T[:n]
        print(f'yield_lag: lag {lag:2d} Myr, phase medians of the ratio to OK22: ' + ' '.join(f'{np.median(r[(tt>=a)&(tt<b)]):6.1f}' for a, b, _ in C.PHASES) + f'; time-median (t>25) {np.median(r[tt>25]):.2f}; fraction >3: {np.mean(r[tt>25]>3):.2f}')
    P.save(fig, 'yield_lag'); plt.close(fig)

@figure
def yield_sources():
    """top: ratio of the pressure per unit star formation to the OK22 yield with the 10 Myr and the 40 Myr rate; bottom: supernova rate of
    the clean layer columns per unit 10 Myr star formation rate, normalised to its median between 40 and 100 Myr"""
    import glob as _g
    rows = []
    for fn in sorted(_g.glob(f'{C.PRFM_DIR}/patch_[0-9][0-9][0-9]_z05.h5')):
        with h5py.File(fn, 'r') as f:
            t = float(f.attrs['time_myr']); acc = np.zeros(6)
            for g in f:
                if not g.startswith('frame_'): continue
                G = f[g]; sep = float(G.attrs['separation_kpc']); r = lambda k: G[k][:].ravel()
                Sig = r('Sigma_gas'); W = r('W_2p'); fi = r('f_intruder'); ok = (Sig > 1) & (W > 0) & ((fi < 0.1) | C.is_merged(sep, t))
                acc += [r('Ptot_2p')[ok].sum(), r('SigSFR_10')[ok].sum(), r('SigSFR_40')[ok].sum(), r('SNrate')[ok].sum(), (r('Sigma_gas_2p') * r('Ptot_2p'))[ok].sum(), r('Sigma_gas_2p')[ok].sum()]
            rows.append((t, *acc))
    R = np.array(rows); T = R[:, 0]; Pb = R[:, 5] / R[:, 6]
    u10 = R[:, 1] / np.maximum(R[:, 2], 1e-30) / 4.81e3 / ok22.ups_tot(Pb); u40 = R[:, 1] / np.maximum(R[:, 3], 1e-30) / 4.81e3 / ok22.ups_tot(Pb)
    sn = R[:, 4] / np.maximum(R[:, 2], 1e-30); sn /= np.median(sn[(T > 40) & (T < 100)])
    sn40 = R[:, 4] / np.maximum(R[:, 3], 1e-30); sn40 /= np.median(sn40[(T > 40) & (T < 100)]); o = T > 10
    fig, ax = plt.subplots(2, 1, figsize=(P.COL1, P.COL1 * 1.05), sharex=True, gridspec_kw=dict(hspace=0.06, height_ratios=[1.3, 1]))
    ax[0].plot(T[o], u10[o], color=P.light, lw=1.3, label='10 Myr rate'); ax[0].plot(T[o], u40[o], color=P.blue, lw=1.5, label='40 Myr rate (OK22)')
    ax[0].axhline(1, color=P.ink, lw=0.7); ax[0].set_yscale('log'); ax[0].set_ylim(0.1, 3e3); ax[0].set_ylabel(r'$(P_{\rm tot,2p}/\Sigma_{\rm SFR})\,/\,\Upsilon_{\rm tot}^{\rm OK22}$')
    P.phases(ax[0], box=True); P.place_legend(ax[0], loc='upper left')
    ax[1].plot(T[o], sn[o], color=P.light, lw=1.3, label='per 10 Myr rate'); ax[1].plot(T[o], sn40[o], color=P.blue, lw=1.5, label='per 40 Myr rate')
    ax[1].axhline(1, color=P.ink, lw=0.7)
    ax[1].set_yscale('log'); ax[1].set_ylim(0.02, 2e3); ax[1].set_ylabel('SN rate per unit SFR (disc = 1)'); ax[1].set_xlabel(r'$t$ [Myr]'); ax[1].set_xlim(10, 226)
    P.place_legend(ax[1], loc='upper left', fontsize=6.5, ncol=2)
    for tp in P.PERICENTRES: ax[1].axvline(tp, color=P.grey, lw=0.6, ls=':' if tp < 200 else '-')
    P.save(fig, 'yield_sources'); plt.close(fig)

@figure
def turb_per_sn():
    """top: turbulent pressure of the two-phase gas per unit supernova rate, summed over the clean layer columns, relative to its median
    between 40 and 100 Myr; bottom: Sigma-weighted median hot volume fraction (T > 2e4 K) of the midplane slab"""
    import glob as _g
    rows = []
    for fn in sorted(_g.glob(f'{C.PRFM_DIR}/patch_[0-9][0-9][0-9]_z05.h5')):
        with h5py.File(fn, 'r') as f:
            t = float(f.attrs['time_myr']); pt = 0.0; sn = 0.0; fv = []; S = []
            for g in f:
                if not g.startswith('frame_'): continue
                G = f[g]; sep = float(G.attrs['separation_kpc']); r = lambda k: G[k][:].ravel()
                Sig = r('Sigma_gas'); W = r('W_2p'); fi = r('f_intruder'); ok = (Sig > 1) & (W > 0) & ((fi < 0.1) | C.is_merged(sep, t))
                pt += r('Pturb_2p')[ok].sum(); sn += r('SNrate')[ok].sum(); fv.append(r('fvol_hot')[ok]); S.append(Sig[ok])
            fv = np.concatenate(fv); S = np.concatenate(S); o = np.argsort(fv); cw = np.cumsum(S[o]) / S.sum()
            rows.append((t, pt, sn, fv[o][np.searchsorted(cw, 0.5)]))
    R = np.array(rows); T = R[:, 0]; ratio = R[:, 1] / np.maximum(R[:, 2], 1e-30); ratio /= np.median(ratio[(T > 40) & (T < 100)]); o = T > 10
    fig, ax = plt.subplots(2, 1, figsize=(P.COL1, P.COL1 * 1.0), sharex=True, gridspec_kw=dict(hspace=0.06, height_ratios=[1.3, 1]))
    ax[0].plot(T[o], ratio[o], color=P.blue, lw=1.5); ax[0].axhline(1, color=P.ink, lw=0.7)
    ax[0].set_yscale('log'); ax[0].set_ylim(0.05, 500); ax[0].set_ylabel(r'$P_{\rm turb,2p}\,/\,$SN rate (disc = 1)'); P.phases(ax[0], box=True)
    ax[1].plot(T[o], R[o, 3], color=P.orange, lw=1.4); ax[1].set_ylim(0, 0.2); ax[1].set_ylabel('hot volume fraction'); ax[1].set_xlabel(r'$t$ [Myr]'); ax[1].set_xlim(10, 226)
    for tp in P.PERICENTRES: ax[1].axvline(tp, color=P.grey, lw=0.6, ls=':' if tp < 200 else '-')
    P.save(fig, 'turb_per_sn'); plt.close(fig)

@figure
def pressure_budget():
    """what holds the layer up against what each source can supply, per snapshot, Sigma_gas-weighted medians over the clean layer
    columns: measured P_tot and W; the pressure PRFM feedback supports at the measured 40 Myr star formation rate, i.e. the P for which
    OK22 eq. 28a returns that rate; and the shear term 0.02 Sigma_gas H S^2 with S the mass-weighted shear rate of the slab"""
    import glob as _g
    P_UNIT = 4.903245584337325e-06; rows = []
    for fn in sorted(_g.glob(f'{C.PRFM_DIR}/patch_[0-9][0-9][0-9]_z05.h5')):
        with h5py.File(fn, 'r') as f:
            t = float(f.attrs['time_myr']); keys = ['Sigma_gas', 'W_2p', 'f_intruder', 'Ptot_2p', 'SigSFR_40', 'H', 'shear']
            acc = {k: [] for k in keys}
            for g in f:
                if not g.startswith('frame_'): continue
                G = f[g]; sep = float(G.attrs['separation_kpc']); r = lambda k: G[k][:].ravel()
                Sig = r('Sigma_gas'); W = r('W_2p'); fi = r('f_intruder'); ok = (Sig > 1) & (W > 0) & ((fi < 0.1) | C.is_merged(sep, t))
                for k in keys: acc[k].append(r(k)[ok])
            for k in keys: acc[k] = np.concatenate(acc[k])
            S = acc['Sigma_gas']
            def wm(x):
                o = np.argsort(x); cw = np.cumsum(S[o]) / S.sum(); return x[o][np.searchsorted(cw, 0.5)]
            pfb = 10 ** ((np.log10(np.maximum(acc['SigSFR_40'], 1e-9)) + 7.43) / 1.18)      # OK22 eq. 28a inverted
            psh = 0.02 * (acc['Sigma_gas'] * 1e6) * acc['H'] * acc['shear'] ** 2 * P_UNIT
            rows.append((t, wm(acc['Ptot_2p']), wm(acc['W_2p']), wm(pfb), wm(psh)))
    R = np.array(rows); T = R[:, 0]; o = T > 25
    fig, ax = P.fig()
    ax.plot(T[o], R[o, 1], color=P.ink, lw=1.7, label=r'$P_{\rm tot,2p}$ measured')
    ax.plot(T[o], R[o, 2], color=P.grey, lw=1.0, ls=':', label=r'$\mathcal{W}_{\rm 2p}$')
    ax.plot(T[o], R[o, 3], color=P.blue, lw=1.4, ls='--', label=r'feedback at the measured SFR')
    ax.plot(T[o], R[o, 4], color=P.orange, lw=1.4, label=r'$0.02\,\Sigma_{\rm gas} H S^2$, flow shear')
    ax.set_yscale('log'); ax.set_xlim(25, 226); ax.set_ylim(3e1, 1e9)
    ax.set_xlabel(r'$t$ [Myr]'); ax.set_ylabel(r'$P/k_{\rm B}$ [K cm$^{-3}$]'); P.phases(ax, box=True, y=0.99)
    P.place_legend(ax, loc='upper left', fontsize=6.3, ncol=2, columnspacing=0.9, handletextpad=0.5)
    for a, b, lab in ((40, 100, 'discs'), (167, 195, 'failure')):
        m = (T >= a) & (T < b); print(f'  {lab:8s} P_tot {np.median(R[m,1]):9.3g}  W {np.median(R[m,2]):9.3g}  feedback-supported {np.median(R[m,3]):9.3g} ({np.median(R[m,3]/R[m,1]):.3f} of P)  shear {np.median(R[m,4]):9.3g} ({np.median(R[m,4]/R[m,1]):.1f} of P)')
    P.save(fig, 'pressure_budget'); plt.close(fig)

@figure
def yield_components():
    """the yield split by component: Ups_X = sum_c P_X,c / sum_c Sigma_SFR,40,c / 4.81e3 for X = thermal, turbulent, Maxwell stress,
    over the clean layer columns, against the OK22 thermal and turbulent yields at the same mean pressure (eqs 26a, 26b)"""
    import glob as _g
    rows = []
    for fn in sorted(_g.glob(f'{C.PRFM_DIR}/patch_[0-9][0-9][0-9]_z05.h5')):
        with h5py.File(fn, 'r') as f:
            t = float(f.attrs['time_myr']); acc = np.zeros(7)
            for g in f:
                if not g.startswith('frame_'): continue
                G = f[g]; sep = float(G.attrs['separation_kpc']); r = lambda k: G[k][:].ravel()
                S = r('Sigma_gas'); W = r('W_2p'); fi = r('f_intruder'); ok = (S > 1) & (W > 0) & ((fi < 0.1) | C.is_merged(sep, t))
                acc += [r('Pth_2p')[ok].sum(), r('Pturb_2p')[ok].sum(), r('Pmag_2p')[ok].sum(), r('Ptot_2p')[ok].sum(), r('SigSFR_40')[ok].sum(), (r('Sigma_gas_2p') * r('Ptot_2p'))[ok].sum(), r('Sigma_gas_2p')[ok].sum()]
            rows.append((t, *acc))
    R = np.array(rows); T = R[:, 0]; Pb = R[:, 6] / R[:, 7]; SF = np.maximum(R[:, 5], 1e-30); o = T > 40
    fig, ax = P.fig()
    ax.plot(T[o], (R[:, 2] / SF / 4.81e3)[o], color=P.blue, lw=1.6, label='turbulent')
    ax.plot(T[o], (R[:, 3] / SF / 4.81e3)[o], color=P.orange, lw=1.4, label='Maxwell stress')
    ax.plot(T[o], (R[:, 1] / SF / 4.81e3)[o], color=P.light, lw=1.4, label='thermal')
    ax.plot(T[o], ok22.ups_turb(Pb[o]), color=P.blue, lw=1.0, ls='--', label='OK22 turbulent')
    ax.plot(T[o], ok22.ups_th(Pb[o]), color=P.grey, lw=1.0, ls=':', label='OK22 thermal')
    ax.set_yscale('log'); ax.set_xlim(40, 226); ax.set_ylim(20, 1e6)
    ax.set_xlabel(r'$t$ [Myr]'); ax.set_ylabel(r'$P_X\,/\,\Sigma_{\rm SFR,40}$ [km s$^{-1}$]'); P.phases(ax, labels=False)
    P.place_legend(ax, loc='upper left', fontsize=6.3, ncol=2, columnspacing=0.9, handletextpad=0.5)
    for a, b, lab in ((40, 100, 'discs'), (167, 195, 'failure')):
        m = (T >= a) & (T < b); ot = np.median(ok22.ups_tot(Pb[m]))
        print(f'  {lab:8s} in units of the OK22 total yield: thermal {np.median(R[m,1]/SF[m]/4.81e3)/ot:6.1f}  turbulent {np.median(R[m,2]/SF[m]/4.81e3)/ot:6.1f}  Maxwell {np.median(R[m,3]/SF[m]/4.81e3)/ot:6.1f}')
    P.save(fig, 'yield_components'); plt.close(fig)

@figure
def yield_regime():
    """where PRFM holds and where it fails, without putting the star formation rate on both axes: the measured 40 Myr rate against the
    rate the OK22 relations return for the measured weight (eq. 28b) and pressure (eq. 28a) of the same snapshot"""
    import glob as _g
    rows = []
    for fn in sorted(_g.glob(f'{C.PRFM_DIR}/patch_[0-9][0-9][0-9]_z05.h5')):
        with h5py.File(fn, 'r') as f:
            t = float(f.attrs['time_myr']); acc = np.zeros(4)
            for g in f:
                if not g.startswith('frame_'): continue
                G = f[g]; sep = float(G.attrs['separation_kpc']); r = lambda k: G[k][:].ravel()
                S = r('Sigma_gas'); W = r('W_2p'); fi = r('f_intruder'); ok = (S > 1) & (W > 0) & ((fi < 0.1) | C.is_merged(sep, t))
                acc += [r('Ptot_2p')[ok].sum(), r('W_2p')[ok].sum(), r('SigSFR_40')[ok].sum(), ok.sum()]
            rows.append((t, *acc))
    R = np.array(rows); T = R[:, 0]; n = R[:, 4]
    Pm = R[:, 1] / n; Wm = R[:, 2] / n; sfr = R[:, 3] / n; o = T > 40
    fig, ax = P.fig()
    for (t0, t1, lab), c in zip(C.PHASES, (P.light, P.blue, P.orange, P.ink)):
        m = o & (T >= t0) & (T < t1)
        if m.sum(): ax.scatter(ok22.sfr_of_W(Wm[m]), sfr[m], s=16, color=c, lw=0, label=lab.replace('\n', ' '))
    g = np.logspace(-6, 0, 10)
    ax.plot(g, g, color=P.ink, lw=0.9); ax.plot(g, g / 10, color=P.grey, lw=0.6, ls=':'); ax.plot(g, g / 100, color=P.grey, lw=0.6, ls=':')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(2e-4, 3e-2); ax.set_ylim(1e-5, 3e-2)
    ax.set_xlabel(r'$\Sigma_{\rm SFR}$ predicted by PRFM from $\mathcal{W}_{\rm 2p}$ [M$_\odot$ yr$^{-1}$ kpc$^{-2}$]')
    ax.set_ylabel(r'$\Sigma_{\rm SFR,40}$ measured [M$_\odot$ yr$^{-1}$ kpc$^{-2}$]')
    P.place_legend(ax)
    for a, b, lab in ((40, 100, 'discs'), (125, 140, '2nd burst'), (169, 195, 'trough'), (195, 226, '3rd burst')):
        m = (T >= a) & (T < b)
        print(f'  {lab:10s} measured / predicted from W = {np.median(sfr[m] / ok22.sfr_of_W(Wm[m])):6.3f}   from P = {np.median(sfr[m] / ok22.sfr_of_Ptot(Pm[m])):6.3f}')
    P.save(fig, 'yield_regime'); plt.close(fig)


@figure
def gas_vs_bound():
    """is the gas gone or is it unbound?  Per snapshot: the gas mass of the clean layer columns, the mass in dense clumps
    (n_H > 100 cm^-3 groups), the part of that clump mass with alpha_tot < 4, and the stellar mass formed in the last 10 Myr."""
    import glob as _g
    rows = []
    for fn in sorted(_g.glob(f'{C.PRFM_DIR}/patch_[0-9][0-9][0-9]_z05.h5')):
        with h5py.File(fn, 'r') as f:
            t = float(f.attrs['time_myr']); mg = 0.0; A = float(f.attrs['patch_kpc']) ** 2
            for g in f:
                if not g.startswith('frame_'): continue
                G = f[g]; sep = float(G.attrs['separation_kpc']); r = lambda k: G[k][:].ravel()
                S = r('Sigma_gas'); W = r('W_2p'); fi = r('f_intruder'); ok = (S > 1) & (W > 0) & ((fi < 0.1) | C.is_merged(sep, t))
                mg += S[ok].sum() * A * 1e6
            rows.append((t, mg))
    G1 = np.array(rows)
    cl = []
    for fn in sorted(_g.glob(f'{C.DATADIR}/clouds/clouds_[0-9]*.npz')):
        d = np.load(fn, allow_pickle=True); t = float(d['t']); R = d['cl100_rows']; col = {n: i for i, n in enumerate(d['names'])}
        if len(R) == 0: cl.append((t, 0.0, 0.0)); continue
        M = R[:, col['M']]; at = R[:, col['alpha_vir']] * (1 + R[:, col['vA']] ** 2 / np.maximum(R[:, col['sig3d']] ** 2, 1e-30))
        cl.append((t, M.sum(), M[at < 4].sum()))
    Cl = np.array(cl)
    s4 = C.read(231, 4, ['StellarFormationTime', 'Masses']); tf = s4['StellarFormationTime'] * C.MYR; ms = s4['Masses'] * 1e10
    my = np.array([ms[(tf > t - 10) & (tf <= t)].sum() for t in Cl[:, 0]])
    fig, ax = P.fig()
    ax.plot(G1[:, 0], G1[:, 1], color=P.grey, lw=1.2, ls=':', label='layer gas')
    ax.plot(Cl[:, 0], Cl[:, 1], color=P.orange, lw=1.5, label=r'clumps, $n_{\rm H} > 100$ cm$^{-3}$')
    ax.plot(Cl[:, 0], Cl[:, 2], color=P.blue, lw=1.5, label=r'of those, $\alpha_{\rm tot} < 4$')
    ax.plot(Cl[:, 0], my, color=P.ink, lw=1.4, label='stars, 10 Myr')
    ax.set_yscale('log'); ax.set_xlim(25, 226); ax.set_ylim(1e3, 2e9)
    ax.set_xlabel(r'$t$ [Myr]'); ax.set_ylabel(r'mass [M$_\odot$]'); P.phases(ax, box=True)
    P.place_legend(ax)
    for a, b, lab in ((40, 100, 'discs'), (169, 195, 'trough')):
        m = (Cl[:, 0] >= a) & (Cl[:, 0] < b); mm = (G1[:, 0] >= a) & (G1[:, 0] < b)
        print(f'  {lab:8s} layer gas {np.median(G1[mm,1]):9.3g}  clumps {np.median(Cl[m,1]):9.3g}  of those bound {np.median(Cl[m,2]):9.3g}  stars/10 Myr {np.median(my[m]):9.3g}')
    P.save(fig, 'gas_vs_bound'); plt.close(fig)

@figure
def alpha_chain():
    """the causal chain in the trough: the effective dispersion of the layer sets the clump dispersions, which set alpha_tot, which
    sets the bound fraction of the dense gas.  Top: sigma_eff of the clean layer columns (Sigma-weighted median) and the median
    sigma_3d of the n_H > 100 clumps.  Bottom: median alpha_tot of those clumps and the fraction of clump mass with alpha_tot < 4."""
    import glob as _g
    lay = []
    for fn in sorted(_g.glob(f'{C.PRFM_DIR}/patch_[0-9][0-9][0-9]_z05.h5')):
        with h5py.File(fn, 'r') as f:
            t = float(f.attrs['time_myr']); Pl = []; Rl = []; Sl = []
            for g in f:
                if not g.startswith('frame_'): continue
                G = f[g]; sep = float(G.attrs['separation_kpc']); r = lambda k: G[k][:].ravel()
                S = r('Sigma_gas'); W = r('W_2p'); fi = r('f_intruder'); ok = (S > 1) & (W > 0) & ((fi < 0.1) | C.is_merged(sep, t))
                Pl.append(r('Ptot_2p')[ok]); Rl.append(r('rho_mid_2p')[ok]); Sl.append(S[ok])
            Pl = np.concatenate(Pl); Rl = np.concatenate(Rl); Sl = np.concatenate(Sl)
            se = np.sqrt(np.maximum(Pl, 0) / np.maximum(Rl * 1e9, 1e-30) / 4.903245584337325e-06)
            o = np.argsort(se); cw = np.cumsum(Sl[o]) / Sl.sum(); lay.append((t, se[o][np.searchsorted(cw, 0.5)]))
    L = np.array(lay)
    cl = []
    for fn in sorted(_g.glob(f'{C.DATADIR}/clouds/clouds_[0-9]*.npz')):
        d = np.load(fn, allow_pickle=True); t = float(d['t']); R = d['cl100_rows']; col = {n: i for i, n in enumerate(d['names'])}
        if len(R) < 5: continue
        M = R[:, col['M']]; sg = R[:, col['sig3d']]
        at = R[:, col['alpha_vir']] * (1 + R[:, col['vA']] ** 2 / np.maximum(sg ** 2, 1e-30))
        cl.append((t, np.median(sg), np.median(at), M[at < 4].sum() / M.sum()))
    Cl = np.array(cl)
    fig, ax = plt.subplots(2, 1, figsize=(P.COL1, P.COL1 * 1.05), sharex=True, gridspec_kw=dict(hspace=0.06, height_ratios=[1, 1]))
    ax[0].plot(L[:, 0], L[:, 1], color=P.orange, lw=1.5, label=r'$\sigma_{\rm eff}$ of the layer')
    ax[0].plot(Cl[:, 0], Cl[:, 1], color=P.blue, lw=1.5, label=r'$\sigma_{\rm 3d}$ of the clumps')
    ax[0].set_yscale('log'); ax[0].set_ylim(0.5, 60); ax[0].set_ylabel(r'$\sigma$ [km s$^{-1}$]'); ax[0].set_xlim(25, 226)
    ax[1].plot(Cl[:, 0], Cl[:, 2], color=P.ink, lw=1.5, label=r'median $\alpha_{\rm tot}$')
    ax[1].axhline(4, color=P.grey, lw=0.8, ls='--')
    a2 = ax[1].twinx(); a2.plot(Cl[:, 0], Cl[:, 3], color=P.blue, lw=1.3, ls=':'); a2.set_ylim(0, 1.05)
    a2.set_ylabel(r'bound fraction of clump mass', color=P.blue); a2.tick_params(axis='y', colors=P.blue, direction='in')
    ax[1].set_yscale('log'); ax[1].set_ylim(0.7, 400); ax[1].set_ylabel(r'$\alpha_{\rm tot}$'); ax[1].set_xlabel(r'$t$ [Myr]'); ax[1].set_xlim(25, 226)
    P.phases(ax[0], box=True); P.place_legend(ax[0])
    for tp in P.PERICENTRES: ax[1].axvline(tp, color=P.grey, lw=0.6, ls=':' if tp < 200 else '-')
    P.place_legend(ax[1])
    P.save(fig, 'alpha_chain'); plt.close(fig)

@figure
def yield_forward():
    """the yield when the pressure at t is paired with the star formation that follows it: x = the forward offset Delta, y = the
    phase-median of [sum_c P_tot,c(t) / sum_c Sigma_SFR,40,c(t + Delta)] / Ups_OK22(Pbar(t)) over the clean layer columns.
    Delta = 0 is the usual backward 40 Myr window."""
    import glob as _g
    rows = []
    for fn in sorted(_g.glob(f'{C.PRFM_DIR}/patch_[0-9][0-9][0-9]_z05.h5')):
        with h5py.File(fn, 'r') as f:
            t = float(f.attrs['time_myr']); acc = np.zeros(4)
            for g in f:
                if not g.startswith('frame_'): continue
                G = f[g]; sep = float(G.attrs['separation_kpc']); r = lambda k: G[k][:].ravel()
                S = r('Sigma_gas'); W = r('W_2p'); fi = r('f_intruder'); ok = (S > 1) & (W > 0) & ((fi < 0.1) | C.is_merged(sep, t))
                acc += [r('Ptot_2p')[ok].sum(), r('SigSFR_40')[ok].sum(), (r('Sigma_gas_2p') * r('Ptot_2p'))[ok].sum(), r('Sigma_gas_2p')[ok].sum()]
            rows.append((t, *acc))
    R = np.array(rows); T = R[:, 0]; Pb = R[:, 3] / R[:, 4]
    D = np.arange(0, 56, 2.0)
    sets = [((40, 100), 'discs, 40 to 100 Myr', P.light), ((107, 125), 'after 2nd passage', P.blue),
            ((125, 165), 'burst and decay', P.green), ((169, 195), 'trough, 169 to 195 Myr', P.orange)]
    fig, ax = P.fig()
    for (a, b), lab, c in sets:
        y = []
        for d in D:
            k = int(round(d / C.SNAP_DT)); n = len(T) - k
            rat = R[:n, 1] / np.maximum(R[k:, 2], 1e-30) / 4.81e3 / ok22.ups_tot(Pb[:n]); tt = T[:n]
            m = (tt >= a) & (tt < b); y.append(np.median(rat[m]) if m.sum() else np.nan)
        ax.plot(D, y, color=c, lw=1.5, label=lab)
        if lab.startswith('trough'): print('  trough: ' + '  '.join(f'{d:.0f}:{v:.1f}' for d, v in zip(D[::5], np.array(y)[::5])))
    ax.axhline(1, color=P.ink, lw=0.8); ax.axhline(3, color=P.grey, lw=0.6, ls=':')
    ax.set_yscale('log'); ax.set_xlim(0, 54); ax.set_ylim(0.5, 1e3)
    ax.set_xlabel(r'star formation measured $\Delta$ Myr after the pressure')
    ax.set_ylabel(r'$(P_{\rm tot,2p}/\Sigma_{\rm SFR})\,/\,\Upsilon_{\rm tot}^{\rm OK22}$')
    P.place_legend(ax)
    P.save(fig, 'yield_forward'); plt.close(fig)

@figure
def weight_sources():
    """the weight of the layer and what provides the gravity: per snapshot the Sigma_gas-weighted median over the clean layer columns of
    W = 1/2 (W_up + W_lo) for all matter and for the gas, the stars and the dark matter separately, and the measured P_tot for comparison"""
    import glob as _g
    rows = []
    for fn in sorted(_g.glob(f'{C.PRFM_DIR}/patch_[0-9][0-9][0-9]_z05.h5')):
        with h5py.File(fn, 'r') as f:
            t = float(f.attrs['time_myr']); acc = {k: [] for k in ('S', 'W', 'P', 'g', 's', 'd')}
            for g in f:
                if not g.startswith('frame_'): continue
                G = f[g]; sep = float(G.attrs['separation_kpc']); r = lambda k: G[k][:].ravel()
                S = r('Sigma_gas'); W = r('W_2p'); fi = r('f_intruder'); ok = (S > 1) & (W > 0) & ((fi < 0.1) | C.is_merged(sep, t))
                half = lambda a, b: 0.5 * (r(a)[ok] + r(b)[ok])
                acc['S'].append(S[ok]); acc['W'].append(W[ok]); acc['P'].append(r('Ptot_2p')[ok])
                acc['g'].append(half('W_gas_up', 'W_gas_lo')); acc['s'].append(half('W_star_up', 'W_star_lo')); acc['d'].append(half('W_dm_up', 'W_dm_lo'))
            A = {k: np.concatenate(v) for k, v in acc.items()}; Sg = A['S']
            def wm(x):
                o = np.argsort(x); cw = np.cumsum(Sg[o]) / Sg.sum(); return x[o][np.searchsorted(cw, 0.5)]
            rows.append((t, wm(A['W']), wm(A['P']), wm(A['g']), wm(A['s']), wm(A['d'])))
    R = np.array(rows); T = R[:, 0]; o = T > 10
    fig, ax = P.fig()
    ax.plot(T[o], R[o, 1], color=P.ink, lw=1.7, label=r'$\mathcal{W}_{\rm 2p}$, all matter')
    ax.plot(T[o], R[o, 2], color=P.grey, lw=1.0, ls=':', label=r'$P_{\rm tot,2p}$')
    ax.plot(T[o], R[o, 4], color=P.orange, lw=1.4, label='from the stars')
    ax.plot(T[o], R[o, 5], color=P.blue, lw=1.4, label='from the dark matter')
    ax.plot(T[o], R[o, 3], color=P.light, lw=1.4, label='from the gas')
    ax.set_yscale('log'); ax.set_xlim(10, 226); ax.set_ylim(3e1, 3e6)
    ax.set_xlabel(r'$t$ [Myr]'); ax.set_ylabel(r'$\mathcal{W}/k_{\rm B}$ [K cm$^{-3}$]'); P.phases(ax, box=True); P.place_legend(ax, ncol=2)
    for a, b, lab in ((40, 100, 'discs'), (110, 125, 'after 2nd'), (169, 195, 'trough'), (195, 226, '3rd burst')):
        m = (T >= a) & (T < b); g = lambda i: np.median(R[m, i])
        print(f'  {lab:10s} W {g(1):9.3g}  gas {g(3)/g(1):.2f}  stars {g(4)/g(1):.2f}  dark matter {g(5)/g(1):.2f}')
    P.save(fig, 'weight_sources'); plt.close(fig)

def _cols(keys, step=4):
    """every clean layer column of every step-th snapshot: dict of arrays plus the time"""
    import glob as _g
    out = {k: [] for k in list(keys) + ['t']}
    for fn in sorted(_g.glob(f'{C.PRFM_DIR}/patch_[0-9][0-9][0-9]_z05.h5')):
        k = int(os.path.basename(fn)[6:9])
        if k % step: continue
        with h5py.File(fn, 'r') as f:
            t = float(f.attrs['time_myr'])
            for g in f:
                if not g.startswith('frame_'): continue
                G = f[g]; sep = float(G.attrs['separation_kpc']); r = lambda kk: G[kk][:].ravel()
                S = r('Sigma_gas'); W = r('W_2p'); fi = r('f_intruder'); ok = (S > 1) & (W > 0) & ((fi < 0.1) | C.is_merged(sep, t))
                for kk in keys: out[kk].append(r(kk)[ok])
                out['t'].append(np.full(ok.sum(), t))
    return {k: np.concatenate(v) for k, v in out.items()}

@figure
def balance_scatter():
    """the vertical balance column by column: P_tot against W for every clean layer column of every 4th snapshot, coloured by phase"""
    d = _cols(['Ptot_2p', 'W_2p'])
    fig, ax = P.fig()
    for (t0, t1, lab), c in zip(C.PHASES, (P.light, P.blue, P.orange, P.ink)):
        m = (d['t'] >= t0) & (d['t'] < t1)
        ax.scatter(d['W_2p'][m], d['Ptot_2p'][m], s=3, color=c, alpha=0.45, lw=0, rasterized=True, label=lab.replace('\n', ' '))
    g = np.logspace(1, 7, 10); ax.plot(g, g, color=P.ink, lw=0.9); ax.plot(g, 2 * g, color=P.grey, lw=0.6, ls=':'); ax.plot(g, g / 2, color=P.grey, lw=0.6, ls=':')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(1e2, 3e6); ax.set_ylim(1e2, 3e6)
    ax.set_xlabel(r'$\mathcal{W}_{\rm 2p}$ [K cm$^{-3}$]'); ax.set_ylabel(r'$P_{\rm tot,2p}$ [K cm$^{-3}$]'); P.place_legend(ax)
    r = d['Ptot_2p'] / d['W_2p']
    print(f'  {len(r)} columns: median P/W {np.median(r):.2f}, 16-84 % {np.percentile(r,16):.2f}-{np.percentile(r,84):.2f}, within a factor 2: {np.mean((r>0.5)&(r<2)):.2f}')
    P.save(fig, 'balance_scatter'); plt.close(fig)

@figure
def sfr_pressure_scatter():
    """the PRFM star formation relation column by column: Sigma_SFR over 40 Myr against the layer P_tot, with the OK22 relation"""
    d = _cols(['Ptot_2p', 'W_2p', 'SigSFR_40'])
    ok = (d['SigSFR_40'] > 0) & (d['Ptot_2p'] > 0)   # P_tot can be negative where the Maxwell stress is
    fig, ax = P.fig()
    for (t0, t1, lab), c in zip(C.PHASES, (P.light, P.blue, P.orange, P.ink)):
        m = ok & (d['t'] >= t0) & (d['t'] < t1)
        ax.scatter(d['Ptot_2p'][m], d['SigSFR_40'][m], s=3, color=c, alpha=0.45, lw=0, rasterized=True, label=lab.replace('\n', ' '))
    g = np.logspace(2, 7, 10); ax.plot(g, ok22.sfr_of_Ptot(g), color=P.ink, lw=1.1, label='OK22 eq. 28a')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(3e2, 3e6); ax.set_ylim(1e-5, 3e-1)
    ax.set_xlabel(r'$P_{\rm tot,2p}$ [K cm$^{-3}$]'); ax.set_ylabel(r'$\Sigma_{\rm SFR,40}$ [M$_\odot$ yr$^{-1}$ kpc$^{-2}$]'); P.place_legend(ax)
    for t0, t1, lab in C.PHASES:
        m = ok & (d['t'] >= t0) & (d['t'] < t1)
        if m.sum() > 20: print(f'  {lab.replace(chr(10)," "):28s} median log(measured/OK22) = {np.median(np.log10(d["SigSFR_40"][m] / ok22.sfr_of_Ptot(d["Ptot_2p"][m]))):+.2f} dex, scatter {np.std(np.log10(d["SigSFR_40"][m] / ok22.sfr_of_Ptot(d["Ptot_2p"][m]))):.2f}')
    P.save(fig, 'sfr_pressure_scatter'); plt.close(fig)

@figure
def pde_vs_weight():
    """the analytic estimate of the weight used by OK22 and by Kruijssen (2012), P_DE, against the weight measured from the
    particle-mesh solve, column by column"""
    d = _cols(['P_DE', 'W_2p'])
    ok = (d['P_DE'] > 0) & (d['W_2p'] > 0)
    fig, ax = P.fig()
    for (t0, t1, lab), c in zip(C.PHASES, (P.light, P.blue, P.orange, P.ink)):
        m = ok & (d['t'] >= t0) & (d['t'] < t1)
        ax.scatter(d['W_2p'][m], d['P_DE'][m], s=3, color=c, alpha=0.45, lw=0, rasterized=True, label=lab.replace('\n', ' '))
    g = np.logspace(1, 7, 10); ax.plot(g, g, color=P.ink, lw=0.9); ax.plot(g, 2 * g, color=P.grey, lw=0.6, ls=':'); ax.plot(g, g / 2, color=P.grey, lw=0.6, ls=':')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(1e2, 3e6); ax.set_ylim(1e2, 3e6)
    ax.set_xlabel(r'$\mathcal{W}_{\rm 2p}$, particle-mesh [K cm$^{-3}$]'); ax.set_ylabel(r'$P_{\rm DE}$, analytic [K cm$^{-3}$]'); P.place_legend(ax)
    r = d['P_DE'][ok] / d['W_2p'][ok]
    print(f'  {ok.sum()} columns: median P_DE/W {np.median(r):.2f}, 16-84 % {np.percentile(r,16):.2f}-{np.percentile(r,84):.2f}')
    for t0, t1, lab in C.PHASES:
        m = ok & (d['t'] >= t0) & (d['t'] < t1)
        if m.sum() > 20: print(f'  {lab.replace(chr(10)," "):28s} median {np.median(d["P_DE"][m]/d["W_2p"][m]):.2f}')
    P.save(fig, 'pde_vs_weight'); plt.close(fig)

@figure
def tdyn_vs_delay():
    """the vertical dynamical time of the layer, t_dyn = 2 H / sigma_eff (OK22 sec. 2), against the delay between a pericentre and the
    burst it produces.  OK22 define their equilibrium as an average over a few t_dyn, so this is the timescale on which their
    quasi-steady assumption is defined."""
    import glob as _g
    rows = []
    for fn in sorted(_g.glob(f'{C.PRFM_DIR}/patch_[0-9][0-9][0-9]_z05.h5')):
        with h5py.File(fn, 'r') as f:
            t = float(f.attrs['time_myr']); H = []; SE = []; S = []
            for g in f:
                if not g.startswith('frame_'): continue
                G = f[g]; sep = float(G.attrs['separation_kpc']); r = lambda k: G[k][:].ravel()
                Sg = r('Sigma_gas'); W = r('W_2p'); fi = r('f_intruder'); ok = (Sg > 1) & (W > 0) & ((fi < 0.1) | C.is_merged(sep, t))
                H.append(r('H')[ok]); SE.append(r('sigma_eff')[ok]); S.append(Sg[ok])
            H = np.concatenate(H); SE = np.concatenate(SE); S = np.concatenate(S)
            td = 2 * H / np.maximum(SE, 1e-30) * 977.8
            o = np.argsort(td); cw = np.cumsum(S[o]) / S.sum()
            rows.append((t, td[o][np.searchsorted(cw, 0.5)], td[o][np.searchsorted(cw, 0.16)], td[o][np.searchsorted(cw, 0.84)]))
    R = np.array(rows); T = R[:, 0]; o = T > 25
    fig, ax = P.fig()
    ax.fill_between(T[o], R[o, 2], R[o, 3], color=P.light, alpha=0.5, lw=0, label='16 to 84 % of the columns')
    ax.plot(T[o], R[o, 1], color=P.blue, lw=1.6, label=r'$t_{\rm dyn} = 2H/\sigma_{\rm eff}$')
    for tp, pk, lab in ((109.5, 126.1, '2nd'), (169.2, 204.4, '3rd')):
        ax.annotate('', xy=(pk, 42), xytext=(tp, 42), arrowprops=dict(arrowstyle='<->', color=P.orange, lw=1.2))
        ax.text(0.5 * (tp + pk), 44, f'{pk-tp:.0f} Myr', ha='center', va='bottom', fontsize=6.5, color=P.orange)
    ax.set_xlim(25, 226); ax.set_ylim(0, 58); ax.set_xlabel(r'$t$ [Myr]'); ax.set_ylabel(r'$t_{\rm dyn}$ [Myr]')
    P.phases(ax, labels=False); P.place_legend(ax)
    print(f'  t_dyn over t>40: median {np.median(R[T>40,1]):.1f} Myr, 16-84 % {np.percentile(R[T>40,1],16):.0f}-{np.percentile(R[T>40,1],84):.0f}')
    P.save(fig, 'tdyn_vs_delay'); plt.close(fig)

@figure
def calibration_compare():
    """the measured yield against three calibrations: OK22 (TIGRESS-classic, solar metallicity, supernovae and FUV only),
    TIGRESS-NCR at solar metallicity, and TIGRESS-NCR at the 0.1 solar metallicity of this run (Kim et al. 2024)"""
    import glob as _g
    rows = []
    for fn in sorted(_g.glob(f'{C.PRFM_DIR}/patch_[0-9][0-9][0-9]_z05.h5')):
        with h5py.File(fn, 'r') as f:
            t = float(f.attrs['time_myr']); acc = np.zeros(5)
            for g in f:
                if not g.startswith('frame_'): continue
                G = f[g]; sep = float(G.attrs['separation_kpc']); r = lambda k: G[k][:].ravel()
                S = r('Sigma_gas'); W = r('W_2p'); fi = r('f_intruder'); ok = (S > 1) & (W > 0) & ((fi < 0.1) | C.is_merged(sep, t))
                acc += [r('Ptot_2p')[ok].sum(), r('SigSFR_40')[ok].sum(), r('W_2p')[ok].sum(), ok.sum(), (r('Sigma_gas_2p') * r('Ptot_2p'))[ok].sum()]
            rows.append((t, *acc))
    R = np.array(rows); T = R[:, 0]; n = R[:, 4]
    ups = R[:, 1] / np.maximum(R[:, 2], 1e-30) / 4.81e3          # measured, ratio of sums
    W4 = (R[:, 3] / n) / 1e4; Pm = R[:, 1] / n
    o = T > 40
    fig, ax = P.fig()
    ax.plot(T[o], ups[o], color=P.ink, lw=1.7, label='measured')
    ax.plot(T[o], ok22.ncr_ups_tot(W4[o], 0.1), color=P.orange, lw=1.4, label=r'TIGRESS-NCR, $Z = 0.1\,Z_\odot$')
    ax.plot(T[o], ok22.ncr_ups_tot(W4[o], 1.0), color=P.blue, lw=1.4, ls='--', label=r'TIGRESS-NCR, $Z = Z_\odot$')
    ax.plot(T[o], ok22.ups_tot(Pm[o]), color=P.grey, lw=1.2, ls=':', label='OK22 (TIGRESS-classic)')
    ax.set_yscale('log'); ax.set_xlim(40, 226); ax.set_ylim(3e2, 3e5)
    ax.set_xlabel(r'$t$ [Myr]'); ax.set_ylabel(r'$\Upsilon_{\rm tot} = P_{\rm tot,2p}/\Sigma_{\rm SFR,40}$ [km s$^{-1}$]')
    P.phases(ax, box=True); P.place_legend(ax)
    for a, b, lab in ((40, 100, 'discs'), (169, 195, 'trough')):
        m = (T >= a) & (T < b)
        print(f'  {lab:8s} measured {np.median(ups[m]):8.0f}  NCR(0.1) {np.median(ok22.ncr_ups_tot(W4[m],0.1)):8.0f}  NCR(1) {np.median(ok22.ncr_ups_tot(W4[m],1.0)):8.0f}  OK22 {np.median(ok22.ups_tot(Pm[m])):8.0f}')
    P.save(fig, 'calibration_compare'); plt.close(fig)

if __name__ == '__main__':
    names = sys.argv[1:] or ['all']
    for n in (FIGS if names == ['all'] else names): FIGS[n]()
