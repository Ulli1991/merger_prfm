"""PRFM per snapshot: equilibrium (P/W), pressure budget, feedback yield vs OK22, and turbulent energy budget
(SN input vs dissipation vs compression work).  Reads prfm/patch_XXX.h5 (all frames), clean patches with Sigma_gas > 1.
out: prfm/prfm_time.npz, prfm/figs/prfm_time.png (+ ~/out)"""
import sys, os, glob, numpy as np, h5py
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C, ok22
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
KB = 1.3807e-16; MSUN = 1.989e33; PC = 3.0857e18; KPC = 1e3 * PC; YR = 3.156e7; MYR = 1e6 * YR
SFR_CGS = MSUN / YR / KPC**2            # Msun/yr/kpc^2 -> g/s/cm^2
suffix = sys.argv[1] if len(sys.argv) > 1 else ''      # '' or '_los'
rows = []
for k in range(C.NSNAP):
    fn = f'{C.DATADIR}/prfm/patch_{k:03d}{suffix}.h5'
    if not os.path.exists(fn): continue
    with h5py.File(fn, 'r') as f:
        t = f.attrs['time_myr']; PU = f.attrs['P_unit']
        for g in f:
            G = f[g]; sep = G.attrs['separation_kpc']
            d = {key: G[key][:].ravel() for key in ['Sigma_gas', 'W', 'W_2p', 'Ptot', 'Ptot_2p', 'Pth_neu', 'Pth_ion', 'Pth_hot', 'Pturb', 'Pturb_2p',
                                                     'Pmag', 'Pmag_2p', 'SigSFR_10', 'SNrate', 'rho_mid_2p', 'H', 'div_v', 'shear', 'f_intruder']}
            sel = (d['Sigma_gas'] > 1) & ((d['f_intruder'] < 0.1) | C.is_merged(sep, t)) & (d['W_2p'] > 0) & (d['Ptot_2p'] > 0)
            if sel.sum() < 3: continue
            d = {key: v[sel] for key, v in d.items()}
            pw = d['Ptot_2p'] / d['W_2p']
            Pt = d['Ptot'].sum()
            fr = dict(warm=(d['Pth_neu'] + d['Pth_ion']).sum() / Pt, hot=d['Pth_hot'].sum() / Pt, turb=d['Pturb'].sum() / Pt, mag=d['Pmag'].sum() / Pt)
            meanP = d['Ptot_2p'].mean(); meanS = d['SigSFR_10'].mean()
            ups = KB * meanP / max(meanS * SFR_CGS, 1e-40) / 1e5                    # km/s
            ups_ok = ok22.ups_tot(meanP)
            sfr_off = np.log10(max(meanS, 1e-12) / ok22.sfr_of_Ptot(meanP))          # dex, binned-mean style
            sig = np.sqrt(np.maximum(d['Pturb_2p'] / PU / np.maximum(d['rho_mid_2p'] * 1e9, 1e-30), 0))   # km/s
            Sg = d['Sigma_gas'] * MSUN / PC**2; H = np.maximum(d['H'], 0.02) * KPC
            e_diss = Sg * (sig * 1e5)**3 / (2 * H) * KPC**2                          # erg/s/kpc^2
            e_sn = d['SNrate'] * 1e51 / MYR                                           # erg/s/kpc^2 (1e51 per SN)
            e_comp = KB * d['Ptot_2p'] * np.maximum(-d['div_v'], 0) * 1e5 / KPC * 2 * H * KPC**2
            e_shear = KB * d['Ptot_2p'] * d['shear'] * 1e5 / KPC * 2 * H * KPC**2
            rows.append((k, t, sep, g, sel.sum(), np.median(pw), np.percentile(pw, 16), np.percentile(pw, 84), np.median(d['Ptot'] / d['W']),
                         fr['warm'], fr['hot'], fr['turb'], fr['mag'], ups, ups_ok, sfr_off, e_sn.mean(), e_diss.mean(), e_comp.mean(), e_shear.mean(),
                         np.median(sig), d['f_intruder'].mean()))
names = 'snap t sep frame n pw pw16 pw84 pw_all f_warm f_hot f_turb f_mag ups ups_ok sfr_off e_sn e_diss e_comp e_shear sig_turb f_intr'.split()
R = {nm: np.array([r[i] for r in rows]) for i, nm in enumerate(names)}
np.savez(f'{C.DATADIR}/prfm/prfm_time{suffix}.npz', **R)
# ---- per-snapshot table (frames pooled by median) and figure
snaps = np.unique(R['snap']); T = np.array([R['t'][R['snap'] == s][0] for s in snaps]); SEP = np.array([R['sep'][R['snap'] == s][0] for s in snaps])
def per(key, f=np.median): return np.array([f(R[key][R['snap'] == s]) for s in snaps])
print('snap   t[Myr]  sep[kpc]  P/W(2p)  P/W(all)  f_warm f_hot f_turb f_mag   Ups/Ups_OK22  SFR off[dex]  E_diss/E_SN  E_comp/E_SN  sig_turb')
for i, s in enumerate(snaps):
    if s % 8: continue
    print(f'{s:4d}  {T[i]:6.1f}  {SEP[i]:6.2f}    {per("pw")[i]:5.2f}    {per("pw_all")[i]:5.2f}    {per("f_warm")[i]:.2f}  {per("f_hot")[i]:.2f}  {per("f_turb")[i]:.2f}  {per("f_mag")[i]:.2f}     '
          f'{per("ups")[i]/per("ups_ok")[i]:8.1f}      {per("sfr_off")[i]:+5.2f}      {per("e_diss")[i]/max(per("e_sn")[i],1e-30):8.1f}    {per("e_comp")[i]/max(per("e_sn")[i],1e-30):8.1f}    {per("sig_turb")[i]:5.1f}')
ink = '#1f2937'; blue = '#2f6fb2'; light = '#9dbfe3'
fig, ax = plt.subplots(5, 1, figsize=(8, 12), sharex=True, gridspec_kw=dict(hspace=0.12))
for a in ax: a.grid(alpha=0.2); a.spines[['top', 'right']].set_visible(False)
ax[0].plot(T, SEP, color=blue, lw=2); ax[0].set_ylabel('separation [kpc]'); ax[0].set_title(f'PRFM per snapshot{suffix}: clean patches, $\\Sigma_{{\\rm gas}}>1$', loc='left', color=ink)
ax[1].fill_between(T, per('pw16'), per('pw84'), color=light, alpha=0.6, lw=0); ax[1].plot(T, per('pw'), color=blue, lw=2, label='two-phase gas'); ax[1].plot(T, per('pw_all'), color=ink, lw=1.2, ls='--', label='all gas')
ax[1].axhline(1, color=ink, lw=0.8); ax[1].set_yscale('log'); ax[1].set_ylim(0.05, 5); ax[1].set_ylabel('$P_{\\rm tot}/\\mathcal{W}$'); ax[1].legend(frameon=False, fontsize=8)
comps = [('f_turb', 'turbulent', '#2f6fb2'), ('f_warm', 'thermal warm/cold', '#7fb3e0'), ('f_hot', 'thermal hot', '#c9dff2'), ('f_mag', 'magnetic', '#4b5563')]
ax[2].stackplot(T, *[per(k) for k, _, _ in comps], labels=[l for _, l, _ in comps], colors=[c for _, _, c in comps], lw=0); ax[2].set_ylim(0, 1); ax[2].set_ylabel('fraction of $P_{\\rm tot}$'); ax[2].legend(frameon=False, fontsize=8, ncol=4, loc='upper left')
ax[3].plot(T, per('ups') / per('ups_ok'), color=blue, lw=2, label='$\\Upsilon_{\\rm eff}/\\Upsilon_{\\rm OK22}$ (pressure per unit SFR)'); ax[3].axhline(1, color=ink, lw=0.8); ax[3].set_yscale('log'); ax[3].set_ylabel('pressure excess'); ax[3].legend(frameon=False, fontsize=8)
ax[4].plot(T, per('e_sn'), color=ink, lw=2, label='SN input ($10^{51}$ erg per SN)'); ax[4].plot(T, per('e_diss'), color=blue, lw=2, label='turbulent dissipation $\\Sigma\\sigma^3/2H$')
ax[4].plot(T, per('e_comp'), color='#7fb3e0', lw=2, label='compression work $P|\\nabla\\cdot v|2H$'); ax[4].plot(T, per('e_shear'), color='#4b5563', lw=1.2, ls='--', label='shear work')
ax[4].set_yscale('log'); ax[4].set_ylabel('erg s$^{-1}$ kpc$^{-2}$'); ax[4].legend(frameon=False, fontsize=8, ncol=2); ax[4].set_xlabel('t [Myr]')
for a in ax: a.axvline(C.T_MERGED, color=ink, lw=0.8, ls=':')
fn = f'{C.DATADIR}/prfm/figs/prfm_time{suffix}.png'; fig.savefig(fn, dpi=140, bbox_inches='tight'); os.system(f'cp {fn} ~/out/'); print('fig', fn)
