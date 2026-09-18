"""Pressure budget by source: P_fb = Upsilon_OK22(P) * Sigma_SFR(10 Myr) (what feedback regulation would supply) versus
P_flow = P_tot - P_fb (what the flow must supply), per snapshot and per phase, with correlations of P_flow against the
measured driving rates.  out: prfm/pressure_sources.npz, prfm/figs/pressure_sources.png (+ ~/out)"""
import sys, os, numpy as np, h5py
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C, ok22
from scipy.stats import spearmanr
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
KB = 1.3807e-16; MSUN = 1.989e33; PC = 3.0857e18; KPC = 1e3 * PC; YR = 3.156e7; SFR_CGS = MSUN / YR / KPC**2
pp = []   # per patch
for k in range(C.NSNAP):
    fn = f'{C.PRFM_DIR}/patch_{k:03d}.h5'
    if not os.path.exists(fn): continue
    with h5py.File(fn, 'r') as f:
        t = f.attrs['time_myr']
        for g in f:
            G = f[g]; sep = G.attrs['separation_kpc']
            d = {key: G[key][:].ravel() for key in ['Sigma_gas', 'W_2p', 'Ptot_2p', 'Pturb_2p', 'Pmag_2p', 'SigSFR_10', 'SigSFR_40', 'H', 'div_v', 'shear', 'f_intruder', 'rho_mid_2p']}
            sel = (d['Sigma_gas'] > 1) & ((d['f_intruder'] < 0.1) | C.is_merged(sep, t)) & (d['W_2p'] > 0) & (d['Ptot_2p'] > 0)
            for i in np.flatnonzero(sel):
                pp.append((k, t, sep) + tuple(d[key][i] for key in d))
names = ['snap', 't', 'sep', 'Sigma_gas', 'W_2p', 'Ptot_2p', 'Pturb_2p', 'Pmag_2p', 'SigSFR_10', 'SigSFR_40', 'H', 'div_v', 'shear', 'f_intruder', 'rho_mid_2p']
R = {nm: np.array([r[i] for r in pp]) for i, nm in enumerate(names)}
P = R['Ptot_2p']; S = R['SigSFR_10']
R['P_fb'] = ok22.ups_tot(P) * 1e5 * S * SFR_CGS / KB          # K cm^-3: pressure feedback would supply at this SFR
R['P_flow'] = P - R['P_fb']
R['e_shear'] = P * R['shear']; R['e_comp'] = P * np.maximum(-R['div_v'], 0); R['inflow'] = R['Sigma_gas'] * np.maximum(-R['div_v'], 0)
np.savez(f'{C.DATADIR}/prfm/pressure_sources.npz', **R)
snaps = np.unique(R['snap']); T = np.array([R['t'][R['snap'] == s][0] for s in snaps]); SEP = np.array([R['sep'][R['snap'] == s][0] for s in snaps])
def per(key, f=np.sum): return np.array([f(R[key][R['snap'] == s]) for s in snaps])
f_fb = np.clip(per('P_fb') / per('Ptot_2p'), 0, 1); pw = per('Ptot_2p') / per('W_2p')
print('area-summed budget per snapshot: P_fb/P_tot = fraction of the pressure that the OK22 yield times the actual SFR accounts for')
print('snap   t    sep   P/W   P_fb/P_tot   P_flow/W')
for i, s in enumerate(snaps):
    if s % 8 == 0: print(f'{s:4d} {T[i]:6.1f} {SEP[i]:5.2f}  {pw[i]:5.2f}   {f_fb[i]:6.2f}     {(1-f_fb[i])*pw[i]:6.2f}')
EP = [*C.PHASES]
print('\nper phase: Spearman of P_flow (patches with P_flow > 0) against candidate drivers')
print('phase                n   f(P_flow>0)  median P_flow/P_tot  | shear work  compression work  inflow rate  Sigma_gas  W  f_intruder  sigma_turb')
for t0, t1, lab in EP:
    m = (R['t'] >= t0) & (R['t'] < t1); q = m & (R['P_flow'] > 0)
    sig = np.sqrt(R['Pturb_2p'] / 4.903e-6 / np.maximum(R['rho_mid_2p'] * 1e9, 1e-30))
    cors = [spearmanr(R['P_flow'][q], v[q])[0] for v in (R['e_shear'], R['e_comp'], R['inflow'], R['Sigma_gas'], R['W_2p'], R['f_intruder'], sig)]
    print(f'{lab:20s} {m.sum():5d}   {q.sum()/m.sum():5.2f}       {np.median((R["P_flow"]/P)[m]):5.2f}          | ' + '  '.join(f'{c:+.2f}' for c in cors))
ink = '#1f2937'; blue = '#2f6fb2'; light = '#9dbfe3'
fig, ax = plt.subplots(3, 1, figsize=(8, 8), sharex=True, gridspec_kw=dict(hspace=0.12))
for a in ax: a.grid(alpha=0.2); a.spines[['top', 'right']].set_visible(False); a.axvline(C.T_MERGED, color=ink, lw=0.8, ls=':')
ax[0].plot(T, SEP, color=blue, lw=2); ax[0].set_ylabel('separation [kpc]'); ax[0].set_title('Pressure by source (clean patches, area-summed)', loc='left', color=ink)
ax[1].plot(T, pw, color=ink, lw=2); ax[1].axhline(1, color=ink, lw=0.8); ax[1].set_yscale('log'); ax[1].set_ylim(0.1, 3); ax[1].set_ylabel('$P_{\\rm tot}/\\mathcal{W}$')
ax[2].fill_between(T, 0, f_fb, color=blue, lw=0, label='feedback: $\\Upsilon_{\\rm OK22}\\,\\Sigma_{\\rm SFR}$'); ax[2].fill_between(T, f_fb, 1, color=light, lw=0, label='flow: $P_{\\rm tot}-\\Upsilon_{\\rm OK22}\\Sigma_{\\rm SFR}$')
ax[2].set_ylim(0, 1); ax[2].set_ylabel('fraction of $P_{\\rm tot}$'); ax[2].legend(frameon=False, fontsize=8, loc='lower left'); ax[2].set_xlabel('t [Myr]')
fn = f'{C.DATADIR}/prfm/figs/pressure_sources.png'; fig.savefig(fn, dpi=140, bbox_inches='tight'); os.system(f'cp {fn} ~/out/'); print('fig', fn)
