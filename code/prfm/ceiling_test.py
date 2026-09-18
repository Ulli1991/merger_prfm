"""Is the largest cluster per 25 Myr a physical ceiling set by the weight, or the extreme value of N draws from the MF?
(a) per 25 Myr bin: M_max vs W_L (SF patches), vs the median midplane density rho_mid,2p of ALL clean patches, vs W_L/sigma_eff^2,
    vs the number of bound clusters N_cl and the total young stellar mass in the bin;
(b) extreme-value prediction: for a power law of index alpha above M_min=300, the expected maximum of N draws is
    M_min N^(1/(alpha-1)) (median: M_min (N/ln2)^(1/(alpha-1))); compare with the observed M_max per bin using alpha from the bin's
    own clusters (and a fixed alpha=1.8) - if it matches, M_max is a size-of-sample statistic, not a mass scale."""
import sys, numpy as np
sys.path.insert(0, '/raven/u/uli/dwarf_merger'); import common as C
from scipy.stats import spearmanr
D = dict(np.load(f'{C.DATADIR}/prfm/patches_clusters.npz')); WL = np.load(f'{C.DATADIR}/prfm/cluster_link.npz')['WL']
E = np.load(f'{C.DATADIR}/clusters/clusters_env.npz')
clean = (D['f_intruder'] < 0.1) | C.is_merged(D['sep'], D['t']); base = clean & (D['Sigma_gas'] > 1) & np.isfinite(WL) & (WL > 0) & (D['t'] > 10)
cleanE = (E['f_intruder'] < 0.1) | C.is_merged(E['sep'], E['t']); cb = cleanE & (E['bound'] > 0) & (E['nucleus'] == 0)
edges = np.arange(10, 235, 25); rows = []
print('bin    M_max    N_cl(indep)  alpha_bin  Myoung_tot   W_L(SF med)  W_L(all med)  rho_mid(all med)  rho_mid(all 90%)  W/sig_eff^2(med)  f_burst | EV pred alpha_bin   EV pred a=1.8')
for t0, t1 in zip(edges[:-1], edges[1:]):
    mc = cb & (E['t'] >= t0) & (E['t'] < t1); mci = mc & (np.mod(E['snap'], 10) == 0)      # independent populations for N
    ma = base & (D['t'] >= t0) & (D['t'] < t1); ms = ma & (D['M_young'] > 500); mai = ma & (np.mod(D['snap'], 10) == 0)
    if mc.sum() < 5 or ma.sum() < 10: continue
    M = E['M'][mci]; Mf = M[M >= 300]; N = len(Mf); al = 1 + N / np.sum(np.log(Mf / 300.)) if N > 5 else np.nan
    Mmax = E['M'][mc].max(); ev = 300 * (N / np.log(2)) ** (1 / (al - 1)) if N > 5 else np.nan; ev18 = 300 * (N / np.log(2)) ** (1 / 0.8)
    Myt = D['M_young'][mai].sum()
    rows.append((t0, Mmax, N, al, Myt, np.median(WL[ms]) if ms.sum() else np.nan, np.median(WL[ma]), np.median(D['rho_mid_2p'][ma]), np.percentile(D['rho_mid_2p'][ma], 90), np.median(WL[ma] / D['sigma_eff'][ma] ** 2), ms.sum() / ma.sum(), ev, ev18))
    r = rows[-1]; print(f'{t0:3.0f}  {r[1]:8.2e}   {r[2]:4d}       {r[3]:5.2f}    {r[4]:8.2e}   {r[5]:9.2e}   {r[6]:9.2e}   {r[7]:9.3e}       {r[8]:9.3e}       {r[9]:9.2e}     {r[10]:.2f}   | {r[11]:9.2e}       {r[12]:9.2e}')
R = np.array(rows); names = ['t0', 'Mmax', 'N', 'alpha', 'Myoung', 'WL_sf', 'WL_all', 'rho_med', 'rho_90', 'Wsig2', 'fburst', 'ev', 'ev18']
print('\nSpearman of log M_max with each quantity over ALL bins (n=%d), and log-log slope:' % len(R))
for i, n in enumerate(names[2:], 2):
    ok = np.isfinite(R[:, i]) & (R[:, i] > 0)
    print(f'  {n:8s}: rho = {spearmanr(np.log10(R[ok, 1]), np.log10(R[ok, i]))[0]:+.2f}   slope = {np.polyfit(np.log10(R[ok, i]), np.log10(R[ok, 1]), 1)[0]:+.2f}')
q = R[:, 0] == 160
print(f'\nquiescent bin (160): M_max {R[q,1][0]:.1e}, W_L(all) {R[q,6][0]:.1e} vs other bins median {np.median(R[~q,6]):.1e}; rho_mid(all) {R[q,7][0]:.2e} vs others {np.median(R[~q,7]):.2e} (range {R[~q,7].min():.2e}-{R[~q,7].max():.2e}); rho_90 {R[q,8][0]:.2e} vs others {np.median(R[~q,8]):.2e}; W/sig^2 {R[q,9][0]:.1e} vs others {np.median(R[~q,9]):.1e}')
print(f'extreme-value check: log10(M_max_obs / EV_pred) per bin with the bin alpha: ' + ' '.join(f'{v:+.2f}' for v in np.log10(R[:, 1] / R[:, 11])) + f'  -> mean {np.nanmean(np.log10(R[:,1]/R[:,11])):+.2f}, scatter {np.nanstd(np.log10(R[:,1]/R[:,11])):.2f} dex')
print(f'                     with alpha = 1.8 for all bins:                      ' + ' '.join(f'{v:+.2f}' for v in np.log10(R[:, 1] / R[:, 12])) + f'  -> mean {np.nanmean(np.log10(R[:,1]/R[:,12])):+.2f}, scatter {np.nanstd(np.log10(R[:,1]/R[:,12])):.2f} dex')
