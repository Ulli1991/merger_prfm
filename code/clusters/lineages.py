"""Cloud lineages by particle-ID overlap, and the lifetime-integrated star formation efficiency.
Needs clouds_KKK.npz with <cat>_member_ids / <cat>_member_off (clouds.py rerun of 2026-09-15).
Link cloud c at k to cloud c' at k+1 if c' receives the largest share of c's members AND c is the largest contributor to c'
(mutual main-branch link).  A lineage = chain of such links.  Per lineage: t_start, t_end, M_max (and its time), the stars formed
from its members in each interval (t_k, t_{k+1}] summed along the chain plus the 10 Myr tail from the last snapshot, the onset state
(properties at the last snapshot BEFORE the first stars form, or the first snapshot if it already forms stars), eps_int = M_star / M_max.
usage: python lineages.py cl10|cl100 [kmax]      out: /ptmp/uli/dwarf_merger/clouds/lineages_<cat>.npz + printed analysis"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common as C
from scipy.stats import spearmanr
cat = sys.argv[1] if len(sys.argv) > 1 else 'cl100'; KMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 221
OUT = f'{C.DATADIR}/clouds'; EP = [*C.PHASES]
snaps = [k for k in range(1, KMAX + 1) if os.path.exists(f'{OUT}/clouds_{k:03d}.npz')]
D = {}
for k in snaps:
    d = np.load(f'{OUT}/clouds_{k:03d}.npz')
    if f'{cat}_member_ids' not in d.files: continue
    ids = d[f'{cat}_member_ids']; off = d[f'{cat}_member_off']; R = d[f'{cat}_rows']; ncl = len(R)
    lab = np.repeat(np.arange(ncl), np.diff(off)) if ncl else np.zeros(0, np.int64)
    o = np.argsort(ids); ids_s = ids[o]; lab_s = lab[o]
    # stars formed from each cloud's members in the next snapshot interval and in the next 10 Myr
    t = float(d['t']); tf = d['star_tform']; sc = d[f'{cat}_star_cl']; sm = d['star_mass']
    tnext = t + C.SNAP_DT
    ms1 = np.zeros(ncl); m1 = (sc >= 0) & (tf <= tnext + 1e-6); np.add.at(ms1, sc[m1], sm[m1])
    D[k] = dict(t=t, R=R, ids=ids_s, lab=lab_s, ncl=ncl, ms1=ms1, ms10=R[:, 16] if ncl else np.zeros(0))
snaps = sorted(D); print(f'{cat}: {len(snaps)} snapshots with member ids')
# ---- links
desc = {}; prog = {}
for k, k2 in zip(snaps[:-1], snaps[1:]):
    a, b = D[k], D[k2]
    if a['ncl'] == 0 or b['ncl'] == 0: continue
    j = np.searchsorted(b['ids'], a['ids']); j = np.minimum(j, len(b['ids']) - 1); hit = b['ids'][j] == a['ids']
    la = a['lab'][hit]; lb = b['lab'][j[hit]]
    pair = la.astype(np.int64) * b['ncl'] + lb; up, cnt = np.unique(pair, return_counts=True); pa = up // b['ncl']; pb = up % b['ncl']
    # main descendant of each a-cloud, main progenitor of each b-cloud
    best_d = {}; best_p = {}
    for x, y, c in zip(pa, pb, cnt):
        if c > best_d.get(x, (-1, 0))[1]: best_d[x] = (y, c)
        if c > best_p.get(y, (-1, 0))[1]: best_p[y] = (x, c)
    for x, (y, c) in best_d.items():
        if best_p.get(y, (-1, 0))[0] == x: desc[(k, x)] = (k2, y); prog[(k2, y)] = (k, x)
# ---- lineages
rows = []; names = ['t_start', 't_end', 'life', 'M_max', 't_Mmax', 'Mstar', 'eps_int', 'eps_onset', 't_onset', 'M_on', 'n_mean_on', 'n_max_on', 'alpha_on', 'sig_on', 'rh_on', 'Sigma_on', 'n_snaps', 'sf_snaps', 'B_on', 'cs_on', 'vA_on', 'sigs_on', 'T_on', 'k_on', 'x_on', 'y_on', 'z_on']
seen = set()
for k in snaps:
    for c in range(D[k]['ncl']):
        if (k, c) in prog or (k, c) in seen: continue
        chain = [(k, c)]
        while chain[-1] in desc: chain.append(desc[chain[-1]])
        seen.update(chain)
        Ms = np.array([D[kk]['R'][cc, 2] for kk, cc in chain]); ts = np.array([D[kk]['t'] for kk, cc in chain])
        ms1 = np.array([D[kk]['ms1'][cc] for kk, cc in chain]); kl, cl_ = chain[-1]
        Mstar = ms1[:-1].sum() + D[kl]['ms10'][cl_]          # intervals along the chain + 10 Myr tail from the last snapshot
        sf = np.flatnonzero(ms1 > 0); i_on = (sf[0] if len(sf) else -1)
        i_state = max(i_on - 1, 0) if i_on >= 0 else -1
        if i_state >= 0:
            kk, cc = chain[i_state]; r = D[kk]['R'][cc]; Sig = 0.5 * r[2] / (np.pi * r[11] ** 2)
            on = (ts[i_on], r[2], r[12], r[13], r[14], r[10], r[11], Sig); ext = (tuple(r[17:22]) if len(r) >= 22 else (np.nan,) * 5) + (kk, r[4], r[5], r[6])
        else: on = (np.nan,) * 8; ext = (np.nan,) * 9
        rows.append((ts[0], ts[-1], ts[-1] - ts[0], Ms.max(), ts[np.argmax(Ms)], Mstar, Mstar / Ms.max(), (Mstar / on[1] if i_state >= 0 else np.nan)) + on + (len(chain), len(sf)) + ext)
L = np.array(rows); np.savez(f'{OUT}/lineages_{cat}.npz', rows=L, names=np.array(names)); col = {n: i for i, n in enumerate(names)}
print(f'{len(L)} lineages; with star formation: {np.sum(L[:, col["Mstar"]] > 0)}; lifetime median {np.median(L[:, col["life"]]):.1f} Myr (SF lineages {np.median(L[L[:, col["Mstar"]] > 0, col["life"]]):.1f})')
# ---- analysis
g = (L[:, col['Mstar']] > 0) & (L[:, col['M_max']] >= 300) & (L[:, col['t_start']] > 5) & (L[:, col['t_end']] < 221)   # complete lineages
print('\nlifetime-integrated efficiency eps_int = M_star,total / M_max per phase (by onset time), SF lineages with M_max >= 300, complete in time')
print('phase                  n    life med [Myr]  M_max med   eps_int med   16-84 %          sig(log eps_int)  sig(log eps_10 for ref)  sig(log Mstar)  sig(log Mmax)  corr(logMmax, log eps)')
E10 = {'before 1st passage': 0.70, 'between passages': 0.81, 'coalescence': 0.85, 'remnant': 0.78}
for t0, t1, lab in EP:
    m = g & (L[:, col['t_onset']] >= t0) & (L[:, col['t_onset']] < t1)
    if m.sum() < 10: print(lab, 'too few', m.sum()); continue
    e = L[m, col['eps_int']]; le = np.log10(e); lm = np.log10(L[m, col['M_max']]); ls = np.log10(L[m, col['Mstar']])
    print(f'{lab:20s} {m.sum():5d}   {np.median(L[m, col["life"]]):6.1f}        {np.median(L[m, col["M_max"]]):8.0f}   {np.median(e):.4f}     {np.percentile(e,16):.4f}-{np.percentile(e,84):.4f}     {np.std(le):.2f}              {E10.get(lab, float('nan')):.2f}                  {np.std(ls):.2f}           {np.std(lm):.2f}          {np.corrcoef(lm, le)[0,1]:+.2f}')
print('\nwhat sets eps_int: Spearman with the PRE-ONSET state (last snapshot before the first stars form)')
print('phase                  n   rho(n_mean)  rho(n_max)  rho(alpha_vir)  rho(sigma_3d)  rho(Sigma)  rho(M)  rho(r_h)  rho(life)')
for t0, t1, lab in EP:
    m = g & (L[:, col['t_onset']] >= t0) & (L[:, col['t_onset']] < t1) & np.isfinite(L[:, col['n_mean_on']])
    if m.sum() < 10: continue
    le = np.log10(L[m, col['eps_int']]); r = lambda key: spearmanr(le, np.log10(L[m, col[key]]))[0]
    print(f'{lab:20s} {m.sum():5d}    {r("n_mean_on"):+.2f}        {r("n_max_on"):+.2f}        {r("alpha_on"):+.2f}          {r("sig_on"):+.2f}         {r("Sigma_on"):+.2f}      {r("M_on"):+.2f}    {r("rh_on"):+.2f}     {r("life"):+.2f}')
print('\nintegrated star mass per lineage (M_star > 500): running index a(>500/1000/3000) per phase, vs bursts 1.65/1.56/1.54/1.39 at >500')
def ri(M, q):
    M = M[M >= q]; n = len(M); return 1 + n / np.sum(np.log(M / q)) if n >= 5 else np.nan
for t0, t1, lab in EP:
    m = g & (L[:, col['t_onset']] >= t0) & (L[:, col['t_onset']] < t1) & (L[:, col['Mstar']] > 500); Ms = L[m, col['Mstar']]
    if m.sum() < 5: continue
    print(f'{lab:20s} n={m.sum():4d}  a = {ri(Ms,500):.2f} / {ri(Ms,1000):.2f} / {ri(Ms,3000):.2f}   M_star max {Ms.max():.2e}   sig(log Mstar) {np.std(np.log10(Ms)):.2f}')
