"""Star cluster mass function and its slope.

For every STEP-th snapshot, take stars younger than AGE_MAX Myr (each star
particle = 4 Msun), link them with a friends-of-friends of linking length
LINK pc, keep groups of at least NMIN members and record mass, position,
half-mass radius and mean age.  The slope alpha of dN/dM ~ M^-alpha is the
maximum-likelihood power-law estimate above M_MIN (Clauset+09):
    alpha = 1 + N / sum ln(M_i / M_MIN),  err = (alpha-1)/sqrt(N)

usage: python cluster_mf.py [--step 10] [--age 10] [--link 5] [--nmin 25] [--mmin 300]
out:   /ptmp/uli/dwarf_merger/clusters/clusters_<tag>.npz  and figs/
"""
import sys, os, glob, argparse, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import numba as nb
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import common as C

ap = argparse.ArgumentParser()
ap.add_argument('--step', type=int, default=1)      # catalogue every STEP-th snapshot
ap.add_argument('--fitstep', type=int, default=10)  # MF fits use snapshots at this spacing (independent populations)
ap.add_argument('--age', type=float, default=10.0)     # Myr, upper age
ap.add_argument('--agemin', type=float, default=0.0)   # Myr, lower age
ap.add_argument('--link', type=float, default=5.0)     # pc
ap.add_argument('--nmin', type=int, default=25)
ap.add_argument('--mmin', type=float, default=300.0)   # Msun, fit lower limit
ap.add_argument('--snaps', type=int, nargs='*', default=None)   # restrict to these snapshots (tests)
ap.add_argument('--suffix', default='')                          # output tag suffix (tests)
a = ap.parse_args()
tag = f'age{a.agemin:g}-{a.age:g}_l{a.link:g}_n{a.nmin}' + a.suffix
NB_MAX = 40000      # groups larger than this get the potential energy from a 4000-member subsample against all members (exact pairwise below)
NSUB = 4000
SOFT_PC = 1.0       # Plummer softening for the pairwise potential
G_PC = 4.30091e-3   # pc (km/s)^2 / Msun
OUT = f'{C.DATADIR}/clusters'; FIG = f'{OUT}/figs'
os.makedirs(FIG, exist_ok=True)

@nb.njit(cache=True)
def _find(parent, i):
    while parent[i] != i:
        parent[i] = parent[parent[i]]
        i = parent[i]
    return i

@nb.njit(cache=True)
def _fof_cells(pos, link, order, cell_start, cell_of, ncell):
    """union-find over a linked-cell list; memory O(N), no pair list."""
    n = pos.shape[0]
    parent = np.arange(n)
    l2 = link * link
    nx, ny, nz = ncell
    for a in range(n):
        p = order[a]
        c = cell_of[p]
        cz = c % nz; cy = (c // nz) % ny; cx = c // (nz * ny)
        for dx in range(-1, 2):
            ix = cx + dx
            if ix < 0 or ix >= nx: continue
            for dy in range(-1, 2):
                iy = cy + dy
                if iy < 0 or iy >= ny: continue
                for dz in range(-1, 2):
                    iz = cz + dz
                    if iz < 0 or iz >= nz: continue
                    cc = (ix * ny + iy) * nz + iz
                    for b in range(cell_start[cc], cell_start[cc + 1]):
                        q = order[b]
                        if q <= p: continue
                        d0 = pos[p, 0] - pos[q, 0]; d1 = pos[p, 1] - pos[q, 1]; d2 = pos[p, 2] - pos[q, 2]
                        if d0 * d0 + d1 * d1 + d2 * d2 < l2:
                            rp = _find(parent, p); rq = _find(parent, q)
                            if rp != rq:
                                parent[rq] = rp
    lab = np.empty(n, np.int64)
    for i in range(n):
        lab[i] = _find(parent, i)
    return lab

@nb.njit(parallel=True, fastmath=True, cache=True)
def _pot_energy(pos_pc, m, soft):
    """total gravitational PE in Msun (km/s)^2 (Plummer-softened pairwise)."""
    n = pos_pc.shape[0]; s2 = soft * soft
    tot = 0.0
    for i in nb.prange(n):
        acc = 0.0
        for j in range(i + 1, n):
            dx = pos_pc[i, 0] - pos_pc[j, 0]; dy = pos_pc[i, 1] - pos_pc[j, 1]; dz = pos_pc[i, 2] - pos_pc[j, 2]
            acc += m[i] * m[j] / np.sqrt(dx * dx + dy * dy + dz * dz + s2)
        tot += acc
    return -G_PC * tot

@nb.njit(parallel=True, fastmath=True, cache=True)
def _pot_sub(pos_pc, m, idx, soft):
    """potential at the subsample members idx from all members (Plummer-softened); PE = 0.5 * M_tot * <phi>_subsample."""
    out = np.zeros(len(idx)); s2 = soft * soft
    for a_ in nb.prange(len(idx)):
        i = idx[a_]; acc = 0.0
        for j in range(pos_pc.shape[0]):
            if j == i: continue
            dx = pos_pc[i, 0] - pos_pc[j, 0]; dy = pos_pc[i, 1] - pos_pc[j, 1]; dz = pos_pc[i, 2] - pos_pc[j, 2]
            acc += m[j] / np.sqrt(dx * dx + dy * dy + dz * dz + s2)
        out[a_] = -G_PC * acc
    return out

@nb.njit(parallel=True, fastmath=True, cache=True)
def _phi_all(pos_pc, m, soft):
    """softened potential at every member from all other members [ (km/s)^2 ]."""
    n = pos_pc.shape[0]; s2 = soft * soft; out = np.zeros(n)
    for i in nb.prange(n):
        acc = 0.0
        for j in range(n):
            if j == i: continue
            dx = pos_pc[i, 0] - pos_pc[j, 0]; dy = pos_pc[i, 1] - pos_pc[j, 1]; dz = pos_pc[i, 2] - pos_pc[j, 2]
            acc += m[j] / np.sqrt(dx * dx + dy * dy + dz * dz + s2)
        out[i] = -G_PC * acc
    return out

@nb.njit(parallel=True, fastmath=True, cache=True)
def _phi_correct(pos_pc, m, phi, rem, soft):
    """subtract from phi the contribution of the removed members rem (indices)."""
    s2 = soft * soft; n = pos_pc.shape[0]
    for i in nb.prange(n):
        acc = 0.0
        for a_ in range(len(rem)):
            j = rem[a_]
            if j == i: continue
            dx = pos_pc[i, 0] - pos_pc[j, 0]; dy = pos_pc[i, 1] - pos_pc[j, 1]; dz = pos_pc[i, 2] - pos_pc[j, 2]
            acc += m[j] / np.sqrt(dx * dx + dy * dy + dz * dz + s2)
        phi[i] += G_PC * acc
    return phi

def self_bound(pos_pc, vel, m, soft, max_iter=50):
    """Subfind-style iterative unbinding: remove members with 0.5 v^2 + phi > 0 (v relative to the mean velocity of the members
    still bound, phi from the members still bound) until no member is removed or fewer than 2 remain.
    Returns the boolean mask of the self-bound remnant."""
    keep = np.ones(len(m), bool); phi = _phi_all(pos_pc, m, soft)
    for _ in range(max_iter):
        if keep.sum() < 2: break
        vc = np.average(vel[keep], axis=0, weights=m[keep])
        e = 0.5 * np.sum((vel - vc) ** 2, 1) + phi
        rem = np.flatnonzero(keep & (e > 0))
        if len(rem) == 0: break
        # remove at most the worst 20 % per iteration (as Subfind does) to avoid overshooting
        if len(rem) > 0.2 * keep.sum(): rem = rem[np.argsort(e[rem])[-int(0.2 * keep.sum()) - 1:]]
        keep[rem] = False; phi = _phi_correct(pos_pc, m, phi, rem, soft)
    return keep

def fof(pos, link):
    """group labels via friends-of-friends with linking length `link`."""
    lo = pos.min(0) - link
    ncell = np.maximum(((pos.max(0) + link - lo) / link).astype(np.int64) + 1, 1)
    ncell = np.minimum(ncell, 2048)
    csize = (pos.max(0) + link - lo) / ncell
    idx = np.minimum(((pos - lo) / csize).astype(np.int64), ncell - 1)
    cell_of = (idx[:, 0] * ncell[1] + idx[:, 1]) * ncell[2] + idx[:, 2]
    order = np.argsort(cell_of, kind='stable')
    cnt = np.bincount(cell_of, minlength=int(np.prod(ncell)))
    cell_start = np.concatenate([[0], np.cumsum(cnt)])
    lab = _fof_cells(pos, link, order, cell_start, cell_of, ncell)
    return np.unique(lab, return_inverse=True)[1]

def find_clusters(k):
    d = np.load(f'{C.DATADIR}/stars/stars_{k:03d}.npz')
    t = float(d['time_myr']); age = t - d['tform_myr']
    # IC disk stars carry tform=0 exactly; only stars formed in the run count
    sel = (age >= a.agemin) & (age < a.age) & (d['tform_myr'] > 0)
    pos = d['pos'][sel].astype(np.float64); m = d['mass'][sel].astype(np.float64) * C.MSUN
    ag = age[sel]
    if sel.sum() < a.nmin:
        return t, []
    lab = fof(pos, a.link * 1e-3)
    cnt = np.bincount(lab)
    good = np.where(cnt >= a.nmin)[0]
    vel = d['vel'][sel].astype(np.float64)
    out = []
    for g in good:
        s_ = lab == g
        mm = m[s_]; pp = pos[s_]; vv = vel[s_]; n_ = int(s_.sum())
        com = np.average(pp, axis=0, weights=mm); vcom = np.average(vv, axis=0, weights=mm)
        r = np.linalg.norm(pp - com, axis=1)
        rh = np.sqrt(np.interp(0.5, np.cumsum(np.sort(mm)) / mm.sum(), np.sort(r) ** 2)) if len(r) > 1 else 0.0
        ke = 0.5 * np.sum(mm * np.sum((vv - vcom) ** 2, axis=1))
        if n_ <= NB_MAX:
            pe = _pot_energy(np.ascontiguousarray((pp - com) * 1e3), mm, SOFT_PC)
        else:
            idx = np.random.default_rng(int(g)).choice(n_, NSUB, replace=False)
            pe = 0.5 * mm.sum() * _pot_sub(np.ascontiguousarray((pp - com) * 1e3), mm, idx, SOFT_PC).mean()
        nucleus = 0.0      # every group is energy-tested and kept (the former N > 40000 / r_h > 20 pc exclusion is gone); r_h and N are stored
        # self-bound remnant (iterative unbinding); for very large groups the per-member potential costs N^2 but there are few of them
        sb = self_bound(np.ascontiguousarray((pp - com) * 1e3), vv, mm, SOFT_PC)
        if sb.sum() >= 2:
            r_sb = np.linalg.norm(pp[sb] - np.average(pp[sb], axis=0, weights=mm[sb]), axis=1)
            rh_sb = np.sqrt(np.interp(0.5, np.cumsum(np.sort(mm[sb])) / mm[sb].sum(), np.sort(r_sb) ** 2)) * 1e3
        else: rh_sb = 0.0
        out.append((k, t, mm.sum(), n_, com[0], com[1], com[2], vcom[0], vcom[1], vcom[2], rh * 1e3,
                    ag[s_].mean(), ag[s_].std(), ke, pe, float(ke + pe < 0) if np.isfinite(pe) else 0.0, nucleus, mm[sb].sum(), int(sb.sum()), rh_sb))
    return t, out

def mle_slope(M, mmin):
    M = M[M >= mmin]
    if len(M) < 5:
        return np.nan, np.nan, len(M)
    al = 1 + len(M) / np.sum(np.log(M / mmin))
    return al, (al - 1) / np.sqrt(len(M)), len(M)

snaps = a.snaps if a.snaps else list(range(0, C.NSNAP, a.step))
rows = []
t0 = time.time()
for k in snaps:
    t, cl = find_clusters(k)
    rows += cl
    if cl:
        ms = np.array([c[2] for c in cl]); bd = np.array([c[15] for c in cl]) > 0
        al, e, n = mle_slope(ms, a.mmin)
        print(f'snap {k:3d} t={t:6.1f} Myr: {len(cl):4d} groups, bound {bd.sum():4d} ({ms[bd].sum()/ms.sum():.2f} of mass), M_max={ms.max():8.0f}, '
              f'alpha(M>{a.mmin:g})={al:.2f}+-{e:.2f} (N={n})  [{time.time()-t0:.0f}s]', flush=True)
    else:
        print(f'snap {k:3d} t={t:6.1f} Myr: no clusters', flush=True)

R = np.array(rows, dtype=np.float64)
names = ['snap', 't', 'M', 'N', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'rh_pc', 'age_mean', 'age_std', 'KE', 'PE', 'bound', 'nucleus', 'M_sb', 'N_sb', 'rh_sb_pc']
np.savez(f'{OUT}/clusters_{tag}.npz', **{n: R[:, i] for i, n in enumerate(names)},
         step=a.step, age_max=a.age, age_min=a.agemin, link_pc=a.link, nmin=a.nmin)
# ---- mass-function fits: independent populations (every fitstep-th snapshot), all groups
fit = (np.mod(R[:, 0], a.fitstep) == 0) & (R[:, 16] == 0)
R = R[fit]
M, T = R[:, 2], R[:, 1]
BD = R[:, 15] > 0

# ---- global slope + bootstrap
al, e, n = mle_slope(M, a.mmin)
boot = [mle_slope(np.random.choice(M, len(M)), a.mmin)[0] for _ in range(500)]
alb, eb, nbd = mle_slope(M[BD], a.mmin)
print(f'\nALL snapshots (every {a.fitstep}th, all groups; {len(M)} groups, {n} above {a.mmin:g} Msun): '
      f'alpha = {al:.2f} +- {e:.2f} (bootstrap {np.nanstd(boot):.2f}), M_max = {M.max():.0f} Msun')
print(f'   bound groups only: alpha = {alb:.2f} +- {eb:.2f} (N={nbd}), bound mass fraction {M[BD].sum()/M.sum():.2f}, bound number fraction {BD.mean():.2f}')
for mm in [100, 200, 300, 500, 1000, 2000]:
    al2, e2, n2 = mle_slope(M, mm); print(f'   M_min={mm:5d}: alpha={al2:.2f}+-{e2:.2f} N={n2}')

# ---- slope in time windows
wins = [(0, 50), (50, 100), (100, 150), (150, 200), (200, 232)]
print('\nslope in time windows:')
wal = []
for lo, hi in wins:
    s = (T >= lo) & (T < hi)
    al2, e2, n2 = mle_slope(M[s], a.mmin); wal.append((lo, hi, al2, e2, n2))
    print(f'   t={lo:3d}-{hi:3d} Myr: alpha={al2:.2f}+-{e2:.2f} N={n2}, M_max={M[s].max() if s.any() else 0:.0f}')

# ---- figures
fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
edges = np.logspace(2, np.log10(M.max()) + 0.1, 25)
h, _ = np.histogram(M, edges); w = np.diff(edges); c = np.sqrt(edges[1:] * edges[:-1])
ok = h > 0
ax[0].errorbar(c[ok], h[ok] / w[ok], yerr=np.sqrt(h[ok]) / w[ok], fmt='o', ms=4, color='k', label='all snapshots')
xx = np.array([a.mmin, M.max()]); A = n * (al - 1) * a.mmin ** (al - 1)
hb, _ = np.histogram(M[BD], edges); okb = hb > 0
ax[0].plot(c[okb], hb[okb] / w[okb], 's', ms=3, mfc='none', color='tab:green', label=rf'bound only, $\alpha={alb:.2f}\pm{eb:.2f}$')
ax[0].plot(xx, A * xx ** -al, 'r-', label=rf'MLE $\alpha={al:.2f}\pm{e:.2f}$ (M>{a.mmin:g})')
ax[0].plot(xx, A * xx ** -2.0 * (xx[0] ** 2 / xx[0] ** al), 'b:', label=r'$\alpha=2$')
ax[0].set(xscale='log', yscale='log', xlabel=r'$M_{\rm cl}$ [M$_\odot$]', ylabel=r'$dN/dM$ [M$_\odot^{-1}$]',
          title=f'age {a.agemin:g}-{a.age:g} Myr, FoF {a.link:g} pc, N>={a.nmin}, no nucleus'); ax[0].legend(fontsize=8)
for lo, hi, al2, e2, n2 in wal:
    s = (T >= lo) & (T < hi)
    if s.sum() > 5:
        ms = np.sort(M[s])[::-1]
        ax[1].plot(ms, np.arange(1, len(ms) + 1), label=f'{lo}-{hi} Myr: α={al2:.2f}±{e2:.2f}')
ax[1].set(xscale='log', yscale='log', xlabel=r'$M_{\rm cl}$ [M$_\odot$]', ylabel=r'$N(>M)$'); ax[1].legend(fontsize=8)
sc = ax[2].scatter(T, M, c=R[:, 10], s=6, cmap='viridis', vmin=0, vmax=20 if a.agemin == 0 else 50)
plt.colorbar(sc, ax=ax[2], label=r'$r_h$ [pc]')
ax[2].set(yscale='log', xlabel='t [Myr]', ylabel=r'$M_{\rm cl}$ [M$_\odot$]')
fig.tight_layout(); fig.savefig(f'{FIG}/cluster_mf_{tag}.png', dpi=130)
print('figure:', f'{FIG}/cluster_mf_{tag}.png')
