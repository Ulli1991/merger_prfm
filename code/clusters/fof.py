"""Friends-of-friends via linked cells with a SPARSE cell table (only occupied cells are stored, so the linking length can be
much smaller than the particle extent).  Cell size = linking length; neighbours are the 27 surrounding cells, found by binary
search in the sorted list of occupied cell keys.  Returns group labels 0..ngroups-1 (singletons included)."""
import numpy as np
import numba as nb

@nb.njit(cache=True)
def _find(parent, i):
    while parent[i] != i:
        parent[i] = parent[parent[i]]
        i = parent[i]
    return i

@nb.njit(cache=True)
def _fof_sparse(pos, link, order, keys_sorted, ukeys, ustart, uend, ny, nz):
    n = pos.shape[0]
    parent = np.arange(n)
    l2 = link * link
    nu = ukeys.shape[0]
    for a in range(n):
        p = order[a]
        key = keys_sorted[a]
        cz = key % nz; cy = (key // nz) % ny; cx = key // (nz * ny)
        for dx in range(-1, 2):
            for dy in range(-1, 2):
                for dz in range(-1, 2):
                    ix = cx + dx; iy = cy + dy; iz = cz + dz
                    if ix < 0 or iy < 0 or iz < 0 or iy >= ny or iz >= nz: continue
                    nkey = (ix * ny + iy) * nz + iz
                    # binary search
                    lo = 0; hi = nu
                    while lo < hi:
                        mid = (lo + hi) // 2
                        if ukeys[mid] < nkey: lo = mid + 1
                        else: hi = mid
                    if lo >= nu or ukeys[lo] != nkey: continue
                    for b in range(ustart[lo], uend[lo]):
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

def fof(pos, link):
    """group labels via friends-of-friends with linking length `link` (same units as pos)."""
    pos = np.ascontiguousarray(pos, dtype=np.float64)
    lo = pos.min(0) - link
    ncell = ((pos.max(0) + link - lo) / link).astype(np.int64) + 2
    idx = ((pos - lo) / link).astype(np.int64)
    key = (idx[:, 0] * ncell[1] + idx[:, 1]) * ncell[2] + idx[:, 2]
    order = np.argsort(key, kind='stable'); ks = key[order]
    ukeys, ustart = np.unique(ks, return_index=True); uend = np.append(ustart[1:], len(ks))
    lab = _fof_sparse(pos, float(link), order, ks, ukeys, ustart.astype(np.int64), uend.astype(np.int64), int(ncell[1]), int(ncell[2]))
    return np.unique(lab, return_inverse=True)[1]
