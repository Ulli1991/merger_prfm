"""Shared helpers for the dwarf-merger movie + PRFM analysis.

Units of the GIZMO run: length kpc, mass 1e10 Msun, velocity km/s, B in Gauss.
Time unit = kpc / (km/s) = 0.9778 Gyr.
"""
import os, numpy as np, h5py
import numba as nb

SNAPDIR = '/ptmp/agurman/snapshots/merger_fid'
DATADIR = '/ptmp/uli/dwarf_merger'
PRFM_DIR = os.environ.get('PRFM_DIR', f'{DATADIR}/prfm')   # patch_XXX*.h5 location; PRFM_DIR=/ptmp/uli/dwarf_merger/prfm_one = single frame after the 2nd passage
NSNAP   = 232
SEP_MERGED = 0.5   # kpc; galaxy separation below which the system is treated as merged (one frame, no intruder cut)
T_MERGED = 106.0    # Myr; first coalescence (sep < SEP_MERGED at snap 109); after this the gas is mixed and the intruder cut is void
# merger phases from the nuclear separation (pericentres at 40.1, 109.5, 169.2, 217.1 Myr; apocentres 73, 140, 190)
PERICENTRES = [40.1, 109.5, 169.2, 217.1]
PHASES = [(25, 40, 'before 1st passage'), (40, 110, '1st to 2nd passage'), (110, 169, '2nd to 3rd passage'), (169, 232, '3rd passage to coalescence')]
def is_merged(sep, t):
    """patches/clusters for which the intruder-fraction cut does not apply"""
    return (np.asarray(sep) < SEP_MERGED) | (np.asarray(t) >= T_MERGED)
SNAP_DT = 0.001 * (3.085678e21 / 1e5) / 3.15576e13   # Myr between snapshots (0.9778)

# ---- units -------------------------------------------------------------
UNIT_L_CM  = 3.085678e21
UNIT_M_G   = 1.989e43
UNIT_V_CMS = 1e5
UNIT_T_S   = UNIT_L_CM / UNIT_V_CMS
MYR        = UNIT_T_S / 3.15576e13          # code time -> Myr  (=977.8)
MSUN       = 1e10                           # code mass -> Msun
PROTONMASS = 1.6726e-24
BOLTZMANN  = 1.3807e-16
GAMMA      = 5. / 3.
XH         = 0.76
G_CGS      = 6.674e-8
PC_CM      = 3.085678e18
KB_KMS     = BOLTZMANN / (1e5**2)           # not used directly

# fixed movie / analysis centre (system barycentre drifts < 0.2 kpc)
CENTER = np.array([149.95, 149.90, 149.75], dtype=np.float64)

def snap_path(i):
    return f'{SNAPDIR}/snap_{i:03d}.hdf5'

def snap_time_myr(i):
    with h5py.File(snap_path(i), 'r') as f:
        return float(f['Header'].attrs['Time']) * MYR

def read(i, ptype, fields):
    """Read a list of datasets for one particle type; returns dict."""
    out = {}
    with h5py.File(snap_path(i), 'r') as f:
        g = f[f'PartType{ptype}']
        for k in fields:
            out[k] = g[k][:]
        out['Time'] = float(f['Header'].attrs['Time'])
    return out

def temperature(u, ne):
    """u in code units (km/s)^2, ne = electron abundance per H -> T in K."""
    mu = 4.0 / (1.0 + 3.0 * XH + 4.0 * XH * ne)
    return (GAMMA - 1.0) * u * UNIT_V_CMS**2 * mu * PROTONMASS / BOLTZMANN

def nH_cgs(rho_code):
    return rho_code * UNIT_M_G / UNIT_L_CM**3 * XH / PROTONMASS

# ---- nested grid definition (movie) ------------------------------------
# (half-size in kpc, n per side).  All centred on CENTER.
LEVELS = [(12.0, 512), (4.0, 512), (1.5, 512)]
VEL_LEVEL = (12.0, 192)            # mass-weighted velocity grid for advection

# ---- numba kernel deposition -------------------------------------------
@nb.njit(inline='always')
def _w_cubic(q):
    # cubic spline, q = r/h with h = full support radius (GIZMO convention)
    if q >= 1.0:
        return 0.0
    if q < 0.5:
        return 1.0 - 6.0 * q * q + 6.0 * q * q * q
    t = 1.0 - q
    return 2.0 * t * t * t

@nb.njit(parallel=True, fastmath=True)
def deposit_kernel(pos, mass, hsml, cen, half, n, rmax_vox, out):
    """Mass-conserving, particle-normalised kernel deposition onto n^3 grid.

    pos:  (N,3) float64 positions; mass: (N,) ; hsml: (N,) support radius.
    Grid spans cen +- half, voxel size dx = 2*half/n.  Particles with
    h < dx are deposited with CIC.  Kernel radius clipped to rmax_vox voxels.
    Work is split over threads by z-slab so no atomics are needed.
    """
    dx = 2.0 * half / n
    nth = nb.get_num_threads()
    N = pos.shape[0]
    for th in nb.prange(nth):
        z0 = (n * th) // nth
        z1 = (n * (th + 1)) // nth
        for p in range(N):
            m = mass[p]
            if m <= 0.0:
                continue
            x = (pos[p, 0] - cen[0] + half) / dx
            y = (pos[p, 1] - cen[1] + half) / dx
            z = (pos[p, 2] - cen[2] + half) / dx
            h = hsml[p] / dx
            if h > rmax_vox:
                h = rmax_vox
            if h < 1.0:
                # CIC
                xc = x - 0.5; yc = y - 0.5; zc = z - 0.5
                i0 = int(np.floor(xc)); j0 = int(np.floor(yc)); k0 = int(np.floor(zc))
                fx = xc - i0; fy = yc - j0; fz = zc - k0
                for dk in range(2):
                    k = k0 + dk
                    if k < z0 or k >= z1:
                        continue
                    wz = fz if dk == 1 else 1.0 - fz
                    for di in range(2):
                        i = i0 + di
                        if i < 0 or i >= n:
                            continue
                        wx = fx if di == 1 else 1.0 - fx
                        for dj in range(2):
                            j = j0 + dj
                            if j < 0 or j >= n:
                                continue
                            wy = fy if dj == 1 else 1.0 - fy
                            out[i, j, k] += m * wx * wy * wz
                continue
            # kernel: bounding box in voxel indices
            ia = int(np.floor(x - h)); ib = int(np.floor(x + h))
            ja = int(np.floor(y - h)); jb = int(np.floor(y + h))
            ka = int(np.floor(z - h)); kb = int(np.floor(z + h))
            if ib < 0 or ia >= n or jb < 0 or ja >= n or kb < z0 or ka >= z1:
                # still need to know whether particle touches this slab at all
                if kb < z0 or ka >= z1:
                    continue
                if ib < 0 or ia >= n or jb < 0 or ja >= n:
                    continue
            # normalisation over the full kernel (all z), so mass is conserved
            wsum = 0.0
            for i in range(max(ia, 0), min(ib, n - 1) + 1):
                ddx = (i + 0.5 - x)
                for j in range(max(ja, 0), min(jb, n - 1) + 1):
                    ddy = (j + 0.5 - y)
                    for k in range(max(ka, 0), min(kb, n - 1) + 1):
                        ddz = (k + 0.5 - z)
                        q = np.sqrt(ddx * ddx + ddy * ddy + ddz * ddz) / h
                        wsum += _w_cubic(q)
            if wsum <= 0.0:
                # kernel smaller than voxel centre spacing: NGP
                i = int(x); j = int(y); k = int(z)
                if 0 <= i < n and 0 <= j < n and z0 <= k < z1:
                    out[i, j, k] += m
                continue
            inv = m / wsum
            for i in range(max(ia, 0), min(ib, n - 1) + 1):
                ddx = (i + 0.5 - x)
                for j in range(max(ja, 0), min(jb, n - 1) + 1):
                    ddy = (j + 0.5 - y)
                    for k in range(max(ka, z0), min(kb, z1 - 1) + 1):
                        ddz = (k + 0.5 - z)
                        q = np.sqrt(ddx * ddx + ddy * ddy + ddz * ddz) / h
                        w = _w_cubic(q)
                        if w > 0.0:
                            out[i, j, k] += inv * w


def deposit_ngp(pos, w, cen, half, n):
    """Fast NGP deposition with bincount (for stars / velocity moments)."""
    dx = 2.0 * half / n
    idx = np.floor((pos - cen + half) / dx).astype(np.int64)
    ok = np.all((idx >= 0) & (idx < n), axis=1)
    idx = idx[ok]
    flat = (idx[:, 0] * n + idx[:, 1]) * n + idx[:, 2]
    g = np.bincount(flat, weights=w[ok], minlength=n**3)
    return g.reshape(n, n, n), ok


def select_box(pos, cen, half, pad=0.0):
    return np.all(np.abs(pos - cen) < half + pad, axis=1)
