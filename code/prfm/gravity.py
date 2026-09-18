"""Isolated-boundary FFT gravity solver (two-level PM) to get the gravitational
acceleration of every gas particle, split by source (gas / stars / DM).

Fine level:  n_f^3 over +-half_f kpc, gas + stars only (DM is too coarse to
             deposit at this resolution; its field is smooth on these scales).
Coarse level: n_c^3 over +-half_c kpc, all components.
g = g_c(all) + [g_f(inner gas+stars) - g_c(inner gas+stars)]
"""
import numpy as np, time
from scipy import fft as sfft
from scipy.ndimage import map_coordinates
import numba as nb

G_KPC = 4.3009e-6          # G in kpc (km/s)^2 / Msun

def _green(n, dx, soft):
    """Green's function on the 2n padded grid (periodic distance), -G/sqrt(r^2+eps^2)."""
    m = 2 * n
    ax = np.arange(m) * dx
    ax = np.where(ax > n * dx, ax - m * dx, ax)          # wrap: distances in [-n dx, n dx)
    x2 = (ax**2)[:, None, None]; y2 = (ax**2)[None, :, None]; z2 = (ax**2)[None, None, :]
    r2 = x2 + y2 + z2
    return -G_KPC / np.sqrt(r2 + soft**2)

class PMLevel:
    def __init__(self, cen, half, n, soft=None, workers=-1):
        self.cen = np.asarray(cen, float); self.half = half; self.n = n
        self.dx = 2 * half / n
        self.soft = soft if soft is not None else 1.5 * self.dx
        self.workers = workers
        t0 = time.time()
        g = _green(n, self.dx, self.soft)
        self.Gk = sfft.rfftn(g, workers=workers)
        del g
        print(f'  PM level n={n} dx={self.dx*1e3:.1f} pc soft={self.soft*1e3:.0f} pc: green {time.time()-t0:.0f}s', flush=True)

    def deposit(self, pos, mass):
        """CIC deposition of masses -> (n,n,n) mass grid; returns grid and inside mask."""
        n = self.n
        x = (pos - self.cen + self.half) / self.dx - 0.5
        inside = np.all((x > 0) & (x < n - 1), axis=1)
        x = x[inside]; w = mass[inside]
        i0 = np.floor(x).astype(np.int64); f = x - i0
        grid = np.zeros(n**3)
        for dx_ in (0, 1):
            wx = f[:, 0] if dx_ else 1 - f[:, 0]
            for dy_ in (0, 1):
                wy = f[:, 1] if dy_ else 1 - f[:, 1]
                for dz_ in (0, 1):
                    wz = f[:, 2] if dz_ else 1 - f[:, 2]
                    idx = ((i0[:, 0] + dx_) * n + (i0[:, 1] + dy_)) * n + (i0[:, 2] + dz_)
                    grid += np.bincount(idx, weights=w * wx * wy * wz, minlength=n**3)
        return grid.reshape(n, n, n), inside

    def potential(self, massgrid):
        n = self.n; m = 2 * n
        pad = np.zeros((m, m, m))
        pad[:n, :n, :n] = massgrid
        pk = sfft.rfftn(pad, workers=self.workers)
        pk *= self.Gk
        phi = sfft.irfftn(pk, s=(m, m, m), workers=self.workers)[:n, :n, :n]
        return np.ascontiguousarray(phi)

    def accel_grid(self, massgrid):
        phi = self.potential(massgrid)
        gx, gy, gz = np.gradient(phi, self.dx)
        return -np.stack([gx, gy, gz]).astype(np.float32)

    def accel_at(self, acc, pos):
        """trilinear interpolation of the acceleration grid at positions (N,3)."""
        x = ((pos - self.cen + self.half) / self.dx - 0.5).T
        out = np.empty((pos.shape[0], 3), np.float32)
        for c in range(3):
            out[:, c] = map_coordinates(acc[c], x, order=1, mode='nearest')
        return out


def two_level_accel(pos_eval, comps, cen, half_f, n_f, half_c, n_c, workers=-1):
    """comps: dict name -> (pos, mass). Returns dict name -> (N_eval,3) accel [ (km/s)^2/kpc ]
    for each component plus 'total'.  DM is coarse-only."""
    t0 = time.time()
    Lf = PMLevel(cen, half_f, n_f, workers=workers)
    Lc = PMLevel(cen, half_c, n_c, workers=workers)
    out = {}
    tot = np.zeros((pos_eval.shape[0], 3), np.float32)
    for name, (p, m) in comps.items():
        gc_grid_all, _ = Lc.deposit(p, m)
        acc = Lc.accel_grid(gc_grid_all)
        a = Lf.accel_at(acc, pos_eval) if False else Lc.accel_at(acc, pos_eval)
        if name != 'dm':
            # fine correction: inner particles only
            gf, inside = Lf.deposit(p, m)
            accf = Lf.accel_grid(gf)
            gci, _ = Lc.deposit(p[inside], m[inside])
            accci = Lc.accel_grid(gci)
            # only apply where the evaluation point lies in the fine box (with margin)
            infine = np.all(np.abs(pos_eval - cen) < half_f - 2 * Lf.dx, axis=1)
            corr = Lf.accel_at(accf, pos_eval[infine]) - Lc.accel_at(accci, pos_eval[infine])
            a[infine] += corr
        out[name] = a
        tot += a
        print(f'  gravity {name}: done ({time.time()-t0:.0f}s)', flush=True)
    out['total'] = tot
    return out


@nb.njit(parallel=True, fastmath=True)
def direct_accel(pts, pos, mass, soft):
    """Brute-force reference (for validation on a few points)."""
    out = np.zeros((pts.shape[0], 3))
    for i in nb.prange(pts.shape[0]):
        ax = 0.0; ay = 0.0; az = 0.0
        for j in range(pos.shape[0]):
            dx = pos[j, 0] - pts[i, 0]; dy = pos[j, 1] - pts[i, 1]; dz = pos[j, 2] - pts[i, 2]
            r2 = dx * dx + dy * dy + dz * dz + soft * soft
            inv = mass[j] / (r2 * np.sqrt(r2))
            ax += dx * inv; ay += dy * inv; az += dz * inv
        out[i, 0] = G_KPC * ax; out[i, 1] = G_KPC * ay; out[i, 2] = G_KPC * az
    return out
