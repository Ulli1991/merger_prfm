"""PRFM patch analysis for one snapshot.

For each galaxy frame (two before coalescence, one after) lay a grid of
PATCH x PATCH kpc columns in the disk plane and measure, per column:
  gas weight W = int rho g_n dz  (direct, from the PM gravity solve, split by
  gas / stars / DM, above and below the density-weighted midplane),
  midplane pressures P_th (neutral/ionised/hot), P_turb, P_mag,
  Sigma_gas, Sigma_*, rho_sd, H, sigma_z, Sigma_SFR (10, 40 Myr), SN rate,
  G0, phase fractions, velocity-gradient diagnostics.
Also the analytic P_DE for comparison.

usage: python patches.py SNAP [SNAP ...]
output: /ptmp/uli/dwarf_merger/prfm/patch_SSS.h5
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, h5py
from scipy.ndimage import map_coordinates
import common as C
from gravity import PMLevel, two_level_accel, direct_accel, G_KPC

PATCH = float(os.environ.get('PRFM_PATCH', 0.5))   # kpc, column side (PRFM_PATCH=0.25, 0.125 for the scale test)
RMAX = 3.0             # kpc, half-extent of patch grid
ZCOL = float(os.environ.get('PRFM_ZCOL', 1.5))   # kpc, column half-height (PRFM_ZCOL=0.5 with PRFM_ZREF=1 = layer weight)
SLAB_MIN = 0.025       # kpc, minimum slab half-thickness for midplane averages
SLAB_FRAC = 0.25       # slab half-thickness = max(SLAB_MIN, SLAB_FRAC * H)
FINE = (4.0, 512)      # gravity fine level: half-size kpc, n
COARSE = (32.0, 512)
T_HOT = 1.0e5
T_2P = 2.0e4          # OK22 'two-phase' gas: T < 2e4 K
SFR_WINDOWS = (10.0, 40.0)
# Msun (km/s)^2 kpc^-3  ->  K cm^-3
P_UNIT = 1.989e33 * 1e10 / (3.085678e21**3) / 1.3807e-16
# kpc^3 code-mass/code-density volume -> nothing; magnetic: erg/cm^3 / k_B
ID_SPLIT = 26_000_000  # PartType2 ids <= split belong to galaxy A
GAS_SPLIT = 14_500_000 # gas ids < split belong to galaxy A (97% pure at t=0)
NZB = 60               # vertical profile bins over +-ZCOL

def frame_axes(nhat):
    nhat = nhat / np.linalg.norm(nhat)
    a = np.array([1.0, 0, 0]) if abs(nhat[0]) < 0.9 else np.array([0, 1.0, 0])
    e1 = np.cross(nhat, a); e1 /= np.linalg.norm(e1)
    e2 = np.cross(nhat, e1)
    return e1, e2, nhat

ONEFRAME = float(os.environ.get('PRFM_ONEFRAME', 'inf'))   # Myr; from this time on a single frame 'M' regardless of the separation
def galaxy_frames(s2pos, s2ids, s2vel, gpos, gvel, gm, tmyr=0.0):
    """Return list of dict(cen, vcen, nhat, label).  Two frames (A, B: own centre, own normal) while the stellar centres are
    more than 0.5 kpc apart and t < ONEFRAME; otherwise one frame M centred midway between the two stellar centres, with the
    mean velocity of all stars within 3 kpc and the normal from the gas angular momentum within 3 kpc."""
    frames = []
    A = s2ids <= ID_SPLIT
    cA = np.median(s2pos[A], 0); cB = np.median(s2pos[~A], 0)
    sep = np.linalg.norm(cA - cB)
    def mk(cen, vc, label, rad):
        d = gpos - cen; r = np.linalg.norm(d, axis=1); sel = r < rad
        L = np.sum(gm[sel, None] * np.cross(d[sel], gvel[sel] - vc), 0)
        return dict(cen=cen, vcen=vc, nhat=L / np.linalg.norm(L), label=label, sep=sep)
    if tmyr >= ONEFRAME:
        c = 0.5 * (cA + cB); vc = s2vel[np.linalg.norm(s2pos - c, axis=1) < 3.0].mean(0)
        frames.append(mk(c, vc, 'M', 3.0)); return frames
    if sep > 0.5:
        for lab, sel, c in (('A', A, cA), ('B', ~A, cB)):
            vc = s2vel[sel][np.linalg.norm(s2pos[sel] - c, axis=1) < 1.5].mean(0)
            frames.append(mk(c, vc, lab, 2.0))
    else:
        c = 0.5 * (cA + cB); vc = s2vel[np.linalg.norm(s2pos - c, axis=1) < 1.5].mean(0)
        frames.append(mk(c, vc, 'M', 3.0))
    return frames

def layer_cut(prof, zmid, zcol, dzb, floor=0.03, rise=3.0):
    """For each patch, walk the vertical mass profile up and down from the midplane bin;
    stop where the mass drops below `floor` x peak or where it rises again by `rise` x
    above the running minimum (a second disk). Returns z_lo, z_hi (kpc, frame coords)."""
    npatch, nzb = prof.shape
    z_lo = np.full(npatch, -zcol); z_hi = np.full(npatch, zcol)
    for p in range(npatch):
        pk = prof[p].max()
        if pk <= 0:
            continue
        b0 = int(np.clip((zmid[p] + zcol) / dzb, 0, nzb - 1))
        for direction, store in ((1, 'hi'), (-1, 'lo')):
            runmin = prof[p, b0]; bmin = b0; b = b0; cut = None
            while 0 <= b + direction < nzb:
                b += direction
                m = prof[p, b]
                if m < runmin:
                    runmin = m; bmin = b
                if m < floor * pk:
                    cut = b; break
                if runmin < 0.5 * pk and m > rise * max(runmin, 1e-3 * pk):
                    cut = bmin; break
            if cut is not None:
                if store == 'hi':
                    z_hi[p] = (cut + 0.5) * dzb - zcol
                else:
                    z_lo[p] = (cut + 0.5) * dzb - zcol
    return z_lo, z_hi

def bsum(idx, w, n):
    return np.bincount(idx, weights=w, minlength=n)[:n]

def analyse(isnap):
    t0 = time.time()
    LOS = bool(os.environ.get('PRFM_LOS'))
    fout = f'{C.DATADIR}/prfm/patch_{isnap:03d}' + ('_los' if LOS else '') + os.environ.get('PRFM_TAG', '') + '.h5'
    if os.path.exists(fout):
        print('exists', fout); return
    g = C.read(isnap, 0, ['Coordinates', 'Velocities', 'Masses', 'Density', 'InternalEnergy', 'ElectronAbundance',
                          'ChemicalAbundancesSG', 'MagneticField', 'EnergyDensityFUV', 'VelocityGradient', 'ParticleIDs'])
    gid = g['ParticleIDs']
    tmyr = g['Time'] * C.MYR
    gpos = g['Coordinates'].astype(np.float64); gvel = g['Velocities'].astype(np.float64)
    gm = g['Masses'].astype(np.float64) * C.MSUN          # Msun
    rho = g['Density'].astype(np.float64)                  # code
    T = C.temperature(g['InternalEnergy'], g['ElectronAbundance'])
    u = g['InternalEnergy'].astype(np.float64)             # (km/s)^2
    xh2 = np.clip(g['ChemicalAbundancesSG'][:, 0], 0, 0.5).astype(np.float64) * 2   # mass fraction of H in H2
    xhii = np.clip(g['ChemicalAbundancesSG'][:, 1], 0, 1).astype(np.float64)
    B = g['MagneticField'].astype(np.float64)              # Gauss
    G0 = g['EnergyDensityFUV'].astype(np.float64)
    vg = g['VelocityGradient'].astype(np.float64).reshape(-1, 3, 3)   # km/s/kpc, d v_i / d x_j (assumed)
    vol = g['Masses'].astype(np.float64) / rho             # kpc^3 (code units)
    s2 = C.read(isnap, 2, ['Coordinates', 'Velocities', 'Masses', 'ParticleIDs'])
    s4 = C.read(isnap, 4, ['Coordinates', 'Masses', 'StellarFormationTime'])
    dm = C.read(isnap, 1, ['Coordinates', 'Masses'])
    s2pos = s2['Coordinates'].astype(np.float64); s2m = s2['Masses'].astype(np.float64) * C.MSUN
    s4pos = s4['Coordinates'].astype(np.float64); s4m = s4['Masses'].astype(np.float64) * C.MSUN
    age = tmyr - s4['StellarFormationTime'].astype(np.float64) * C.MYR
    dmpos = dm['Coordinates'].astype(np.float64); dmm = dm['Masses'].astype(np.float64) * C.MSUN
    print(f'snap {isnap} t={tmyr:.1f} Myr read ({time.time()-t0:.0f}s)', flush=True)

    # ---- gravity at gas particle positions, split by source
    spos = np.vstack([s2pos, s4pos]); sm = np.concatenate([s2m, s4m])
    comps = {'gas': (gpos, gm), 'star': (spos, sm), 'dm': (dmpos, dmm)}
    acc = two_level_accel(gpos, comps, C.CENTER, FINE[0], FINE[1], COARSE[0], COARSE[1])
    print(f'  gravity done ({time.time()-t0:.0f}s)', flush=True)
    # validation on 64 random gas particles inside the fine box (all particles, direct sum)
    rng = np.random.default_rng(isnap)
    infine = np.all(np.abs(gpos - C.CENTER) < 2.5, axis=1)
    chk = rng.choice(np.flatnonzero(infine), 64, replace=False)
    allpos = np.vstack([gpos, spos, dmpos]); allm = np.concatenate([gm, sm, dmm])
    aref = direct_accel(gpos[chk], allpos, allm, 0.02)
    err = np.linalg.norm(acc['total'][chk] - aref, axis=1) / np.linalg.norm(aref, axis=1)
    print(f'  gravity check: median rel err {np.median(err):.3f}, 90% {np.percentile(err,90):.3f} ({time.time()-t0:.0f}s)', flush=True)
    del allpos, allm
    # stellar + DM density grids for rho_sd
    Lf = PMLevel(C.CENTER, FINE[0], FINE[1]); Lc = PMLevel(C.CENTER, COARSE[0], COARSE[1])
    star_grid, _ = Lf.deposit(spos, sm); star_grid /= Lf.dx**3
    dm_grid, _ = Lc.deposit(dmpos, dmm); dm_grid /= Lc.dx**3
    def dens_at(p):
        xf = ((p - C.CENTER + FINE[0]) / Lf.dx - 0.5).T
        xc = ((p - C.CENTER + COARSE[0]) / Lc.dx - 0.5).T
        return (map_coordinates(star_grid, xf, order=1, mode='nearest'),
                map_coordinates(dm_grid, xc, order=1, mode='nearest'))

    # SN catalogue for SN rate per patch
    cat = np.load(f'{C.DATADIR}/stars/catalog.npz')
    snsel = (cat['tsn_myr'] > tmyr - 10) & (cat['tsn_myr'] <= tmyr)
    xsn = cat['xsn'][snsel]

    frames = galaxy_frames(s2pos, s2['ParticleIDs'], s2['Velocities'].astype(np.float64), gpos, gvel, gm, tmyr)
    if LOS:   # line-of-sight test: same centres, other normals (box axes and tilts about e1 towards e2)
        extra = []
        for fr in frames:
            e1, e2, n = frame_axes(fr['nhat'])
            for lab, nh in (('x', [1.0, 0, 0]), ('y', [0, 1.0, 0]), ('z', [0, 0, 1.0])):
                extra.append(dict(fr, nhat=np.array(nh), label=fr['label'] + '_' + lab))
            for th in (30, 60, 90):
                r = np.radians(th)
                extra.append(dict(fr, nhat=np.cos(r) * n + np.sin(r) * e2, label=f"{fr['label']}_t{th}"))
        frames += extra
        print(f'  LOS mode: {len(frames)} frames: ' + ' '.join(fr['label'] for fr in frames), flush=True)
    nb_ = int(2 * RMAX / PATCH); npatch = nb_ * nb_
    ZREF = {}; PATCH_REF = 0.5; NB_REF = int(2 * RMAX / PATCH_REF)   # PRFM_ZREF=1: centre each column on the gas midplane of the full-column run (same LOS mode, 0.5 kpc grid)
    if os.environ.get('PRFM_ZREF'):
        with h5py.File(f'{C.DATADIR}/prfm/patch_{isnap:03d}' + ('_los' if LOS else '') + os.environ.get('PRFM_REF_TAG', '') + '.h5', 'r') as fref:
            PATCH_REF = float(fref.attrs['patch_kpc']); NB_REF = int(fref.attrs['nb'])
            for gname in fref:
                if gname.startswith('frame_'): ZREF[gname.replace('frame_', '')] = fref[gname]['zmid'][:].ravel()
        print(f'  column reference midplanes loaded for frames {list(ZREF)} (reference grid {PATCH_REF} kpc, nb {NB_REF})', flush=True)
    with h5py.File(fout + '.tmp', 'w') as f:
        f.attrs['time_myr'] = tmyr; f.attrs['patch_kpc'] = PATCH; f.attrs['rmax_kpc'] = RMAX
        f.attrs['zcol_kpc'] = ZCOL; f.attrs['nb'] = nb_; f.attrs['P_unit'] = P_UNIT
        f.attrs['grav_check_median_relerr'] = float(np.median(err))
        for fr in frames:
            e1, e2, n = frame_axes(fr['nhat'])
            grp = f.create_group(f'frame_{fr["label"]}')
            grp.attrs['center'] = fr['cen']; grp.attrs['vcen'] = fr['vcen']; grp.attrs['nhat'] = n
            grp.attrs['e1'] = e1; grp.attrs['e2'] = e2; grp.attrs['separation_kpc'] = fr['sep']
            zref = ZREF.get(fr['label'], None)
            def to_frame(p):
                d = p - fr['cen']
                x, y, z = d @ e1, d @ e2, d @ n
                if zref is not None:            # z measured from the reference midplane of the 0.5 kpc reference column the point falls in
                    ix = np.clip(np.floor((x + RMAX) / PATCH_REF).astype(np.int64), 0, NB_REF - 1)
                    iy = np.clip(np.floor((y + RMAX) / PATCH_REF).astype(np.int64), 0, NB_REF - 1)
                    z = z - zref[ix * NB_REF + iy]
                return x, y, z
            def patch_index(x, y, z, zc=ZCOL):
                ix = np.floor((x + RMAX) / PATCH).astype(np.int64); iy = np.floor((y + RMAX) / PATCH).astype(np.int64)
                ok = (ix >= 0) & (ix < nb_) & (iy >= 0) & (iy < nb_) & (np.abs(z) < zc)
                return np.where(ok, ix * nb_ + iy, -1)
            A = PATCH**2
            # ---- gas
            gx, gy, gz = to_frame(gpos)
            pid = patch_index(gx, gy, gz); ing = pid >= 0
            pidg = pid[ing]; mg = gm[ing]; zg = gz[ing]
            vn = (gvel[ing] - fr['vcen']) @ n
            Sig = bsum(pidg, mg, npatch) / A                          # Msun/kpc^2
            zmid = bsum(pidg, mg * zg, npatch) / np.maximum(Sig * A, 1e-30)
            Hg = np.sqrt(bsum(pidg, mg * (zg - zmid[pidg])**2, npatch) / np.maximum(Sig * A, 1e-30))
            vbar = bsum(pidg, mg * vn, npatch) / np.maximum(Sig * A, 1e-30)
            # intruder gas fraction (gas from the other galaxy) in the column
            gidc = gid[ing]
            own = (gidc < GAS_SPLIT) if fr['label'] == 'A' else (gidc >= GAS_SPLIT) if fr['label'] == 'B' else (gidc < GAS_SPLIT)
            f_intr = 1.0 - bsum(pidg[own], mg[own], npatch) / np.maximum(Sig * A, 1e-30)
            # vertical mass profile per patch and layer cut around the midplane
            dzb = 2 * ZCOL / NZB
            zb = np.clip(((zg + ZCOL) / dzb).astype(np.int64), 0, NZB - 1)
            prof = np.bincount(pidg * NZB + zb, weights=mg, minlength=npatch * NZB)[:npatch * NZB].reshape(npatch, NZB)
            z_lo, z_hi = layer_cut(prof, zmid, ZCOL, dzb)
            inlayer = (zg > z_lo[pidg]) & (zg < z_hi[pidg])
            Sig_layer = bsum(pidg[inlayer], mg[inlayer], npatch) / A
            sigz_col = np.sqrt(bsum(pidg, mg * (vn - vbar[pidg])**2, npatch) / np.maximum(Sig * A, 1e-30))
            # weight, split by source and by side
            dz = zg - zmid[pidg]
            up = dz > 0
            Wt = {}
            col2p = T[ing] < T_2P
            for comp in ('gas', 'star', 'dm', 'total'):
                gn = acc[comp][ing] @ n                                # (km/s)^2/kpc
                Wt[comp + '_up'] = bsum(pidg[up], -(mg * gn)[up], npatch) / A * P_UNIT
                Wt[comp + '_lo'] = bsum(pidg[~up], (mg * gn)[~up], npatch) / A * P_UNIT
                if comp == 'total':
                    Wt['2p_up'] = bsum(pidg[up & col2p], -(mg * gn)[up & col2p], npatch) / A * P_UNIT
                    Wt['2p_lo'] = bsum(pidg[~up & col2p], (mg * gn)[~up & col2p], npatch) / A * P_UNIT
                    Wt['layer_up'] = bsum(pidg[up & inlayer], -(mg * gn)[up & inlayer], npatch) / A * P_UNIT
                    Wt['layer_lo'] = bsum(pidg[~up & inlayer], (mg * gn)[~up & inlayer], npatch) / A * P_UNIT
                    l2 = inlayer & col2p
                    Wt['layer2p_up'] = bsum(pidg[up & l2], -(mg * gn)[up & l2], npatch) / A * P_UNIT
                    Wt['layer2p_lo'] = bsum(pidg[~up & l2], (mg * gn)[~up & l2], npatch) / A * P_UNIT
            # midplane slab
            hs = np.maximum(SLAB_MIN, SLAB_FRAC * Hg)
            inslab = np.abs(dz) < hs[pidg]
            ps = pidg[inslab]; ms = mg[inslab]; Vs = 2 * hs * A            # kpc^3
            Ts = T[ing][inslab]; us = u[ing][inslab]; xi = xhii[ing][inslab]
            hot = Ts > T_HOT; ion = (xi > 0.5) & ~hot; neu = ~ion & ~hot
            Pth = bsum(ps, ms * (C.GAMMA - 1) * us, npatch) / Vs * P_UNIT
            Pth_neu = bsum(ps[neu], (ms * (C.GAMMA - 1) * us)[neu], npatch) / Vs * P_UNIT
            Pth_ion = bsum(ps[ion], (ms * (C.GAMMA - 1) * us)[ion], npatch) / Vs * P_UNIT
            Pth_hot = bsum(ps[hot], (ms * (C.GAMMA - 1) * us)[hot], npatch) / Vs * P_UNIT
            vns = vn[inslab]
            vbar_s = bsum(ps, ms * vns, npatch) / np.maximum(bsum(ps, ms, npatch), 1e-30)
            Pturb = bsum(ps, ms * (vns - vbar_s[ps])**2, npatch) / Vs * P_UNIT
            Pturb_neu = bsum(ps[~hot], (ms * (vns - vbar_s[ps])**2)[~hot], npatch) / Vs * P_UNIT
            # OK22 "two-phase" (T < 2e4 K) midplane quantities
            twop = Ts < T_2P
            vbar_2p = bsum(ps[twop], (ms * vns)[twop], npatch) / np.maximum(bsum(ps[twop], ms[twop], npatch), 1e-30)
            Pth_2p = bsum(ps[twop], (ms * (C.GAMMA - 1) * us)[twop], npatch) / Vs * P_UNIT
            Pturb_2p = bsum(ps[twop], (ms * (vns - vbar_2p[ps])**2)[twop], npatch) / Vs * P_UNIT
            rho_mid_2p = bsum(ps[twop], ms[twop], npatch) / Vs
            # magnetic: Maxwell stress zz = (B^2 - 2 Bn^2)/8pi  [erg/cm^3]; volume weighted
            Bs = B[ing][inslab]; Bn = Bs @ n
            Pi = (np.sum(Bs**2, 1) - 2 * Bn**2) / (8 * np.pi) / C.BOLTZMANN
            Pmagz = (np.sum(Bs**2, 1) / (8 * np.pi)) / C.BOLTZMANN
            vs = vol[ing][inslab]
            Pmag = bsum(ps, vs * Pi, npatch) / Vs
            Pmag_2p = bsum(ps[twop], (vs * Pi)[twop], npatch) / Vs
            Pmag_tot = bsum(ps, vs * Pmagz, npatch) / Vs
            rho_mid = bsum(ps, ms, npatch) / Vs                          # Msun/kpc^3
            fvol_hot = bsum(ps[hot], vs[hot], npatch) / Vs
            fvol_fill = bsum(ps, vs, npatch) / Vs
            G0_mid = bsum(ps, ms * G0[ing][inslab], npatch) / np.maximum(bsum(ps, ms, npatch), 1e-30)
            fH2 = bsum(pidg, mg * xh2[ing], npatch) / np.maximum(Sig * A, 1e-30)
            fHII = bsum(pidg, mg * xhii[ing], npatch) / np.maximum(Sig * A, 1e-30)
            fcold = bsum(pidg[T[ing] < 500], mg[T[ing] < 500], npatch) / np.maximum(Sig * A, 1e-30)
            # velocity-gradient diagnostics (mass-weighted in slab)
            vgs = vg[ing][inslab]
            div = np.trace(vgs, axis1=1, axis2=2)
            S = 0.5 * (vgs + np.swapaxes(vgs, 1, 2)) - div[:, None, None] / 3 * np.eye(3)
            shear = np.sqrt(np.sum(S**2, axis=(1, 2)))
            Om = 0.5 * (vgs - np.swapaxes(vgs, 1, 2)); vort = np.sqrt(2 * np.sum(Om**2, axis=(1, 2)))
            msum = np.maximum(bsum(ps, ms, npatch), 1e-30)
            div_m = bsum(ps, ms * div, npatch) / msum
            shear_m = bsum(ps, ms * shear, npatch) / msum
            vort_m = bsum(ps, ms * vort, npatch) / msum
            # ---- stars & SN
            sx, sy, sz = to_frame(spos)
            sp = patch_index(sx, sy, sz); ins = sp >= 0
            Sig_star = bsum(sp[ins], sm[ins], npatch) / A
            nx, ny, nz = to_frame(s4pos)
            np_ = patch_index(nx, ny, nz); inn = np_ >= 0
            SigSFR = {}
            nlay = inn.copy(); nlay[inn] = (nz[inn] > z_lo[np_[inn]]) & (nz[inn] < z_hi[np_[inn]])
            for w in SFR_WINDOWS:
                yz = inn & (age < w) & (age >= 0)
                SigSFR[w] = bsum(np_[yz], s4m[yz], npatch) / A / (w * 1e6)      # Msun/yr/kpc^2
                yl = yz & nlay
                SigSFR[f'{w}_layer'] = bsum(np_[yl], s4m[yl], npatch) / A / (w * 1e6)
            snx, sny, snz = to_frame(xsn)
            snp = patch_index(snx, sny, snz); insn = snp >= 0
            SNrate = bsum(snp[insn], np.ones(insn.sum()), npatch) / A / 10.0    # per Myr per kpc^2
            # ---- rho_sd at the patch midplane
            ix, iy = np.divmod(np.arange(npatch), nb_)
            zr = zref[np.clip(((ix + 0.5) * PATCH / PATCH_REF).astype(np.int64), 0, NB_REF - 1) * NB_REF + np.clip(((iy + 0.5) * PATCH / PATCH_REF).astype(np.int64), 0, NB_REF - 1)] if zref is not None else 0.0
            pc = fr['cen'] + ((ix + 0.5) * PATCH - RMAX)[:, None] * e1 + ((iy + 0.5) * PATCH - RMAX)[:, None] * e2 + (zmid + zr)[:, None] * n
            rs, rd = dens_at(pc)
            rho_sd = rs + rd
            # analytic P_DE (OK22 eq. 8 style) with sigma_eff from measured midplane P_tot / rho_mid
            Ptot = Pth + Pturb + Pmag
            sig_eff = np.sqrt(np.maximum(Ptot / P_UNIT, 0) / np.maximum(rho_mid, 1e-30))    # km/s
            P_DE = (np.pi * G_KPC * Sig**2 / 2 + Sig * np.sqrt(2 * G_KPC * np.maximum(rho_sd, 0)) * sig_eff) * P_UNIT
            # ---- write
            out = dict(Sigma_gas=Sig / 1e6, zmid=zmid, H=Hg, sigma_z_col=sigz_col, vbar_z=vbar,
                       Pth=Pth, Pth_neu=Pth_neu, Pth_ion=Pth_ion, Pth_hot=Pth_hot, Pturb=Pturb, Pturb_neu=Pturb_neu,
                       Pmag=Pmag, Pmag_tot=Pmag_tot, Ptot=Ptot, rho_mid=rho_mid / 1e9,
                       Pth_2p=Pth_2p, Pturb_2p=Pturb_2p, Pmag_2p=Pmag_2p, Ptot_2p=Pth_2p + Pturb_2p + Pmag_2p,
                       rho_mid_2p=rho_mid_2p / 1e9, Sigma_gas_2p=bsum(pidg[col2p], mg[col2p], npatch) / A / 1e6, fvol_hot=fvol_hot, fvol_fill=fvol_fill,
                       G0_mid=G0_mid, fH2=fH2, fHII=fHII, fcold=fcold, div_v=div_m, shear=shear_m, vort=vort_m,
                       Sigma_star=Sig_star / 1e6, rho_star_mid=rs / 1e9, rho_dm_mid=rd / 1e9, rho_sd=rho_sd / 1e9,
                       SigSFR_10=SigSFR[10.0], SigSFR_40=SigSFR[40.0], SNrate=SNrate, P_DE=P_DE, sigma_eff=sig_eff,
                       SigSFR_10_layer=SigSFR['10.0_layer'], SigSFR_40_layer=SigSFR['40.0_layer'],
                       f_intruder=f_intr, Sigma_gas_layer=Sig_layer / 1e6, z_lo=z_lo, z_hi=z_hi,
                       slab_half=hs, Npart_slab=bsum(ps, np.ones(len(ps)), npatch), Npart_col=bsum(pidg, np.ones(len(pidg)), npatch))
            for k_, v_ in Wt.items():
                out['W_' + k_] = v_
            out['W'] = 0.5 * (Wt['total_up'] + Wt['total_lo'])
            out['W_2p'] = 0.5 * (Wt['2p_up'] + Wt['2p_lo'])
            out['W_layer'] = 0.5 * (Wt['layer_up'] + Wt['layer_lo'])
            out['W_layer2p'] = 0.5 * (Wt['layer2p_up'] + Wt['layer2p_lo'])
            grp.create_dataset('zprofile', data=prof.reshape(nb_, nb_, NZB).astype(np.float32))
            for k_, v_ in out.items():
                grp.create_dataset(k_, data=np.asarray(v_, np.float64).reshape(nb_, nb_))
            print(f'  frame {fr["label"]}: sep={fr["sep"]:.2f} kpc; patches with Sigma>1: {np.sum(Sig/1e6>1)}; '
                  f'median W/Ptot = {np.median((out["W"]/np.maximum(Ptot,1e-30))[Sig/1e6>1]):.2f} ({time.time()-t0:.0f}s)', flush=True)
    os.replace(fout + '.tmp', fout)
    print(f'snap {isnap} done in {time.time()-t0:.0f}s', flush=True)

if __name__ == '__main__':
    for a in sys.argv[1:]:
        analyse(int(a))
