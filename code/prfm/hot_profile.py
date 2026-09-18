"""Volume fraction of gas above 2e4 K as a function of height above the midplane, in the galaxy frame, for a few snapshots;
also the volume-weighted temperature distribution within |z| < 50 pc.  Volume of a particle = m / rho."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C
from patches import galaxy_frames, ID_SPLIT
for k in [int(x) for x in sys.argv[1:]] or [60, 150, 180]:
    g = C.read(k, 0, ['Coordinates', 'Velocities', 'Masses', 'Density', 'InternalEnergy', 'ElectronAbundance'])
    s2 = C.read(k, 2, ['Coordinates', 'ParticleIDs', 'Velocities'])
    t = g['Time'] * C.MYR
    pos = g['Coordinates'].astype(np.float64); m = g['Masses'].astype(np.float64); rho = g['Density'].astype(np.float64)
    T = C.temperature(g['InternalEnergy'], g['ElectronAbundance']); vol = m / rho
    fr = galaxy_frames(s2['Coordinates'], s2['ParticleIDs'], s2['Velocities'].astype(np.float64), pos, g['Velocities'].astype(np.float64), m, t)[0]
    d = pos - fr['cen']; z = d @ fr['nhat']; R = np.linalg.norm(d - np.outer(z, fr['nhat']), axis=1)
    sel = R < 2.0
    z = z[sel]; TT = T[sel]; vv = vol[sel]; rr = rho[sel]
    print(f'\nsnap {k}, t = {t:.0f} Myr, frame {fr["label"]}: gas within R < 2 kpc')
    print('  |z| range [pc]   volume [kpc^3]   f_vol(T>2e4)   f_vol(T>1e5)   f_vol(T>1e6)   f_mass(T>2e4)')
    edges = [0, 25, 50, 100, 200, 400, 800, 1600, 3200]
    for a, b in zip(edges[:-1], edges[1:]):
        mm = (np.abs(z) * 1e3 >= a) & (np.abs(z) * 1e3 < b)
        if mm.sum() < 10: continue
        V = vv[mm].sum(); h = TT[mm] > 2e4
        print(f'  {a:5d}-{b:5d}    {V:12.4g}   {vv[mm][h].sum()/V:10.3f}   {vv[mm][TT[mm]>1e5].sum()/V:10.3f}   {vv[mm][TT[mm]>1e6].sum()/V:10.3f}   {m[sel][mm][h].sum()/m[sel][mm].sum():10.3f}')
    mm = np.abs(z) * 1e3 < 50
    print('  within 50 pc: volume-weighted temperature percentiles [K]:', ' '.join(f'{p}%={np.percentile(np.repeat(TT[mm], np.maximum((vv[mm]/vv[mm].min()).astype(int), 1)[:0] ) if False else TT[mm], p):.3g}' for p in (16, 50, 84, 99)))
    o = np.argsort(TT[mm]); cw = np.cumsum(vv[mm][o]) / vv[mm].sum()
    print('  within 50 pc: volume-weighted T median = %.3g K, 84%% = %.3g K, 99%% = %.3g K; n_H median = %.3g cm^-3' % (TT[mm][o][np.searchsorted(cw, 0.5)], TT[mm][o][np.searchsorted(cw, 0.84)], TT[mm][o][np.searchsorted(cw, 0.99)], np.median(C.nH_cgs(rr[mm]))))
