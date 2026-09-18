import sys; sys.path.insert(0, '/raven/u/uli/dwarf_merger')
import numpy as np, common as C
for k in [30, 60, 90, 120]:
    g = C.read(k, 0, ['Masses', 'Density', 'InternalEnergy', 'ElectronAbundance'])
    T = C.temperature(g['InternalEnergy'], g['ElectronAbundance']); u = g['InternalEnergy'].astype(np.float64)
    cs = np.sqrt((C.GAMMA - 1) * u)            # isothermal sound speed km/s
    nH = C.nH_cgs(g['Density']); m = g['Masses']
    for lab, s in [('T<500K', T < 500), ('T<100K', T < 100), ('nH>10', nH > 10), ('nH>100', nH > 100), ('nH>1000', nH > 1000)]:
        if s.sum() == 0: print(f'snap {k} {lab}: none'); continue
        w = m[s] / m[s].sum(); o = np.argsort(cs[s]); cw = np.cumsum(w[o])
        p = np.interp([0.16, 0.5, 0.84], cw, cs[s][o])
        print(f'snap {k} {lab:7s}: frac of gas mass {m[s].sum()/m.sum():.3f}  cs mass-wtd mean {np.sum(w*cs[s]):.2f}  16/50/84 % {p[0]:.2f} {p[1]:.2f} {p[2]:.2f} km/s   T median {np.interp(0.5, cw, T[s][o]):.0f} K', flush=True)
