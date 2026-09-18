"""Bound cluster mass over all stellar mass formed in the same age window, per snapshot (every STEP-th):
Gamma_young (age<10 Myr) and Gamma_old (20-50 Myr), whole system and clean-patch subset not needed.
usage: python gamma_time.py [STEP]"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, common as C
step = int(sys.argv[1]) if len(sys.argv) > 1 else 10
cats = {'young(0-10)': np.load(f'{C.DATADIR}/clusters/clusters_age0-10_l5_n25.npz'), 'old(20-50)': np.load(f'{C.DATADIR}/clusters/clusters_age20-50_l5_n25.npz')}
win = {'young(0-10)': (0, 10), 'old(20-50)': (20, 50)}
print('snap    t   ' + '   '.join(f'{k}: M_formed  M_grouped  M_bound  Gamma_all Gamma_bound' for k in cats))
res = []
for k in range(0, C.NSNAP, step):
    st = np.load(f'{C.DATADIR}/stars/stars_{k:03d}.npz'); t = float(st['time_myr']); age = t - st['tform_myr']
    row = [k, t]; line = f'{k:4d} {t:6.1f}'
    for name, c in cats.items():
        lo, hi = win[name]
        Mf = np.sum(st['mass'][(age >= lo) & (age < hi) & (st['tform_myr'] > 0)]) * C.MSUN
        s = (c['snap'] == k) & (c['nucleus'] == 0)
        Mg = c['M'][s].sum(); Mb = c['M'][s & (c['bound'] > 0)].sum()
        row += [Mf, Mg, Mb]; line += f'   {Mf:9.0f} {Mg:9.0f} {Mb:9.0f}  {Mg/max(Mf,1):.2f}  {Mb/max(Mf,1):.2f}'
    print(line); res.append(row)
res = np.array(res); np.savez(f'{C.DATADIR}/clusters/gamma_time.npz', data=res)
for j, name in enumerate(cats):
    Mf, Mg, Mb = res[:, 2 + 3 * j], res[:, 3 + 3 * j], res[:, 4 + 3 * j]
    for lo, hi in [(0, 50), (50, 100), (100, 150), (150, 232)]:
        e = (res[:, 1] >= lo) & (res[:, 1] < hi)
        print(f'{name} {lo:3d}-{hi:3d} Myr: Gamma_bound = {Mb[e].sum()/max(Mf[e].sum(),1):.2f}  (grouped {Mg[e].sum()/max(Mf[e].sum(),1):.2f})')
