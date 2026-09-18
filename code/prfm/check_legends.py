"""For every figure in paper_figs: check that the legend lies inside the axes and that no plotted line or point falls
inside the legend's bounding box.  Prints one line per axes with a legend."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import paperstyle as P
_orig_save = P.save
BAD = []
def check(fig_, name, outdirs=None):
    fig_.canvas.draw(); rend = fig_.canvas.get_renderer()
    for ia, ax in enumerate(fig_.axes):
        leg = ax.get_legend()
        if leg is None: continue
        lb = leg.get_window_extent(rend)
        ab = ax.get_window_extent(rend)
        out = not (lb.x0 >= ab.x0 - 1 and lb.x1 <= ab.x1 + 1 and lb.y0 >= ab.y0 - 1 and lb.y1 <= ab.y1 + 1)
        pad = 1.0
        box = [lb.x0 - pad, lb.x1 + pad, lb.y0 - pad, lb.y1 + pad]
        hit = 0; tot = 0
        for ln in ax.get_lines():
            if not ln.get_visible() or ln.get_label().startswith('_'): continue
            xy = ln.get_xydata()
            if len(xy) == 0: continue
            d = ax.transData.transform(xy); d = d[np.isfinite(d).all(1)]
            if len(xy) == 2 and (np.allclose(xy[:, 0], xy[0, 0]) or np.allclose(xy[:, 1], xy[0, 1])): continue   # axhline/axvline guides
            tot += len(d); hit += int(((d[:, 0] > box[0]) & (d[:, 0] < box[1]) & (d[:, 1] > box[2]) & (d[:, 1] < box[3])).sum())
        for co in ax.collections:
            off = co.get_offsets()
            if off is None or len(off) == 0: continue
            d = co.get_offset_transform().transform(off); d = d[np.isfinite(d).all(1)]
            tot += len(d); hit += int(((d[:, 0] > box[0]) & (d[:, 0] < box[1]) & (d[:, 1] > box[2]) & (d[:, 1] < box[3])).sum())
        tag = 'OUTSIDE-AXES' if out else ('OVERLAP' if hit else 'ok')
        if tag != 'ok': BAD.append((name, ia, tag, hit, tot))
        print(f'  {name:22s} axes {ia}: {tag:12s} points in legend box {hit:5d} of {tot}')
    _orig_save(fig_, name, *( [outdirs] if outdirs else [] ))
P.save = check
import paper_figs as F
names = sys.argv[1:] or list(F.FIGS)
for n in names:
    try: F.FIGS[n]()
    except Exception as e: print(f'  {n:22s} FAILED: {type(e).__name__}: {e}')
    plt.close('all')
print('\nproblems:', len(BAD))
for b in BAD: print('  ', b)
