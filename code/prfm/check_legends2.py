"""Report, per figure, where the legend ended up, whether the axis had to be stretched, and how much empty space that left."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import paperstyle as P
_orig = P.save
def check(fig_, name, outdirs=None):
    fig_.canvas.draw(); rend = fig_.canvas.get_renderer()
    for ia, ax in enumerate(fig_.axes):
        leg = ax.get_legend()
        if leg is None: continue
        lb = leg.get_window_extent(rend); ab = ax.get_window_extent(rend)
        fx = (lb.x0 - ab.x0) / ab.width; fy = (lb.y0 - ab.y0) / ab.height
        fw = lb.width / ab.width; fh = lb.height / ab.height
        # how much of the vertical range is empty above the highest datum and below the lowest
        ys = []
        for ln in ax.get_lines():
            d = ln.get_ydata()
            if len(d) > 2: ys.append(np.asarray(d, float))
        for co in ax.collections:
            off = co.get_offsets()
            if off is not None and len(off): ys.append(np.asarray(off)[:, 1])
        lo, hi = ax.get_ylim(); log = ax.get_yscale() == 'log'
        gap = ''
        if ys:
            y = np.concatenate(ys); y = y[np.isfinite(y)]
            if log: y = y[y > 0]
            if len(y):
                if log: top = (np.log10(hi) - np.log10(y.max())) / (np.log10(hi) - np.log10(lo)); bot = (np.log10(y.min()) - np.log10(lo)) / (np.log10(hi) - np.log10(lo))
                else: top = (hi - y.max()) / (hi - lo); bot = (y.min() - lo) / (hi - lo)
                gap = f'empty top {top:.2f} bottom {bot:.2f}'
        print(f'  {name:22s} ax{ia}: legend at x {fx:.2f} y {fy:.2f} (w {fw:.2f} h {fh:.2f})  fontsize {leg.get_texts()[0].get_fontsize():.1f}  {gap}')
    _orig(fig_, name, *( [outdirs] if outdirs else [] ))
P.save = check
import paper_figs as F
for n in (sys.argv[1:] or list(F.FIGS)):
    try: F.FIGS[n]()
    except Exception as e: print(f'  {n:22s} FAILED {e}')
    plt.close('all')
