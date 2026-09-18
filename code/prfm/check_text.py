import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import paperstyle as P
_o=P.save
def chk(fig_, name, outdirs=None):
    fig_.canvas.draw(); rend=fig_.canvas.get_renderer()
    for ia,ax in enumerate(fig_.axes):
        leg=ax.get_legend()
        if leg is None: continue
        lb=leg.get_window_extent(rend)
        for tx in ax.texts:
            if not tx.get_text().strip(): continue
            tb=tx.get_window_extent(rend)
            bb=tx.get_bbox_patch()
            if bb is not None:
                try:
                    p=bb.get_window_extent(rend); tb=tb.union([tb,p])
                except Exception: pass
            pad=2.0
            if tb.x1+pad>lb.x0 and tb.x0-pad<lb.x1 and tb.y1+pad>lb.y0 and tb.y0-pad<lb.y1:
                print(f'  {name:22s} ax{ia}: legend touches text {tx.get_text()[:26]!r}')
    _o(fig_, name, *([outdirs] if outdirs else []))
P.save=chk
import paper_figs as F
for n in (sys.argv[1:] or list(F.FIGS)):
    try: F.FIGS[n]()
    except Exception as e: print(f'  {n} FAILED {e}')
    plt.close('all')
print('done')
