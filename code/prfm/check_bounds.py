"""any text, legend or label of an axes that crosses its own panel edge"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import paperstyle as P
_o=P.save; BAD=[]
def chk(fig_, name, outdirs=None):
    fig_.canvas.draw(); rend=fig_.canvas.get_renderer()
    for ia,ax in enumerate(fig_.axes):
        ab=ax.get_window_extent(rend)
        for kind,obj in [('text',t) for t in ax.texts]+([('legend',ax.get_legend())] if ax.get_legend() else []):
            try: b=obj.get_window_extent(rend)
            except Exception: continue
            if kind=='text' and not obj.get_text().strip(): continue
            if b.x0 < ab.x0-1 or b.x1 > ab.x1+1 or b.y0 < ab.y0-1 or b.y1 > ab.y1+1:
                lbl = obj.get_text()[:24] if kind=='text' else 'legend'
                print(f'  {name:22s} ax{ia}: {kind} {lbl!r} crosses the panel edge'); BAD.append(name)
    _o(fig_, name, *([outdirs] if outdirs else []))
P.save=chk
import paper_figs as F
for n in (sys.argv[1:] or list(F.FIGS)):
    try: F.FIGS[n]()
    except Exception as e: print(f'  {n} FAILED {e}')
    plt.close('all')
print('out-of-bounds items:', len(BAD))
