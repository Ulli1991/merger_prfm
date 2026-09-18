"""Paper (MNRAS-like) matplotlib style shared by the figure scripts: serif fonts, ticks inside on all four sides, minor ticks,
no grid, single-column width 3.32 in (84 mm)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
COL1 = 3.32; COL2 = 7.0            # inches
ink = '#1f2937'; blue = '#2f6fb2'; light = '#9dbfe3'; grey = '#6b7280'; orange = '#c2410c'; green = '#2e7d4f'
import common as C
PHASES = [(25, 40, 'before 1st\npassage'), (40, 110, '1st to 2nd\npassage'), (110, 169, '2nd to 3rd\npassage'), (169, 226, '3rd passage to\ncoalescence')]
PERICENTRES = C.PERICENTRES

def use():
    plt.rcParams.update({
        'font.family': 'serif', 'font.serif': ['STIXGeneral', 'Times New Roman', 'DejaVu Serif'], 'mathtext.fontset': 'stix',
        'font.size': 8.5, 'axes.labelsize': 9, 'axes.titlesize': 9, 'legend.fontsize': 7.5, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
        'axes.linewidth': 0.7, 'lines.linewidth': 1.2,
        'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True,
        'xtick.minor.visible': True, 'ytick.minor.visible': True,
        'xtick.major.size': 3.5, 'xtick.minor.size': 2, 'ytick.major.size': 3.5, 'ytick.minor.size': 2,
        'xtick.major.width': 0.7, 'ytick.major.width': 0.7, 'xtick.minor.width': 0.5, 'ytick.minor.width': 0.5,
        'axes.grid': False, 'legend.frameon': False, 'legend.handlelength': 1.8,
        'savefig.dpi': 300, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02,
        'axes.spines.top': True, 'axes.spines.right': True, 'axes.prop_cycle': plt.cycler(color=[blue, orange, ink, green, light, grey]),
    })

def fig(width=COL1, aspect=0.72, **kw):
    return plt.subplots(figsize=(width, width * aspect), **kw)

def phases(ax, labels=True, y=0.97, fontsize=7, box=False):
    """phase boundaries as thin dotted lines, coalescence as a solid line; labels along the top inside the axes.
    A label is placed at the centre of the part of its phase that is actually visible, nudged inward if it would
    cross the panel edge, and dropped if the visible span is too narrow to hold it."""
    for tp in PERICENTRES: ax.axvline(tp, color=grey, lw=0.6, ls=':' if tp < 200 else '-')
    if not labels: return
    x0, x1 = ax.get_xlim(); txt = []
    for t0, t1, lab in PHASES:
        a, b = max(t0, min(x0, x1)), min(t1, max(x0, x1))
        if b - a <= 0.02 * abs(x1 - x0): continue
        t = ax.text(0.5 * (a + b), y, lab, transform=ax.get_xaxis_transform(), ha='center', va='top',
                    fontsize=fontsize, color=grey,
                    bbox=dict(facecolor='white', edgecolor='none', alpha=0.85, pad=1.5) if box else None)
        txt.append((t, a, b))
    fig = ax.figure; fig.canvas.draw(); rend = fig.canvas.get_renderer(); ab = ax.get_window_extent(rend)
    for t, a, b in txt:
        tb = t.get_window_extent(rend)
        span = abs(ax.transData.transform((b, 0))[0] - ax.transData.transform((a, 0))[0])
        if tb.width > span * 1.02 and tb.width > 0.9 * ab.width:
            t.remove(); continue                                   # no room for it at all
        if tb.x0 < ab.x0 + 1:                                      # hanging off the left edge
            t.set_ha('left'); t.set_x(ax.transData.inverted().transform((ab.x0 + 2, 0))[0])
        elif tb.x1 > ab.x1 - 1:                                    # off the right edge
            t.set_ha('right'); t.set_x(ax.transData.inverted().transform((ab.x1 - 2, 0))[0])

def _data_in_box(ax, box, rend=None):
    """number of plotted points (lines and scatter) inside a display-coordinate box [x0, x1, y0, y1]; any text of the
    axes (phase labels, annotations) that intersects the box counts heavily so the legend never lands on it"""
    import numpy as np
    n = 0
    rend = rend or ax.figure.canvas.get_renderer()
    for tx in ax.texts:
        if not tx.get_visible() or not tx.get_text().strip(): continue
        try: tb = tx.get_window_extent(rend)
        except Exception: continue
        if tb.x1 > box[0] and tb.x0 < box[1] and tb.y1 > box[2] and tb.y0 < box[3]: n += 1000
    for ln in ax.get_lines():
        if not ln.get_visible(): continue
        xy = ln.get_xydata()
        if len(xy) == 0: continue
        if len(xy) == 2 and (np.allclose(xy[:, 0], xy[0, 0]) or np.allclose(xy[:, 1], xy[0, 1])): continue   # axhline/axvline guides
        d = ax.transData.transform(xy); d = d[np.isfinite(d).all(1)]
        n += int(((d[:, 0] > box[0]) & (d[:, 0] < box[1]) & (d[:, 1] > box[2]) & (d[:, 1] < box[3])).sum())
    for co in ax.collections:
        off = co.get_offsets()
        if off is None or len(off) == 0: continue
        d = co.get_offset_transform().transform(np.asarray(off)); d = d[np.isfinite(d).all(1)]
        n += int(((d[:, 0] > box[0]) & (d[:, 0] < box[1]) & (d[:, 1] > box[2]) & (d[:, 1] < box[3])).sum())
    return n

def _grow(ax, corner, frac=0.10):
    """extend the axis range away from `corner` so that the legend there gets free space"""
    import numpy as np
    lo, hi = ax.get_ylim(); log = ax.get_yscale() == 'log'
    if log:
        l, h = np.log10(lo), np.log10(hi); d = (h - l) * frac
        ax.set_ylim(10 ** (l - d), hi) if 'lower' in corner else ax.set_ylim(lo, 10 ** (h + d))
    else:
        d = (hi - lo) * frac
        ax.set_ylim(lo - d, hi) if 'lower' in corner else ax.set_ylim(lo, hi + d)

def place_legend(ax, loc=None, pad=1.5, grow=True, min_fontsize=6.0, **kw):
    """Put the legend inside the axes where nothing is plotted under it.  Tries the corners, keeping `loc` first if
    given; if every corner is covered it extends the y range on that side by exactly the amount the legend needs, and
    only if that fails does it shrink the legend.  Never places the legend outside the axes."""
    import numpy as np
    fig = ax.figure
    corners = ['upper left', 'upper right', 'lower left', 'lower right']
    order = ([loc] if loc in corners else []) + [c for c in corners if c != loc]
    if loc and loc not in corners: order = [loc] + order

    def measure(leg):
        fig.canvas.draw(); rend = fig.canvas.get_renderer()
        lb = leg.get_window_extent(rend); ab = ax.get_window_extent(rend)
        inside = lb.x0 >= ab.x0 - 1 and lb.x1 <= ab.x1 + 1 and lb.y0 >= ab.y0 - 1 and lb.y1 <= ab.y1 + 1
        hit = _data_in_box(ax, [lb.x0 - pad, lb.x1 + pad, lb.y0 - pad, lb.y1 + pad], rend) if inside else 10 ** 9
        return lb, ab, inside, hit

    def fit(leg, corner):
        """extend the y range on the legend's side by exactly what the legend needs over its own x span"""
        for _ in range(3):
            lb, ab, inside, hit = measure(leg)
            if hit == 0: return True
            if not inside: return False
            hf = (lb.height + 2 * pad) / ab.height
            x0, x1 = (ax.transData.inverted().transform([(lb.x0 - pad, lb.y0), (lb.x1 + pad, lb.y0)]))[:, 0]
            xs, ys = [], []
            for ln in ax.get_lines():
                xy = ln.get_xydata()
                if len(xy) < 3: continue
                xs.append(xy[:, 0]); ys.append(xy[:, 1])
            for co in ax.collections:
                off = co.get_offsets()
                if off is not None and len(off): off = np.asarray(off); xs.append(off[:, 0]); ys.append(off[:, 1])
            if not xs: return False
            X = np.concatenate(xs); Y = np.concatenate(ys)
            m = np.isfinite(X) & np.isfinite(Y) & (X >= min(x0, x1)) & (X <= max(x0, x1))
            lo, hi = ax.get_ylim(); log = ax.get_yscale() == 'log'
            if log: m &= (Y > 0)
            if not m.any(): return False
            up = 'upper' in corner
            yd = Y[m].max() if up else Y[m].min()
            if log:
                l, h, d = np.log10(lo), np.log10(hi), np.log10(yd)
                if up: ax.set_ylim(lo, 10 ** (l + (d - l) / max(1 - hf, 0.2)))
                else:  ax.set_ylim(10 ** (h - (h - d) / max(1 - hf, 0.2)), hi)
            else:
                if up: ax.set_ylim(lo, lo + (yd - lo) / max(1 - hf, 0.2))
                else:  ax.set_ylim(hi - (hi - yd) / max(1 - hf, 0.2), hi)
        return measure(leg)[3] == 0

    base = kw.pop('fontsize', plt.rcParams['legend.fontsize']); ncol0 = kw.pop('ncol', 1)
    ylim0 = ax.get_ylim(); best = None                        # (hits, corner, fontsize, ncol, ylim)
    def note(c, fs, ncol):
        nonlocal best
        h = measure(ax.get_legend())[3]
        if best is None or h < best[0]: best = (h, c, fs, ncol, ax.get_ylim())
        return h
    for fs in [base] + ([f for f in (base - 0.5, base - 1.0, base - 1.5, base - 2.0) if f >= min_fontsize]):
        for ncol in ([ncol0] if ncol0 > 1 else [1, 2]):
            for c in order:                                   # first pass: a corner that is already clear
                leg = ax.legend(loc=c, fontsize=fs, ncol=ncol, **kw)
                if note(c, fs, ncol) == 0: return leg
            if not grow: continue
            for c in order:                                   # second pass: make room on that side
                ax.set_ylim(*ylim0)
                leg = ax.legend(loc=c, fontsize=fs, ncol=ncol, **kw)
                if fit(leg, c): return leg
                note(c, fs, ncol)
            ax.set_ylim(*ylim0)
    h, c, fs, ncol, yl = best                                 # nothing was clear: the least bad of everything tried
    ax.set_ylim(*yl)
    leg = ax.legend(loc=c, fontsize=fs, ncol=ncol, **kw)
    for _ in range(12):                                       # keep making room until the last few points are out
        if measure(leg)[3] == 0: break
        _grow(ax, c, 0.08)
    return leg

def save(fig_, name, outdirs=('/ptmp/uli/dwarf_merger/prfm/figs/paper', '/u/uli/merger_prfm/figs')):
    import os
    for d in outdirs:
        os.makedirs(d, exist_ok=True); fig_.savefig(f'{d}/{name}.png'); fig_.savefig(f'{d}/{name}.pdf')
    print('saved', name, 'to', outdirs)
