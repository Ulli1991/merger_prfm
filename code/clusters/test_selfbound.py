import sys, numpy as np, time; sys.path.insert(0,'/raven/u/uli/dwarf_merger/clusters')
P='/raven/u/uli/dwarf_merger/clusters/cluster_mf.py'
src=open(P).read().split('# ---- ')[0].replace("a = ap.parse_args()","a = ap.parse_args(['--step','10'])")
ns={'__file__': P, '__name__': 'cluster_mf_test'}; exec(compile(src, P, 'exec'), ns)
for k in (60, 130):
    t0=time.time(); t, out = ns['find_clusters'](k); R=np.array(out); M=R[:,2]; Msb=R[:,17]; bd=R[:,15]
    print(f'snap {k} t={t:.0f}: {len(out)} groups in {time.time()-t0:.0f}s;  mass in groups {M.sum():.3g}, E<0 bound frac {M[bd>0].sum()/M.sum():.2f}, self-bound frac {Msb.sum()/M.sum():.2f}', flush=True)
    for i in np.argsort(M)[-4:]: print(f'   M={M[i]:8.0f} N={int(R[i,3])} rh={R[i,10]:5.1f}pc bound={int(bd[i])}  M_sb={Msb[i]:8.0f} ({Msb[i]/M[i]:.2f}) rh_sb={R[i,19]:5.1f}pc', flush=True)
