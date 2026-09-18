import sys, glob, numpy as np
sys.path.insert(0,'/raven/u/uli/dwarf_merger'); sys.path.insert(0,'/raven/u/uli/dwarf_merger/clusters')
import common as C
from fof import fof
EP=[*C.PHASES]
files=[f for f in sorted(glob.glob(f'{C.DATADIR}/clouds/clouds_[0-9]*.npz')) if int(f[-7:-4])%10==0]
ev={'cl10':[],'cl100':[]}; sizes={'cl10':[],'cl100':[]}
for fn in files:
    d=np.load(fn); k=int(d['snap']); t=float(d['t']); sm=d['star_mass']
    if len(sm)<10: continue
    pos=d['star_pos'].astype(np.float64); lab=fof(pos,50e-3); cnt=np.bincount(lab)
    for cat in ev:
        sc=d[f'{cat}_star_cl']; R=d[f'{cat}_rows']
        sf=R[R[:,16]>500]; sizes[cat].append(np.column_stack([np.full(len(sf),t),sf[:,2],sf[:,11],sf[:,16]]))
        for g in np.flatnonzero(cnt>=25):
            m=lab==g; M=sm[m].sum(); par=sc[m]; ok=par>=0
            if ok.sum()==0: ev[cat].append((t,M,0,0.0,0.0,np.nan)); continue
            u,inv=np.unique(par[ok],return_inverse=True); mp=np.bincount(inv,weights=sm[m][ok]); dom=u[np.argmax(mp)]
            ev[cat].append((t,M,len(u),mp.max()/M,sm[m][ok].sum()/M,R[dom,11]))
for cat in ev:
    E=np.array(ev[cat]); S=np.vstack(sizes[cat])
    print(f'\n{cat}: 50 pc events (>=25 stars) -> parent clouds; SF clouds (M_star>500): r_h')
    print('phase                 n_ev  f(1 cloud) f(<=2) f(<=3)  n_cl med (16-84)  dominant frac med  f_in_cloud | r_h of dominant cloud med [pc] | SF clouds r_h med 16-84 [pc]  M med')
    for t0,t1,lab in EP:
        m=(E[:,0]>=t0)&(E[:,0]<t1)&(E[:,2]>0); s=(S[:,0]>=t0)&(S[:,0]<t1)
        print(f'{lab:20s} {m.sum():5d}   {np.mean(E[m,2]==1):.2f}      {np.mean(E[m,2]<=2):.2f}   {np.mean(E[m,2]<=3):.2f}    {np.median(E[m,2]):3.0f} ({np.percentile(E[m,2],16):.0f}-{np.percentile(E[m,2],84):.0f})        {np.median(E[m,3]):.2f}            {np.median(E[m,4]):.2f}    |      {np.nanmedian(E[m,5]):5.1f}               |   {np.median(S[s,2]):5.1f}  {np.percentile(S[s,2],16):.1f}-{np.percentile(S[s,2],84):.1f}     {np.median(S[s,1]):.0f}')
        mb=m&(E[:,1]>1e4)
        if mb.sum()>=5: print(f'      events > 1e4 Msun: n={mb.sum()}, f(1 cloud) {np.mean(E[mb,2]==1):.2f}, n_cl med {np.median(E[mb,2]):.0f}, dominant frac med {np.median(E[mb,3]):.2f}, dominant r_h med {np.nanmedian(E[mb,5]):.1f} pc')
