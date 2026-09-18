"""Single-panel figures for the story page, one per step (prfm/figs/story_*.png)."""
import numpy as np, sys, os, glob, h5py
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C, ok22
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.optimize import minimize
ink='#1f2937'; blue='#2f6fb2'; light='#9dbfe3'; grey='#6b7280'; orange='#c2410c'; FIG='/ptmp/uli/dwarf_merger/prfm/figs'
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': True, 'grid.alpha': 0.2, 'legend.frameon': False})
def new(): fig,ax=plt.subplots(figsize=(7.2,4.0)); return fig,ax
def phases(ax):
    for t0,t1,lab in [*C.PHASES]:
        ax.axvspan(t0,t1,color='k',alpha=0.03 if lab in('between passages','remnant') else 0.0,lw=0); ax.text(0.5*(t0+t1),ax.get_ylim()[1],lab,ha='center',va='top',fontsize=8,color=grey)
    ax.axvline(C.T_MERGED,color=ink,lw=0.8,ls=':')
# 0 separation
L=np.load('/ptmp/uli/dwarf_merger/prfm/layer_verdict_z05.npz')
fig,ax=new(); ax.plot(L['TF'],L['sepF'],color=blue,lw=2); ax.set_ylim(0,4.5); phases(ax); ax.set_xlabel('t [Myr]'); ax.set_ylabel('$d_{\\rm nuclei}$  [kpc]'); fig.savefig(f'{FIG}/story_0_sep.png',dpi=150,bbox_inches='tight'); plt.close(fig)
# 1 equilibrium
fig,ax=new(); ax.plot(L['TF'],L['mwF'],color=light,lw=2,label='full column, 3 kpc'); ax.plot(L['TL'],L['mwL'],color=blue,lw=2.2,label='gas layer, |z| < 0.5 kpc'); ax.axhline(1,color=ink,lw=0.8); ax.set_yscale('log'); ax.set_ylim(0.15,3); phases(ax)
ax.set_xlabel('t [Myr]'); ax.set_ylabel('$P_{\\rm tot}\\,/\\,\\mathcal{W}$'); ax.legend(loc='lower left'); fig.savefig(f'{FIG}/story_1_equilibrium.png',dpi=150,bbox_inches='tight'); plt.close(fig)
# 2 feedback share
P=dict(np.load('/ptmp/uli/dwarf_merger/prfm/pressure_sources.npz')); T=np.unique(P['snap']); per=lambda f: np.array([f(P['snap']==s) for s in T])
tt=per(lambda m:P['t'][m][0]); ffb=per(lambda m:np.clip(np.sum(P['P_fb'][m])/np.sum(P['Ptot_2p'][m]),0,1)); ok=tt>10
fig,ax=new(); ax.fill_between(tt[ok],0,ffb[ok],color=blue,lw=0,label='supplied by feedback: $\\Upsilon_{\\rm tot}(P)\\,\\Sigma_{\\rm SFR}$'); ax.fill_between(tt[ok],ffb[ok],1,color=light,lw=0,label='supplied by the flow'); ax.set_ylim(0,1); phases(ax); ax.set_xlabel('t [Myr]'); ax.set_ylabel('$P_{\\rm feedback}\\,/\\,P_{\\rm tot}$'); ax.legend(loc='lower left'); fig.savefig(f'{FIG}/story_2_feedback.png',dpi=150,bbox_inches='tight'); plt.close(fig)
# 3 shear closure (quiescent cells): P_turb,res vs Sigma H S_p^2
rows=[]
for fn in sorted(glob.glob(f'{C.DATADIR}/prfm/turb_*.npz')):
    d=dict(np.load(fn)); k=int(d['snap'][0])
    with h5py.File(f'{C.DATADIR}/prfm/patch_{k:03d}.h5','r') as f:
        for fr in np.unique(d['frame']):
            g=f[fr]; m=d['frame']==fr; p=d['patch'][m].astype(int); H=g['H'][:].ravel()[p]
            for j,i in enumerate(np.flatnonzero(m)): rows.append([d['t'][i],d['sep'][i],d['Sigma_gas'][i],d['Pturb_2p'][i],d['f_intruder'][i],d['SigSFR_10'][i],d['sig_tot'][i],d['sig_res'][i],d['shear_patch'][i],H[j]])
R=np.array(rows); t,sep,Sig,Pt,fi,sfr,st,sr,Sp,H=R.T; ok=((fi<0.1)|C.is_merged(sep,t))&(Pt>0)&(st>0)&(t>10)&(sfr==0)&(Sp>0)
X=Sig*1e6*H*Sp**2*4.903245584337325e-06; Pres=Pt*(sr/st)**2
fig,ax=new(); ax.scatter(X[ok],Pres[ok],s=8,c=blue,alpha=0.45,label='quiescent patches (no star formation in 10 Myr)'); xx=np.logspace(3,8,50); ax.plot(xx,0.02*xx,color=ink,lw=2,label='$P_{\\rm turb} = \\varepsilon_S\\,\\Sigma H S^2$,  shear efficiency $\\varepsilon_S$ = 0.02')
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(1e3,2e6); ax.set_ylim(1e1,3e4); ax.set_xlabel('$\\Sigma\\,H\\,S^2$  [K cm$^{-3}$]'); ax.set_ylabel('$P_{\\rm turb}$  [K cm$^{-3}$]'); ax.legend(loc='upper left'); fig.savefig(f'{FIG}/story_3_shear.png',dpi=150,bbox_inches='tight'); plt.close(fig)
# 4 dynamo: field strength from Pmag_tot = B^2/8pi (NOT Pmag_2p, which is the Maxwell stress and can be negative); prfm/bfield_time.npz
Rb=np.load(f'{C.DATADIR}/prfm/bfield_time.npz')['rows']   # columns: snap, t, <B>_mass [muG], median B, median Pmag_tot/(Pth+Pturb), median stress share
fig,ax=new(); ax.plot(Rb[:,1],Rb[:,2],color=blue,lw=2.2,label='$\\langle B\\rangle$ over clean patches, mass-weighted [$\\mu$G]'); ax.plot(Rb[:,1],Rb[:,4],color=ink,lw=1.5,ls='--',label='$P_{\\rm mag}\\,/\\,(P_{\\rm th}+P_{\\rm turb})$, median patch')
ax.set_yscale('log'); ax.set_ylim(3e-3,30); phases(ax); ax.set_xlabel('t [Myr]'); ax.set_ylabel('$B$ [$\\mu$G]   and   $P_{\\rm mag}/(P_{\\rm th}+P_{\\rm turb})$'); ax.legend(loc='lower right'); fig.savefig(f'{FIG}/story_4_dynamo.png',dpi=150,bbox_inches='tight'); plt.close(fig)
# 5 cluster MF per phase
E=np.load('/ptmp/uli/dwarf_merger/clusters/clusters_env.npz'); cleanE=(E['f_intruder']<0.1)|C.is_merged(E['sep'],E['t']); cb=cleanE&(E['bound']>0)&(E['nucleus']==0)
def schechter(M,mmin=300.):
    M=M[M>=mmin]; n=len(M); x=np.logspace(np.log10(mmin),8,2000)
    def nll(p):
        a,lMc=p; Mc=10**lMc
        if a<0 or a>3.5 or lMc<np.log10(mmin) or lMc>8: return 1e30
        Z=np.trapezoid(x**-a*np.exp(-x/Mc),x); return -(np.sum(-a*np.log(M)-M/Mc)-n*np.log(Z))
    best=min((minimize(nll,[a0,l0],method='Nelder-Mead') for a0 in (1.3,2.0) for l0 in (3.5,4.5,5.5)),key=lambda r:r.fun); return best.x[0],10**best.x[1],n
EP=[*C.PHASES]; cols=['#9dbfe3','#5b93c9','#2f6fb2','#1b3f6b']
fig,ax=new(); edges=np.logspace(2.5,5.6,14)
for (t0,t1,lab),c in zip(EP,cols):
    M=E['M'][cb&(E['t']>=t0)&(E['t']<t1)]; n,_=np.histogram(M,edges); dl=np.diff(np.log10(edges)); x=np.sqrt(edges[:-1]*edges[1:]); y=n/dl/len(M); a,Mc,nn=schechter(M); k=n>0
    ax.plot(x[k],y[k],'o',color=c,ms=5,label=f'{lab}:  α = {a:.2f},  $M_c$ = {Mc:.1e}'); xx=np.logspace(2.5,5.6,100); yy=xx**(1-a)*np.exp(-xx/Mc); yy*=y[k][0]/np.interp(x[k][0],xx,yy); ax.plot(xx,yy,color=c,lw=1.5)
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_ylim(1e-4,3); ax.set_xlabel('$M_{\\rm cluster}$  [M$_\\odot$]'); ax.set_ylabel('dN / dlog M'); ax.legend(fontsize=8.5); fig.savefig(f'{FIG}/story_5_mf.png',dpi=150,bbox_inches='tight'); plt.close(fig)
# 6 reservoir: M_max vs M_young per patch
D=dict(np.load('/ptmp/uli/dwarf_merger/prfm/patches_clusters.npz')); clean=(D['f_intruder']<0.1)|C.is_merged(D['sep'],D['t']); m=clean&(D['M_young']>500)&(D['M_max']>0)&(D['Sigma_gas']>1)
fig,ax=new(); ax.scatter(D['M_young'][m],D['M_max'][m],s=7,c=blue,alpha=0.4,label=f'{m.sum()} patches with a bound cluster'); xx=np.logspace(2.7,6.5,10); ax.plot(xx,0.5*xx,color=ink,lw=2,label='$M_{\\max} = 0.5\\,M_{\\rm young}$'); ax.plot(xx,xx,color=grey,lw=1,ls=':',label='$M_{\\max} = M_{\\rm young}$')
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(5e2,3e6); ax.set_ylim(1e2,3e5); ax.set_xlabel('$M_{\\rm young}$  [M$_\\odot$]'); ax.set_ylabel('$M_{\\max}$  [M$_\\odot$]'); ax.legend(loc='upper left'); fig.savefig(f'{FIG}/story_6_reservoir.png',dpi=150,bbox_inches='tight'); plt.close(fig)
# 7 ceiling
Rb=np.load('/ptmp/uli/dwarf_merger/prfm/cluster_link_theory.npz')['bins']; act=Rb[:,-1]>0.09e-3
fig,ax=new(); W=np.logspace(3.5,5.8,50); ax.plot(W,0.5*ok22.sfr_of_W(W)*0.25*1e7,color=ink,lw=2,label='PRFM ceiling  $M_{\\max} = 0.5\\,\\Sigma_{\\rm SFR}(\\mathcal{W})\\,A\\,\\tau$  (nothing fitted)')
ax.scatter(Rb[act,4],Rb[act,1],s=60,c=blue,zorder=3,label='largest bound cluster per 25 Myr interval'); ax.scatter(Rb[~act,4],Rb[~act,1],s=90,facecolors='none',edgecolors=orange,lw=2,zorder=3,label='quiescent interval 160–185 Myr')
for r in Rb: ax.annotate(f'{r[0]:.0f}',(r[4],r[1]),xytext=(5,-3),textcoords='offset points',fontsize=8,color=grey)
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(5e3,8e5); ax.set_ylim(5e2,5e5); ax.set_xlabel('$\\mathcal{W}$  [K cm$^{-3}$]'); ax.set_ylabel('$M_{\\max}$  [M$_\\odot$]'); ax.legend(loc='upper left',fontsize=8.5); fig.savefig(f'{FIG}/story_7_ceiling.png',dpi=150,bbox_inches='tight'); plt.close(fig)
for f in sorted(glob.glob(f'{FIG}/story_*.png')): print(os.path.basename(f), os.path.getsize(f)//1024,'KB'); os.system(f'cp {f} ~/out/')
