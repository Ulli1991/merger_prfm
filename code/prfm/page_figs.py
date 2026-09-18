import numpy as np, sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C, ok22
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.optimize import minimize
ink='#1f2937'; blue='#2f6fb2'; light='#9dbfe3'; grey='#6b7280'; orange='#c2410c'
FIG='/ptmp/uli/dwarf_merger/prfm/figs'
plt.rcParams.update({'font.size': 9.5, 'axes.titlesize': 10})
# ---------- Fig A: the ceiling
R=np.load('/ptmp/uli/dwarf_merger/prfm/cluster_link_theory.npz')['bins']; act=R[:,-1]>0.09e-3
fig,ax=plt.subplots(figsize=(6.4,4.6)); ax.grid(alpha=0.2); ax.spines[['top','right']].set_visible(False); ax.set_xscale('log'); ax.set_yscale('log')
W=np.logspace(3.5,5.8,50); ax.plot(W,0.5*ok22.sfr_of_W(W)*0.25*1e7,color=ink,lw=2,label='PRFM ceiling  $M_{\\max} = 0.5\\,\\Sigma_{\\rm SFR}(\\mathcal{W})\\,A\\,\\tau$  (no free parameter)')
ax.scatter(R[act,4],R[act,1],s=60,c=blue,zorder=3,label='largest bound cluster, 25 Myr bins with star formation')
ax.scatter(R[~act,4],R[~act,1],s=90,facecolors='none',edgecolors=orange,lw=2,zorder=3,label='quiescent interval 160–185 Myr')
for r in R: ax.annotate(f'{r[0]:.0f}',(r[4],r[1]),xytext=(5,-3),textcoords='offset points',fontsize=8,color=grey)
ax.set_xlabel('layer weight of the heaviest star-forming patches (90th percentile) [K cm$^{-3}$]'); ax.set_ylabel('largest bound cluster [M$_\\odot$]')
ax.set_xlim(5e3,8e5); ax.set_ylim(5e2,5e5); ax.legend(frameon=False,fontsize=8,loc='upper left')
fig.savefig(f'{FIG}/ceiling.png',dpi=150,bbox_inches='tight'); plt.close(fig)
# ---------- Fig B: cluster MF per phase with Schechter fits
E=np.load('/ptmp/uli/dwarf_merger/clusters/clusters_env.npz'); cleanE=(E['f_intruder']<0.1)|C.is_merged(E['sep'],E['t']); cb=cleanE&(E['bound']>0)&(E['nucleus']==0)
def schechter(M,mmin=300.):
    M=M[M>=mmin]; n=len(M); x=np.logspace(np.log10(mmin),8,2000)
    def nll(p):
        a,lMc=p; Mc=10**lMc
        if a<0 or a>3.5 or lMc<np.log10(mmin) or lMc>8: return 1e30
        Z=np.trapezoid(x**-a*np.exp(-x/Mc),x); return -(np.sum(-a*np.log(M)-M/Mc)-n*np.log(Z))
    best=min((minimize(nll,[a0,l0],method='Nelder-Mead') for a0 in (1.3,2.0) for l0 in (3.5,4.5,5.5)),key=lambda r:r.fun); return best.x[0],10**best.x[1],n
EP=[*C.PHASES]; cols=['#9dbfe3','#5b93c9','#2f6fb2','#1b3f6b']
fig,ax=plt.subplots(figsize=(6.4,4.6)); ax.grid(alpha=0.2); ax.spines[['top','right']].set_visible(False); ax.set_xscale('log'); ax.set_yscale('log')
edges=np.logspace(2.5,5.6,14)
for (t0,t1,lab),c in zip(EP,cols):
    M=E['M'][cb&(E['t']>=t0)&(E['t']<t1)]; n,_=np.histogram(M,edges); dl=np.diff(np.log10(edges)); x=np.sqrt(edges[:-1]*edges[1:]); y=n/dl/len(M)
    a,Mc,nn=schechter(M); ok=n>0; ax.plot(x[ok],y[ok],'o',color=c,ms=5,label=f'{lab}:  $\\alpha$ = {a:.2f},  $M_c$ = {Mc:.1e} M$_\\odot$')
    xx=np.logspace(2.5,5.6,100); yy=xx**(1-a)*np.exp(-xx/Mc); yy*=y[ok][0]/np.interp(x[ok][0],xx,yy); ax.plot(xx,yy,color=c,lw=1.5)
ax.set_xlabel('bound cluster mass [M$_\\odot$]'); ax.set_ylabel('dN / dlog M  (normalised)'); ax.set_ylim(1e-4,3); ax.legend(frameon=False,fontsize=8)
fig.savefig(f'{FIG}/mf_phases.png',dpi=150,bbox_inches='tight'); plt.close(fig)
# ---------- Fig C: pressure sources and dynamo
P=dict(np.load('/ptmp/uli/dwarf_merger/prfm/pressure_sources.npz')); T=np.unique(P['snap']); KB=1.3807e-16
def per(f): return np.array([f(P['snap']==s) for s in T])
tt=per(lambda m:P['t'][m][0]); sep=per(lambda m:P['sep'][m][0])
B=per(lambda m:np.sum(np.sqrt(8*np.pi*P['Pmag_2p'][m]*KB)*1e6*P['Sigma_gas'][m])/P['Sigma_gas'][m].sum())
ratio=per(lambda m:np.median(P['Pmag_2p'][m]/np.maximum(P['Pturb_2p'][m],1e-30)))
ffb=per(lambda m:np.clip(np.sum(P['P_fb'][m])/np.sum(P['Ptot_2p'][m]),0,1))
sig=per(lambda m:np.median(np.sqrt(P['Pturb_2p'][m]/4.903e-6/np.maximum(P['rho_mid_2p'][m]*1e9,1e-30))))
ok=tt>10
fig,ax=plt.subplots(4,1,figsize=(7,9),sharex=True,gridspec_kw=dict(hspace=0.12))
for a in ax: a.grid(alpha=0.2); a.spines[['top','right']].set_visible(False); a.axvline(C.T_MERGED,color=ink,lw=0.8,ls=':')
ax[0].plot(tt,sep,color=blue,lw=2); ax[0].set_ylabel('nuclear separation [kpc]')
ax[1].fill_between(tt[ok],0,ffb[ok],color=blue,lw=0,label='feedback: $\\Upsilon_{\\rm OK22}\\,\\Sigma_{\\rm SFR}$'); ax[1].fill_between(tt[ok],ffb[ok],1,color=light,lw=0,label='flow: shear-driven turbulence + field'); ax[1].set_ylim(0,1); ax[1].set_ylabel('share of $P_{\\rm tot}$'); ax[1].legend(frameon=False,fontsize=8,loc='lower left')
ax[2].plot(tt[ok],sig[ok],color=blue,lw=2); ax[2].set_ylabel('$\\sigma_z$ [km s$^{-1}$]'); ax[2].set_ylim(0,14)
ax[3].plot(tt[ok],B[ok],color=blue,lw=2,label='$\\langle B\\rangle$ [$\\mu$G], mass-weighted'); ax[3].plot(tt[ok],ratio[ok],color=ink,lw=1.5,ls='--',label='$P_{\\rm mag}/P_{\\rm turb}$, median patch'); ax[3].set_yscale('log'); ax[3].set_ylim(3e-3,30); ax[3].set_ylabel('dynamo'); ax[3].legend(frameon=False,fontsize=8,loc='lower right'); ax[3].set_xlabel('t [Myr]')
fig.savefig(f'{FIG}/sources_dynamo.png',dpi=150,bbox_inches='tight'); plt.close(fig)
for f in ['ceiling','mf_phases','sources_dynamo']: os.system(f'cp {FIG}/{f}.png ~/out/'); print(f, os.path.getsize(f'{FIG}/{f}.png')//1024,'KB')
