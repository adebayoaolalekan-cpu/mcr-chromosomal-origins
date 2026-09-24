import json, numpy as np, math, collections
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch
meta=json.load(open('metadata.json')); clade=json.load(open('clades.json'))
sites=json.load(open('catalytic_sites.json'))
def readfa(p):
    d={};h=None
    for l in open(p):
        l=l.strip()
        if l.startswith('>'): h=l[1:]; d[h]=''
        else: d[h]+=l
    return d
aln=readfa('aln_full.faa'); raw=readfa('all_pet.faa')
MCR1='MCR-1.1|NG_050417.1'
a=aln[MCR1]; pos2col={}; p=0
for c,ch in enumerate(a):
    if ch!='-': p+=1; pos2col[p]=c
ncol=len(a)
grp={'MCR (n=11)':[k for k in aln if k.startswith('MCR-')],
     'EptA/MCR clade,\nchromosomal (n=65)':[k for k in aln if clade.get(k)=='EptA/MCR' and not k.startswith('MCR-')],
     'EptB clade (n=36)':[k for k in aln if clade.get(k)=='EptB'],
     'CptA/EptC clade (n=49)':[k for k in aln if clade.get(k)=='CptA/EptC']}
def entropy_profile(keys):
    out=[]
    for c in range(ncol):
        col=[aln[k][c] for k in keys if aln[k][c]!='-']
        if len(col)<max(3,0.4*len(keys)): out.append(np.nan); continue
        cnt=collections.Counter(col); n=len(col)
        H=-sum((v/n)*math.log2(v/n) for v in cnt.values())
        out.append(1-H/math.log2(20))
    return np.array(out)
prof={g:entropy_profile(k) for g,k in grp.items()}
KD=dict(A=1.8,R=-4.5,N=-3.5,D=-3.5,C=2.5,Q=-3.5,E=-3.5,G=-0.4,H=-3.2,I=4.5,L=3.8,K=-3.9,M=1.9,F=2.8,P=-1.6,S=-0.8,T=-0.7,W=-0.9,Y=-1.3,V=4.2)
s1=raw[MCR1]
w=19; hyd=[sum(KD.get(c,0) for c in s1[i:i+w])/w for i in range(len(s1)-w+1)]
fig=plt.figure(figsize=(7.2,5.4))
gs=fig.add_gridspec(3,1,height_ratios=[0.9,0.9,2.0],hspace=0.55)
axA=fig.add_subplot(gs[0]); 
axA.add_patch(FancyBboxPatch((1,0.25),185,0.5,boxstyle='round,pad=0.02',fc='#cfe3f5',ec='#1b5e9e',lw=0.8))
axA.add_patch(FancyBboxPatch((190,0.25),len(s1)-190,0.5,boxstyle='round,pad=0.02',fc='#f7dfd0',ec='#b3541e',lw=0.8))
axA.text(93,0.5,'N-terminal transmembrane domain',ha='center',va='center',fontsize=6.5)
axA.text(365,0.5,'periplasmic catalytic (sulfatase-like) domain',ha='center',va='center',fontsize=6.5)
for i,s in enumerate(sorted(sites,key=lambda z:int(z['site'][1:]))):
    pos=int(s['site'][1:])
    col='#c1272d' if s['site'] in ('E246','T285','H390','H395','D465','H466','H478') else '#6a3d9a'
    top=0.95 if i%2==0 else 1.16
    axA.plot([pos,pos],[0.78,top],color=col,lw=0.9)
    axA.text(pos,top+0.03,s['site'],rotation=90,fontsize=4.4,ha='center',va='bottom',color=col)
axA.set_xlim(0,len(s1)+6); axA.set_ylim(0,2.0); axA.axis('off')
axA.text(0,1.85,'a  Domain architecture of MCR-1.1 (541 aa) with catalytic (red) and disulfide-forming (purple) residues',fontsize=7)
axB=fig.add_subplot(gs[1])
axB.plot(range(1+w//2,len(hyd)+1+w//2),hyd,color='#1b5e9e',lw=0.8)
axB.axhline(1.6,color='#999999',lw=0.6,ls='--')
axB.fill_between(range(1+w//2,len(hyd)+1+w//2),1.6,hyd,where=(np.array(hyd)>1.6),color='#cfe3f5')
axB.set_xlim(0,len(s1)+6); axB.set_ylabel('Kyte–Doolittle\nhydropathy',fontsize=6)
axB.tick_params(labelsize=6); axB.set_title('b  Hydropathy profile (window = 19)',fontsize=7,loc='left')
for s in ['top','right']: axB.spines[s].set_visible(False)
axC=fig.add_subplot(gs[2])
cols={'MCR (n=11)':'#c1272d','EptA/MCR clade,\nchromosomal (n=65)':'#1b5e9e','EptB clade (n=36)':'#2e7d32','CptA/EptC clade (n=49)':'#6a3d9a'}
xs=np.arange(1,len(s1)+1)
for g,pr in prof.items():
    v=np.array([pr[pos2col[i]] for i in xs])
    k=15; sm=np.convolve(np.nan_to_num(v),np.ones(k)/k,mode='same')
    axC.plot(xs,sm,color=cols[g],lw=0.9,label=g)
for s in sites:
    pos=int(s['site'][1:])
    axC.axvline(pos,color='#dddddd',lw=0.5,zorder=0)
axC.set_xlim(0,len(s1)+6); axC.set_ylim(0,1.0)
axC.set_xlabel('position in MCR-1.1',fontsize=7); axC.set_ylabel('relative conservation\n(1 − normalised entropy, 15-residue mean)',fontsize=6)
axC.legend(fontsize=5.5,frameon=False,ncol=2,loc='upper left')
axC.tick_params(labelsize=6)
axC.set_title('c  Positional conservation within each clade, mapped to MCR-1.1 coordinates',fontsize=7,loc='left')
for s in ['top','right']: axC.spines[s].set_visible(False)
plt.savefig('figures/Figure3_domains.png',dpi=600,bbox_inches='tight')
plt.savefig('figures/Figure3_domains.pdf',bbox_inches='tight')
print('fig3 ok')
