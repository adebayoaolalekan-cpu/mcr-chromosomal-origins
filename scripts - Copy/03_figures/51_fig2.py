import json, re, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
def readfa(p):
    d={};h=None
    for l in open(p):
        l=l.strip()
        if l.startswith('>'): h=l[1:].split()[0]; d[h]=''
        else: d[h]+=l
    return d
aln=readfa('v2_aln_full.faa')
def ident(a,b):
    m=s=0
    for x,y in zip(aln[a],aln[b]):
        if x=='-' or y=='-': continue
        s+=1; m+=(x==y)
    return 100.0*m/s
mcr=sorted([k for k in aln if k.startswith('MCR-')],key=lambda x:int(re.search(r'MCR-(\d+)',x).group(1)))
M=np.array([[100.0 if a==b else ident(a,b) for b in mcr] for a in mcr])
NEW=[('MCR-1.1',98.7,62.4,'Moraxella sp. MSG13C03'),('MCR-2.1',81.0,64.1,'Moraxella sp. MSG13C03'),
     ('MCR-3.1',95.0,66.9,'Aeromonas media'),('MCR-4.1',81.5,80.2,'Shewanella baltica'),
     ('MCR-5.1',45.1,45.1,'Cupriavidus sp.'),('MCR-6.1',82.3,65.1,'Moraxella sp. MSG13C03'),
     ('MCR-7.1',78.3,78.3,'Aeromonas jandaei'),('MCR-8.1',43.8,47.9,'Enterobacterales'),
     ('MCR-9.1',87.4,99.0,'Buttiauxella cochleicola'),('MCR-10.1',84.6,84.6,'Buttiauxella cochleicola'),
     ('MCR-12.1',40.6,41.0,'Pseudomonas sp. RC25-4')]
fig=plt.figure(figsize=(7.4,3.8))
gs=fig.add_gridspec(1,2,width_ratios=[1.02,1.12],wspace=0.46)
ax=fig.add_subplot(gs[0,0])
im=ax.imshow(M,cmap='YlGnBu',vmin=30,vmax=100)
lab=[m.split('|')[0].replace('MCR-','') for m in mcr]
ax.set_xticks(range(len(lab))); ax.set_xticklabels(lab,rotation=90,fontsize=6.5)
ax.set_yticks(range(len(lab))); ax.set_yticklabels(lab,fontsize=6.5)
for i in range(len(lab)):
    for j in range(len(lab)):
        if i!=j: ax.text(j,i,'%.0f'%M[i,j],ha='center',va='center',fontsize=4.4,color=('white' if M[i,j]>70 else '#222222'))
ax.set_title('a  Pairwise amino acid identity between\n     the MCR reference alleles (%)',fontsize=7.5,loc='left')
cb=fig.colorbar(im,ax=ax,fraction=0.045); cb.ax.tick_params(labelsize=6)
ax2=fig.add_subplot(gs[0,1])
ys=np.arange(len(NEW))[::-1]
new=[v[1] for v in NEW]; old=[v[2] for v in NEW]
ax2.barh(ys,old,color='#cfd8dc',height=0.70,label='clustered dataset')
ax2.barh(ys,new,color='#1b5e9e',height=0.42,label='unclustered dataset')
ax2.set_yticks(ys); ax2.set_yticklabels([v[0] for v in NEW],fontsize=6.5)
for y,v,o in zip(ys,NEW,[v[3] for v in NEW]):
    ax2.text(max(v[1],v[2])+1.4,y,'%.1f%%  %s'%(v[1],o),va='center',fontsize=5.4,style='italic')
ax2.set_xlim(0,148); ax2.set_xticks([0,20,40,60,80,100]); ax2.tick_params(labelsize=6.5)
ax2.set_xlabel('identity to the closest full-length chromosomal relative (%)',fontsize=6.8)
ax2.set_title('b  Closest chromosomal relative of each family,\n     before and after removal of the clustering step',fontsize=7.5,loc='left')
for s in ['top','right']: ax2.spines[s].set_visible(False)
ax2.legend(fontsize=6.2,frameon=False,loc='upper left',bbox_to_anchor=(0.0,-0.12),ncol=2)
ax2.annotate('99.0% against a 400-residue\npartial sequence',xy=(99,ys[8]),xytext=(104,ys[8]+1.15),fontsize=5.0,color='#8a1c1c',
             arrowprops=dict(arrowstyle='-',color='#8a1c1c',lw=0.5))
plt.savefig('figures/Figure2_identity.png',dpi=600,bbox_inches='tight')
plt.savefig('figures/Figure2_identity.pdf',bbox_inches='tight')
print('fig2 ok')
