import json, re, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
meta=json.load(open('metadata.json'))
def readfa(p):
    d={};h=None
    for l in open(p):
        l=l.strip()
        if l.startswith('>'): h=l[1:]; d[h]=''
        else: d[h]+=l
    return d
aln=readfa('aln_full.faa')
def ident(a,b):
    m=s=0
    for x,y in zip(aln[a],aln[b]):
        if x=='-' or y=='-': continue
        s+=1; m+=(x==y)
    return 100.0*m/s
mcr=sorted([k for k in aln if k.startswith('MCR-')],key=lambda x:int(re.search(r'MCR-(\d+)',x).group(1)))
chrom=[k for k in aln if not k.startswith('MCR-')]
M=np.array([[100.0 if a==b else ident(a,b) for b in mcr] for a in mcr])
best=[]
for a in mcr:
    sc=sorted(((ident(a,c),c) for c in chrom),reverse=True)[0]
    best.append((a,sc[0],sc[1]))
fig=plt.figure(figsize=(7.2,3.5))
gs=fig.add_gridspec(1,2,width_ratios=[1.05,1.0],wspace=0.42)
ax=fig.add_subplot(gs[0,0])
im=ax.imshow(M,cmap='YlGnBu',vmin=30,vmax=100)
lab=[m.split('|')[0].replace('MCR-','') for m in mcr]
ax.set_xticks(range(len(lab))); ax.set_xticklabels(lab,rotation=90,fontsize=6)
ax.set_yticks(range(len(lab))); ax.set_yticklabels(lab,fontsize=6)
for i in range(len(lab)):
    for j in range(len(lab)):
        if i!=j: ax.text(j,i,'%.0f'%M[i,j],ha='center',va='center',fontsize=4.0,color=('white' if M[i,j]>70 else '#222222'))
ax.set_title('a  Pairwise amino acid identity\n     between MCR reference alleles (%)',fontsize=7,loc='left')
cb=fig.colorbar(im,ax=ax,fraction=0.045); cb.ax.tick_params(labelsize=5.5)
ax2=fig.add_subplot(gs[0,1])
ys=np.arange(len(best))[::-1]
vals=[b[1] for b in best]
names=[b[0].split('|')[0] for b in best]
orgs=[meta[b[2]]['organism'] for b in best]
cols=['#c1272d' if v>=75 else ('#e08214' if v>=60 else '#4575b4') for v in vals]
ax2.barh(ys,vals,color=cols,height=0.62)
ax2.set_yticks(ys); ax2.set_yticklabels(names,fontsize=6)
for y,v,o in zip(ys,vals,orgs):
    ax2.text(v+1.2,y,'%.1f%%  %s'%(v,o),va='center',fontsize=5.2)
ax2.set_xlim(0,132); ax2.set_xticks([0,20,40,60,80,100]); ax2.tick_params(labelsize=6)
ax2.set_xlabel('identity to closest chromosomal homologue (%)',fontsize=6.5)
ax2.set_title('b  Closest chromosomal relative of each MCR family',fontsize=7,loc='left')
for s in ['top','right']: ax2.spines[s].set_visible(False)
plt.savefig('figures/Figure2_identity.png',dpi=600,bbox_inches='tight')
plt.savefig('figures/Figure2_identity.pdf',bbox_inches='tight')
print('fig2 ok')
