import json
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mp
S=json.load(open('/home/claude/pet/analysis/synteny.json'))
def colour(name):
    n=name.lower()
    if 'mcr' in n or 'eptA' in name or 'PEtN' in name: return '#c1272d'
    if any(k in n for k in ['tnp','transposase','mobile','insertion','is5','is1182']): return '#e08214'
    if any(k in n for k in ['mob','mbe','nik','tra','rep','par','conjug','pi protein','excl']): return '#4575b4'
    return '#c8c8c8'
rows=list(S.items())
H=2.9
fig,ax=plt.subplots(figsize=(7.4,1.45*len(rows)+1.5))
for i,(name,feats) in enumerate(rows):
    y=(len(rows)-i)*H
    xs=[f[0] for f in feats]+[f[1] for f in feats]
    lo,hi=min(xs),max(xs); span=hi-lo
    ax.plot([0,1],[y,y],color='#555555',lw=0.6,zorder=0)
    ax.text(0.0,y+1.18,name,ha='left',va='bottom',fontsize=6.4,fontweight='bold')
    ax.text(1.004,y,'%.1f kb'%(span/1000.0),fontsize=5.2,va='center')
    k=0
    for s,e,st,nm in feats:
        x0=(s-lo)/span; w=(e-s)/span
        if 'transposon' in nm or nm=='mobile element':
            ax.add_patch(plt.Rectangle((x0,y-0.42),w,0.84,fc='#fdf0e0',ec='#e08214',lw=0.7,ls='--',zorder=1))
            ax.text(x0+w/2,y+0.46,nm,ha='center',fontsize=4.6,color='#b35806')
            continue
        c=colour(nm)
        hl=min(w*0.4,0.014)
        if st=='-': ax.arrow(x0+w,y,-w,0,width=0.22,head_width=0.34,head_length=hl,fc=c,ec='none',length_includes_head=True,zorder=2)
        else: ax.arrow(x0,y,w,0,width=0.22,head_width=0.34,head_length=hl,fc=c,ec='none',length_includes_head=True,zorder=2)
        big=(c=='#c1272d')
        off=0.40 if k%2==0 else -0.36
        if big:
            yy=0.42 if not any(('transposon' in f[3] or f[3]=='mobile element') and f[0]<=s and e<=f[1] for f in feats) else 0.90
            ax.text(x0+w/2,y+yy,nm,ha='center',va='bottom',fontsize=6.0,color='#c1272d',fontweight='bold')
        else:
            ax.text(x0+w/2,y+off,nm,ha='center',va=('bottom' if off>0 else 'top'),rotation=0,
                    fontsize=5.0,color='#444444')
        k+=1
ax.set_xlim(-0.01,1.09); ax.set_ylim(-0.4,(len(rows)+0.7)*H); ax.axis('off')
leg=[mp.Patch(fc='#c1272d',label='phosphoethanolamine transferase gene'),
     mp.Patch(fc='#e08214',label='transposase / insertion sequence'),
     mp.Patch(fc='#4575b4',label='plasmid replication, mobilisation, conjugation'),
     mp.Patch(fc='#c8c8c8',label='other coding sequence')]
ax.legend(handles=leg,fontsize=6.4,frameon=False,loc='upper left',ncol=2,bbox_to_anchor=(0.02,0.015))
plt.savefig('/home/claude/pet/analysis/rev/figures/Figure4_synteny.png',dpi=600,bbox_inches='tight')
plt.savefig('/home/claude/pet/analysis/rev/figures/Figure4_synteny.pdf',bbox_inches='tight')
print('ok')
