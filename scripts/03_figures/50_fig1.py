import dendropy, json, collections, re
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
meta=json.load(open('/home/claude/pet/analysis/metadata.json'))
clade=json.load(open('/home/claude/pet/analysis/clades.json'))
res=json.load(open('out/reservoir_v2.json'))
KEEP={r['acc'] for r in res[:24]}
COL={'EptA/MCR':'#1b5e9e','EptB':'#2e7d32','CptA/EptC':'#6a3d9a'}
MCRC='#c1272d'
t=dendropy.Tree.get(path='ml.treefile',schema='newick',preserve_underscores=True)
t.reroot_at_midpoint(update_bipartitions=False,suppress_unifurcations=True)
def keeptip(lab):
    return lab.startswith('MCR-') or lab.split('|')[0] in KEEP
# collapse maximal clades with no kept tip and >=2 leaves
collapsed=[]
def mark(nd):
    lv=[l.taxon.label for l in nd.leaf_iter()]
    if not nd.is_leaf() and len(lv)>=2 and not any(keeptip(x) for x in lv):
        collapsed.append(nd); return True
    for c in list(nd.child_node_iter()): mark(c)
    return False
mark(t.seed_node)
cset={id(n) for n in collapsed}
def label_for(nd):
    orgs=[meta[l.taxon.label]['organism'].split()[0] for l in nd.leaf_iter() if l.taxon.label in meta]
    c=collections.Counter(orgs).most_common(2)
    nm=c[0][0] if c else 'mixed'
    if len(c)>1 and c[1][1]>=c[0][1]: nm=nm+' and others'
    elif len(set(orgs))>2: nm=nm+' and others'
    return '%s (%d sequences)'%(nm,len(list(nd.leaf_iter())))
rows=[]
def walk(nd):
    if id(nd) in cset: rows.append(('clade',nd)); return
    if nd.is_leaf(): rows.append(('tip',nd)); return
    for c in nd.child_node_iter(): walk(c)
walk(t.seed_node)
y={}
for i,(k,nd) in enumerate(rows): y[id(nd)]=i
def ypos(n):
    if id(n) in y: return y[id(n)]
    v=[ypos(c) for c in n.child_node_iter()]; y[id(n)]=sum(v)/len(v); return y[id(n)]
ypos(t.seed_node)
x={}
def xpos(n,acc=0.0):
    x[id(n)]=acc
    if id(n) in cset: return
    for c in n.child_node_iter(): xpos(c,acc+(c.edge.length or 0))
xpos(t.seed_node)
def depth(nd):
    return max((c.edge.length or 0)+depth(c) for c in nd.child_node_iter()) if not nd.is_leaf() else 0.0
fig,ax=plt.subplots(figsize=(7.6,len(rows)*0.135+1.6))
def draw(n):
    if n.parent_node is not None:
        ax.plot([x[id(n.parent_node)],x[id(n)]],[y[id(n)],y[id(n)]],color='#2b2b2b',lw=0.8)
    if id(n) in cset: return
    ch=list(n.child_node_iter())
    if ch:
        ax.plot([x[id(n)],x[id(n)]],[y[id(ch[0])],y[id(ch[-1])]],color='#2b2b2b',lw=0.8)
        for c in ch: draw(c)
draw(t.seed_node)
xmax=max(x[id(nd)]+ (depth(nd) if id(nd) in cset else 0) for k,nd in rows)
for k,nd in rows:
    yy=y[id(nd)]
    if k=='clade':
        d=depth(nd); x0=x[id(nd)]
        ax.add_patch(Polygon([[x0,yy],[x0+d,yy-0.42],[x0+d,yy+0.42]],closed=True,facecolor='#c9d4de',edgecolor='#5b6b78',lw=0.5))
        ax.plot([x0+d,xmax*1.02],[yy,yy],color='#dddddd',lw=0.3,ls=':')
        ax.text(xmax*1.03,yy,label_for(nd),fontsize=5.6,va='center',color='#444444',style='italic')
    else:
        lab=nd.taxon.label; acc,org=lab.split('|',1)
        ismcr=lab.startswith('MCR-')
        col=MCRC if ismcr else COL.get(clade.get(lab,''),'#777777')
        txt=(acc) if ismcr else (acc+'  '+org.replace('_',' '))
        ax.plot([x[id(nd)],xmax*1.02],[yy,yy],color='#dddddd',lw=0.3,ls=':')
        ax.text(xmax*1.03,yy,txt,fontsize=5.9,va='center',color=col,fontweight=('bold' if ismcr else 'normal'))
for n in t.preorder_node_iter():
    if n.is_leaf() or id(n) in cset or n.parent_node is None or not n.label: continue
    if id(n) not in x or id(n) not in y: continue
    m=re.match(r'([\d.]+)/([\d.]+)',n.label)
    if not m: continue
    sh,bb=float(m.group(1)),float(m.group(2))
    if sh>=80 and bb>=95:
        ax.plot(x[id(n)],y[id(n)],'o',ms=2.8,color='#000000',mec='none')
    elif sh>=70 or bb>=85:
        ax.plot(x[id(n)],y[id(n)],'o',ms=1.8,color='#8a8a8a',mec='none')
for name,col in COL.items():
    idxs=[y[id(nd)] for k,nd in rows for l in nd.leaf_iter() if clade.get(l.taxon.label)==name]
    if not idxs: continue
    xb=xmax*2.02
    lo,hi=min(idxs),max(idxs)
    ax.plot([xb,xb],[lo-0.4,hi+0.4],color=col,lw=2.6,solid_capstyle='butt')
    if hi-lo>=3:
        ax.text(xb*1.012,(lo+hi)/2,name,rotation=90,va='center',ha='left',fontsize=7.5,color=col)
    else:
        ax.text(xb*1.012,(lo+hi)/2,name,va='center',ha='left',fontsize=7.5,color=col)
ax.set_ylim(-4,len(rows)+1); ax.set_xlim(-0.02,xmax*2.55); ax.axis('off')
ax.plot([0,0.2],[-1.6,-1.6],color='black',lw=1.2)
ax.text(0.1,-2.9,'0.2 substitutions per site',ha='left',fontsize=7)
handles=[Line2D([],[],color=MCRC,lw=3,label='mobilised MCR reference alleles'),
         Line2D([],[],color='#1b5e9e',lw=3,label='chromosomal members of the EptA/MCR clade'),
         Line2D([],[],color='#2e7d32',lw=3,label='EptB clade'),
         Line2D([],[],color='#6a3d9a',lw=3,label='CptA/EptC and OpgE-like clade'),
         Line2D([],[],marker='o',color='black',lw=0,label='SH-aLRT ≥ 80 and ultrafast bootstrap ≥ 95',ms=4),
         Line2D([],[],marker='o',color='#8a8a8a',lw=0,label='SH-aLRT ≥ 70 or ultrafast bootstrap ≥ 85',ms=3),
         Line2D([],[],marker='>',color='#5b6b78',lw=0,label='collapsed chromosomal subtree',ms=5)]
ax.legend(handles=handles,loc='upper left',bbox_to_anchor=(0.02,-0.01),fontsize=7,frameon=False,ncol=2,columnspacing=1.2)
plt.savefig('figures/Figure1_phylogeny.png',dpi=600,bbox_inches='tight')
plt.savefig('figures/Figure1_phylogeny.pdf',bbox_inches='tight')
print('rows drawn:',len(rows),' collapsed clades:',len(collapsed))
