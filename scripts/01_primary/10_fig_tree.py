import dendropy, json, re, collections
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
meta=json.load(open('metadata.json')); clade=json.load(open('clades.json'))
COL={'EptA/MCR':'#1b5e9e','EptB':'#2e7d32','CptA/EptC':'#6a3d9a'}
MCRC='#c1272d'
DONORS={'WP_490288700.1','WP_284034866.1','WP_220053229.1','WP_499669208.1','WP_508945341.1','WP_499756849.1','WP_261396260.1','WP_232921000.1','WP_227311958.1','WP_454751270.1','WP_507779059.1','WP_463893905.1'}
taxa=dendropy.TaxonNamespace()
t=dendropy.Tree.get(path='tree_ml.nwk',schema='newick',taxon_namespace=taxa,preserve_underscores=True)
b=dendropy.TreeList.get(path='boot_trees.nwk',schema='newick',taxon_namespace=taxa,preserve_underscores=True)
t.encode_bipartitions(); cnt=collections.Counter()
for x in b:
    x.encode_bipartitions()
    for bp in x.bipartition_encoding: cnt[bp.split_bitmask]+=1
sup={id(nd):cnt.get(nd.bipartition.split_bitmask,0) for nd in t.internal_nodes() if nd.bipartition is not None}
t.reroot_at_midpoint(update_bipartitions=False,suppress_unifurcations=True)
leaves=list(t.leaf_node_iter()); y={}
for i,l in enumerate(leaves): y[id(l)]=i
def ypos(n):
    if n.is_leaf(): return y[id(n)]
    v=[ypos(c) for c in n.child_node_iter()]; y[id(n)]=sum(v)/len(v); return y[id(n)]
ypos(t.seed_node)
x={}
def xpos(n,acc=0.0):
    x[id(n)]=acc
    for c in n.child_node_iter(): xpos(c,acc+(c.edge.length or 0))
xpos(t.seed_node)
fig,ax=plt.subplots(figsize=(7.4,14.0))
for n in t.preorder_node_iter():
    if n.parent_node is not None:
        ax.plot([x[id(n.parent_node)],x[id(n)]],[y[id(n)],y[id(n)]],color='#2b2b2b',lw=0.65)
    ch=list(n.child_node_iter())
    if ch: ax.plot([x[id(n)],x[id(n)]],[y[id(ch[0])],y[id(ch[-1])]],color='#2b2b2b',lw=0.65)
xmax=max(x.values())
for l in leaves:
    lab=l.taxon.label; acc,org=lab.split('|',1)
    ismcr=lab.startswith('MCR-')
    col=MCRC if ismcr else COL.get(clade.get(lab,''),'#777777')
    txt=(acc+'  ('+meta[lab]['acc']+')') if ismcr else (acc+'  '+org.replace('_',' '))
    if (not ismcr) and acc in DONORS: txt+='  ◀'
    ax.plot([x[id(l)],xmax*1.015],[y[id(l)],y[id(l)]],color='#d8d8d8',lw=0.3,ls=':')
    ax.text(xmax*1.025,y[id(l)],txt,fontsize=4.1,va='center',color=col,
            fontweight=('bold' if ismcr or acc in DONORS else 'normal'))
for n in t.internal_nodes():
    s=sup.get(id(n))
    if s is not None and s>=70 and n.parent_node is not None:
        ax.plot(x[id(n)],y[id(n)],'o',ms=(1.7 if s<90 else 2.7),color=('#7a7a7a' if s<90 else '#000000'),mec='none')
# clade brackets
for name,col in COL.items():
    idxs=[y[id(l)] for l in leaves if clade.get(l.taxon.label)==name]
    if not idxs: continue
    xb=xmax*1.40
    ax.plot([xb,xb],[min(idxs),max(idxs)],color=col,lw=2.4,solid_capstyle='butt')
    ax.text(xb*1.012,(min(idxs)+max(idxs))/2,name,rotation=90,va='center',ha='left',fontsize=7,color=col)
ax.set_ylim(-4,len(leaves)+1); ax.set_xlim(-0.02,xmax*1.55); ax.axis('off')
ax.plot([0,0.2],[-2.0,-2.0],color='black',lw=1.2)
ax.text(0.1,-3.4,'0.2 substitutions per site',ha='center',fontsize=6.5)
handles=[Line2D([],[],color=MCRC,lw=3,label='mobilised MCR reference alleles'),
         Line2D([],[],color='#1b5e9e',lw=3,label='chromosomal members, EptA/MCR clade'),
         Line2D([],[],color='#2e7d32',lw=3,label='EptB clade'),
         Line2D([],[],color='#6a3d9a',lw=3,label='CptA/EptC clade'),
         Line2D([],[],marker='o',color='black',lw=0,label='bootstrap ≥ 90%',ms=4),
         Line2D([],[],marker='o',color='#7a7a7a',lw=0,label='bootstrap 70–89%',ms=3),
         Line2D([],[],marker='<',color='black',lw=0,label='closest chromosomal relative of an MCR family',ms=4)]
ax.legend(handles=handles,loc='upper left',bbox_to_anchor=(0.0,0.14),fontsize=6,frameon=False)
plt.savefig('figures/Figure1_phylogeny.png',dpi=600,bbox_inches='tight')
plt.savefig('figures/Figure1_phylogeny.pdf',bbox_inches='tight')
print('ok')
