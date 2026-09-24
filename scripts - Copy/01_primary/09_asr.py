import numpy as np, json, re, collections, dendropy
from scipy.stats import gamma as gdist
import pyvolve.empirical_matrices as em
AA='ARNDCQEGHILKMFPSTWYV'
# pyvolve matrices use alphabetical order ACDEFGHIKLMNPQRSTVWY
PV='ACDEFGHIKLMNPQRSTVWY'
R=np.array(em.lg_matrix,dtype=float); pi=np.array(em.lg_freqs,dtype=float)
idx=[PV.index(a) for a in AA]
R=R[np.ix_(idx,idx)]; pi=pi[idx]; pi=pi/pi.sum()
R=(R+R.T)/2
Q=R*pi[None,:]
np.fill_diagonal(Q,0); np.fill_diagonal(Q,-Q.sum(1))
mu=-(np.diag(Q)*pi).sum(); Q=Q/mu
ev,V=np.linalg.eig(Q); Vi=np.linalg.inv(V)
def P(t):
    return np.real(V@np.diag(np.exp(ev*t))@Vi)
ALPHA=1.385; K=4
qs=[(2*i+1)/(2*K) for i in range(K)]
rates=np.array([gdist.ppf(q,a=ALPHA,scale=1/ALPHA) for q in qs]); rates=rates/rates.mean()
def readfa(p):
    d={};h=None
    for l in open(p):
        l=l.strip()
        if l.startswith('>'): h=l[1:]; d[h]=''
        else: d[h]+=l
    return d
aln=readfa('aln_full.faa')
taxa=dendropy.TaxonNamespace()
tree=dendropy.Tree.get(path='tree_ml.nwk',schema='newick',taxon_namespace=taxa,preserve_underscores=True)
tree.reroot_at_midpoint(update_bipartitions=False, suppress_unifurcations=True)
tree.is_rooted=True
ncol=len(next(iter(aln.values())))
tips=[n for n in tree.leaf_node_iter()]
for n in tips: assert n.taxon.label in aln, n.taxon.label
def tip_vec(lab):
    s=aln[lab]; M=np.ones((ncol,20))
    for c,ch in enumerate(s):
        if ch in AA:
            M[c]=0; M[c][AA.index(ch)]=1
    return M
TIP={n.taxon.label:tip_vec(n.taxon.label) for n in tips}
nodes=list(tree.postorder_node_iter())
def down_pass(rate):
    D={}
    for n in nodes:
        if n.is_leaf(): D[id(n)]=TIP[n.taxon.label].copy()
        else:
            v=np.ones((ncol,20))
            for ch in n.child_node_iter():
                t=(ch.edge.length or 1e-8)*rate
                v*= D[id(ch)]@P(t).T
            D[id(n)]=v
    return D
def up_pass(D,rate):
    U={}
    root=tree.seed_node
    U[id(root)]=np.tile(pi,(ncol,1))
    for n in tree.preorder_node_iter():
        if n is root: continue
        p=n.parent_node
        v=U[id(p)].copy()
        for sib in p.child_node_iter():
            if sib is n: continue
            t=(sib.edge.length or 1e-8)*rate
            v*= D[id(sib)]@P(t).T
        t=(n.edge.length or 1e-8)*rate
        U[id(n)]=v@P(t)
    return U
post=collections.defaultdict(lambda: np.zeros((ncol,20)))
wsum=np.zeros(ncol)
for rate in rates:
    D=down_pass(rate); U=up_pass(D,rate)
    sitelk=(D[id(tree.seed_node)]*np.tile(pi,(ncol,1))).sum(1)
    w=sitelk/K
    wsum+=w
    for n in nodes:
        if n.is_leaf(): continue
        m=D[id(n)]*U[id(n)]
        m=m/m.sum(1,keepdims=True)
        post[id(n)]+= m*w[:,None]
for k in post: post[k]=post[k]/wsum[:,None]
# node selection
def mrca(labels): return tree.mrca(taxon_labels=labels)
targets={}
pairs={'mcr-1/2/6 ancestor':['MCR-1.1|NG_050417.1','MCR-2.1|NG_051171.1','MCR-6.1|NG_055781.1'],
       'mcr-1/2/6 + Moraxella donor':['MCR-1.1|NG_050417.1','MCR-2.1|NG_051171.1','MCR-6.1|NG_055781.1','WP_490288700.1|Faucicola_osloensis','WP_284034866.1|Moraxellaceae'],
       'mcr-3/7 ancestor':['MCR-3.1|NG_055505.1','MCR-7.1|NG_056413.1'],
       'mcr-9/10 ancestor':['MCR-9.1|MK070339','MCR-10.1|NG_066767.1'],
       'mcr-4 + Shewanella donor':['MCR-4.1|NG_057470.1','WP_220053229.1|Shewanella'],
       'mcr-5 + Cupriavidus donor':['MCR-5.1|NG_055658.1','WP_227311958.1|Cupriavidus_sp_MP-37','WP_454751270.1|Cupriavidus'],
       'EptA/MCR clade ancestor':[l.taxon.label for l in tree.leaf_node_iter() if l.taxon.label.startswith('MCR-')]}
sites=json.load(open('catalytic_sites.json'))
print('%-30s %-6s %s'%('node','ntips','ancestral residues at catalytic/structural sites (posterior)'))
out={}
for name,labs in pairs.items():
    nd=mrca(labs); m=post[id(nd)]
    ntips=len(list(nd.leaf_iter()))
    txt=[]
    for s in sites:
        c=s['col']; j=int(np.argmax(m[c])); txt.append('%s:%s(%.2f)'%(s['site'],AA[j],m[c][j]))
    print('%-30s %-6d %s'%(name,ntips,' '.join(txt[:8])))
    seq=''.join(AA[int(np.argmax(m[c]))] for c in range(ncol) if m[c].max()>0)
    mean_pp=float(np.mean(m.max(1)))
    out[name]=dict(ntips=ntips,mean_posterior=round(mean_pp,3),
                   sites={s['site']:[AA[int(np.argmax(m[s['col']]))],round(float(m[s['col']].max()),3)] for s in sites})
    print('    mean marginal posterior over %d columns: %.3f'%(ncol,mean_pp))
json.dump(out,open('asr_nodes.json','w'),indent=1)
# ancestral vs extant MCR-1: differences in catalytic domain
nd=mrca(pairs['mcr-1/2/6 + Moraxella donor']); m=post[id(nd)]
anc=[AA[int(np.argmax(m[c]))] for c in range(ncol)]
m1=aln['MCR-1.1|NG_050417.1']
diff=sum(1 for c in range(ncol) if m1[c]!='-' and anc[c]!=m1[c])
tot=sum(1 for c in range(ncol) if m1[c]!='-')
print('\nMCR-1.1 differs from the reconstructed pre-mobilisation ancestor at %d of %d aligned positions (%.1f%%)'%(diff,tot,100*diff/tot))
