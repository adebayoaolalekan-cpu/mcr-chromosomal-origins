import json, statistics, dendropy, collections
res=json.load(open('out/reservoir_v2.json'))
L=[r['length'] for r in res]
print('chromosomal/related n=%d  median length %d  mean %.1f'%(len(L),statistics.median(L),statistics.mean(L)))
low=[r for r in res if r['tm']<3]
print('with fewer than 3 predicted segments: %d  lengths %s'%(len(low),sorted(r['length'] for r in low)))
print('with tm==1: %d  %s'%(sum(1 for r in res if r['tm']==1),[ (r['acc'],r['length']) for r in res if r['tm']==1]))
cl=json.load(open('out/closest_v2.json'))
for k in ['MCR-8.1|NG_061399.1','MCR-9.1|MK070339','MCR-12.1|NG_245195.1']:
    v=cl[k]
    print(k.split('|')[0],'full-only %s %.1f | incl-partial %s %.1f'%(v['best_full_only']['hit'].split('|')[0],v['best_full_only']['local'],v['best_incl_partial']['hit'].split('|')[0],v['best_incl_partial']['local']))
sp=json.load(open('out/cptA_split.json'))
t=dendropy.Tree.get(path='v2.treefile',schema='newick',preserve_underscores=True)
tips={l.taxon.label:l for l in t.leaf_node_iter()}
op=tips['OpgE_ECOLI_P75785']; nd=op.parent_node
rows=[]
while nd is not None:
    lv=[l.taxon.label for l in nd.leaf_iter()]
    n_op=sum(1 for x in lv if sp.get(x)=='OpgE-like')
    n_cp=sum(1 for x in lv if sp.get(x) and sp[x]!='OpgE-like')
    rows.append((len(lv),n_op,n_cp,nd.label))
    nd=nd.parent_node
print('\nclades containing OpgE, from the tip outwards (total, OpgE-like, CptA/EptC-like, support):')
for r in rows[:9]: print('   ',r)
print('\nannotated in the dataset: OpgE-like %d, CptA/EptC-like %d'%(sum(1 for v in sp.values() if v=='OpgE-like'),sum(1 for v in sp.values() if v!='OpgE-like')))
