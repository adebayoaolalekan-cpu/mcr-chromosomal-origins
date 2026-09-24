import dendropy, json, collections, re
meta=json.load(open('metadata.json'))
taxa=dendropy.TaxonNamespace()
t=dendropy.Tree.get(path='tree_ml.nwk',schema='newick',taxon_namespace=taxa,preserve_underscores=True)
t.reroot_at_midpoint(update_bipartitions=True,suppress_unifurcations=True)
EPTA='NP_418538.2|Escherichia_coli_str_K-12_substr_MG1655'
EPTB='NP_418002.2|Escherichia_coli_str_K-12_substr_MG1655'
EPTC='NP_418390.1|Escherichia_coli_str_K-12_substr_MG1655'
mcr=[l.taxon.label for l in t.leaf_node_iter() if l.taxon.label.startswith('MCR-')]
def grow(anchor, excl):
    nd=t.find_node_with_taxon_label(anchor)
    best=nd
    while nd.parent_node is not None:
        nd=nd.parent_node
        labs=[l.taxon.label for l in nd.leaf_iter()]
        if any(e in labs for e in excl): break
        best=nd
    return [l.taxon.label for l in best.leaf_iter()], best
cA,ndA=grow(EPTA,[EPTB,EPTC])
cB,ndB=grow(EPTB,[EPTA,EPTC])
cC,ndC=grow(EPTC,[EPTA,EPTB])
print('EptA/MCR clade: %d tips (%d MCR, %d chromosomal)'%(len(cA),sum(1 for x in cA if x.startswith('MCR-')),sum(1 for x in cA if not x.startswith('MCR-'))))
print('EptB clade    : %d tips'%len(cB))
print('CptA/EptC clade: %d tips'%len(cC))
assigned=set(cA)|set(cB)|set(cC)
rest=[l.taxon.label for l in t.leaf_node_iter() if l.taxon.label not in assigned]
print('unassigned    : %d tips'%len(rest))
print('  MCR in EptA clade:',sorted(x.split("|")[0] for x in cA if x.startswith('MCR-')))
print('  MCR outside     :',sorted(x.split("|")[0] for x in mcr if x not in cA))
cl={}
for x in cA: cl[x]='EptA/MCR'
for x in cB: cl[x]='EptB'
for x in cC: cl[x]='CptA/EptC'
for x in rest: cl[x]='unassigned'
json.dump(cl,open('clades.json','w'),indent=1)
# genus composition of EptA/MCR clade
def org(l): return meta[l]['organism']
gen=collections.Counter(org(x).split()[0] for x in cA if not x.startswith('MCR-'))
print('\nEptA/MCR clade chromosomal members span %d genus-level labels'%len(gen))
print(sorted(gen))
