import dendropy, json, re, collections
meta=json.load(open('metadata.json'))
def fam(lab):
    if lab.startswith('MCR-'): return 'MCR'
    d=meta[lab]['desc'].lower()
    if 'kdo(2)' in d: return 'EptB'
    if 'cpta' in d or 'eptc' in d: return 'CptA/EptC'
    if 'epta' in d: return 'EptA'
    return 'PET-other'
taxa=dendropy.TaxonNamespace()
ml=dendropy.Tree.get(path='tree_ml.nwk',schema='newick',taxon_namespace=taxa,preserve_underscores=True)
boots=dendropy.TreeList.get(path='boot_trees.nwk',schema='newick',taxon_namespace=taxa,preserve_underscores=True)
ml.encode_bipartitions()
counts=collections.Counter()
for t in boots:
    t.encode_bipartitions()
    for b in t.bipartition_encoding: counts[b.split_bitmask]+=1
for nd in ml.internal_nodes():
    if nd.bipartition is not None: nd.label=str(counts.get(nd.bipartition.split_bitmask,0))
ml.reroot_at_midpoint(update_bipartitions=True, suppress_unifurcations=True)
ml.is_rooted=True
ml.write(path='tree_ml_boot.nwk',schema='newick',suppress_rooting=True)
print(collections.Counter(fam(t.label) for t in taxa))
mcr=[t.label for t in taxa if t.label.startswith('MCR-')]
mr=ml.mrca(taxon_labels=mcr)
clade=[l.taxon.label for l in mr.leaf_iter()]
print('\n=== smallest clade with all MCR: %d tips (bootstrap %s) ==='%(len(clade),mr.label))
for c in sorted(clade, key=lambda x:(fam(x),x)): print('  %-12s %s'%(fam(c),c))
# sister group of each mcr family
print('\n=== per-family placement ===')
for m in sorted(mcr, key=lambda x:int(re.search(r'MCR-(\d+)',x).group(1))):
    nd=ml.find_node_with_taxon_label(m)
    p=nd.parent_node; sib=[l.taxon.label for l in p.leaf_iter() if l.taxon.label!=m]
    # climb until sister contains a non-MCR
    while p is not None and all(s.startswith('MCR-') for s in sib) and p.parent_node is not None:
        p=p.parent_node; sib=[l.taxon.label for l in p.leaf_iter() if not l.taxon.label.startswith('MCR-')]
    print('%-14s nearest non-MCR relatives (node bs=%s): %s'%(m,p.label,', '.join(sib[:6]) if sib else 'none'))
