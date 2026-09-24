import json, dendropy
def rd(p):
    d={};n=None;b=[]
    for l in open(p):
        if l[0]=='>':
            if n:d[n]=''.join(b)
            n=l[1:].strip().split()[0];b=[]
        else:b.append(l.strip())
    d[n]=''.join(b);return d
aln=rd('aln_full.faa')
sites=json.load(open('/home/claude/pet/analysis/catalytic_sites.json'))
cols={s['site']:s['col'] for s in sites}
INVAR=['E246','T285','H390','H395','D465','H466']
VAR=['C281','C291','C356','C364','C414','C422','K333','H478','S284','N329']
t=dendropy.Tree.get(path='asr1.treefile',schema='newick',preserve_underscores=True)
t.reroot_at_midpoint(update_bipartitions=False)
mcrset={k for k in aln if k.startswith('MCR-')}
hdr=None; state={}
for line in open('asr1.state'):
    if line.startswith('#'): continue
    f=line.rstrip('\n').split('\t')
    if hdr is None: hdr=f; continue
    node,site,st=f[0],int(f[1])-1,f[2]
    if site in cols.values():
        state.setdefault(node,{})[site]=(st,max(float(x) for x in f[3:]))
# maximal mcr-only clades
groups=[]
for nd in t.postorder_node_iter():
    lv={l.taxon.label for l in nd.leaf_iter()}
    if lv and lv<=mcrset: groups.append((len(lv),nd))
groups.sort(key=lambda x:-x[0])
used=set(); lineages=[]
for sz,nd in groups:
    lv={l.taxon.label for l in nd.leaf_iter()}
    if lv & used: continue
    used|=lv; lineages.append((sorted(x.split('|')[0] for x in lv), nd))
print('independent mobilised lineages in the ML topology:',len(lineages))
out={}
for names,nd in lineages:
    stem=nd.parent_node
    sis={l.taxon.label for l in stem.leaf_iter()} - {l.taxon.label for l in nd.leaf_iter()}
    slab=stem.label
    row={}
    for s in INVAR+VAR:
        st,pp=state.get(slab,{}).get(cols[s],('?',0))
        row[s]=[st,round(pp,3)]
    out['+'.join(names)]={'stem_node':slab,'stem_support':slab,'sister_n':len(sis),
        'sister_examples':[x.split('|')[0] for x in sorted(sis)][:4],'states':row}
print('\nreconstructed states at the stem node of each mobilised lineage (the node it shares with chromosomal relatives)')
print('%-26s %-14s %s'%('lineage','stem support',' '.join('%-6s'%s for s in INVAR+VAR)))
for k,v in out.items():
    print('%-26s %-14s %s'%(k[:26],str(v['stem_node'])[:14],' '.join('%s%-5s'%(v['states'][s][0],('%.2f'%v['states'][s][1]).lstrip('0')) for s in INVAR+VAR)))
json.dump(out,open('out/asr_stem.json','w'),indent=1)
