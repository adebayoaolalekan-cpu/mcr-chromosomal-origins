import json, dendropy, collections
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
ref='MCR-1.1|NG_050417.1'
cols={s['site']:s['col'] for s in sites}
INVAR=['E246','T285','H390','H395','D465','H466']
VAR=[s for s in cols if s not in INVAR]

t=dendropy.Tree.get(path='asr1.treefile',schema='newick',preserve_underscores=True)
# outgroup-free: use midpoint root for reporting mcr-subtending nodes
mcr=[k for k in aln if k.startswith('MCR-')]
lab={}
for nd in t.preorder_node_iter():
    if nd.label: lab[nd.label]=nd
# read states
hdr=None; state={}
for line in open('asr1.state'):
    if line.startswith('#'): continue
    f=line.rstrip('\n').split('\t')
    if hdr is None: hdr=f; continue
    node,site,st=f[0],int(f[1]),f[2]
    if site-1 in cols.values():
        probs=[float(x) for x in f[3:]]
        state.setdefault(node,{})[site-1]=(st,max(probs))
print('nodes with states:',len(state))
aas=hdr[3:]
# identify the parent node of each mcr tip and the node subtending each mcr lineage
tn={l.taxon.label:l for l in t.leaf_node_iter()}
res={}
for m in mcr:
    p=tn[m].parent_node
    nlab=p.label
    if nlab not in state: 
        print('no state for',m,nlab); continue
    row={}
    for s,c in cols.items():
        st,pp=state[nlab].get(c,('?',0))
        row[s]={'state':st,'pp':round(pp,3),'mcr1_residue':aln[ref][c]}
    res[m.split('|')[0]]={'node':nlab,'states':row}
json.dump(res,open('out/asr_primary.json','w'),indent=1)
print('\nancestral states at the node subtending each mobilised family (MCR-1.1 numbering)')
print('%-10s %-8s %s'%('family','node',' '.join('%-6s'%s for s in INVAR+VAR)))
for f,v in res.items():
    print('%-10s %-8s %s'%(f,v['node'],' '.join('%s%-5s'%(v['states'][s]['state'],('%.2f'%v['states'][s]['pp']).lstrip('0')) for s in INVAR+VAR)))
