import json, dendropy
def rd(p):
    d={};n=None;b=[]
    for l in open(p):
        if l[0]=='>':
            if n:d[n]=''.join(b)
            n=l[1:].strip().split()[0];b=[]
        else:b.append(l.strip())
    d[n]=''.join(b);return d
alnA=rd('/home/claude/pet/analysis/aln_full.faa'); alnB=rd('v2_aln_full.faa')
sites=json.load(open('/home/claude/pet/analysis/catalytic_sites.json'))
ref='MCR-1.1|NG_050417.1'
def col2res(seq,c): return sum(1 for x in seq[:c+1] if x!='-')-1
def res2col(seq,ri):
    c=0
    for i,ch in enumerate(seq):
        if ch!='-':
            if c==ri: return i
            c+=1
mapB={}
for s in sites:
    ri=col2res(alnA[ref],s['col']); mapB[s['site']]=res2col(alnB[ref],ri)
VAR=['C281','C291','C356','C364','C414','C422','K333','H478','S284','N329']
t=dendropy.Tree.get(path='asr2.treefile',schema='newick',preserve_underscores=True)
mcr={k for k in alnB if k.startswith('MCR-')}
state={}
hdr=None
for line in open('asr2.state'):
    if line.startswith('#'): continue
    f=line.rstrip('\n').split('\t')
    if hdr is None: hdr=f; continue
    c=int(f[1])-1
    if c in mapB.values(): state.setdefault(f[0],{})[c]=(f[2],max(float(x) for x in f[3:]))
groups=[]
for nd in t.postorder_node_iter():
    lv={l.taxon.label for l in nd.leaf_iter()}
    if lv and lv<=mcr: groups.append((len(lv),nd))
groups.sort(key=lambda x:-x[0])
used=set(); lin=[]
for sz,nd in groups:
    lv={l.taxon.label for l in nd.leaf_iter()}
    if lv&used: continue
    used|=lv; lin.append((sorted(x.split('|')[0] for x in lv),nd))
print('lineages:',len(lin))
prev=json.load(open('out/asr_stem.json'))
out={}
print('%-26s %s'%('lineage',' '.join('%-6s'%s for s in VAR)))
for names,nd in lin:
    stem=nd.parent_node; lab=stem.label if stem is not None else None
    row={}
    for s in VAR:
        st,pp=state.get(lab,{}).get(mapB[s],('?',0)); row[s]=[st,round(pp,3)]
    out['+'.join(names)]={'stem_node':lab,'states':row}
    print('%-26s %s'%('+'.join(names)[:26],' '.join('%s%-5s'%(row[s][0],('%.2f'%row[s][1]).lstrip('0')) for s in VAR)))
json.dump(out,open('out/asr_stem_rooted.json','w'),indent=1)
# agreement with the primary analysis for the eight cysteine/K333/H478 positions
KEY=['C281','C291','C356','C364','C414','C422','K333','H478']
agree=tot=0
for k,v in prev.items():
    kk=k if k in out else None
    if kk is None:
        for cand in out:
            if set(cand.split('+'))<=set(k.split('+')) or set(k.split('+'))<=set(cand.split('+')): kk=cand; break
    if kk is None: continue
    for s in KEY:
        tot+=1; agree+= (v['states'][s][0]==out[kk]['states'][s][0])
print('\nagreement at the eight informative active-site positions between the two rootings: %d of %d'%(agree,tot))
