import json
def rd(p):
    d={};n=None;b=[]
    for l in open(p):
        if l[0]=='>':
            if n:d[n]=''.join(b)
            n=l[1:].strip().split()[0];b=[]
        else:b.append(l.strip())
    d[n]=''.join(b);return d
aln=rd('v2_aln_full.faa'); raw=rd('v2_dataset.faa')
meta=json.load(open('/home/claude/pet/analysis/metadata.json'))
sites=json.load(open('/home/claude/pet/analysis/catalytic_sites.json'))
mcr=[k for k in aln if k.startswith('MCR-')]
outs={'OpgB_ECOLI_P39401','AslA_ECOLI_P25549','YidJ_ECOLI_P31447','YdeN_ECOLI_P77318'}
chrom=[k for k in aln if k not in mcr and k not in outs]

# map catalytic columns from old aln_full (808 cols) onto v2_aln_full via MCR-1.1 residue index
old=rd('/home/claude/pet/analysis/aln_full.faa')
ref='MCR-1.1|NG_050417.1'
def col2res(seq,col):
    if seq[col]=='-': return None
    return sum(1 for c in seq[:col+1] if c!='-')-1
def res2col(seq,ri):
    c=0
    for i,ch in enumerate(seq):
        if ch!='-':
            if c==ri: return i
            c+=1
INV=[]; VAR=[]
newcol={}
for s in sites:
    ri=col2res(old[ref],s['col'])
    nc=res2col(aln[ref],ri)
    newcol[s['site']]=(nc,s['mcr'])
    (INV if s['mcr_pct']==100.0 and s['chrom_pct']==100.0 else VAR).append(s['site'])
print('invariant across all:',INV)
print('variable positions:',VAR)

rows=[]
for c in chrom:
    best=(0,None)
    for m in mcr:
        A,B=aln[c],aln[m]; ov=idn=0
        for x,y in zip(A,B):
            if x=='-' or y=='-': continue
            ov+=1
            if x==y: idn+=1
        v=100*idn/ov if ov else 0
        if v>best[0]: best=(v,m)
    hits=sum(1 for st,(nc,aa) in newcol.items() if aln[c][nc]==aa)
    md=meta.get(c,{})
    rows.append({'acc':c.split('|')[0],'organism':md.get('organism',c.split('|')[1] if '|' in c else 'UniProtKB record'),
       'length':len(raw[c]),'tm':md.get('tm'),'catalytic_16':hits,
       'closest_mcr':best[1].split('|')[0] if best[1] else None,'identity':round(best[0],1)})
rows.sort(key=lambda r:-r['identity'])
json.dump(rows,open('out/reservoir_v2.json','w'),indent=1)
print('\ntop 20 candidates')
print('%-16s %-34s %5s %4s %5s %-9s %s'%('accession','organism','len','TM','16sites','closest','id%'))
for r in rows[:20]:
    print('%-16s %-34s %5s %4s %5s %-9s %.1f'%(r['acc'],str(r['organism'])[:34],r['length'],r['tm'],r['catalytic_16'],r['closest_mcr'],r['identity']))
n70=sum(1 for r in rows if r['identity']>=70); n60=sum(1 for r in rows if r['identity']>=60)
print('\n>=70%%: %d   >=60%%: %d   total chromosomal/other: %d'%(n70,n60,len(rows)))
