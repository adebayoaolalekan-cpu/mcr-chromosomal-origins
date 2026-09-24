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
old=rd('/home/claude/pet/analysis/aln_full.faa')
meta=json.load(open('/home/claude/pet/analysis/metadata.json'))
sites=json.load(open('/home/claude/pet/analysis/catalytic_sites.json'))
ref='MCR-1.1|NG_050417.1'
mcr=[k for k in aln if k.startswith('MCR-')]
outs={'OpgB_ECOLI_P39401','AslA_ECOLI_P25549','YidJ_ECOLI_P31447','YdeN_ECOLI_P77318'}
chrom=[k for k in aln if k not in mcr and k not in outs]

def col2res(seq,col):
    return sum(1 for c in seq[:col+1] if c!='-')-1 if seq[col]!='-' else None
def res2col(seq,ri):
    c=0
    for i,ch in enumerate(seq):
        if ch!='-':
            if c==ri: return i
            c+=1
newcol={}
for s in sites:
    ri=col2res(old[ref],s['col']); nc=res2col(aln[ref],ri)
    newcol[s['site']]=(nc, aln[ref][nc])
print('MCR-1.1 residues at the 16 positions:', {k:v[1] for k,v in newcol.items()})

KD={'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2,'X':0,'B':-3.5,'Z':-3.5,'U':2.5}
def tmcount(seq,w=19,thr=1.6):
    h=[sum(KD.get(c,0) for c in seq[i:i+w])/w for i in range(len(seq)-w+1)]
    segs=[];i=0
    while i<len(h):
        if h[i]>thr:
            j=i
            while j<len(h) and h[j]>thr: j+=1
            segs.append((i+1,j+w-1)); i=j
        else: i+=1
    segs=[x for x in segs if x[1]-x[0]+1>=w]
    merged=[]
    for a,b in segs:
        if merged and a<=merged[-1][1]: merged[-1][1]=max(merged[-1][1],b)
        else: merged.append([a,b])
    return len(merged)

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
    SET14=[k for k in newcol if k not in ('S284','N329')]
    hits=sum(1 for st in SET14 if aln[c][newcol[st][0]]==newcol[st][1])
    md=meta.get(c,{})
    org=md.get('organism')
    if not org:
        parts=c.split('_')
        org=' '.join(parts[1:3]) if len(parts)>2 else c
    rows.append({'id':c,'acc':(c.split('|')[0] if '|' in c else c.split('_')[0]),'organism':org,'source':'RefSeq' if '|' in c else 'UniProtKB',
       'length':len(raw[c]),'tm':tmcount(raw[c]),
       'catalytic_14':hits,'closest_mcr':best[1].split('|')[0],'identity':round(best[0],1)})
rows.sort(key=lambda r:-r['identity'])
json.dump(rows,open('out/reservoir_v2.json','w'),indent=1)
print('\n%-18s %-30s %-9s %4s %3s %4s %-9s %s'%('accession','organism','source','len','TM','16','closest','id%'))
for r in rows[:22]:
    print('%-18s %-30s %-9s %4d %3d %4d %-9s %.1f'%(r['acc'],str(r['organism'])[:32],r['source'],r['length'],r['tm'],r['catalytic_14'],r['closest_mcr'],r['identity']))
print()
for t in (90,80,70,60,50):
    print('>=%d%%: %d'%(t,sum(1 for r in rows if r['identity']>=t)))
print('complete 14/14 active-site set:',sum(1 for r in rows if r['catalytic_14']==14),'of',len(rows))
