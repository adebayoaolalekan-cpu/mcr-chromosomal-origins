import json, itertools
def rd(p):
    d={};n=None;b=[]
    for l in open(p):
        if l[0]=='>':
            if n:d[n]=''.join(b)
            n=l[1:].strip().split()[0];b=[]
        else:b.append(l.strip())
    d[n]=''.join(b);return d
aln=rd('v2_aln_full.faa')
raw=rd('v2_dataset.faa')
mcr=[k for k in aln if k.startswith('MCR-')]
out_names={'OpgB_ECOLI_P39401','AslA_ECOLI_P25549','YidJ_ECOLI_P31447','YdeN_ECOLI_P77318'}
chrom=[k for k in aln if k not in mcr and k not in out_names]
print('mcr',len(mcr),'chromosomal/other',len(chrom))
PARTIAL={k for k in raw if len(raw[k])<450}
print('partials:',PARTIAL)

def ident(a,b):
    A,B=aln[a],aln[b]
    ov=idn=0
    for x,y in zip(A,B):
        if x=='-' or y=='-': continue
        ov+=1
        if x==y: idn+=1
    # global: denominator = columns where at least one is a residue
    tot=sum(1 for x,y in zip(A,B) if x!='-' or y!='-')
    return (100*idn/ov if ov else 0, 100*idn/tot if tot else 0, ov)

res={}
for m in sorted(mcr):
    rows=[]
    for c in chrom:
        lo,gl,ov=ident(m,c)
        rows.append((lo,gl,ov,c))
    rows.sort(reverse=True)
    best_all=rows[0]
    rows_f=[r for r in rows if r[3] not in PARTIAL]
    best_f=rows_f[0]
    res[m]={'best_incl_partial':{'hit':best_all[3],'local':round(best_all[0],1),'global':round(best_all[1],1),'ovl':best_all[2]},
            'best_full_only':{'hit':best_f[3],'local':round(best_f[0],1),'global':round(best_f[1],1),'ovl':best_f[2]},
            'top5':[{'hit':r[3],'local':round(r[0],1),'global':round(r[1],1)} for r in rows_f[:5]]}
    print('%-14s %-42s local=%5.1f glob=%5.1f | (incl partial: %-30s %5.1f)'%(m,best_f[3][:42],best_f[0],best_f[1],best_all[3][:30],best_all[0]))
json.dump(res,open('out/closest_v2.json','w'),indent=1)
