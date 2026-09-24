import json, itertools, collections, re
meta=json.load(open('metadata.json'))
def readfa(p):
    d={};h=None
    for l in open(p):
        l=l.strip()
        if l.startswith('>'): h=l[1:]; d[h]=''
        else: d[h]+=l
    return d
aln=readfa('aln_full.faa')
labs=list(aln)
mcr=[l for l in labs if l.startswith('MCR-')]
chrom=[l for l in labs if not l.startswith('MCR-')]
def ident(a,b):
    A,B=aln[a],aln[b]; m=s=0
    for x,y in zip(A,B):
        if x=='-' and y=='-': continue
        if x=='-' or y=='-': continue
        s+=1; m+= (x==y)
    return 100.0*m/s if s else 0.0
rows=[]
print('%-10s %-8s %-46s %-34s %s'%('MCR','best%','closest chromosomal homologue','organism','2nd best'))
best_map={}
for m in sorted(mcr,key=lambda x:int(re.search(r'MCR-(\d+)',x).group(1))):
    sc=sorted(((ident(m,c),c) for c in chrom),reverse=True)
    b=sc[0]; s2=sc[1]
    org=meta[b[1]]['organism']
    print('%-10s %-8.1f %-46s %-34s %s (%.1f%%)'%(m.split('|')[0],b[0],b[1].split('|')[0]+' '+meta[b[1]]['desc'][:28],org,s2[1].split('|')[0],s2[0]))
    best_map[m]=[(round(x,1),y,meta[y]['organism'],meta[y]['desc']) for x,y in sc[:5]]
json.dump(best_map,open('mcr_closest_chromosomal.json','w'),indent=1)
# MCR vs MCR identity
print('\nPairwise identity among MCR family reference proteins (%):')
order=sorted(mcr,key=lambda x:int(re.search(r'MCR-(\d+)',x).group(1)))
print('        '+' '.join('%6s'%o.split('|')[0].replace('MCR-','') for o in order))
mat={}
for a in order:
    row=[]
    for b in order:
        v=100.0 if a==b else ident(a,b); row.append(v); mat[(a,b)]=v
    print('%-7s '%a.split('|')[0].replace('MCR-','')+' '.join('%6.1f'%v for v in row))
vals=[mat[(a,b)] for a,b in itertools.combinations(order,2)]
print('\nBetween-family identity: min %.1f max %.1f mean %.1f'%(min(vals),max(vals),sum(vals)/len(vals)))
# EptA E. coli vs MCR
ec='NP_418538.2|Escherichia_coli_str_K-12_substr_MG1655'
print('\nIdentity to E. coli K-12 EptA (NP_418538.2):')
for a in order: print('  %-8s %.1f'%(a.split('|')[0],ident(a,ec)))
