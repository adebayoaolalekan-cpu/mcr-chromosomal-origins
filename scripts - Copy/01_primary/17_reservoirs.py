import json, re, csv
meta=json.load(open('metadata.json')); clade=json.load(open('clades.json'))
def readfa(p):
    d={};h=None
    for l in open(p):
        l=l.strip()
        if l.startswith('>'): h=l[1:]; d[h]=''
        else: d[h]+=l
    return d
aln=readfa('aln_full.faa')
def ident(a,b):
    m=s=0
    for x,y in zip(aln[a],aln[b]):
        if x=='-' or y=='-': continue
        s+=1; m+=(x==y)
    return 100.0*m/s
mcr=[k for k in aln if k.startswith('MCR-')]
chrom=[k for k in aln if not k.startswith('MCR-')]
rows=[]
for c in chrom:
    best=max((ident(c,m),m) for m in mcr)
    rows.append(dict(acc=meta[c]['acc'],organism=meta[c]['organism'],desc=meta[c]['desc'],
                     clade=clade.get(c,''),length=meta[c]['length'],
                     best_mcr=best[1].split('|')[0],best_identity=round(best[0],1),
                     epta_identity=round(ident(c,'NP_418538.2|Escherichia_coli_str_K-12_substr_MG1655'),1)))
rows.sort(key=lambda r:-r['best_identity'])
with open('TableS1_chromosomal_homologues.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print('Top chromosomal candidates (highest identity to any MCR family):')
print('%-16s %-34s %-12s %-7s %-10s'%('accession','organism','closest MCR','id%','clade'))
for r in rows[:18]:
    print('%-16s %-34s %-12s %-7.1f %-10s'%(r['acc'],r['organism'][:33],r['best_mcr'],r['best_identity'],r['clade']))
n60=sum(1 for r in rows if r['best_identity']>=60); n70=sum(1 for r in rows if r['best_identity']>=70)
print('\nchromosomal proteins >=60%% identity to an MCR family: %d; >=70%%: %d'%(n60,n70))
ent=[r for r in rows if r['clade']=='EptA/MCR']
print('EptA/MCR clade chromosomal members: %d; median best identity %.1f%%'%(len(ent),sorted(r['best_identity'] for r in ent)[len(ent)//2]))
