import json, collections
meta=json.load(open('metadata.json')); clade=json.load(open('clades.json'))
def readfa(p):
    d={};h=None
    for l in open(p):
        l=l.strip()
        if l.startswith('>'): h=l[1:]; d[h]=''
        else: d[h]+=l
    return d
aln=readfa('aln_full.faa'); seqs=readfa('all_pet.faa')
MCR1='MCR-1.1|NG_050417.1'
s=seqs[MCR1]
sites={'E246':246,'C281':281,'S284':284,'T285':285,'C291':291,'N329':329,'K333':333,'C356':356,'C364':364,'H390':390,'H395':395,'C414':414,'C422':422,'D465':465,'H466':466,'H478':478}
bad=[k for k,v in sites.items() if s[v-1]!=k[0]]
print('residue mismatches in MCR-1.1:',bad)
a=aln[MCR1]; pos2col={}; p=0
for c,ch in enumerate(a):
    if ch!='-': p+=1; pos2col[p]=c
mcr=[k for k in aln if k.startswith('MCR-')]
chrom=[k for k in aln if not k.startswith('MCR-')]
epta=[k for k in chrom if clade.get(k)=='EptA/MCR']
print('n mcr=%d n chrom=%d n eptA-clade chrom=%d'%(len(mcr),len(chrom),len(epta)))
def prof(col,g):
    cnt=collections.Counter(aln[k][col] for k in g); tot=sum(cnt.values())
    top,n=cnt.most_common(1)[0]; return top,100.0*n/tot
rows=[]
print('%-6s %-14s %-14s %-14s'%('site','MCR(11)','EptA chrom(65)','all chrom(150)'))
for k in sites:
    col=pos2col[sites[k]]
    m=prof(col,mcr); e=prof(col,epta); c=prof(col,chrom)
    print('%-6s %-14s %-14s %-14s'%(k,'%s %.0f%%'%m,'%s %.0f%%'%e,'%s %.0f%%'%c))
    rows.append(dict(site=k,col=col,mcr=m[0],mcr_pct=round(m[1],1),epta=e[0],epta_pct=round(e[1],1),chrom=c[0],chrom_pct=round(c[1],1)))
json.dump(rows,open('catalytic_sites.json','w'),indent=1)
