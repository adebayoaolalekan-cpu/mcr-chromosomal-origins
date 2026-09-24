import json, math, collections, re
meta=json.load(open('metadata.json'))
def readfa(p):
    d={};h=None
    for l in open(p):
        l=l.strip()
        if l.startswith('>'): h=l[1:]; d[h]=''
        else: d[h]+=l
    return d
aln=readfa('aln_full.faa'); seqs=readfa('all_pet.faa')
MCR1='MCR-1.1|NG_050417.1'; EPTA='NP_418538.2|Escherichia_coli_str_K-12_substr_MG1655'
s=seqs[MCR1]
print('MCR-1.1 length', len(s))
# residues reported as catalytic / metal-binding / disulfide in MCR-1 structures
sites={'E246':246,'T285':285,'H395':395,'D465':465,'H466':466,'H478':478,
       'C281':281,'C291':291,'C356':356,'C364':364,'C414':414,'C422':422,
       'K333':333,'H478b':478,'N329':329,'H390':390,'S284':284}
print('residue identity check in MCR-1.1 (1-based):')
for k,v in sites.items():
    print('   %-6s -> %s'%(k,s[v-1]))
# map MCR-1 positions to alignment columns
a=aln[MCR1]; pos2col={}; p=0
for c,ch in enumerate(a):
    if ch!='-': p+=1; pos2col[p]=c
def col_profile(col, group):
    cnt=collections.Counter(aln[k][col] for k in group)
    tot=sum(cnt.values()); top,n=cnt.most_common(1)[0]
    return top, 100.0*n/tot, cnt
mcr=[k for k in aln if k.startswith('MCR-')]
chrom=[k for k in aln if not k.startswith('MCR-')]
eptA_clade=[k for k in chrom if 'epta' in meta[k]['desc'].lower() or '--lipid a' in meta[k]['desc'].lower()]
print('\nconservation at catalytic/metal sites (MCR-1 numbering):')
print('%-6s %-4s %-22s %-22s %s'%('site','col','MCR set (n=%d)'%len(mcr),'all chromosomal (n=%d)'%len(chrom),'EptA-like (n=%d)'%len(eptA_clade)))
rows=[]
for k in ['E246','C281','S284','T285','C291','N329','K333','C356','C364','H390','H395','C414','C422','D465','H466','H478']:
    v=sites[k]; col=pos2col[v]
    a1=col_profile(col,mcr); a2=col_profile(col,chrom); a3=col_profile(col,eptA_clade)
    print('%-6s %-4d %-22s %-22s %s'%(k,col,'%s %.0f%%'%(a1[0],a1[1]),'%s %.0f%%'%(a2[0],a2[1]),'%s %.0f%%'%(a3[0],a3[1])))
    rows.append(dict(site=k,col=col,mcr=a1[0],mcr_pct=round(a1[1],1),chrom=a2[0],chrom_pct=round(a2[1],1),epta=a3[0],epta_pct=round(a3[1],1)))
json.dump(rows,open('catalytic_sites.json','w'),indent=1)
# transmembrane helix prediction (Kyte-Doolittle)
KD=dict(A=1.8,R=-4.5,N=-3.5,D=-3.5,C=2.5,Q=-3.5,E=-3.5,G=-0.4,H=-3.2,I=4.5,L=3.8,K=-3.9,M=1.9,F=2.8,P=-1.6,S=-0.8,T=-0.7,W=-0.9,Y=-1.3,V=4.2)
def tm_segments(seq,w=19,thr=1.6):
    h=[sum(KD.get(c,0) for c in seq[i:i+w])/w for i in range(len(seq)-w+1)]
    segs=[];i=0
    while i<len(h):
        if h[i]>thr:
            j=i
            while j<len(h) and h[j]>thr: j+=1
            segs.append((i+1,j+w-1)); i=j
        else: i+=1
    return [s for s in segs if s[1]-s[0]+1>=w]
tmc=collections.defaultdict(list)
for k,sq in seqs.items():
    n=len(tm_segments(sq))
    grp='mobilised MCR' if k.startswith('MCR-') else 'chromosomal'
    tmc[grp].append(n)
    meta[k]['tm']=n
for g,v in tmc.items():
    print('\npredicted TM segments, %s: median %d (range %d-%d)'%(g,sorted(v)[len(v)//2],min(v),max(v)))
# domain boundary: last TM of MCR-1 and EptA
print('MCR-1.1 TM segments:',tm_segments(seqs[MCR1]))
print('EptA E. coli TM segments:',tm_segments(seqs[EPTA]))
json.dump(meta,open('metadata.json','w'),indent=1)
