import warnings, json, re; warnings.filterwarnings('ignore')
from Bio.Seq import Seq
from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds
from pyfamsa import Aligner, Sequence
def readfa(p,up=True):
    d={};h=None
    for l in open(p):
        l=l.rstrip()
        if l.startswith('>'): h=l[1:].split()[0]; d[h]=''
        else: d[h]+=(l.upper() if up else l)
    return d
mcr=readfa('/home/claude/pet/data/mcr_ref_nt.fasta')
don=readfa('/home/claude/pet/data/loci/donor_cds.fna')
def gc(s): return 100*sum(c in 'GC' for c in s)/len(s)
def gc3(s):
    t=s[2::3]; return 100*sum(c in 'GC' for c in t)/len(t)
pairs=[('MCR-1.1','Faucicola_osloensis_CP190772.1_2453447-2455138_minus'),
       ('MCR-2.1','Faucicola_osloensis_CP190772.1_2453447-2455138_minus'),
       ('MCR-6.1','Faucicola_osloensis_CP190772.1_2453447-2455138_minus'),
       ('MCR-1.1','Moraxella_sp_DAMCJP010000009.1_20006-21652_minus'),
       ('MCR-2.1','Moraxella_sp_DAMCJP010000009.1_20006-21652_minus'),
       ('MCR-6.1','Moraxella_sp_DAMCJP010000009.1_20006-21652_minus'),
       ('MCR-4.1','Shewanella_sp_CP080411.1_2613104-2614729_minus'),
       ('MCR-12.1','Pseudomonas_sp_CM194505.1_700636-702288_plus')]
key={k.split('|')[0]:k for k in mcr}
res={}
print('%-9s %-28s %-7s %-7s %-7s %-7s'%('mcr','chromosomal homologue','dN','dS','dN/dS','aa id%'))
for m,d in pairs:
    a=mcr[key[m]]; b=don[d]
    pa=str(Seq(a).translate(table=11)).rstrip('*'); pb=str(Seq(b).translate(table=11)).rstrip('*')
    al=Aligner(guide_tree='upgma').align([Sequence(b'A',pa.encode()),Sequence(b'B',pb.encode())])
    A={r.id.decode():r.sequence.decode() for r in al}
    ia=ib=0; ca=[]; cb=[]; ident=0; n=0
    for x,y in zip(A['A'],A['B']):
        if x!='-' and y!='-':
            ca.append(a[ia*3:ia*3+3]); cb.append(b[ib*3:ib*3+3]); ident+= (x==y); n+=1
        if x!='-': ia+=1
        if y!='-': ib+=1
    dn,ds=cal_dn_ds(CodonSeq(''.join(ca)),CodonSeq(''.join(cb)),method='NG86')
    w=dn/ds if ds else float('nan')
    print('%-9s %-28s %-7.4f %-7.4f %-7.4f %-7.1f'%(m,d.split('_')[0],dn,ds,w,100*ident/n))
    res[m+'|'+d.split('_')[0]]=dict(dn=round(dn,4),ds=round(ds,4),w=round(w,4),aa_id=round(100*ident/n,1))
print('\nGC content of coding sequences:')
for m in ['MCR-1.1','MCR-2.1','MCR-4.1','MCR-6.1','MCR-12.1']:
    s=mcr[key[m]]; print('  %-9s GC=%.1f%% GC3=%.1f%%'%(m,gc(s),gc3(s)))
for d,s in don.items():
    print('  %-28s GC=%.1f%% GC3=%.1f%%'%(d.split('_')[0],gc(s),gc3(s)))
json.dump(res,open('dnds_donor.json','w'),indent=1)
