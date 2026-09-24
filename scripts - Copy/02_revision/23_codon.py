import os, re, json
from Bio.Seq import Seq
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
seqs={}
for k,v in mcr.items(): seqs[k.split('|')[0].replace('.','_')]=v
name={'Faucicola':'Faucicola_osloensis','Shewanella':'Shewanella_sp','Pseudomonas':'Pseudomonas_sp','Moraxella':'Moraxella_sp'}
for k,v in don.items(): seqs[name[k.split('_')[0]]]=v
prot={k:str(Seq(v).translate(table=11)).rstrip('*') for k,v in seqs.items()}
order=sorted(seqs)
al=Aligner(guide_tree='upgma').align([Sequence(k.encode(),prot[k].encode()) for k in order])
P={r.id.decode():r.sequence.decode() for r in al}
def bt(k):
    s=seqs[k]; out=[];i=0
    for ch in P[k]:
        if ch=='-': out.append('---')
        else: out.append(s[i*3:i*3+3]); i+=1
    return ''.join(out)
C={k:bt(k) for k in order}
L=len(C[order[0]])
# drop columns with any gap and any stop codon
keep=[]
for i in range(0,L,3):
    cods=[C[k][i:i+3] for k in order]
    if any('-' in c for c in cods): continue
    if any(str(Seq(c).translate(table=11))=='*' for c in cods): continue
    keep.append(i)
print('codon columns retained: %d of %d'%(len(keep),L//3))
clean={k:''.join(C[k][i:i+3] for i in keep) for k in order}
os.makedirs('rev/codeml',exist_ok=True)
with open('rev/codeml/codon.phy','w') as f:
    f.write(' %d %d\n'%(len(order),len(keep)*3))
    for k in order: f.write('%-24s %s\n'%(k[:24],clean[k]))
with open('rev/codeml/codon.fasta','w') as f:
    for k in order: f.write('>%s\n%s\n'%(k,clean[k]))
json.dump(order,open('rev/codeml/order.json','w'))
print('taxa:',order)
