from pyfamsa import Aligner, Sequence
from pytrimal import Alignment, AutomaticTrimmer
import re
seqs=[]; name=None; buf=[]
for line in open('v2_dataset.faa'):
    if line.startswith('>'):
        if name: seqs.append((name,''.join(buf)))
        name=line[1:].strip().split()[0]; buf=[]
    else: buf.append(line.strip())
seqs.append((name,''.join(buf)))
print('n=',len(seqs))
al=Aligner(threads=1)
msa=al.align([Sequence(n.encode(),s.encode()) for n,s in seqs])
with open('v2_aln_full.faa','w') as fh:
    for s in msa: fh.write('>%s\n%s\n'%(s.id.decode(),s.sequence.decode()))
a=Alignment.load('v2_aln_full.faa')
t=AutomaticTrimmer(method='automated1').trim(a)
t.dump('v2_aln_automated1.faa',format='fasta')
print('full cols',len(msa[0].sequence),'trimmed cols',len(t.sequences[0]))
