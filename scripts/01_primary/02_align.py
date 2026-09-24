from pyfamsa import Aligner, Sequence
from pytrimal import Alignment, AutomaticTrimmer
seqs=[]
h=None
for l in open('all_pet.faa'):
    l=l.strip()
    if l.startswith('>'): h=l[1:]; seqs.append([h,''])
    else: seqs[-1][1]+=l
records=[Sequence(h.encode(), s.encode()) for h,s in seqs]
al=Aligner(guide_tree="upgma", threads=0).align(records)
with open('aln_full.faa','w') as o:
    for r in al:
        o.write('>'+r.id.decode()+'\n'+r.sequence.decode()+'\n')
print('aligned columns:', len(al[0].sequence))
a=Alignment.load('aln_full.faa')
for method in ['automated1','gappyout','strict']:
    t=AutomaticTrimmer(method=method).trim(a)
    print(method, '-> columns', len(t.sequences[0]))
    t.dump(f'aln_{method}.faa', format='fasta')
