import json, os
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
s=json.load(open('fam/summary.json'))
for f in s:
    lines=open('fam/%s.phy'%f).read().split('\n')[1:]
    recs=[]
    for ln in lines:
        if not ln.strip(): continue
        nm,sq=ln.split()[0], ln.split()[1]
        recs.append(SeqRecord(Seq(sq), id=nm))
    aln=MultipleSeqAlignment(recs)
    dm=DistanceCalculator('identity').get_distance(aln)
    tr=DistanceTreeConstructor().nj(dm)
    tr.root_at_midpoint()
    for c in tr.find_clades():
        if c.name and c.name.startswith('Inner'): c.name=None
        c.confidence=None
    Phylo.write(tr,'fam/%s.nwk'%f,'newick',format_branch_length='%0.6f')
    txt=open('fam/%s.nwk'%f).read().replace("'","").strip()
    open('fam/%s.nwk'%f,'w').write(txt+'\n')
    print(f,'tree ok')
