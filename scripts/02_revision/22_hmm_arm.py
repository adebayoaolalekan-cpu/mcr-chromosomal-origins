import pyhmmer, json, collections, re, sys
def nm(x): return x if isinstance(x,str) else x.decode()
from pyfamsa import Aligner, Sequence
from Bio.Seq import Seq
def readfa(p,up=False):
    d={};h=None
    for l in open(p):
        l=l.rstrip()
        if l.startswith('>'): h=l[1:]; d[h]=''
        else: d[h]+=(l.upper() if up else l)
    return d
chrom={k.split()[0]:v for k,v in readfa('/home/claude/pet/data/pet_chromosomal.faa').items()}
cat={}
for f in ['/home/claude/pet/data/abricate_ncbi.fa','/home/claude/pet/data/abricate_resfinder.fa']:
    for h,s in readfa(f,up=True).items():
        p=h.split('~~~'); m=re.match(r'(mcr-\d+)\.\d+',p[1])
        if m: cat[p[1].split('_')[0]]=(m.group(1),str(Seq(s).translate(table=11)).rstrip('*'))
fam=collections.defaultdict(list)
for a,(f,p) in cat.items(): fam[f].append((a,p))
alpha=pyhmmer.easel.Alphabet.amino(); bg=pyhmmer.plan7.Background(alpha); bl=pyhmmer.plan7.Builder(alpha)
hmms=[]; cut={}
for f,mem in sorted(fam.items()):
    if len(mem)<3: continue
    al=Aligner(guide_tree='upgma').align([Sequence(a.encode(),p.encode()) for a,p in mem])
    msa=pyhmmer.easel.TextMSA(name=f.encode(),sequences=[pyhmmer.easel.TextSequence(name=r.id if isinstance(r.id,bytes) else r.id.encode(),sequence=r.sequence.decode() if isinstance(r.sequence,bytes) else r.sequence) for r in al]).digitize(alpha)
    hmm,_,_=bl.build_msa(msa,bg); hmms.append(hmm)
    q=[pyhmmer.easel.TextSequence(name=a.encode(),sequence=p).digitize(alpha) for a,p in mem]
    sc=[h.score for top in pyhmmer.hmmer.hmmsearch([hmm],q,E=10.0,cpus=1) for h in top]
    cut[f]=min(sc)
    print('profile %-7s n=%-3d cutoff %.0f'%(f,len(mem),cut[f]),flush=True)
cq=[pyhmmer.easel.TextSequence(name=k.encode(),sequence=v).digitize(alpha) for k,v in chrom.items()]
best={}
for hmm in hmms:
    f=nm(hmm.name)
    for top in pyhmmer.hmmer.hmmsearch([hmm],cq,E=10.0,cpus=1):
        for h in top:
            a=nm(h.name)
            if a not in best or h.score>best[a][1]: best[a]=(f,h.score)
    print('scored against',f,flush=True)
det=[(a,f,round(s,1),round(cut[f])) for a,(f,s) in best.items() if s>=cut[f]]
print('\ncandidates reaching a family cutoff: %d of %d'%(len(det),len(chrom)),flush=True)
top10=sorted(best.items(),key=lambda x:-x[1][1])[:12]
print('\n%-16s %-8s %9s %9s %s'%('accession','profile','score','cutoff','call'))
for a,(f,s) in top10:
    print('%-16s %-8s %9.1f %9d %s'%(a,f,s,round(cut[f]),'detected' if s>=cut[f] else 'missed'))
json.dump({'cutoffs':{k:round(v,1) for k,v in cut.items()},'best':{a:[f,round(s,1)] for a,(f,s) in best.items()},'detected':det},open('rev/out/detection_hmm.json','w'),indent=1)
