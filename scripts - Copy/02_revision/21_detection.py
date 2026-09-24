import json, re, collections
from pyfamsa import Aligner, Sequence
from Bio.Seq import Seq
import pyhmmer
def readfa(p,up=False):
    d={};h=None
    for l in open(p):
        l=l.rstrip()
        if l.startswith('>'): h=l[1:]; d[h]=''
        else: d[h]+=(l.upper() if up else l)
    return d
ncbi=readfa('/home/claude/pet/data/abricate_ncbi.fa',up=True)
res=readfa('/home/claude/pet/data/abricate_resfinder.fa',up=True)
cat={}
for src in (ncbi,res):
    for h,s in src.items():
        p=h.split('~~~')
        m=re.match(r'(mcr-\d+)\.(\d+)',p[1])
        if m: cat[p[1]]=(m.group(1),str(Seq(s).translate(table=11)).rstrip('*'))
fam=collections.defaultdict(list)
for allele,(f,pr) in cat.items(): fam[f].append((allele,pr))
print('catalogued alleles: %d in %d families'%(len(cat),len(fam)))
chrom={k.split()[0]:v for k,v in readfa('/home/claude/pet/data/pet_chromosomal.faa').items()}
# identities from the existing multiple alignment (mcr family references vs candidates)
aln=readfa('/home/claude/pet/analysis/aln_full.faa')
def ident(a,b):
    m=s=0
    for x,y in zip(aln[a],aln[b]):
        if x=='-' or y=='-': continue
        s+=1; m+=(x==y)
    return 100.0*m/s
mcr=[k for k in aln if k.startswith('MCR-')]
cand=[k for k in aln if not k.startswith('MCR-')]
best={c:max((ident(c,m),m.split('|')[0]) for m in mcr) for c in cand}
vals=sorted(((v[0],c,v[1]) for c,v in best.items()),reverse=True)
print('\nhighest protein identity of any chromosomal candidate to an mcr family reference: %.1f%%'%vals[0][0])
for thr in (90,80,70,60,50):
    print('  candidates at >= %d%%: %d of %d'%(thr,sum(1 for v in vals if v[0]>=thr),len(vals)))
# within-family allele spread, to show the family reference is representative
spread={}
for f,mem in fam.items():
    if len(mem)<2: continue
    ref=[p for a,p in mem if a.endswith('.1')]
    if not ref: continue
    r=ref[0]
    al=Aligner(guide_tree='upgma').align([Sequence(a.encode(),p.encode()) for a,p in mem])
    A={x.id.decode():x.sequence.decode() for x in al}
    rid=[a for a,p in mem if p==r][0]
    ids=[]
    for a,p in mem:
        if a==rid: continue
        m=s=0
        for x,y in zip(A[rid],A[a]):
            if x=='-' or y=='-': continue
            s+=1; m+=(x==y)
        ids.append(100*m/s)
    spread[f]=(round(min(ids),1),round(max(ids),1),len(mem))
print('\nwithin-family identity to the .1 reference allele (min, max, n):')
print({k:v for k,v in sorted(spread.items())})
# profile arm
alpha=pyhmmer.easel.Alphabet.amino(); bg=pyhmmer.plan7.Background(alpha); builder=pyhmmer.plan7.Builder(alpha)
hmms=[]; cut={}
for f,mem in fam.items():
    if len(mem)<3: continue
    al=Aligner(guide_tree='upgma').align([Sequence(a.encode(),p.encode()) for a,p in mem])
    msa=pyhmmer.easel.TextMSA(name=f.encode(),sequences=[pyhmmer.easel.TextSequence(name=r.id.decode().encode(),sequence=r.sequence.decode()) for r in al]).digitize(alpha)
    hmm,_,_=builder.build_msa(msa,bg); hmms.append(hmm)
    q=[pyhmmer.easel.TextSequence(name=a.encode(),sequence=p).digitize(alpha) for a,p in mem]
    sc=[h.score for top in pyhmmer.hmmer.hmmsearch([hmm],q,E=1e10) for h in top]
    cut[f]=min(sc)
print('\nfamily profiles: %d'%len(hmms))
print('trusted-cutoff surrogate per family:',{k:round(v) for k,v in sorted(cut.items())})
cq=[pyhmmer.easel.TextSequence(name=k.encode(),sequence=v).digitize(alpha) for k,v in chrom.items()]
hits=collections.defaultdict(dict)
for top in pyhmmer.hmmer.hmmsearch(hmms,cq,E=1000.0,cpus=2):
    f=top.query.name.decode()
    for h in top: hits[h.name.decode()][f]=h.score
above=[(a,f,round(s,1)) for a,d in hits.items() for f,s in d.items() if f in cut and s>=cut[f]]
print('\ncandidates reaching a family trusted-cutoff surrogate: %d of %d'%(len(set(a for a,_,_ in above)),len(chrom)))
for a in sorted(above,key=lambda x:-x[2])[:8]: print('   ',a)
mx={}
for a,d in hits.items():
    f,s=max(d.items(),key=lambda x:x[1]); mx[a]=(f,round(s,1),round(cut.get(f,0)))
top=sorted(mx.items(),key=lambda x:-x[1][1])[:10]
print('\nhighest-scoring candidates against any family profile (score, cutoff):')
for a,(f,s,c) in top: print('   %-16s %-8s %7.1f  cutoff %d'%(a,f,s,c))
json.dump({'best_identity':{c:[round(v[0],1),v[1]] for c,v in best.items()},
           'cutoffs':{k:round(v,1) for k,v in cut.items()},
           'above_cutoff':above,
           'max_scores':{k:[v[0],v[1],v[2]] for k,v in mx.items()}},open('rev/out/detection.json','w'),indent=1)
