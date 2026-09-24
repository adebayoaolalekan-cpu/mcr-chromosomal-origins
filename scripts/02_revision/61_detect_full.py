import json, re, collections
from pyfamsa import Aligner, Sequence
from Bio.Seq import Seq
import pyhmmer

def readfa(p, up=False):
    d={};h=None
    for l in open(p):
        l=l.rstrip()
        if l.startswith('>'): h=l[1:].split()[0]; d[h]=''
        else: d[h]+=(l.upper() if up else l)
    return d

ncbi=readfa('/home/claude/pet/data/abricate_ncbi.fa',up=True)
resf=readfa('/home/claude/pet/data/abricate_resfinder.fa',up=True)
cat={}
for src in (ncbi,resf):
    for h,s in src.items():
        p=h.split('~~~')
        if len(p)<2: continue
        m=re.match(r'(mcr-\d+)\.(\d+)',p[1])
        if m: cat[p[1]]=(m.group(1), str(Seq(s).translate(table=11)).rstrip('*'))
fam=collections.defaultdict(list)
for allele,(f,pr) in cat.items(): fam[f].append((allele,pr))
print('catalogued alleles: %d in %d families'%(len(cat),len(fam)))
print('alleles per family:',{k:len(v) for k,v in sorted(fam.items(), key=lambda x:int(x[0].split('-')[1]))})

v2=readfa('v2_dataset.faa')
up=readfa('uniprot/all_uniprot.faa')
outs={'OpgB_ECOLI_P39401','AslA_ECOLI_P25549','YidJ_ECOLI_P31447','YdeN_ECOLI_P77318'}
query={k:v for k,v in v2.items() if not k.startswith('MCR-') and k not in outs}
for k in ('A0A220T317_Moraxella_sp_MSG47C17_mcr6like','A0A2R3NW71_Aeromonas_hydrophila_mcr5like'):
    if k in up: query[k]=up[k]
print('query set: %d proteins'%len(query))

alpha=pyhmmer.easel.Alphabet.amino()
bg=pyhmmer.plan7.Background(alpha)
builder=pyhmmer.plan7.Builder(alpha)
hmms=[]; cut={}
for f,mem in fam.items():
    if len(mem)<3: continue
    al=Aligner(guide_tree='upgma').align([Sequence(a.encode(),p.encode()) for a,p in mem])
    msa=pyhmmer.easel.TextMSA(name=f.encode(),
        sequences=[pyhmmer.easel.TextSequence(name=r.id.decode().encode(),sequence=r.sequence.decode()) for r in al]).digitize(alpha)
    hmm,_,_=builder.build_msa(msa,bg); hmms.append(hmm)
    q=[pyhmmer.easel.TextSequence(name=a.encode(),sequence=p).digitize(alpha) for a,p in mem]
    sc=[h.score for top in pyhmmer.hmmer.hmmsearch([hmm],q,E=1e6) for h in top]
    cut[f]=min(sc)
print('profiles built for %d families: %s'%(len(hmms),sorted(cut)))
print('thresholds:',{k:round(v,1) for k,v in sorted(cut.items())})

cq=[pyhmmer.easel.TextSequence(name=k.encode(),sequence=v).digitize(alpha) for k,v in query.items()]
hits=collections.defaultdict(dict)
for top in pyhmmer.hmmer.hmmsearch(hmms,cq,E=1000.0,cpus=2):
    qn=top.query.name; f=qn if isinstance(qn,str) else qn.decode()
    
    for h in top:
        hn=h.name; hits[hn if isinstance(hn,str) else hn.decode()][f]=h.score
best={}
for a,d in hits.items():
    f,s=max(d.items(),key=lambda x:x[1]); best[a]=[f,round(s,1)]
above=[(a,f,round(s,1)) for a,d in hits.items() for f,s in d.items() if f in cut and s>=cut[f]]
det=sorted(set(a for a,_,_ in above))
print('\nscored: %d of %d query proteins'%(len(hits),len(query)))
print('reaching a family threshold: %d'%len(det))
for a in det: print('   detected:',a,best[a])
rank=sorted(((v[1],a,v[0]) for a,v in best.items() if a not in det),reverse=True)[:6]
print('\nhighest scoring proteins that were NOT detected:')
for s,a,f in rank: print('   %-45s %-7s %8.1f   threshold %.1f'%(a,f,s,cut[f]))
json.dump({'n_query':len(query),'n_scored':len(hits),'cutoffs':{k:round(v,1) for k,v in cut.items()},
           'detected':det,'best':best,
           'top_undetected':[[a,f,s,round(cut[f],1)] for s,a,f in rank]},
          open('out/detection_full.json','w'),indent=1)
print('\nwritten out/detection_full.json')
