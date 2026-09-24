import re, os, subprocess, json
from pyfamsa import Aligner, Sequence
from Bio.Seq import Seq

recs={}
name=None; buf=[]
for line in open('/home/claude/pet/data/abricate_ncbi.fa'):
    if line.startswith('>'):
        if name: recs[name]=''.join(buf)
        name=line[1:].strip(); buf=[]
    else: buf.append(line.strip())
recs[name]=''.join(buf)

fams={}
for h,s in recs.items():
    m=re.search(r'~~~(mcr-\d+)\.(\S+?)~~~', h, re.I)
    if not m: continue
    f=m.group(1).lower(); allele=m.group(1)+'.'+m.group(2)
    s=s.upper()
    if len(s)%3: continue
    aa=str(Seq(s).translate())
    if aa.count('*')!=1 or not aa.endswith('*'): continue
    fams.setdefault(f,{})[allele]=(s,aa[:-1])

out={}
for f,d in sorted(fams.items(), key=lambda x:int(x[0].split('-')[1])):
    # dereplicate identical protein+nt
    seen={}; keep={}
    for a,(nt,aa) in d.items():
        if nt in seen: continue
        seen[nt]=a; keep[a]=(nt,aa)
    if len(keep)<4: 
        print(f,'skip, only',len(keep),'unique nt'); continue
    al=Aligner(threads=1)
    msa=al.align([Sequence(a.encode(),aa.encode()) for a,(nt,aa) in keep.items()])
    prot={s.id.decode():s.sequence.decode() for s in msa}
    order=list(prot)
    L=len(prot[order[0]])
    cod={}
    for a in order:
        nt=keep[a][0][:-3]; p=0; row=[]
        for c in prot[a]:
            if c=='-': row.append('---')
            else: row.append(nt[p*3:p*3+3]); p+=1
        cod[a]=''.join(row)
    # strip gappy columns (any gap)
    keepcols=[i for i in range(L) if all(cod[a][i*3:i*3+3]!='---' for a in order)]
    phy=os.path.join('fam',f+'.phy')
    with open(phy,'w') as fh:
        fh.write(' %d %d\n'%(len(order), len(keepcols)*3))
        for a in order:
            seq=''.join(cod[a][i*3:i*3+3] for i in keepcols)
            fh.write('%s  %s\n'%(a.replace('.','_').replace('-','_')[:30].ljust(30), seq))
    out[f]={'n':len(order),'codons':len(keepcols),'alleles':order}
    print(f,'n=%d codons=%d'%(len(order),len(keepcols)))
json.dump(out,open('fam/summary.json','w'),indent=1)
