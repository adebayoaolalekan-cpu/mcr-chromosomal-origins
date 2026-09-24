import re, collections
def readfa(p):
    d={};h=None
    for l in open(p):
        l=l.rstrip()
        if l.startswith('>'): h=l[1:]; d[h]=''
        else: d[h]+=l
    return d
chrom=readfa('/home/claude/pet/data/pet_chromosomal.faa')
mcr=readfa('/home/claude/pet/data/mcr_ref_aa.fasta')
def org(h):
    m=re.search(r'\[([^\]]+)\]',h); return m.group(1) if m else 'unknown'
def label(acc,h):
    o=org(h).replace(' ','_').replace('.','').replace('(','').replace(')','').replace('/','-')
    return f"{acc}|{o}"
meta={}
out=open('all_pet.faa','w')
for h,s in chrom.items():
    acc=h.split()[0]
    lab=label(acc,h)
    meta[lab]=dict(acc=acc,organism=org(h),desc=h.split(' ',1)[1] if ' ' in h else '',cls='chromosomal',length=len(s))
    out.write(f">{lab}\n{s}\n")
for h,s in mcr.items():
    name,acc=h.split('|')
    lab=f"{name}|{acc}"
    meta[lab]=dict(acc=acc,organism='plasmid-borne',desc=name,cls='mobilised',length=len(s))
    out.write(f">{lab}\n{s}\n")
out.close()
import json
json.dump(meta,open('metadata.json','w'),indent=1)
print('total sequences:',len(meta))
print('chromosomal:',sum(1 for v in meta.values() if v['cls']=='chromosomal'))
print('mobilised:',sum(1 for v in meta.values() if v['cls']=='mobilised'))
