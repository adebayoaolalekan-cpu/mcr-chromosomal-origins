import re, itertools, collections, json
from Bio.Seq import Seq
def readfa(p):
    d={};h=None
    for l in open(p):
        l=l.rstrip()
        if l.startswith('>'): h=l[1:]; d[h]=''
        else: d[h]+=l.upper()
    return d
ncbi=readfa('/home/claude/pet/data/abricate_ncbi.fa'); res=readfa('/home/claude/pet/data/abricate_resfinder.fa')
fam=collections.defaultdict(dict)
for h,s in ncbi.items():
    p=h.split('~~~')
    m=re.match(r'mcr-(\d+)\.(\d+)$',p[1])
    if m: fam['mcr-'+m.group(1)][p[1]]=(s,p[2])
for h,s in res.items():
    p=h.split('~~~')
    m=re.match(r'mcr-(\d+)\.(\d+)_',p[1])
    if m and m.group(1)=='9': fam['mcr-9'][p[1].split('_')[0]]=(s,p[2])
def gc(s): return 100.0*sum(c in 'GC' for c in s)/len(s)
def gc3(s): 
    t=s[2::3]; return 100.0*sum(c in 'GC' for c in t)/len(t)
print('%-8s %-4s %-7s %-7s %-8s'%('family','n','GC%','GC3%','len'))
out={}
for f in sorted(fam,key=lambda x:int(x.split('-')[1])):
    v=fam[f]; ref=v.get(f+'.1') or list(v.values())[0]
    print('%-8s %-4d %-7.1f %-7.1f %-8d'%(f,len(v),gc(ref[0]),gc3(ref[0]),len(ref[0])))
    out[f]=dict(n=len(v),gc=round(gc(ref[0]),1),gc3=round(gc3(ref[0]),1),length=len(ref[0]),acc=ref[1])
json.dump(out,open('mcr_gc.json','w'),indent=1)
# within-family dN/dS (Nei-Gojobori, pairwise) for families with >=4 variants
from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds
print('\nwithin-family pairwise dN/dS (Nei-Gojobori, mean over pairs):')
print('%-8s %-5s %-9s %-9s %-9s %s'%('family','pairs','mean dN','mean dS','mean w','max w'))
sel={}
for f in sorted(fam,key=lambda x:int(x.split('-')[1])):
    v=fam[f]
    if len(v)<4: continue
    names=sorted(v); L=len(v[names[0]][0])
    names=[n for n in names if len(v[n][0])==L]
    if len(names)<4: 
        Lc=collections.Counter(len(v[n][0]) for n in sorted(v)); L=Lc.most_common(1)[0][0]
        names=[n for n in sorted(v) if len(v[n][0])==L]
    dn_l=[];ds_l=[];w=[]
    for a,b in itertools.combinations(names,2):
        try:
            dn,ds=cal_dn_ds(CodonSeq(v[a][0]),CodonSeq(v[b][0]),method='NG86')
        except Exception: continue
        if dn is None or ds is None: continue
        dn_l.append(dn); ds_l.append(ds)
        if ds>0.002: w.append(dn/ds)
    if not w: continue
    sel[f]=dict(pairs=len(w),dn=round(sum(dn_l)/len(dn_l),4),ds=round(sum(ds_l)/len(ds_l),4),w=round(sum(w)/len(w),3),wmax=round(max(w),3),n_seqs=len(names))
    print('%-8s %-5d %-9.4f %-9.4f %-9.3f %.3f'%(f,len(w),sel[f]['dn'],sel[f]['ds'],sel[f]['w'],sel[f]['wmax']))
json.dump(sel,open('mcr_dnds.json','w'),indent=1)
