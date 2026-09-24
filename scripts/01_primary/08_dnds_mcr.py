import re, itertools, json, collections, warnings
warnings.filterwarnings('ignore')
from Bio.Seq import Seq
from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds
from pyfamsa import Aligner, Sequence
def readfa(p,up=False):
    d={};h=None
    for l in open(p):
        l=l.rstrip()
        if l.startswith('>'): h=l[1:]; d[h]=''
        else: d[h]+=(l.upper() if up else l)
    return d
nt=readfa('/home/claude/pet/data/mcr_ref_nt.fasta',up=True)
names=sorted(nt,key=lambda x:int(re.search(r'MCR-(\d+)',x).group(1)))
prot={k:str(Seq(v).translate(table=11)).rstrip('*') for k,v in nt.items()}
al=Aligner(guide_tree="upgma").align([Sequence(k.encode(),prot[k].encode()) for k in names])
pal={r.id.decode():r.sequence.decode() for r in al}
def backtranslate(name):
    p=pal[name]; s=nt[name]; out=[];i=0
    for ch in p:
        if ch=='-': out.append('---')
        else: out.append(s[i*3:i*3+3]); i+=1
    return ''.join(out)
cod={n:backtranslate(n) for n in names}
print('codon alignment columns:',len(cod[names[0]])//3)
res={}
print('\nPairwise dN/dS between MCR family reference alleles (NG86):')
print('%-8s %-8s %-8s %-8s %-7s'%('A','B','dN','dS','dN/dS'))
for a,b in itertools.combinations(names,2):
    A=''.join(x for x,y in zip([cod[a][i:i+3] for i in range(0,len(cod[a]),3)],[cod[b][i:i+3] for i in range(0,len(cod[b]),3)]) if '-' not in x and '-' not in y)
    B=''.join(y for x,y in zip([cod[a][i:i+3] for i in range(0,len(cod[a]),3)],[cod[b][i:i+3] for i in range(0,len(cod[b]),3)]) if '-' not in x and '-' not in y)
    try: dn,ds=cal_dn_ds(CodonSeq(A),CodonSeq(B),method='NG86')
    except Exception as e: continue
    if ds and ds>0:
        w=dn/ds
        res[(a.split('|')[0],b.split('|')[0])]=(round(dn,4),round(ds,4),round(w,4))
        print('%-8s %-8s %-8.4f %-8.4f %-7.4f'%(a.split('|')[0],b.split('|')[0],dn,ds,w))
json.dump({f'{k[0]}|{k[1]}':v for k,v in res.items()},open('dnds_between_mcr.json','w'),indent=1)
ws=[v[2] for v in res.values()]
print('\nbetween-family dN/dS: n=%d, median %.3f, range %.3f-%.3f'%(len(ws),sorted(ws)[len(ws)//2],min(ws),max(ws)))
close=[(k,v) for k,v in res.items() if v[1]<1.0]
print('\nclosely related pairs (dS<1): ')
for k,v in sorted(close,key=lambda x:x[1][1]): print('  %-9s %-9s dN=%.4f dS=%.4f w=%.3f'%(k[0],k[1],*v))
