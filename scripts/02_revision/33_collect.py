import json, re, os, glob, math, statistics
R={}
R['topology']=json.load(open('out/topology_test.json'))
R['closest']=json.load(open('out/closest_v2.json'))
R['reservoir']=json.load(open('out/reservoir_v2.json'))
R['asr_stem']=json.load(open('out/asr_stem.json'))
det=json.load(open('out/detection_hmm.json'))
R['detection']={'cutoffs':det['cutoffs'],'n_detected':len(det['detected']),
   'best':sorted(((v[1],k,v[0]) for k,v in det['best'].items()),reverse=True)[:5]}
pw=json.load(open('codeml/pairwise_ml.json'))
ds=[x['ds'] for x in pw]
R['saturation']={'n_pairs':len(ds),'n_dS_gt1':sum(1 for x in ds if x>1),
   'median_dS':round(statistics.median(ds),2),'max_dS':round(max(ds),2),'min_dS':round(min(ds),3),
   'usable':[{'a':x['a'],'b':x['b'],'dN':round(x['dn'],4),'dS':round(x['ds'],4),'w':round(x['w'],4)} for x in pw if x['ds']<=1]}
fam={}
for d in sorted(glob.glob('fam/mcr-*_M0')):
    f=os.path.basename(d).split('_')[0]
    txt=open(os.path.join(d,'out.txt')).read()
    w=float(re.search(r'omega \(dN/dS\) =\s+([\d.]+)',txt).group(1))
    dn=float(re.search(r'tree length for dN:\s+([\d.]+)',txt).group(1))
    dsv=float(re.search(r'tree length for dS:\s+([\d.]+)',txt).group(1))
    n=int(open('fam/%s.phy'%f).readline().split()[0]); L=int(open('fam/%s.phy'%f).readline().split()[1])
    fam[f]={'n_alleles':n,'omega':None if dsv==0 else round(w,3),'dN_tree':dn,'dS_tree':dsv}
R['within_family']=fam
def lnL(p):
    m=re.search(r'lnL\(ntime:\s*\d+\s+np:\s*\d+\):\s*(-?[\d.]+)',open(p).read()); return float(m.group(1))
sm={}
for f in ['mcr-1','mcr-3']:
    try:
        d={M:lnL('fam/%s_M%s/out.txt'%(f,M)) for M in ['1','2','7','8']}
        D12=2*(d['2']-d['1']); D78=2*(d['8']-d['7'])
        sm[f]={'lnL':d,'M1a_M2a':{'2dL':round(D12,3),'p':round(math.exp(-D12/2),4)},
               'M7_M8':{'2dL':round(D78,3),'p':round(math.exp(-D78/2),4)}}
    except Exception as e: sm[f]=str(e)
R['site_models']=sm
R['dataset']={'n_v2':191,'n_primary':161,'n_chromosomal_refseq':150,'n_uniprot':28,'n_outgroup':4,
  'n_mcr':11,'genera':44,'records_above_genus':22,'taxon_queries':45,
  'aln_primary_full':808,'aln_primary_trim':313,'aln_v2_full':1038,'aln_v2_trim':444}
R['identical_pairs']=[['A0A220T317 (Moraxella sp. MSG47C17)','MCR-6.1'],['A0A2R3NW71 (Aeromonas hydrophila)','MCR-5.1']]
json.dump(R,open('out/results.json','w'),indent=1)
print(json.dumps({k:(v if not isinstance(v,(list,dict)) else '...') for k,v in R.items()},indent=1))
print('within family:',json.dumps(R['within_family'],indent=1))
print('site models:',json.dumps(R['site_models'],indent=1))
print('saturation:',R['saturation']['n_dS_gt1'],'of',R['saturation']['n_pairs'],'median',R['saturation']['median_dS'])
