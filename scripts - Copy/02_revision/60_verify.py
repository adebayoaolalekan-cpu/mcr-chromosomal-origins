import json, collections
det=json.load(open('out/detection_hmm.json'))
print('detection set size:',len(det['best']))
print('families with profiles:',sorted(det['cutoffs']))
res=json.load(open('out/reservoir_v2.json'))
ids={r['acc'] for r in res}
print('reservoir n:',len(res),'A0A220T317 present:', 'A0A220T317' in ids,' A0A2R3NW71 present:','A0A2R3NW71' in ids)
full=[r for r in res if r['length']>=450]
for t in (90,80,70,60):
    print('  >=%d%% full-length: %d  (incl partial %d)'%(t,sum(1 for r in full if r['identity']>=t),sum(1 for r in res if r['identity']>=t)))
print('genera in top20:',collections.Counter(r['organism'].split()[0] for r in full[:20]).most_common())
print('TM min/max over the whole set:',min(r['tm'] for r in res),max(r['tm'] for r in res))
tm=collections.Counter(r['tm'] for r in res)
print('TM distribution:',dict(sorted(tm.items())))
best=det['best']
have=[r for r in res if r['acc'] in best]
print('overlap of detection set with reservoir table:',len(have))
if have:
    hi=sorted(have,key=lambda r:-r['identity'])[:6]
    for r in hi: print('   %-16s %5.1f%%  best profile %s'%(r['acc'],r['identity'],best[r['acc']]))
R=json.load(open('out/results.json'))
print('median dS:',R['saturation']['median_dS'])
print('detection best five:',R['detection']['best'])
