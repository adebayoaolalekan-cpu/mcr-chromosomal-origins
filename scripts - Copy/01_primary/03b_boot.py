import veryfasttree, random, time, os, sys, contextlib
def readfa(p):
    d=[];h=None
    for l in open(p):
        l=l.strip()
        if l.startswith('>'): h=l[1:]; d.append([h,''])
        else: d[-1][1]+=l
    return d
aln=readfa('aln_automated1.faa')
ncol=len(aln[0][1]); random.seed(42)
t0=time.time()
with open('boot_trees.nwk','w') as out:
    for rep in range(100):
        idx=[random.randrange(ncol) for _ in range(ncol)]
        txt=''.join('>'+h+'\n'+''.join(s[i] for i in idx)+'\n' for h,s in aln)
        bt=veryfasttree.run(alignment=txt, lg=True, threads=4, quiet=True, nopr=True, nosupport=True)
        out.write(bt if bt.endswith('\n') else bt+'\n')
print('bootstraps done %.1fs'%(time.time()-t0))
