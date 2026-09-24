import veryfasttree, random, time
t0=time.time()
tree = veryfasttree.run(alignment='aln_automated1.faa', lg=True, gamma=True, threads=4, quiet=True, nopr=True)
open('tree_ml.nwk','w').write(tree if tree.endswith('\n') else tree+'\n')
print('ML tree done in %.1fs'%(time.time()-t0))
# bootstrap: resample alignment columns
def readfa(p):
    d=[];h=None
    for l in open(p):
        l=l.strip()
        if l.startswith('>'): h=l[1:]; d.append([h,''])
        else: d[-1][1]+=l
    return d
aln=readfa('aln_automated1.faa')
ncol=len(aln[0][1]); random.seed(42)
with open('boot_trees.nwk','w') as out:
    for rep in range(100):
        idx=[random.randrange(ncol) for _ in range(ncol)]
        txt=''.join('>'+h+'\n'+''.join(s[i] for i in idx)+'\n' for h,s in aln)
        bt=veryfasttree.run(alignment=txt, lg=True, threads=4, quiet=True, nopr=True, nosupport=True)
        out.write(bt if bt.endswith('\n') else bt+'\n')
        if (rep+1)%20==0: print('bootstrap',rep+1,'%.1fs'%(time.time()-t0), flush=True)
print('done %.1fs'%(time.time()-t0))
