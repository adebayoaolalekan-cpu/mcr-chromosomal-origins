import json, numpy as np, re
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
btw=json.load(open('dnds_between_mcr.json'))
don=json.load(open('dnds_donor.json'))
gc=json.load(open('mcr_gc.json'))
fig,axes=plt.subplots(1,3,figsize=(7.4,2.7))
ax=axes[0]
w=[v[2] for v in btw.values()]
ax.hist(w,bins=np.arange(0,0.45,0.05),color='#4575b4',edgecolor='white')
ax.axvline(1.0,color='#c1272d',ls='--',lw=0.8)
ax.set_xlabel('dN/dS between MCR families',fontsize=6.5); ax.set_ylabel('number of pairs',fontsize=6.5)
ax.tick_params(labelsize=6); ax.set_title('a  Between-family selection',fontsize=7,loc='left')
for s in ['top','right']: ax.spines[s].set_visible(False)
ax=axes[1]
pairs=[(k,v) for k,v in don.items() if v['ds']>0]
names=[k.replace('|',' vs\n') for k,_ in pairs]
vals=[v['w'] for _,v in pairs]
ax.barh(range(len(vals)),vals,color='#e08214',height=0.6)
ax.set_yticks(range(len(vals))); ax.set_yticklabels(names,fontsize=4.8)
ax.set_xlabel('dN/dS, mobilised gene vs\nchromosomal progenitor',fontsize=6.5)
ax.tick_params(labelsize=6); ax.set_xlim(0,0.42)
ax.set_title('b  Mobilised vs chromosomal',fontsize=7,loc='left')
for s in ['top','right']: ax.spines[s].set_visible(False)
ax=axes[2]
fams=sorted(gc,key=lambda x:int(x.split('-')[1]))
x=np.arange(len(fams))
ax.bar(x-0.2,[gc[f]['gc'] for f in fams],width=0.4,label='GC',color='#1b5e9e')
ax.bar(x+0.2,[gc[f]['gc3'] for f in fams],width=0.4,label='GC3',color='#9ecae1')
ax.axhline(50.8,color='#c1272d',ls='--',lw=0.8)
ax.text(len(fams)-0.4,51.6,'E. coli genome GC',fontsize=4.8,color='#c1272d',ha='right',style='italic')
ax.set_xticks(x); ax.set_xticklabels([f.replace('mcr-','') for f in fams],fontsize=6)
ax.set_xlabel('mcr family',fontsize=6.5); ax.set_ylabel('% G+C',fontsize=6.5)
ax.set_ylim(30,85); ax.tick_params(labelsize=6); ax.legend(fontsize=5.5,frameon=False)
ax.set_title('c  Base composition of mcr genes',fontsize=7,loc='left')
for s in ['top','right']: ax.spines[s].set_visible(False)
plt.tight_layout()
plt.savefig('figures/Figure5_selection.png',dpi=600,bbox_inches='tight')
plt.savefig('figures/Figure5_selection.pdf',bbox_inches='tight')
print('fig5 ok')
