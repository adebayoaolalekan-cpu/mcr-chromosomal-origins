import json, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=json.load(open('out/results.json'))
pw=json.load(open('codeml/pairwise_ml.json'))
gc=json.load(open('/home/claude/pet/analysis/mcr_gc.json'))
wf=R['within_family']
fig=plt.figure(figsize=(7.4,2.9))
gs=fig.add_gridspec(1,3,width_ratios=[1,1,1.05],wspace=0.42)

ax=fig.add_subplot(gs[0,0])
ds=[x['ds'] for x in pw]
ax.hist(ds,bins=np.arange(0,95,5),color='#4575b4',edgecolor='white',lw=0.4)
ax.axvline(1.0,color='#c1272d',lw=1.2,ls='--')
ax.text(3.0,ax.get_ylim()[1]*0.92,'dS = 1',fontsize=6,color='#c1272d')
ax.set_xlabel('synonymous substitutions per site (dS)',fontsize=6.6)
ax.set_ylabel('between-family comparisons',fontsize=6.6)
ax.tick_params(labelsize=6)
ax.set_title('a  Synonymous sites are saturated\n     between families',fontsize=7.5,loc='left')
ax.text(0.40,0.62,'%d of %d pairs\nabove dS = 1\nmedian %.1f'%(R['saturation']['n_dS_gt1'],R['saturation']['n_pairs'],R['saturation']['median_dS']),
        transform=ax.transAxes,fontsize=6.2)
for s in ['top','right']: ax.spines[s].set_visible(False)

ax=fig.add_subplot(gs[0,1])
fams=[f for f in ['mcr-1','mcr-2','mcr-3','mcr-4','mcr-8','mcr-10'] if wf[f]['omega'] is not None]
xs=[wf[f]['dS_tree'] for f in fams]; ys=[wf[f]['omega'] for f in fams]
ax.scatter(xs,ys,s=[max(18,wf[f]['n_alleles']*2.2) for f in fams],color='#1b5e9e',alpha=0.85,edgecolor='white',lw=0.6)
OFF={'mcr-1':(7,3),'mcr-10':(7,4),'mcr-8':(8,-2),'mcr-4':(6,-10),'mcr-2':(6,4),'mcr-3':(8,2)}
for f,x,y in zip(fams,xs,ys):
    ax.annotate(f,(x,y),textcoords='offset points',xytext=OFF.get(f,(5,4)),fontsize=6,style='italic')
ax.axhline(1.0,color='#999999',lw=0.8,ls=':')
ax.set_xlim(0,1.0); ax.set_ylim(0,1.05)
ax.set_xlabel('synonymous tree length within the family',fontsize=6.6)
ax.set_ylabel('dN/dS (model M0)',fontsize=6.6)
ax.tick_params(labelsize=6)
ax.set_title('b  Purifying selection within\n     every family',fontsize=7.5,loc='left')
ax.text(0.02,0.92,'mcr-5 omitted: dS = 0',transform=ax.transAxes,fontsize=6,color='#666666',style='italic')
for s in ['top','right']: ax.spines[s].set_visible(False)

ax=fig.add_subplot(gs[0,2])
order=['mcr-1','mcr-2','mcr-3','mcr-4','mcr-5','mcr-6','mcr-7','mcr-8','mcr-9','mcr-10','mcr-12']
order=[f for f in order if f in gc]
x=np.arange(len(order))
ax.bar(x-0.19,[gc[f]['gc'] for f in order],width=0.36,color='#1b5e9e',label='overall G+C')
ax.bar(x+0.19,[gc[f]['gc3'] for f in order],width=0.36,color='#9ecae1',label='third position')
ax.axhline(50.8,color='#c1272d',lw=0.9,ls='--')
ax.text(len(order)-0.4,51.8,'E. coli',fontsize=5.6,color='#c1272d',ha='right',style='italic')
ax.set_xticks(x); ax.set_xticklabels([f.replace('mcr-','') for f in order],fontsize=6)
ax.set_xlabel('mcr family',fontsize=6.6,style='italic')
ax.set_ylabel('G+C content (%)',fontsize=6.6)
ax.set_ylim(30,82); ax.tick_params(labelsize=6)
ax.legend(fontsize=5.8,frameon=False,loc='upper left',ncol=2)
ax.set_title('c  Base composition varies\n     across the families',fontsize=7.5,loc='left')
for s in ['top','right']: ax.spines[s].set_visible(False)

plt.savefig('figures/Figure5_selection.png',dpi=600,bbox_inches='tight')
plt.savefig('figures/Figure5_selection.pdf',bbox_inches='tight')
print('fig5 ok')
