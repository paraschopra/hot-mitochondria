"""Per-condition matrix-cytosol ΔTm (revised 25 conditions from 16_revised.py) with gene-bootstrap 95% CIs, for figure panel B."""
import pandas as pd, numpy as np
rng=np.random.default_rng(7)
t = pd.read_csv('data/meltome_labeled.csv', low_memory=False)
t = t[(t.Organism!='Human body fluids')]
t['Tm']=t['Melting point [°C]']
meta=pd.read_csv('data/meltome_datasets.csv'); meta['dataset']=meta['Dataset ID'].str.replace(':','_')
t=t.merge(meta[['dataset','Dataset source']], on='dataset', how='left')
euk=['Homo sapiens','Mus musculus','Danio rerio','Drosophila melanogaster','Caenorhabditis elegans','Saccharomyces cerevisiae']
t=t[t.Organism.isin(euk)]
t['gene']=np.where(t.Organism=='Homo sapiens', t.sym.str.upper(), np.where(t.one2one, t.human_syms, None))
t=t[t.gene.notna()]
t['grp']=np.where(t.mc_loc=='Matrix','M', np.where(t.cyto_ref,'C', None)); t=t[t.grp.notna()]
t['cond']=t['Dataset source']+'|'+t.Organism+'|'+t['Strain/tissue'].astype(str)+'|'+t['Cells/lysate']+'|'+t['Lysis conditions']
g=t.groupby(['cond','dataset','gene','grp']).agg(Tm=('Tm','mean')).reset_index()
genes=sorted(g.gene.unique()); gi={x:i for i,x in enumerate(genes)}; g['gi']=g.gene.map(gi)
def wmedian(x,w):
    o=np.argsort(x); x=x[o]; w=w[o]; c=np.cumsum(w); return x[np.searchsorted(c, c[-1]/2)]
blocks=[]
for (c,ds),d in g.groupby(['cond','dataset']):
    d=d.dropna(subset=['Tm']); blocks.append((c, d.gi.values, d.Tm.values, (d.grp=='M').values))
def est(w):
    res={}
    for c,gidx,tm,ism in blocks:
        ww=w[gidx]
        if ww[ism].sum()<10 or ww[~ism].sum()<30: continue
        res.setdefault(c,[]).append(wmedian(tm[ism],ww[ism]) - wmedian(tm[~ism],ww[~ism]))
    return {c:np.mean(v) for c,v in res.items()}
E=pd.Series(est(np.ones(len(genes))))
bs=pd.DataFrame([est(rng.poisson(1.0,len(genes)).astype(float)) for _ in range(500)])
n=g.groupby(['cond','grp']).gene.nunique().unstack()
out=pd.DataFrame({'d':E, 'lo':bs.quantile(.025), 'hi':bs.quantile(.975), 'n_matrix':n['M'], 'n_cyto':n['C']}).loc[E.index]
out.index.name='cond'; out.to_csv('results/conditions_revised.csv'); print(out.round(2).to_string())
