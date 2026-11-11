"""Inner-membrane vs matrix: is any mitochondrial set closer to the putative heat source more stable than the cytosol?
Same estimator as 16_revised.py (per-dataset median difference vs cytosolic reference, lab-aware conditions,
two-level gene + condition bootstrap), split by lysate vs intact cells. Answers the objection that heating is
confined to the inner membrane (Fahimi, Lynch & Matta 2026)."""
import pandas as pd, numpy as np
rng=np.random.default_rng(11)
t = pd.read_csv('data/meltome_labeled.csv', low_memory=False)
t = t[(t.Organism!='Human body fluids')]
t['Tm']=t['Melting point [°C]']
meta=pd.read_csv('data/meltome_datasets.csv'); meta['dataset']=meta['Dataset ID'].str.replace(':','_')
t=t.merge(meta[['dataset','Dataset source']], on='dataset', how='left')
euk=['Homo sapiens','Mus musculus','Danio rerio','Drosophila melanogaster','Caenorhabditis elegans','Saccharomyces cerevisiae']
t=t[t.Organism.isin(euk)]
t['gene']=np.where(t.Organism=='Homo sapiens', t.sym.str.upper(), np.where(t.one2one, t.human_syms, None))
t=t[t.gene.notna()]
ox=t.pathway.fillna('').str.contains('OXPHOS')
sets={'matrix (all)':t.mc_loc=='Matrix',
      'matrix, soluble (non-OXPHOS)':(t.mc_loc=='Matrix')&~ox,
      'inner membrane (all)':t.mc_loc=='MIM',
      'inner membrane, non-OXPHOS':(t.mc_loc=='MIM')&~ox,
      'OXPHOS subunits':ox&t.mc_loc.isin(['Matrix','MIM','mito_other'])}
t['cond']=t['Dataset source']+'|'+t.Organism+'|'+t['Strain/tissue'].astype(str)+'|'+t['Cells/lysate']+'|'+t['Lysis conditions']
def wmedian(x,w):
    o=np.argsort(x); x=x[o]; w=w[o]; c=np.cumsum(w); return x[np.searchsorted(c, c[-1]/2)]
out=open('results/subcompartments.txt','w')
def P(*a):
    s=' '.join(str(x) for x in a); print(s); out.write(s+'\n')
P('Median Tm difference vs cytosolic reference by mitochondrial set; mean over conditions, two-level bootstrap 95% CI (500 reps)')
rows=[]
for name,sel in sets.items():
    d=t[sel|t.cyto_ref].copy(); d['grp']=np.where(sel[d.index],'M','C')
    g=d.groupby(['cond','dataset','gene','grp']).agg(Tm=('Tm','mean')).reset_index().dropna(subset=['Tm'])
    genes=sorted(g.gene.unique()); gi={x:i for i,x in enumerate(genes)}; g['gi']=g.gene.map(gi)
    blocks=[(c,dd.gi.values,dd.Tm.values,(dd.grp=='M').values) for (c,ds),dd in g.groupby(['cond','dataset'])]
    def est(w):
        res={}
        for c,gidx,tm,ism in blocks:
            ww=w[gidx]
            if ww[ism].sum()<10 or ww[~ism].sum()<30: continue
            res.setdefault(c,[]).append(wmedian(tm[ism],ww[ism]) - wmedian(tm[~ism],ww[~ism]))
        return {c:np.mean(v) for c,v in res.items()}
    E=pd.Series(est(np.ones(len(genes))))
    strata={'lysate':[c for c in E.index if '|lysate|' in c],'intact cells':[c for c in E.index if '|cells|' in c]}
    boot={k:[] for k in strata}
    for b in range(500):
        eb=est(rng.poisson(1.0,len(genes)).astype(float))
        for k,cs in strata.items():
            cs2=rng.choice(cs,len(cs)); boot[k].append(np.mean([eb[c] for c in cs2 if c in eb]))
    n=g[g.grp=='M'].gene.nunique()
    for k,cs in strata.items():
        m=np.mean([E[c] for c in cs]); lo,hi=np.percentile(boot[k],[2.5,97.5])
        rows.append((name,k,m,lo,hi,len(cs),n))
        P(f'  {name:30s} [{k:12s}] k={len(cs):2d} genes={n:4d}  dTm={m:+.2f} °C  95% CI [{lo:+.2f}, {hi:+.2f}]')
pd.DataFrame(rows,columns=['set','stratum','d','lo','hi','k','n_genes']).to_csv('results/subcompartments.csv',index=False)
out.close()
