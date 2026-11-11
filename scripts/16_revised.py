"""Revised analysis after independent review: per-dataset estimates, lab-aware conditions, gene-cluster bootstrap,
censoring-aware ranks, required-slope inference with calibration uncertainty."""
import pandas as pd, numpy as np
from scipy import stats
rng=np.random.default_rng(42)
t = pd.read_csv('data/meltome_labeled.csv', low_memory=False)
t = t[(t.Organism!='Human body fluids')]
t['Tm']=t['Melting point [°C]']; t['nonmelter']=t['Protein stability class']=='non-melter'
meta=pd.read_csv('data/meltome_datasets.csv'); meta['dataset']=meta['Dataset ID'].str.replace(':','_')
t=t.merge(meta[['dataset','Dataset source','Melting range [°C]']], on='dataset', how='left')
euk=['Homo sapiens','Mus musculus','Danio rerio','Drosophila melanogaster','Caenorhabditis elegans','Saccharomyces cerevisiae']
t=t[t.Organism.isin(euk)]
# gene key shared across species: human symbol (nonhuman only if strict 1:1)
t['gene']=np.where(t.Organism=='Homo sapiens', t.sym.str.upper(), np.where(t.one2one, t.human_syms, None))
t=t[t.gene.notna()]
t['grp']=np.where(t.mc_loc=='Matrix','M', np.where(t.cyto_ref,'C', None)); t=t[t.grp.notna()]
t['cond']=t['Dataset source']+'|'+t.Organism+'|'+t['Strain/tissue'].astype(str)+'|'+t['Cells/lysate']+'|'+t['Lysis conditions']
out=open('results/revised_results.txt','w')
def P(*a):
    s=' '.join(str(x) for x in a); print(s); out.write(s+'\n')
# per dataset, per gene Tm (mean over protein isoforms in that dataset)
g=t.groupby(['cond','dataset','gene','grp']).agg(Tm=('Tm','mean'), nm=('nonmelter','mean')).reset_index()
conds=sorted(g.cond.unique()); genes=sorted(g.gene.unique()); gi={x:i for i,x in enumerate(genes)}
g['gi']=g.gene.map(gi)
def wmedian(x,w):
    o=np.argsort(x); x=x[o]; w=w[o]; c=np.cumsum(w); return x[np.searchsorted(c, c[-1]/2)]
blocks=[]
for c in conds:
    for ds,d in g[g.cond==c].groupby('dataset'):
        d=d.dropna(subset=['Tm'])
        blocks.append((c, ds, d.gi.values, d.Tm.values, (d.grp=='M').values))
def cond_estimates(weights):
    res={}
    for c,ds,gidx,tm,ism in blocks:
        w=weights[gidx]
        if w[ism].sum()<10 or w[~ism].sum()<30: continue
        res.setdefault(c,[]).append(wmedian(tm[ism],w[ism]) - wmedian(tm[~ism],w[~ism]))
    return {c:np.mean(v) for c,v in res.items()}
est=cond_estimates(np.ones(len(genes)))
E=pd.Series(est); P('Per-condition matrix-cytosol median Tm difference (replicate datasets averaged within lab x condition):'); P(E.round(2).to_string())
strata={'all':list(E.index),'lysate':[c for c in E.index if '|lysate|' in c],'intact cells':[c for c in E.index if '|cells|' in c],
        'non-mammal':[c for c in E.index if not ('Homo' in c or 'Mus' in c)]}
B=1000; boot={k:[] for k in strata}
for b in range(B):
    w=rng.poisson(1.0,len(genes)).astype(float)      # gene-cluster bootstrap shared across all datasets/species
    eb=cond_estimates(w)
    for k,cs in strata.items():
        cs2=rng.choice(cs,len(cs))                    # resample conditions too
        boot[k].append(np.mean([eb[c] for c in cs2 if c in eb]))
summ={}
for k,cs in strata.items():
    m=np.mean([E[c] for c in cs]); lo,hi=np.percentile(boot[k],[2.5,97.5]); summ[k]=(m,lo,hi,len(cs))
    P(f'  [{k:12s}] k={len(cs):2d} mean dTm={m:+.2f} °C  two-level bootstrap 95% CI [{lo:+.2f}, {hi:+.2f}]')
# ---- calibration with uncertainty: species bootstrap of slope ----
lys=pd.read_csv('data/meltome_labeled.csv', low_memory=False)
lys=lys[(lys['Cells/lysate']=='lysate')&(lys['Lysis conditions']=='0.5% NP40')&(lys.Organism!='Human body fluids')]
cal=lys.groupby('Organism').agg(OGT=('OGT [°C]','first'), medTm=('Melting point [°C]','median')).reset_index()
eu=cal[cal.Organism.isin(euk+['Arabidopsis thaliana'])]
def bslope(df):
    s=[]
    for _ in range(4000):
        i=rng.integers(0,len(df),len(df)); d=df.iloc[i]
        if d.OGT.nunique()<3: continue
        s.append(stats.linregress(d.OGT,d.medTm).slope)
    return np.array(s)
sa=bslope(cal); se=bslope(eu)
P(f'\nCalibration slope (species bootstrap): all 13 species median {np.median(sa):.2f} [{np.percentile(sa,2.5):.2f}, {np.percentile(sa,97.5):.2f}]; eukaryotes median {np.median(se):.2f} [{np.percentile(se,2.5):.2f}, {np.percentile(se,97.5):.2f}]')
for k in ['lysate','all','intact cells']:
    hi=summ[k][2]
    for dT in [10,15]:
        req=hi/dT
        P(f'  [{k}] hot matrix +{dT} °C is compatible only if coupling slope < {req:.3f}; P(slope<req): all-species {np.mean(sa<req):.3f}, eukaryotes {np.mean(se<req):.3f}')
# ---- censoring-aware: probability of superiority, non-melters = +inf (Hanley-McNeil SE) ----
P('\nProbability of superiority P(matrix Tm > cytosolic Tm), non-melters treated as most stable:')
ys=[];ses=[]
for c,d in g.groupby('cond'):
    v=d.Tm.where(d.nm<0.5, np.inf); d=d.assign(v=v).dropna(subset=['v'])
    a=d[d.grp=='M'].groupby('gene').v.mean(); bb=d[d.grp=='C'].groupby('gene').v.mean()
    if len(a)<10: continue
    u=stats.mannwhitneyu(a,bb).statistic/(len(a)*len(bb)); q1=u/(2-u); q2=2*u*u/(1+u)
    s=np.sqrt((u*(1-u)+(len(a)-1)*(q1-u*u)+(len(bb)-1)*(q2-u*u))/(len(a)*len(bb))); ys.append(u); ses.append(s)
ys=np.array(ys); ses=np.array(ses); w=1/ses**2; Q=np.sum(w*(ys-np.sum(w*ys)/w.sum())**2); tau2=max(0,(Q-(len(ys)-1))/(w.sum()-np.sum(w**2)/w.sum())); ws=1/(ses**2+tau2)
m=np.sum(ws*ys)/ws.sum(); s=np.sqrt(1/ws.sum()); P(f'  pooled P(sup)={m:.3f} [{m-1.96*s:.3f}, {m+1.96*s:.3f}] over {len(ys)} conditions (0.5 = no difference)')
out.close()
