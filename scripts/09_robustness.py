import pandas as pd, numpy as np
from scipy import stats
rng=np.random.default_rng(1)
t = pd.read_csv('data/meltome_labeled.csv', low_memory=False)
t = t[t.Organism!='Human body fluids']
t['Tm']=t['Melting point [°C]']
t['cond'] = t['Organism']+'|'+t['Strain/tissue'].astype(str)+'|'+t['Cells/lysate']+'|'+t['Lysis conditions']
go=t['Gene Ontology (cellular component)'].fillna('')
t['go_matrix']=go.str.contains('mitochondrial matrix')
t['oxphos']=t.pathway.str.contains('OXPHOS')
t['mitoribo']=t.pathway.str.contains('Mitochondrial ribosome')
t['tca']=t.pathway.str.contains('TCA cycle')
t['mtdna']=t.pathway.str.contains('mtDNA')
t['chap']=t.pathway.str.contains('Chaperones')
t['nonmelter']=t['Protein stability class']=='non-melter'
out=open('results/robustness.txt','w')
def P(*a):
    s=' '.join(str(x) for x in a); print(s); out.write(s+'\n')
def dl(y,se):
    w=1/se**2; ybar=np.sum(w*y)/w.sum(); Q=np.sum(w*(y-ybar)**2); k=len(y)
    tau2=max(0,(Q-(k-1))/(w.sum()-np.sum(w**2)/w.sum())) if k>1 else 0; ws=1/(se**2+tau2); est=np.sum(ws*y)/ws.sum(); return est, np.sqrt(1/ws.sum())
def per_cond(selector, metric='Tm', min_n=10):
    ys=[];ses=[]
    for cond,g in t.groupby('cond'):
        gg=g.groupby('Protein ID').agg(Tm=('Tm','mean'), sel=(selector,'first'), cy=('cyto_ref','first'), nm=('nonmelter','mean'))
        if metric=='rank':  # non-melters ranked as most stable
            gg['v']=gg.Tm.where(gg.nm<0.5, 1e3)
            gg=gg.dropna(subset=['v'])
            # convert to percentile within cytosolic reference: report median percentile shift in °C-equivalent via cyto quantiles
            a=gg[gg.sel==True].v; b=gg[gg.cy==True].v
            if len(a)<min_n or len(b)<50: continue
            # Hodges-Lehmann-free: probability of superiority
            u=stats.mannwhitneyu(a,b).statistic/(len(a)*len(b)); ys.append(u); ses.append(np.sqrt(u*(1-u)*(len(a)+len(b)+1)/(12*len(a)*len(b)))*np.sqrt(12)); continue
        gg=gg.dropna(subset=['Tm'])
        a=gg[gg.sel==True].Tm.values; b=gg[gg.cy==True].Tm.values
        if len(a)<min_n or len(b)<50: continue
        d=[np.median(rng.choice(a,len(a)))-np.median(rng.choice(b,len(b))) for _ in range(500)]
        ys.append(np.median(a)-np.median(b)); ses.append(np.std(d))
    if not ys: return None
    e,s=dl(np.array(ys),np.array(ses)); return e,s,len(ys)
t['is_matrix']=t.mc_loc=='Matrix'; t['is_mim']=t.mc_loc=='MIM'; t['is_mito_any']=t.mc_loc.isin(['Matrix','MIM','mito_other'])
P('Pooled (random-effects) median-Tm difference vs cytosolic reference, all eukaryote conditions')
robrows=[]
for name in ['is_matrix','is_mim','is_mito_any','go_matrix','oxphos','mitoribo','tca','mtdna','chap']:
    r=per_cond(name)
    if r: robrows.append((name,r[0],r[0]-1.96*r[1],r[0]+1.96*r[1],r[2]))
    if r: P(f'  {name:12s} k={r[2]:2d}  dTm={r[0]:+.2f} °C  95%CI [{r[0]-1.96*r[1]:+.2f}, {r[0]+1.96*r[1]:+.2f}]')
pd.DataFrame(robrows,columns=['set','est','lo','hi','k']).to_csv('results/robustness_pooled.csv',index=False)
P('\nProbability of superiority (matrix more stable than cytosolic), non-melters treated as most stable; 0.5 = no difference')
r=per_cond('is_matrix', metric='rank'); P(f'  matrix: P(sup)={r[0]:.3f} ± {1.96*r[1]:.3f} (k={r[2]})')
# non-melter fractions
nm = t.groupby(['cond']).apply(lambda g: pd.Series({'nm_matrix':g[g.is_matrix].nonmelter.mean(),'nm_cyto':g[g.cyto_ref].nonmelter.mean()}))
P('\nNon-melter fraction (matrix vs cytosolic):'); P(nm.dropna().round(3).to_string())
out.close()
