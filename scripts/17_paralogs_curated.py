import pandas as pd, numpy as np
from scipy import stats
rng=np.random.default_rng(3)
# curated same-family pairs (mito member, cytosolic member); each protein appears once
pairs=[('IDH2','IDH1'),('GOT2','GOT1'),('MDH2','MDH1'),('ACO2','ACO1'),('SHMT2','SHMT1'),('ME2','ME1'),('ALDH2','ALDH1A1'),('CKMT1A','CKB'),('CKMT2','CKM'),
('TXN2','TXN'),('TXNRD2','TXNRD1'),('PRDX3','PRDX1'),('ACAT1','ACAT2'),('HMGCS2','HMGCS1'),('PCK2','PCK1'),('GPT2','GPT'),('ACSS1','ACSS2'),('NME4','NME1'),
('HSPA9','HSPA8'),('TRAP1','HSP90AA1'),('TUFM','EEF1A1'),('GFM1','EEF2'),('MTHFD2','MTHFD1'),('ACOT2','ACOT1'),('CA5B','CA2'),('ARG2','ARG1'),('PPIF','PPIA'),
('DNAJA3','DNAJA1'),('TOP1MT','TOP1'),('MTIF2','EIF5B'),('AK3','AK1'),('PITRM1','IDE'),('TFB2M','DIMT1'),('MSRB2','MSRB1'),('PPA2','PPA1'),('DGUOK','DCK'),
('DARS2','DARS1'),('LARS2','LARS1'),('IARS2','IARS1'),('VARS2','VARS1'),('FARS2','FARSA'),('YARS2','YARS1'),('WARS2','WARS1'),('HARS2','HARS1'),('NARS2','NARS1'),
('SARS2','SARS1'),('TARS2','TARS1'),('RARS2','RARS1'),('EARS2','EPRS1'),('CARS2','CARS1'),('MARS2','MARS1'),('AARS2','AARS1')]
t=pd.read_csv('data/meltome_labeled.csv', low_memory=False); t=t[t.Organism=='Homo sapiens']
meta=pd.read_csv('data/meltome_datasets.csv'); meta['dataset']=meta['Dataset ID'].str.replace(':','_')
t=t.merge(meta[['dataset','Dataset source']], on='dataset')
t['S']=t.sym.str.upper(); t['Tm']=t['Melting point [°C]']
t['cond']=t['Dataset source']+'|'+t['Strain/tissue'].astype(str)+'|'+t['Cells/lysate']+'|'+t['Lysis conditions']
loc=t.drop_duplicates('S').set_index('S').mc_loc
valid=[(a,b) for a,b in pairs if loc.get(a) in ('Matrix','MIM','mito_other') and loc.get(b)=='nonmito']
bad=[(a,b,loc.get(a),loc.get(b)) for a,b in pairs if (a,b) not in valid]
out=open('results/paralogs_curated.txt','w')
def P(s): print(s); out.write(s+'\n')
P(f'curated pairs {len(pairs)}, passing localization check {len(valid)}; excluded: {bad}')
for state in ['lysate','cells',None]:
    h=t if state is None else t[t['Cells/lysate']==state]
    per=h.groupby(['cond','dataset','S']).Tm.mean().groupby(['cond','S']).mean()  # replicate datasets averaged within condition
    rec=[]
    for c in per.index.get_level_values(0).unique():
        m=per.loc[c]
        for a,b in valid:
            if a in m and b in m and pd.notna(m[a]) and pd.notna(m[b]): rec.append((a,b,c,m[a]-m[b]))
    R=pd.DataFrame(rec,columns=['mito','cyto','cond','d']); pr=R.groupby(['mito','cyto']).d.mean()
    bs=[np.mean(rng.choice(pr.values,len(pr))) for _ in range(5000)]
    P(f'[{state or "all"}] n_pairs={len(pr)} mean ΔTm(mito-cyto)={pr.mean():+.2f} °C, bootstrap 95% CI [{np.percentile(bs,2.5):+.2f}, {np.percentile(bs,97.5):+.2f}], median {pr.median():+.2f}, '
      f'frac mito more stable {np.mean(pr>0):.2f}, Wilcoxon p={stats.wilcoxon(pr).pvalue:.2g}; pairs with ΔTm>+5: {int((pr>5).sum())}')
    if state is None: pr.round(2).sort_values().to_csv('results/paralogs_curated_pairs.csv')
