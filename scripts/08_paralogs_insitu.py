import pandas as pd, numpy as np
from scipy import stats
t = pd.read_csv('data/meltome_labeled.csv', low_memory=False)
t['Tm'] = t['Melting point [°C]']
out = open('results/CDE_results.txt','w')
def P(*a):
    s=' '.join(str(x) for x in a); print(s); out.write(s+'\n')
hum = t[t.Organism=='Homo sapiens'].copy(); hum['S']=hum.sym.str.upper()
# ---------- C. matrix/cytosol paralog pairs (human) ----------
pairs = [('IDH2','IDH1'),('GOT2','GOT1'),('MDH2','MDH1'),('ACO2','ACO1'),('SHMT2','SHMT1'),('ME2','ME1'),('ME3','ME1'),('ALDH2','ALDH1A1'),
 ('CKMT1A','CKB'),('CKMT2','CKM'),('TXN2','TXN'),('TXNRD2','TXNRD1'),('PRDX3','PRDX1'),('ACAT1','ACAT2'),('HMGCS2','HMGCS1'),('PCK2','PCK1'),
 ('GPT2','GPT'),('ACSS1','ACSS2'),('NME4','NME1'),('HSPA9','HSPA8'),('TRAP1','HSP90AA1'),('TRAP1','HSP90AB1'),('TUFM','EEF1A1'),('GFM1','EEF2'),
 ('DARS2','DARS1'),('LARS2','LARS1'),('IARS2','IARS1'),('VARS2','VARS1'),('FARS2','FARSA'),('YARS2','YARS1'),('WARS2','WARS1'),('HARS2','HARS1'),
 ('NARS2','NARS1'),('SARS2','SARS1'),('TARS2','TARS1'),('RARS2','RARS1'),('EARS2','EPRS1'),('MTHFD2','MTHFD1'),('MTHFD1L','MTHFD1'),('ACSF3','ACSF2'),
 ('SOD2','SOD1'),('GLRX5','GLRX3'),('PPA2','PPA1'),('ADK','ADK'),('AK2','AK1'),('AK3','AK1'),('AK4','AK1'),('GLUD1','GLUD1'),('NFS1','NFS1'),
 ('ALDH4A1','ALDH1A1'),('ALDH18A1','ALDH1A1'),('ECHS1','ECH1'),('HADHA','EHHADH'),('ACADM','ACOX1'),('SUCLG1','SUCLG1'),('CS','CS'),
 ('PDHA1','PDHA1'),('DLD','DLD'),('POLG','POLA1'),('TOP1MT','TOP1'),('TFB2M','DIMT1'),('MRPL12','RPLP1'),('PITRM1','IDE'),('LONP1','LONP1'),
 ('CLPX','CLPX'),('ATAD3A','ATAD3A'),('HSPE1','HSPE1'),('PPIF','PPIA'),('GRPEL1','GRPEL1'),('DNAJA3','DNAJA1'),('GCDH','ACADM'),('GATM','GATM'),
 ('PDK1','PDK1'),('NIT2','NIT1'),('ACOT2','ACOT1'),('ACOT13','ACOT13'),('HAGH','HAGH'),('GSR','GSR'),('GPX4','GPX1'),('MSRB2','MSRB1'),('ARG2','ARG1'),
 ('CA5B','CA2'),('GLS','GLS2'),('ENDOG','DNASE1'),('PNPT1','EXOSC10'),('SUPV3L1','SKIV2L'),('MTIF2','EIF5B'),('TSFM','EEF1B2'),('MTG1','GNL2'),
 ('DTYMK','DTYMK'),('TK2','DCK'),('DGUOK','DCK'),('RRM2B','RRM2'),('NUDT2','NUDT2'),('PNKD','HAGH'),('ETHE1','HAGH'),('LACTB2','HAGH'),('HIBCH','ECH1'),('AUH','ECH1')]
pairs = [(a,b) for a,b in pairs if a!=b]
for state in ['lysate','cells']:
    h = hum[hum['Cells/lysate']==state]
    # per-dataset z-free approach: within each dataset, difference for each pair, then average across datasets
    recs=[]
    for ds,g in h.groupby('dataset'):
        m = g.groupby('S').Tm.mean()
        for a,b in pairs:
            if a in m and b in m and pd.notna(m[a]) and pd.notna(m[b]): recs.append((a,b,ds,m[a]-m[b]))
    R = pd.DataFrame(recs, columns=['mito','cyto','ds','d'])
    if len(R)==0: continue
    pr = R.groupby(['mito','cyto']).d.agg(['mean','count']).reset_index()
    w = stats.wilcoxon(pr['mean'])
    P(f'C. paralog pairs [{state}]: n_pairs={len(pr)}, median dTm(mito-cyto)={pr["mean"].median():.2f}, mean={pr["mean"].mean():.2f} ± {pr["mean"].std()/np.sqrt(len(pr)):.2f} (SEM), Wilcoxon p={w.pvalue:.3g}, frac mito more stable={np.mean(pr["mean"]>0):.2f}')
    P(pr.sort_values('mean').round(2).to_string(index=False))
    pr.to_csv(f'results/C_paralogs_{state}.csv', index=False)

# ---------- D. in situ: intact cells vs lysate, same cell line & lysis ----------
P('\nD. Per-protein shift Tm(intact cells) - Tm(lysate); matrix vs cytosolic reference')
for line,lysis in [('K562','0.5% NP40'),('K562','PBS'),('Jurkat','PBS/NP40')]:
    if line=='Jurkat':
        c = hum[(hum['Strain/tissue']=='Jurkat')&(hum['Cells/lysate']=='cells')]; l = hum[(hum['Strain/tissue']=='Jurkat')&(hum['Cells/lysate']=='lysate')]
    else:
        c = hum[(hum['Strain/tissue']==line)&(hum['Cells/lysate']=='cells')&(hum['Lysis conditions']==lysis)]
        l = hum[(hum['Strain/tissue']==line)&(hum['Cells/lysate']=='lysate')&(hum['Lysis conditions']==lysis)]
    if line=='K562' and lysis=='PBS': c = c[c.dataset.isin(['ma_0028','ma_0029'])]
    cm = c.groupby('Protein ID').agg(Tm=('Tm','mean'), mc=('mc_loc','first'), cy=('cyto_ref','first')); lm = l.groupby('Protein ID').Tm.mean()
    j = cm.join(lm.rename('Tl'), how='inner').dropna(); j['sh']=j.Tm-j.Tl
    a = j[j.mc=='Matrix'].sh; b = j[j.cy].sh
    P(f'  {line} {lysis}: n_matrix={len(a)}, n_cyto={len(b)}, median shift matrix={a.median():.2f}, cyto={b.median():.2f}, relative={a.median()-b.median():.2f} °C, MWU p={stats.mannwhitneyu(a,b).pvalue:.2g}')
    # restricted to proteins melting in the 40-50 window in lysate (where respiration would still be active)
    jj = j[(j.Tl>=40)&(j.Tl<=50)]; a2=jj[jj.mc=='Matrix'].sh; b2=jj[jj.cy].sh
    P(f'     lysate Tm 40-50 only: n_matrix={len(a2)}, relative shift={a2.median()-b2.median():.2f} °C')

# ---------- E. fraction of matrix proteins melted at ambient+10 °C ----------
P('\nE. Fraction of proteins with Tm below (growth temp + ΔT)')
amb = {'Homo sapiens':37,'Mus musculus':37,'Danio rerio':28,'Drosophila melanogaster':25,'Caenorhabditis elegans':20,'Saccharomyces cerevisiae':30}
for org,T0 in amb.items():
    g = t[(t.Organism==org)]
    for state in ['lysate','cells']:
        gs = g[g['Cells/lysate']==state]
        if len(gs)==0: continue
        pm = gs.groupby('Protein ID').agg(Tm=('Tm','mean'), mc=('mc_loc','first'), cy=('cyto_ref','first'), cls=('Protein stability class', lambda x: x.mode().iloc[0]))
        # treat non-melters as not melted (conservative toward hot-mito)
        for dT in [0,5,10,13]:
            thr=T0+dT
            fm = ((pm[pm.mc=='Matrix'].Tm<thr)).mean(); fc=((pm[pm.cy].Tm<thr)).mean()
            P(f'  {org:26s} {state:6s} T={thr:4.0f}°C (ambient+{dT:2d}): matrix {fm*100:5.1f}%  cytosolic {fc*100:5.1f}%  (n_matrix={int((pm.mc=="Matrix").sum())})')
out.close()
