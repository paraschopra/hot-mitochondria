import pandas as pd, numpy as np, re
from scipy import stats
t = pd.read_csv('data/meltome_S4_long.csv')
meta = pd.read_csv('data/meltome_datasets.csv'); meta['dataset']=meta['Dataset ID'].str.replace(':','_')
t = t.merge(meta[['dataset','Organism','Strain/tissue','OGT [°C]','Cells/lysate','Lysis conditions']], on='dataset')
t['acc'] = t['Protein ID'].str.split('_').str[0].str.replace(r'-\d+$','',regex=True)
u = pd.read_csv('data/uniprot_annot.tsv', sep='\t').rename(columns={'Entry':'acc'})
t = t.merge(u, on='acc', how='left')
go = t['Gene Ontology (cellular component)'].fillna(''); loc = t['Subcellular location [CC]'].fillna('')
t['annot'] = (go!='') | (loc!='')
t['mito'] = go.str.contains('mitochondri') | loc.str.contains('Mitochondrion')
t['matrix'] = go.str.contains('mitochondrial matrix') | loc.str.contains('Mitochondrion matrix')
t['memb'] = go.str.contains('membrane') & ~go.str.contains('mitochondri')
t['cyto_only'] = (go.str.contains('cytosol|cytoplasm|nucleus|nucleoplasm')|loc.str.contains('Cytoplasm|Nucleus')) & ~t.mito & ~go.str.contains('extracellular|endoplasmic|Golgi|lysosom|vacuol|chloroplast|plastid|peroxisom')
t.to_csv('data/meltome_annotated.csv', index=False)
rows=[]
for ds,g in t.groupby('dataset'):
    g = g[g.annot]
    if len(g)<50: continue
    tm = g['Melting point [°C]']
    m = g[g.mito]; c = g[g.cyto_only]
    mx = g[g.matrix]
    r = dict(dataset=ds, org=g.Organism.iloc[0], tissue=g['Strain/tissue'].iloc[0], state=g['Cells/lysate'].iloc[0], lysis=g['Lysis conditions'].iloc[0], OGT=g['OGT [°C]'].iloc[0],
             n=len(g), n_mito=len(m), n_matrix=len(mx), n_cyto=len(c),
             tm_mito=m['Melting point [°C]'].median(), tm_cyto=c['Melting point [°C]'].median(), tm_matrix=mx['Melting point [°C]'].median(),
             auc_mito=m.AUC.median(), auc_cyto=c.AUC.median(), auc_matrix=mx.AUC.median(),
             nonmelt_mito=(m['Protein stability class']=='non-melter').mean(), nonmelt_cyto=(c['Protein stability class']=='non-melter').mean())
    r['dTm']=r['tm_mito']-r['tm_cyto']; r['dTm_matrix']=r['tm_matrix']-r['tm_cyto']; r['dAUC']=r['auc_mito']-r['auc_cyto']
    a=m['Melting point [°C]'].dropna(); b=c['Melting point [°C]'].dropna()
    r['p_MW'] = stats.mannwhitneyu(a,b).pvalue if len(a)>5 and len(b)>5 else np.nan
    rows.append(r)
R = pd.DataFrame(rows)
pd.set_option('display.width',300)
print(R[['dataset','org','tissue','state','lysis','n','n_mito','n_matrix','n_cyto','tm_mito','tm_cyto','dTm','dTm_matrix','dAUC','nonmelt_mito','nonmelt_cyto','p_MW']].round(3).to_string())
R.to_csv('results/within_species_mito_vs_cyto.csv', index=False)
