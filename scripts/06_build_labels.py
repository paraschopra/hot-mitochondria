"""Assign uniform mito labels: MitoCarta3 (human) propagated by Alliance orthology; cytosolic reference from GO."""
import pandas as pd, numpy as np, gzip
t = pd.read_csv('data/meltome_annotated.csv', low_memory=False); t['Organism']=t.Organism.str.strip()
mc = pd.read_csv('data/mitocarta3_human.csv')
mc_sym = dict(zip(mc.Symbol.str.upper(), mc['MitoCarta3.0_SubMitoLocalization']))
mc_path = dict(zip(mc.Symbol.str.upper(), mc['MitoCarta3.0_MitoPathways'].fillna('')))
# gene symbol per protein: UniProt gene_primary else suffix of Protein ID
suffix = t['Protein ID'].str.split('_',n=1).str[1]
t['sym'] = t['Gene Names (primary)'].fillna(suffix).astype(str)
nous = ~t['Protein ID'].astype(str).str.contains('_') & t['Gene Names (primary)'].isna()
t.loc[nous,'sym'] = t.loc[nous,'Protein ID'].astype(str)
# fill human GO/length by gene symbol (EMBL datasets use symbols, not accessions)
hs_ = pd.read_csv('data/uniprot_human_sp.tsv', sep='\t')
gmap={}; lmap={}
for _,r in hs_.iterrows():
    names=[str(r['Gene Names (primary)'])]+str(r['Gene Names (synonym)']).split()
    for n in names:
        n=n.upper()
        if n not in gmap or str(r['Gene Names (primary)']).upper()==n: gmap[n]=r['Gene Ontology (cellular component)']; lmap[n]=r['Length']
ish = (t.Organism=='Homo sapiens') & t['Gene Ontology (cellular component)'].isna()
t.loc[ish,'Gene Ontology (cellular component)'] = t.loc[ish,'sym'].str.upper().map(gmap)
t.loc[ish,'Length'] = t.loc[ish,'sym'].str.upper().map(lmap)
orth = pd.read_csv('data/alliance_orth.tsv.gz', sep='\t', comment='#')
orth = orth[(orth.Gene1SpeciesName=='Homo sapiens')]
sp_map = {'Mus musculus':'Mus musculus','Danio rerio':'Danio rerio','Drosophila melanogaster':'Drosophila melanogaster','Caenorhabditis elegans':'Caenorhabditis elegans','Saccharomyces cerevisiae':'Saccharomyces cerevisiae'}
h2 = {}
for sp in sp_map:
    o = orth[orth.Gene2SpeciesName.str.startswith(sp.split()[0])]
    d = {}
    for s2, s1 in zip(o.Gene2Symbol.astype(str), o.Gene1Symbol.astype(str)):
        d.setdefault(str(s2).upper(), set()).add(str(s1).upper())
    h2[sp] = d
rev = {}
for sp in sp_map:
    o = orth[orth.Gene2SpeciesName.str.startswith(sp.split()[0])]
    r = {}
    for s2, s1 in zip(o.Gene2Symbol.astype(str), o.Gene1Symbol.astype(str)):
        r.setdefault(str(s1).upper(), set()).add(str(s2).upper())
    rev[sp] = r
def human_syms(org, sym):
    s = str(sym).upper()
    if org=='Homo sapiens': return {s}
    if org in h2: return h2[org].get(s, set())
    return set()
lab=[]; hs=[]
for org, sym in zip(t.Organism, t.sym):
    H = human_syms(org, sym); hs.append(';'.join(sorted(H)))
    locs = {mc_sym[h] for h in H if h in mc_sym}
    if not H: lab.append('no_ortholog')
    elif any('Matrix' in l for l in locs): lab.append('Matrix')
    elif 'MIM' in locs: lab.append('MIM')
    elif locs: lab.append('mito_other')
    else: lab.append('nonmito')
t['human_syms']=hs; t['mc_loc']=lab
# strict: one human ortholog, and that human gene has exactly one ortholog in this species (1:1)
strict=[]
for org, sym, H, l in zip(t.Organism, t.sym, hs, lab):
    if org=='Homo sapiens': strict.append(True); continue
    H=[h for h in H.split(';') if h]
    strict.append(org in rev and len(H)==1 and len(rev[org].get(H[0],set()))==1)
t['one2one']=strict
go = t['Gene Ontology (cellular component)'].fillna('')
secretory = go.str.contains('extracellular|endoplasmic reticulum|Golgi|lysosom|vacuol|chloroplast|plastid|peroxisom|endosom|plasma membrane')
t['cyto_ref'] = (t.mc_loc=='nonmito') & go.str.contains('cytosol|cytoplasm|nucleus|nucleoplasm') & ~secretory & ~go.str.contains('mitochondri')
t['pathway'] = [ ';'.join(sorted({mc_path.get(h,'') for h in str(H).split(';') if h})) for H in hs]
t.to_csv('data/meltome_labeled.csv', index=False)
print(pd.crosstab(t.Organism, t.mc_loc))
print(t.groupby('Organism').cyto_ref.sum())
