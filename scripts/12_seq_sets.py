import pandas as pd, numpy as np
rng=np.random.default_rng(0)
t = pd.read_csv('data/meltome_labeled.csv', low_memory=False, usecols=['Organism','sym','acc','mc_loc','cyto_ref','human_syms','Sequence'])
hs = pd.read_csv('data/uniprot_human_sp_seq.tsv', sep='\t').dropna(subset=['Gene Names (primary)'])
hseq = dict(zip(hs['Gene Names (primary)'].str.upper(), hs.Sequence))
mc = pd.read_csv('data/mitocarta3_human.csv'); mloc = dict(zip(mc.Symbol.str.upper(), mc['MitoCarta3.0_SubMitoLocalization']))
H = t[t.Organism=='Homo sapiens'].copy(); H['S']=H.sym.str.upper()
cy = sorted(set(H[H.cyto_ref].S) & set(hseq)); mx = sorted(set(H[H.mc_loc=='Matrix'].S) & set(hseq))
# orthologs measured in other species
rec=[]
for org in ['Caenorhabditis elegans','Drosophila melanogaster','Saccharomyces cerevisiae','Danio rerio']:
    g=t[(t.Organism==org)].dropna(subset=['Sequence']).drop_duplicates('acc')
    for hsy,seq,acc in zip(g.human_syms, g.Sequence, g.acc):
        Hs=[h for h in str(hsy).split(';') if h]
        if len(Hs)==1 and Hs[0] in hseq: rec.append((org, acc, Hs[0], seq))
O=pd.DataFrame(rec, columns=['org','id','human','seq'])
O['grp']=O.human.map(lambda h: 'matrix' if mloc.get(h)=='Matrix' else ('cyto' if h in set(cy) else 'other'))
O=O[O.grp!='other']
# keep human genes with >=1 ortholog preferentially, cap sizes
orth_h = set(O.human)
mx_sel = [g for g in mx if g in orth_h] ; mx_sel += list(rng.choice([g for g in mx if g not in orth_h], size=min(150,len([g for g in mx if g not in orth_h])), replace=False))
cy_o = [g for g in cy if g in orth_h]; cy_sel = list(rng.choice(cy_o, size=min(500,len(cy_o)), replace=False)) + list(rng.choice([g for g in cy if g not in orth_h], 150, replace=False))
O = O[O.human.isin(set(mx_sel)|set(cy_sel))]
rows=[('Homo sapiens', g, g, 'matrix', hseq[g]) for g in mx_sel]+[('Homo sapiens', g, g, 'cyto', hseq[g]) for g in cy_sel]
rows += [(r.org, r.id, r.human, r.grp, r.seq) for r in O.itertuples()]
S=pd.DataFrame(rows, columns=['org','id','human','grp','seq'])
S['seq']=S.seq.str[40:1040]   # drop N-terminal 40 aa (presequence) for all; cap length
S=S[S.seq.str.len()>=50]
S.to_csv('data/seq_set.csv', index=False)
print(S.groupby(['org','grp']).size())
