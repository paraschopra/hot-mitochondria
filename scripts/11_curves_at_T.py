import json, pandas as pd, numpy as np
lab = pd.read_csv('data/meltome_labeled.csv', low_memory=False, usecols=['sym','Organism','mc_loc','cyto_ref'])
lab = lab[lab.Organism=='Homo sapiens'].drop_duplicates('sym'); lab['S']=lab.sym.str.upper(); L = lab.drop_duplicates('S').set_index('S')
d = json.load(open('data/full_dataset.json'))
lines = ['HEK293T','U937','HL60','HAOEC','HaCaT','pTcells','HepG2']
rec=[]
for x in d:
    if x['runName'] not in lines: continue
    pid=x['proteinId']; sym=pid.split('_')[-1].upper() if '_' in pid else pid.upper()
    if sym not in L.index: continue
    mb = x['meltingBehaviour']; T=np.array([m['temperature'] for m in mb]); F=np.array([m['fold_change'] for m in mb])
    o=np.argsort(T); T=T[o]; F=F[o]
    def at(t0): return np.interp(t0, T, F) if T.min()<=t0<=T.max() else np.nan
    rec.append((x['runName'], sym, L.loc[sym,'mc_loc'], L.loc[sym,'cyto_ref'], at(41), at(44), at(47), at(50)))
R = pd.DataFrame(rec, columns=['line','sym','mc','cyto','f41','f44','f47','f50'])
R = R.groupby(['sym','mc','cyto'])[['f41','f44','f47','f50']].median().reset_index()
with open('results/curves_at_T.txt','w') as out:
    for grp,sel in [('matrix',R.mc=='Matrix'),('MIM',R.mc=='MIM'),('cytosolic',R.cyto==True)]:
        s=R[sel]; line=f'{grp:9s} n={len(s):5d}  median soluble fraction after 3 min at 41/44/47/50 °C: ' + ' / '.join(f'{s[c].median():.2f}' for c in ['f41','f44','f47','f50']) + f'   | frac proteins <50% soluble at 47: {(s.f47<0.5).mean():.2f}, at 50: {(s.f50<0.5).mean():.2f}'
        print(line); out.write(line+'\n')
