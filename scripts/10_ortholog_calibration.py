"""Ortholog-paired calibration: for the same protein family, how much does Tm change per °C of growth temperature?"""
import pandas as pd, numpy as np
from scipy import stats
t = pd.read_csv('data/meltome_labeled.csv', low_memory=False)
t['Tm']=t['Melting point [°C]']
lys = t[(t['Cells/lysate']=='lysate')&(t['Lysis conditions']=='0.5% NP40')]
hum = lys[lys.Organism=='Homo sapiens'].groupby(lys.sym.str.upper()).Tm.median()
ogt = {'Danio rerio':28,'Drosophila melanogaster':28,'Caenorhabditis elegans':20,'Saccharomyces cerevisiae':30,'Mus musculus':37}
out=open('results/ortholog_calibration.txt','w')
rows=[]
for org,T in ogt.items():
    g = lys[lys.Organism==org].dropna(subset=['Tm'])
    recs=[]
    for hs,tm,mc in zip(g.human_syms, g.Tm, g.mc_loc):
        H=[h for h in str(hs).split(';') if h in hum.index and pd.notna(hum[h])]
        if len(H)==1: recs.append((hum[H[0]]-tm, mc))
    R=pd.DataFrame(recs,columns=['d','mc'])
    dall=R.d.median(); dmat=R[R.mc=='Matrix'].d.median(); dnon=R[R.mc=='nonmito'].d.median()
    rows.append((org,37-T,len(R),dall,dmat,dnon))
    s=f'{org:26s} ΔOGT={37-T:2d}  n_orth={len(R):4d}  median ΔTm(human-species): all={dall:+.2f}  matrix={dmat:+.2f}  nonmito={dnon:+.2f}'; print(s); out.write(s+'\n')
R=pd.DataFrame(rows,columns=['org','dT','n','d','dmat','dnon'])
# regression through origin of ortholog ΔTm on ΔOGT
b=np.sum(R.dT*R.d)/np.sum(R.dT**2); s=f'ortholog-paired slope (through origin) = {b:.3f} °C Tm per °C'; print(s); out.write(s+'\n')
r=stats.linregress(R.dT,R.d); s=f'OLS slope={r.slope:.3f}±{r.stderr:.3f}, intercept={r.intercept:.2f}'; print(s); out.write(s+'\n')
