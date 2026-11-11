"""Sequence 'thermometer' (TemStaPro: ProtT5 embeddings + MLPs trained on prokaryotic growth temperature).
Calibrate score-per-°C from human-vs-ortholog cytosolic pairs, then read the within-species matrix-minus-cytosolic
offset in every species (tests a universal hot matrix, not a human-only one)."""
import pandas as pd, numpy as np
from scipy import stats
rng=np.random.default_rng(5)
R=pd.read_csv('results/temstapro_scores.csv')
h=R[R.org=='Homo sapiens'].set_index('id').score
ogt={'Homo sapiens':37,'Caenorhabditis elegans':20,'Drosophila melanogaster':28,'Saccharomyces cerevisiae':30,'Danio rerio':28}  # Meltome OGTs
out=open('results/temstapro_summary.txt','w')
def P(s): print(s); out.write(s+'\n')
pairs={}
for org,T in ogt.items():
    if org=='Homo sapiens': continue
    o=R[(R.org==org)&(R.grp=='cyto')]; o=o[o.human.isin(h.index)]
    pairs[org]=(37-T, (h.loc[o.human].values-o.score.values))
def slope(sample=False):
    x=[];y=[]
    for org,(dT,d) in pairs.items():
        dd=rng.choice(d,len(d)) if sample else d; x.append(dT); y.append(dd.mean())
    x=np.array(x);y=np.array(y); return np.sum(x*y)/np.sum(x*x)
k=slope(); ks=np.array([slope(True) for _ in range(2000)])
for org,(dT,d) in pairs.items(): P(f'positive control {org:26s} ΔOGT={dT:2d}: human-minus-ortholog score (cytosolic) {d.mean():+.3f} ± {d.std()/np.sqrt(len(d)):.3f}, n={len(d)}')
P(f'calibration k = {k:.4f} score/°C (bootstrap 95% CI {np.percentile(ks,2.5):.4f}..{np.percentile(ks,97.5):.4f})')
rows=[]
for org in ogt:
    a=R[(R.org==org)&(R.grp=='matrix')].score.values; b=R[(R.org==org)&(R.grp=='cyto')].score.values
    bs=[rng.choice(a,len(a)).mean()-rng.choice(b,len(b)).mean() for _ in range(2000)]
    d=a.mean()-b.mean(); lo,hi=np.percentile(bs,[2.5,97.5]); rows.append((org,d,lo,hi))
    P(f'{org:26s} matrix-minus-cytosolic score {d:+.3f} [{lo:+.3f},{hi:+.3f}] -> implied {d/k:+.1f} °C [{lo/k:+.1f},{hi/k:+.1f}]  (hot matrix predicts +10..+15)')
D=np.array([r[1] for r in rows]); P(f'mean over 5 species: {D.mean():+.3f} -> {D.mean()/k:+.1f} °C; species with positive offset: {(D>0).sum()}/5')
