import pandas as pd, numpy as np
from scipy import stats
rng = np.random.default_rng(0)
t = pd.read_csv('data/meltome_labeled.csv', low_memory=False)
t['Tm'] = t['Melting point [°C]']; t['len'] = t['Length']
t['cond'] = t['Organism'] + '|' + t['Strain/tissue'].astype(str) + '|' + t['Cells/lysate'] + '|' + t['Lysis conditions']
out = open('results/main_results.txt','w')
def P(*a):
    s=' '.join(str(x) for x in a); print(s); out.write(s+'\n')

# ---------- A. calibration: proteome median Tm vs growth temperature ----------
lys = t[(t['Cells/lysate']=='lysate') & (t['Lysis conditions']=='0.5% NP40') & (t.Organism!='Human body fluids')]
cal = lys.groupby('Organism').agg(OGT=('OGT [°C]','first'), medTm=('Tm','median'), n=('Tm','count')).reset_index()
fit = stats.linregress(cal.OGT, cal.medTm)
euk = cal[cal.Organism.isin(['Saccharomyces cerevisiae','Arabidopsis thaliana','Caenorhabditis elegans','Drosophila melanogaster','Danio rerio','Mus musculus','Homo sapiens'])]
fit_e = stats.linregress(euk.OGT, euk.medTm)
P('A. Calibration (lysate, NP40):'); P(cal.round(2).to_string())
P(f'  all species: slope={fit.slope:.3f} ± {fit.stderr:.3f} °C Tm per °C OGT, r={fit.rvalue:.3f}')
P(f'  eukaryotes : slope={fit_e.slope:.3f} ± {fit_e.stderr:.3f}, r={fit_e.rvalue:.3f}')
slope = fit.slope

# ---------- B. matrix vs cytosolic reference, per condition ----------
def boot_med_diff(a,b,n=2000):
    a=np.asarray(a); b=np.asarray(b)
    d=[np.median(rng.choice(a,len(a)))-np.median(rng.choice(b,len(b))) for _ in range(n)]
    return np.percentile(d,[2.5,97.5]), np.std(d)
rows=[]
for cond,g in t[t.Organism!='Human body fluids'].groupby('cond'):
    # average replicates per protein within a condition
    gg = g.groupby('Protein ID').agg(Tm=('Tm','mean'), mc=('mc_loc','first'), cyto=('cyto_ref','first'), L=('len','first')).dropna(subset=['Tm'])
    a = gg[gg.mc=='Matrix'].Tm; b = gg[gg.cyto].Tm; m_all = gg[gg.mc.isin(['Matrix','MIM','mito_other'])].Tm
    if len(a)<20 or len(b)<50: continue
    ci, se = boot_med_diff(a,b)
    # length-adjusted (OLS on log length) difference
    h = gg[(gg.mc=='Matrix')|gg.cyto].dropna(subset=['L']).copy(); h['M']=(h.mc=='Matrix').astype(float); h['lL']=np.log(h.L)
    X = np.column_stack([np.ones(len(h)), h.M, h.lL]); beta,*_ = np.linalg.lstsq(X, h.Tm, rcond=None)
    rows.append(dict(cond=cond, n_matrix=len(a), n_cyto=len(b), med_matrix=a.median(), med_cyto=b.median(),
        d_matrix=a.median()-b.median(), ci_lo=ci[0], ci_hi=ci[1], se=se, d_allmito=m_all.median()-b.median(), d_lenadj=beta[1],
        frac_matrix_below50=(a<50).mean(), frac_cyto_below50=(b<50).mean(), p=stats.mannwhitneyu(a,b).pvalue))
B = pd.DataFrame(rows)
pd.set_option('display.width',300)
P('\nB. Matrix vs cytosolic Tm per condition (replicates averaged):'); P(B.round(3).to_string())
# random-effects (DerSimonian-Laird) pooled estimate
def dl(y,se):
    w=1/se**2; ybar=np.sum(w*y)/w.sum(); Q=np.sum(w*(y-ybar)**2); k=len(y)
    tau2=max(0,(Q-(k-1))/(w.sum()-np.sum(w**2)/w.sum())); ws=1/(se**2+tau2); est=np.sum(ws*y)/ws.sum(); s=np.sqrt(1/ws.sum())
    return est, s, tau2
for name, sub in [('all conditions',B), ('lysates',B[B.cond.str.contains('lysate')]), ('intact cells',B[B.cond.str.contains('cells') & ~B.cond.str.contains('lysate')]),
                  ('human',B[B.cond.str.startswith('Homo')]), ('non-mammal eukaryotes',B[~B.cond.str.startswith(('Homo','Mus'))])]:
    if len(sub)==0: continue
    est,s,tau2 = dl(sub.d_matrix.values, sub.se.values)
    P(f'  pooled [{name}] k={len(sub)}: dTm={est:.2f} °C (95% CI {est-1.96*s:.2f}..{est+1.96*s:.2f}), tau={np.sqrt(tau2):.2f}; implied ΔT={est/slope:.1f} °C (upper {(est+1.96*s)/slope:.1f}); max single-condition dTm={sub.d_matrix.max():.2f} -> ΔT={sub.d_matrix.max()/slope:.1f}')
B.to_csv('results/B_matrix_vs_cyto.csv', index=False)
out.close()
