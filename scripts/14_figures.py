import pandas as pd, numpy as np, matplotlib, sys; matplotlib.use('Agg')
GRID = 'grid' in sys.argv   # 2x2 portrait layout for the PDF summary
import matplotlib.pyplot as plt
from scipy import stats
BLUE='#2a78d6'; ORANGE='#eb6834'; GRAY='#8a8985'; INK='#0b0b0b'; INK2='#52514e'
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':INK2,'axes.labelcolor':INK,'xtick.color':INK2,'ytick.color':INK2})
t = pd.read_csv('data/meltome_labeled.csv', low_memory=False); t['Tm']=t['Melting point [°C]']
lys = t[(t['Cells/lysate']=='lysate')&(t['Lysis conditions']=='0.5% NP40')&(t.Organism!='Human body fluids')]
cal = lys.groupby('Organism').agg(OGT=('OGT [°C]','first'), medTm=('Tm','median')).reset_index()
fit = stats.linregress(cal.OGT, cal.medTm)
B = pd.read_csv('results/B_matrix_vs_cyto.csv')
if GRID:
    fig, ax = plt.subplots(2,2, figsize=(10,10.5), gridspec_kw={'height_ratios':[1.3,1]}); ax=ax.flatten()
else:
    fig, ax = plt.subplots(1,4, figsize=(17,5.6), gridspec_kw={'width_ratios':[1,1.35,1,0.9]})
# A calibration
a=ax[0]; a.scatter(cal.OGT, cal.medTm, s=28, color=BLUE, zorder=3)
x=np.linspace(10,75,10); a.plot(x, fit.intercept+fit.slope*x, color=GRAY, lw=1.5)
for _,r in cal.iterrows():
    a.annotate(r.Organism.split()[0][0]+'. '+r.Organism.split()[1], (r.OGT, r.medTm), fontsize=6.5, color=INK2, xytext=(3,-3), textcoords='offset points')
a.set_xlabel('Growth temperature (°C)'); a.set_ylabel('Median proteome Tm (°C)')
a.set_title(f'A  Proteome Tm tracks operating temperature\nslope {fit.slope:.2f} °C/°C (r={fit.rvalue:.2f}, 13 species)', loc='left', fontsize=9.5)
# B forest: revised 25 conditions, grouped lysate / intact cells, each group sorted by value
B=pd.read_csv('results/conditions_revised.csv')
sp={'Homo sapiens':'human','Mus musculus':'mouse','Danio rerio':'zebrafish','Drosophila melanogaster':'fly','Caenorhabditis elegans':'worm','Saccharomyces cerevisiae':'yeast'}
def lab_(c):
    lab_src,org,cell,state,buf=c.split('|'); cell=cell.replace('SII cells','S2 cells').replace('colon cancer spheroids','colon spheroids')
    return f"{cell} · {sp[org]} · {lab_src}"+(' · PBS' if (state=='lysate' and buf=='PBS') else '')
B['state']=np.where(B.cond.str.contains('|lysate|',regex=False),'lysate','cells')
pooled={'lysate':(-0.17,-1.00,0.58),'cells':(1.26,0.72,1.77)}
rows=[]
for st,title in [('cells','Intact cells (13)'),('lysate','Lysates: intrinsic stability (12)')]:
    d=B[B.state==st].sort_values('d',ascending=False)
    rows.append(('head',title,None,None,None))
    for _,r in d.iterrows(): rows.append(('pt',lab_(r.cond),r.d,r.lo,r.hi))
    rows.append(('pool','Pooled',*pooled[st]))
a=ax[1]; y=np.arange(len(rows))[::-1]
a.axvspan(0.55*10, 0.81*15, color=ORANGE, alpha=0.15, lw=0); a.text(0.55*10+0.2, y[0]+0.3, 'prediction if matrix\nruns 10–15 °C hotter', color=ORANGE, fontsize=7.5, va='top')
a.axvline(0,color=GRAY,lw=1)
labels=[]
for yi,(k,n,e,lo,hi) in zip(y,rows):
    if k=='head': labels.append(n); continue
    if k=='pool': a.errorbar(e, yi, xerr=[[e-lo],[hi-e]], fmt='D', ms=6, color=INK, elinewidth=1.8)
    else: a.errorbar(e, yi, xerr=[[e-lo],[hi-e]], fmt='o', ms=3.5, color=BLUE, elinewidth=1.1)
    labels.append(n)
a.set_yticks(y); a.set_yticklabels(labels, fontsize=6.2)
for tl,(k,*_) in zip(a.get_yticklabels(),rows):
    if k in ('head','pool'): tl.set_fontweight('bold'); tl.set_color(INK)
a.set_xlim(-3.5,12.5); a.set_xlabel('Matrix − cytosolic median Tm (°C)')
a.set_title('B  Matrix proteins are not more stable\n25 conditions, sorted within each group', loc='left', fontsize=9.5)
# C paralogs + pathways
RB=pd.read_csv('results/robustness_pooled.csv').set_index('set')
names={'is_matrix':'Matrix (MitoCarta)','go_matrix':'Matrix (GO)','is_mim':'Inner membrane','oxphos':'OXPHOS','tca':'TCA cycle','mitoribo':'Mitoribosome','mtdna':'mtDNA maintenance'}
rob=[(v,RB.loc[k,'est'],RB.loc[k,'lo'],RB.loc[k,'hi']) for k,v in names.items()]
rob.append(('Curated paralogs (n=47)', -0.45, -1.51, 0.60))
a=ax[2]; yy=np.arange(len(rob))[::-1]
a.axvspan(0.55*10, 0.81*15, color=ORANGE, alpha=0.15, lw=0); a.axvline(0,color=GRAY,lw=1)
for yi,(n,e,lo,hi) in zip(yy,rob): a.errorbar(e, yi, xerr=[[e-lo],[hi-e]], fmt='s', ms=4.5, color=BLUE, elinewidth=1.4)
a.set_yticks(yy); a.set_yticklabels([r[0] for r in rob], fontsize=7.5); a.set_xlim(-3.5,12.5)
a.set_xlabel('ΔTm vs cytosolic reference (°C)')
a.set_title('C  No sub-compartment, pathway or\nparalog set reaches the predicted band', loc='left', fontsize=9.5)
a=ax[3]; seq=[('Human',-1.7,-4.2,1.0),('Zebrafish',-6.2,-10.9,-1.3),('Fly (S2)',-2.9,-8.6,2.8),('Worm',-1.3,-4.6,2.2),('Yeast',-4.5,-9.2,0.1)]
yy=np.arange(len(seq))[::-1]; a.axvspan(10,15,color=ORANGE,alpha=0.15,lw=0); a.axvline(0,color=GRAY,lw=1)
for yi,(n,e,lo,hi) in zip(yy,seq): a.errorbar(e, yi, xerr=[[e-lo],[hi-e]], fmt='D', ms=4.5, color=BLUE, elinewidth=1.4)
a.set_yticks(yy); a.set_yticklabels([s[0] for s in seq], fontsize=8); a.set_xlim(-12,17)
a.set_xlabel('Implied matrix − cytosol\noperating temperature (°C)')
a.set_title('D  Sequence thermometer\n(TemStaPro, 5 species)', loc='left', fontsize=9.5)
plt.tight_layout(); plt.savefig('figures/fig1_grid.png' if GRID else 'figures/fig1_main.png', dpi=200); print('ok')
