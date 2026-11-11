import torch, pandas as pd, numpy as np, re, sys, time, glob, os
sys.path.insert(0,'tools/TemStaPro')
from MLP import MLP_C2H2
from transformers import T5Tokenizer, T5EncoderModel
dev = 'mps'
tok = T5Tokenizer.from_pretrained('Rostlab/prot_t5_xl_half_uniref50-enc', do_lower_case=False)
model = T5EncoderModel.from_pretrained('Rostlab/prot_t5_xl_half_uniref50-enc', torch_dtype=torch.float32).to(dev).eval()
S = pd.read_csv('data/seq_set.csv')
embf = 'data/seq_emb.npy'
if os.path.exists(embf): E=np.load(embf)
else:
    E=np.zeros((len(S),1024),dtype=np.float32); t0=time.time()
    order = np.argsort(S.seq.str.len().values)
    bs=8
    for i in range(0,len(order),bs):
        idx=order[i:i+bs]
        seqs=[' '.join(re.sub(r'[UZOB]','X',s)) for s in S.seq.values[idx]]
        enc=tok(seqs, add_special_tokens=True, padding='longest', return_tensors='pt').to(dev)
        with torch.no_grad(): h=model(**enc).last_hidden_state
        for k,j in enumerate(idx):
            L=len(S.seq.values[j]); E[j]=h[k,:L].mean(0).float().cpu().numpy()
        if i % 200 == 0: print(i, f'{time.time()-t0:.0f}s', flush=True)
    np.save(embf,E)
print('nan rows', np.isnan(E).any(1).sum())
X=torch.tensor(E)
res={}
for thr in [40,45,50,55,60,65]:
    ps=[]
    for f in sorted(glob.glob(f'tools/TemStaPro/models/mean_major_imbal-{thr}_s*.pt')):
        m=MLP_C2H2(1024,256,128); ck=torch.load(f,map_location='cpu',weights_only=False); sd=ck.get('state_dict',ck); sd={k.split('model.',1)[-1] if k.startswith('model.model.') else k:v for k,v in sd.items()}; sd={('model.'+k if not k.startswith('model.') else k):v for k,v in sd.items()}; m.load_state_dict(sd); m.eval()
        with torch.no_grad(): ps.append(m(X).squeeze(1).numpy())
    res[f'p{thr}']=np.mean(ps,0)
R=pd.concat([S.drop(columns='seq'), pd.DataFrame(res)],axis=1)
R['score']=R[[f'p{t}' for t in [40,45,50,55,60,65]]].sum(1)   # ~ number of thresholds exceeded
R.to_csv('results/temstapro_scores.csv', index=False)
print(R.groupby(['org','grp'])[['p40','p45','p50','score']].mean().round(3))
