import json, pandas as pd
d = json.load(open('data/full_dataset.json'))
rows = [(x['proteinId'], x.get('uniprotAccession'), x['runName'], x.get('meltingPoint'), x.get('quantNormMeltingPoint')) for x in d]
df = pd.DataFrame(rows, columns=['proteinId','acc','run','tm','tm_qn'])
df.to_csv('data/meltome_flat.csv', index=False)
print(df.groupby('run').agg(n=('tm','size'), n_tm=('tm','count'), med=('tm','median')).to_string())
