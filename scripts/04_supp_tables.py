import pandas as pd
s3 = pd.read_excel('data/meltome_supp/S3.xlsx', sheet_name='Meltome data set')
s3.to_csv('data/meltome_datasets.csv', index=False)
pd.set_option('display.width',250); pd.set_option('display.max_colwidth',40)
print(s3[['Dataset ID','Organism','Strain/tissue','OGT [°C]','Cells/lysate','Replicate/experiment','Lysis conditions','Number of fitted melting curves']].to_string())
x = pd.ExcelFile('data/meltome_supp/S4.xlsx')
out=[]
for s in x.sheet_names[1:]:
    d = x.parse(s); d['dataset']=s; out.append(d)
t = pd.concat(out); t.to_csv('data/meltome_S4_long.csv', index=False); print(t.shape); print(t.head())
