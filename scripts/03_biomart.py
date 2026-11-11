import requests, pandas as pd, io
URL='https://www.ensembl.org/biomart/martservice'
def q(dataset, attrs):
    x = f'''<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" datasetConfigVersion="0.6"><Dataset name="{dataset}" interface="default">''' + ''.join(f'<Attribute name="{a}"/>' for a in attrs) + '</Dataset></Query>'
    r = requests.post(URL, data={'query':x}, timeout=600); r.raise_for_status()
    return pd.read_csv(io.StringIO(r.text), sep='\t')
# human->orthologs
for sp in ['mmusculus','dmelanogaster','celegans','scerevisiae','drerio']:
    df = q('hsapiens_gene_ensembl', ['ensembl_gene_id',f'{sp}_homolog_ensembl_gene',f'{sp}_homolog_orthology_type',f'{sp}_homolog_perc_id'])
    df.to_csv(f'data/orth_human_{sp}.tsv', sep='\t', index=False); print(sp, df.shape)
for ds in ['hsapiens','mmusculus','dmelanogaster','celegans','scerevisiae','drerio']:
    df = q(f'{ds}_gene_ensembl', ['ensembl_gene_id','uniprotswissprot'])
    df.to_csv(f'data/uniprot_map_{ds}.tsv', sep='\t', index=False); print(ds, df.shape)
