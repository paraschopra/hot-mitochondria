import pandas as pd, requests, io, time, concurrent.futures as cf
df = pd.read_csv('data/meltome_flat.csv')
accs = sorted(set(df.acc.dropna()))
print(len(accs))
fields = 'accession,organism_id,gene_primary,length,go_c,cc_subcellular_location,sequence,mass'
def fetch(batch):
    q = ' OR '.join(f'accession:{a}' for a in batch)
    for t in range(5):
        try:
            r = requests.post('https://rest.uniprot.org/uniprotkb/search', data={'query':q,'fields':fields,'format':'tsv','size':500}, timeout=120)
            if r.ok: return r.text
        except Exception as e: pass
        time.sleep(3)
    r = requests.get('https://rest.uniprot.org/uniprotkb/search', params={'query':q,'fields':fields,'format':'tsv','size':500}, timeout=120)
    return r.text
batches = [accs[i:i+150] for i in range(0,len(accs),150)]
out=[]
with cf.ThreadPoolExecutor(6) as ex:
    for i,t in enumerate(ex.map(fetch,batches)):
        try: out.append(pd.read_csv(io.StringIO(t), sep='\t'))
        except Exception as e: print('fail',i,t[:200])
u = pd.concat(out).drop_duplicates('Entry')
u.to_csv('data/uniprot_annot.tsv', sep='\t', index=False)
print(u.shape)
