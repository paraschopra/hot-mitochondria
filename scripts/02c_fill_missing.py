import pandas as pd, requests, time, io, re
t = pd.read_csv('data/meltome_S4_long.csv')
acc = t['Protein ID'].str.split('_').str[0].str.replace(r'-\d+$','',regex=True)
u = pd.read_csv('data/uniprot_annot.tsv', sep='\t')
have=set(u.Entry)
accs = sorted(set(a for a in acc.dropna() if re.match(r'^[A-Z0-9]{6,10}$', a)) - have)
print('missing', len(accs))
fields='accession,organism_id,gene_primary,length,go_c,cc_subcellular_location,sequence'
out=[u]
for i in range(0, len(accs), 4000):
    chunk = accs[i:i+4000]
    r = requests.post('https://rest.uniprot.org/idmapping/run', data={'from':'UniProtKB_AC-ID','to':'UniProtKB','ids':','.join(chunk)}); job=r.json()['jobId']
    while True:
        rr = requests.get(f'https://rest.uniprot.org/idmapping/status/{job}', allow_redirects=False)
        if rr.status_code in (303,302): break
        s = rr.json()
        if s.get('jobStatus')=='FINISHED' or 'results' in s or 'failedIds' in s: break
        time.sleep(5)
    time.sleep(3)
    r = requests.get(f'https://rest.uniprot.org/idmapping/uniprotkb/results/stream/{job}', params={'fields':fields,'format':'tsv'}, timeout=1200)
    d = pd.read_csv(io.StringIO(r.text), sep='\t') if r.ok and r.text.strip() else pd.DataFrame()
    print(i, len(chunk), d.shape, flush=True); out.append(d)
u = pd.concat(out).drop_duplicates('Entry'); u.to_csv('data/uniprot_annot.tsv', sep='\t', index=False); print(u.shape)
