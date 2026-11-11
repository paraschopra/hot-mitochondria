import pandas as pd, requests, time, io, re
t = pd.read_csv('data/meltome_S4_long.csv')
acc = t['Protein ID'].str.split('_').str[0].str.replace(r'-\d+$','',regex=True)
accs = sorted(set(acc.dropna()))
accs = [a for a in accs if re.match(r'^[A-Z0-9]{6,10}$', a)]
print('n accs', len(accs))
fields='accession,organism_id,gene_primary,length,go_c,cc_subcellular_location,sequence'
out=[]
for i in range(0, len(accs), 20000):
    chunk = accs[i:i+20000]
    r = requests.post('https://rest.uniprot.org/idmapping/run', data={'from':'UniProtKB_AC-ID','to':'UniProtKB','ids':','.join(chunk)}); r.raise_for_status()
    job = r.json()['jobId']; print('job', job)
    while True:
        rr = requests.get(f'https://rest.uniprot.org/idmapping/status/{job}', allow_redirects=False)
        if rr.status_code in (303,302): break
        s = rr.json(); print(s.get('jobStatus'))
        if s.get('jobStatus')=='FINISHED' or 'results' in s: break
        time.sleep(10)
    time.sleep(5)
    r = requests.get(f'https://rest.uniprot.org/idmapping/uniprotkb/results/stream/{job}', params={'fields':fields,'format':'tsv'}, timeout=1200)
    r.raise_for_status()
    d = pd.read_csv(io.StringIO(r.text), sep='\t'); print(d.shape); out.append(d)
u = pd.concat(out).drop_duplicates('Entry'); u.to_csv('data/uniprot_annot.tsv', sep='\t', index=False); print(u.shape)
print(u['Organism (ID)'].value_counts().head(20))
