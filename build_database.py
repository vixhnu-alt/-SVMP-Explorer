"""Download real, accession-linked UniProt sequences and feature annotations.
Run on a computer with internet: python build_database.py
Never substitutes fabricated FASTA data for failed downloads.
"""
from pathlib import Path
import json, csv, time
import requests

BASE = Path(__file__).resolve().parent
DB = BASE / 'database'
refs = json.loads((DB/'references.json').read_text())
session = requests.Session()
session.headers['User-Agent'] = 'SVMP-Education-Project/1.0 (UniProt REST client)'

def get_entry(acc):
    r = session.get(f'https://rest.uniprot.org/uniprotkb/{acc}.json', timeout=25)
    r.raise_for_status()
    return r.json()

def extract_domain_features(entry):
    output = []
    for feature in entry.get('features', []):
        typ = feature.get('type','')
        desc = feature.get('description','') or ''
        start = feature.get('location',{}).get('start',{}).get('value')
        end = feature.get('location',{}).get('end',{}).get('value')
        if not isinstance(start,int) or not isinstance(end,int):
            continue
        if typ in {'Domain','Region','Motif','Chain','Propeptide','Signal peptide','Active site','Binding site'}:
            output.append({'type':typ,'description':desc,'start':start,'end':end,
                           'evidence':'UniProt entry annotation'})
    return output

def main():
    downloaded, failures = [], []
    for ref in refs:
        acc = ref['accession']
        try:
            e = get_entry(acc)
            seq = e.get('sequence',{}).get('value','').upper()
            if not seq or not set(seq) <= set('ACDEFGHIKLMNPQRSTVWYBXZJUO'):
                raise ValueError('missing or invalid protein sequence')
            family = ' '.join(e.get('comments',[])[i].get('texts',[{}])[0].get('value','')
                              for i in range(len(e.get('comments',[])))
                              if e['comments'][i].get('commentType')=='SIMILARITY' and e['comments'][i].get('texts'))
            name = e.get('proteinDescription',{}).get('recommendedName',{}).get('fullName',{}).get('value','')
            org = e.get('organism',{}).get('scientificName','')
            # Warn rather than silently change published class label.
            if 'P-I' in family or 'P-II' in family or 'P-III' in family:
                import re
                m = re.search(r'\bP-III\b|\bP-II\b|\bP-I\b', family)
                if m and m.group()!=ref['class']:
                    raise ValueError(f'UniProt class conflict: {m.group()} vs {ref["class"]}')
            status = e.get('entryType','')
            # Keep candidate status transparent; no synthetic C-domain boundaries.
            downloaded.append({**ref,'name':name,'species':org,'sequence':seq,
                               'length':len(seq),'review_status':status,'features':extract_domain_features(e),
                               'uniprot_url':f'https://www.uniprot.org/uniprotkb/{acc}',
                               'family_annotation':family})
            print('OK',acc,ref['class'],len(seq),org)
        except Exception as ex:
            failures.append({'accession':acc,'error':str(ex)})
            print('FAILED',acc,str(ex))
        time.sleep(.15)
    if downloaded:
        (DB/'reference_records.json').write_text(json.dumps(downloaded,indent=2))
        with (DB/'reference_sequences.fasta').open('w') as out:
            for r in downloaded:
                out.write(f'>{r["accession"]}|{r["class"]}|{r["species"]}\n')
                for i in range(0,len(r['sequence']),70):out.write(r['sequence'][i:i+70]+'\n')
    (DB/'download_status.json').write_text(json.dumps({'downloaded':len(downloaded),'failed':failures},indent=2))
    print(f'\nDownloaded {len(downloaded)}/{len(refs)} real reference sequences')
    if len(downloaded)<9:print('WARNING: reference coverage is limited. App will flag weak predictions.')

if __name__=='__main__':main()
