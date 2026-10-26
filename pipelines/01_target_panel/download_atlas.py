"""Download public normal-expression references to the user-authorized D-drive project."""
import hashlib, json, time, zipfile, re, html
from pathlib import Path
from urllib.parse import urljoin
import requests

ROOT=Path(r'D:\Downloads\Virtual_Biotech_Study_20261025')
PAGES=['https://www.proteinatlas.org/humanproteome/tissue/data',
       'https://www.proteinatlas.org/humanproteome/single+cell/single+cell+type/data']
NAMES=['rna_tissue_consensus.tsv.zip','rna_single_cell_type.tsv.zip','rna_single_cell_type_cell_types.tsv.zip']
manifest=[]
for page in PAGES:
    response=requests.get(page,timeout=60); response.raise_for_status()
    (ROOT/'literature'/('hpa_'+('tissue' if '/tissue/' in page else 'single_cell')+'_documentation.html')).write_text(response.text,encoding='utf-8')
    for href, label in re.findall(r'<a\b[^>]*href=[\"\x27]([^\"\x27]+)[\"\x27][^>]*>(.*?)</a>',response.text,flags=re.S|re.I):
        name=html.unescape(re.sub('<[^>]+>','',label)).strip()
        if name not in NAMES or any(x['file']==name for x in manifest): continue
        url=urljoin(page,html.unescape(href)); dest=ROOT/'data'/name
        start=time.monotonic()
        if not dest.exists():
            print('DOWNLOAD',name,url,flush=True)
            with requests.get(url,timeout=(30,120),stream=True) as r:
                r.raise_for_status()
                with open(dest.with_suffix('.part'),'wb') as out:
                    for chunk in r.iter_content(2**20): out.write(chunk)
            dest.with_suffix('.part').replace(dest)
        with zipfile.ZipFile(dest) as z:
            if z.testzip(): raise ValueError('Corrupt zip '+name)
            print('ARCHIVE',name,[(i.filename,i.file_size) for i in z.infolist()],flush=True)
        manifest.append({'file':name,'url':url,'documentation':page,'retrieved':'2026-10-25',
                         'bytes':dest.stat().st_size,'sha256':hashlib.sha256(dest.read_bytes()).hexdigest()})
        print('DONE',name,round(time.monotonic()-start,2),flush=True)
(ROOT/'data'/'hpa_download_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
if set(NAMES)!={m['file'] for m in manifest}: raise ValueError('Missing atlas datasets')
