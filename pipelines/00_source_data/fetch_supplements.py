"""Retrieve primary-publication supplementary assets, validating file type."""
from pathlib import Path
import concurrent.futures,hashlib,json,re,urllib.request,urllib.parse
ROOT=Path(r'D:\CD8_Immunotherapy_20261023'); DEST=ROOT/'data'/'raw'/'publication_supplements'; DEST.mkdir(exist_ok=True)
def fetch(url):
    try:
        r=urllib.request.urlopen(url,timeout=35); b=r.read()
        return url,b,r.geturl()
    except Exception as e: return url,None,type(e).__name__+': '+str(e)
pages=[
 'https://www.nature.com/articles/s41591-019-0522-3',
 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6689255/fullTextXML',
 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6641984/fullTextXML',
 'https://www.omar-lab.com/sctrial/_modules/sctrial/datasets.html']
results=list(concurrent.futures.ThreadPoolExecutor(4).map(fetch,pages)); links=[]; log=[]
for i,(url,b,actual) in enumerate(results):
    log.append({'url':url,'result':actual,'bytes':len(b) if b else 0})
    if b:
        (ROOT/'data'/'metadata'/f'supplement_source_{i}.txt').write_bytes(b)
        h=b.decode('utf-8',errors='replace')
        for s in re.findall(r'https?[^\s"<>]+?(?:xlsx|xls)',h): links.append(s.replace('&amp;','&'))
        for s in re.findall(r'(?:href|xlink:href)=["\']([^"\']+\.(?:xlsx|xls))["\']',h):
            links.append(urllib.parse.urljoin(actual,s))
links.extend(['https://pmc.ncbi.nlm.nih.gov/articles/instance/6689255/bin/NIHMS1531727-supplement-2.xlsx'])
links=sorted(set(s for s in links if any(t in s for t in ['522','6689255','1531727','2018.10.038','S0092867418313941','S009286741831394'])) )
for url,b,actual in concurrent.futures.ThreadPoolExecutor(3).map(fetch,links):
    ok=bool(b and b[:2]==b'PK')
    name=urllib.parse.unquote(url.split('/')[-1])
    if ok: (DEST/name).write_bytes(b)
    log.append({'url':url,'result':actual,'bytes':len(b) if b else 0,'valid_xlsx':ok})
    print('ASSET',ok,name,len(b) if b else 0,flush=True)
(ROOT/'data'/'metadata'/'supplement_download_log.json').write_text(json.dumps(log,indent=2),encoding='utf-8')
print('DISCOVERED_URLS',links,flush=True)
