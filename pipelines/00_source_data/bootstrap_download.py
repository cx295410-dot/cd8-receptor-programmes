"""Download public GEO processed data to a user-authorized D: project only."""
from pathlib import Path
import concurrent.futures, gzip, hashlib, json, os, re, time, urllib.request

ROOT = Path(r'D:\CD8_Immunotherapy_20261023')
assert ROOT.drive.upper() == 'D:' and ROOT.is_dir()
os.environ['TMP'] = os.environ['TEMP'] = str(ROOT/'cache')
MANIFEST = ROOT/'data'/'metadata'/'download_manifest.json'

def download(item):
    url, dest = item
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name+'.part')
    for attempt in range(1,4):
        try:
            if dest.exists() and dest.stat().st_size:
                print('EXISTS',dest.name,flush=True)
                break
            req=urllib.request.Request(url,headers={'User-Agent':'AcademicResearch/1.0'})
            with urllib.request.urlopen(req,timeout=60) as response, part.open('wb') as out:
                expected=int(response.headers.get('Content-Length',0))
                n=0; t=time.monotonic()
                while block:=response.read(1024*1024):
                    out.write(block); n+=len(block)
                    if time.monotonic()-t>25:
                        print('PROGRESS',dest.name,round(n/2**20,1),'MB',flush=True); t=time.monotonic()
                if expected and n!=expected: raise IOError(f'Incomplete download {n}/{expected}')
            part.replace(dest)
            print('DOWNLOADED',dest.name,dest.stat().st_size,flush=True)
            break
        except Exception as exc:
            print('RETRY',attempt,dest.name,type(exc).__name__,str(exc)[:180],flush=True)
            if attempt==3: raise
            time.sleep(2*attempt)
    h=hashlib.sha256()
    with dest.open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''): h.update(b)
    return {'url':url,'path':str(dest),'bytes':dest.stat().st_size,'sha256':h.hexdigest(),'retrieved_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())}

def main():
    tasks=[]
    for accession,prefix in [('GSE123813','GSE123nnn'),('GSE120575','GSE120nnn')]:
        base=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{accession}/suppl/'
        page=urllib.request.urlopen(base,timeout=60).read().decode()
        (ROOT/'data'/'metadata'/f'{accession}_ftp_index.html').write_text(page,encoding='utf-8')
        names=re.findall(r'href="([^"/]+\.gz)"',page)
        for name in names:
            tasks.append((base+name,ROOT/'data'/'raw'/accession/name))
    rows=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        for row in pool.map(download,tasks):
            rows.append(row); MANIFEST.write_text(json.dumps(rows,ensure_ascii=False,indent=2),encoding='utf-8')
    print('MANIFEST',len(rows),'files',sum(r['bytes'] for r in rows),flush=True)
    for row in rows:
        path=Path(row['path'])
        if 'metadata' in path.name or 'patient_ID' in path.name or '_tcr.' in path.name:
            print('\nHEADER',path.name,flush=True)
            with gzip.open(path,'rt',encoding='utf-8') as f:
                for i in range(6):
                    s=f.readline()
                    if not s: break
                    print(s[:2500].rstrip(),flush=True)

if __name__=='__main__': main()
