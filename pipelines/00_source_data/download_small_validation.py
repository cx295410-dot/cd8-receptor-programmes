from pathlib import Path
import urllib.request,json,hashlib,gzip,time,os
os.environ['PYTHONIOENCODING']='utf-8'
import sys
sys.stdout.reconfigure(encoding='utf-8',errors='replace')
ROOT=Path('D:/CD8_Immunotherapy_20261023/analysis_extensions/oncology_letters_extension')
OUT=ROOT/'new_data'; OUT.mkdir(exist_ok=True,parents=True)
files=[('GSE176021','GSE176021_CD8_annotations.rds.gz'),('GSE176021','filelist.txt'),('GSE190265','GSE190265_TPM_France3.csv.gz'),('GSE190265','GSE190265_samples_info_France3.csv.gz')]
records=[]
for acc,name in files:
 base='https://ftp.ncbi.nlm.nih.gov/geo/series/'+acc[:-3]+'nnn/'+acc+'/suppl/'
 dest=OUT/(acc+'_'+name if name=='filelist.txt' else name)
 if not dest.exists():
  part=dest.with_suffix(dest.suffix+'.part'); size=0; tick=time.monotonic()
  with urllib.request.urlopen(urllib.request.Request(base+name,headers={'User-Agent':'Mozilla/5.0'}),timeout=60) as r,part.open('wb') as w:
   while b:=r.read(1024*1024):
    w.write(b); size+=len(b)
    if time.monotonic()-tick>20: print(name,size,'bytes',flush=True);tick=time.monotonic()
   expected=r.headers.get('Content-Length')
  if expected: assert size==int(expected)
  part.replace(dest)
 if name.endswith('.gz'):
  with gzip.open(dest,'rb') as f:
   while f.read(8*1024*1024): pass
 rec=dict(accession=acc,name=name,url=base+name,path=str(dest),bytes=dest.stat().st_size,sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),gzip_crc_verified=name.endswith('.gz'))
 records.append(rec); print('VERIFIED',name,rec['bytes'],flush=True)
(ROOT/'new_data_manifest.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
sys.path.insert(0,'D:/CD8_Immunotherapy_20261023/analysis_extensions/packages')
import rdata,pandas as pd
d=rdata.read_rds(OUT/'GSE176021_CD8_annotations.rds.gz')
print('ANNOTATION_TYPE',type(d),flush=True)
if isinstance(d,pd.DataFrame):
 print('SHAPE',d.shape,'COLS',list(d.columns),flush=True); print(d.head().to_string(),flush=True)
 d.to_csv(OUT/'GSE176021_CD8_annotations.csv.gz',index=True)
 (ROOT/'GSE176021_annotation_schema.json').write_text(json.dumps(dict(shape=list(d.shape),columns=list(d.columns),index_name=d.index.name),indent=2),encoding='utf-8')
else: print(str(d)[:2000],flush=True)
for name in ['GSE190265_TPM_France3.csv.gz','GSE190265_samples_info_France3.csv.gz']:
 d=pd.read_csv(OUT/name,sep=';',index_col=0); print(name,d.shape,'first columns',list(d.columns)[:8],'first IDs',list(d.index[:4]),flush=True)
print((OUT/'GSE176021_filelist.txt').read_text()[:5000],flush=True)
