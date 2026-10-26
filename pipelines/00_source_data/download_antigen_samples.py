"""Download public processed tumor GEX/VDJ for all Table-8 antigen patients."""
from pathlib import Path
import urllib.request,concurrent.futures,json,hashlib,gzip,re,time,tarfile,sys,shutil,io
import pandas as pd
sys.stdout.reconfigure(encoding='utf-8',errors='replace')
ROOT=Path('D:/CD8_Immunotherapy_20261023');OUT=ROOT/'analysis_extensions/GSE176021_antigen_validation';RAW=OUT/'raw';RAW.mkdir(exist_ok=True)
table=pd.read_excel(OUT/'supplements/41586_2021_3752_MOESM8_ESM.xlsx',sheet_name='Supplementary Table 8',header=2)
table=table[table.Patient.astype(str).str.fullmatch(r'(?:MD01|MD043|NY016)-\d{3}')].copy()
table['CDR3_AA_clean']=table['TCR clonotype (AA)'].str.replace('*','',regex=False).str.strip()
assert table['CDR3_AA_clean'].str.fullmatch('[A-Z]+(?:-[A-Z]+)*').all()
table['compound_dual_TRB_entry']=table['CDR3_AA_clean'].str.contains('-',regex=False)
table.to_csv(OUT/'functional_antigen_clonotypes_original_table8.csv',index=False)
patients=sorted(table.Patient.unique());print('ANTIGEN_PATIENTS',patients,flush=True)
items=[]
filelist=ROOT/'analysis_extensions/oncology_letters_extension/new_data/GSE176021_filelist.txt'
for line in filelist.read_text().splitlines():
 cols=line.split('\t')
 if len(cols)<5 or cols[0]!='File':continue
 # Three official VDJ filenames contain vdt/tdj typos; retain exact filenames.
 name=cols[1];m=re.match(r'(GSM\d+)_((?:MD01|MD043|NY016)-\d{3})_tumor_(\d+)(\.(?:vdj|vdt|tdj))?\.tar\.gz$',name)
 if not m or m.group(2) not in patients:continue
 gsm,patient,lane,vdj=m.groups()
 items.append(dict(name=name,bytes=int(cols[3]),gsm=gsm,patient=patient,lane=int(lane),kind='VDJ' if vdj else 'GEX',url='https://ftp.ncbi.nlm.nih.gov/geo/samples/'+gsm[:-3]+'nnn/'+gsm+'/suppl/'+name))
assert items and set(i['patient'] for i in items)==set(patients)
plan=dict(antigen_source='Nature 2021 Supplementary Table 8, original functional MANAFEST/ViraFEST and author TCR cloning labels; preserves assay origin and VDJdb column',patients=patients,tissue='tumor only; excludes normal, lymph node and metastatic tumor',expected_files=len(items),expected_compressed_bytes=sum(i['bytes'] for i in items),E_access_performed=False,items=items)
(OUT/'download_plan.json').write_text(json.dumps(plan,indent=2))
assert shutil.disk_usage(RAW).free>plan['expected_compressed_bytes']*2+4*1024**3
print('DOWNLOAD_PLAN',len(items),'files',plan['expected_compressed_bytes'],'bytes',flush=True)
def download(item):
 dest=RAW/item['name'];record=dest.with_suffix(dest.suffix+'.verified.json')
 if record.exists() and dest.stat().st_size==item['bytes']:
  r=json.loads(record.read_text());assert hashlib.sha256(dest.read_bytes()).hexdigest()==r['sha256'];return r
 last=None
 for attempt in range(3):
  part=dest.with_suffix(dest.suffix+'.part')
  try:
   if not part.exists() or part.stat().st_size!=item['bytes']:
    req=urllib.request.Request(item['url'],headers={'User-Agent':'Mozilla/5.0 research reproducibility'})
    with urllib.request.urlopen(req,timeout=90) as response,part.open('wb') as w:
     while block:=response.read(1024*1024):w.write(block)
   size=part.stat().st_size;sha=hashlib.sha256()
   with part.open('rb') as f:
    while block:=f.read(4*1024*1024):sha.update(block)
   assert size==item['bytes'],(item['name'],size,item['bytes'])
   # Iterate actual member streams to verify gzip CRC and tar readability without extraction.
   layers=1
   if item['kind']=='VDJ':
    payload=part.read_bytes();layers=0
    while payload.startswith(b'\x1f\x8b'):
     payload=gzip.decompress(payload);layers+=1
     assert layers<=3
    archive=tarfile.open(fileobj=io.BytesIO(payload),mode='r:')
   else:archive=tarfile.open(part,mode='r|gz')
   with archive:
    members=[]
    for member in archive:
     members.append(member.name)
     if member.isfile():
      f=archive.extractfile(member)
      while f.read(4*1024*1024):pass
   with gzip.open(part,'rb') as f:
    while f.read(4*1024*1024):pass
   part.replace(dest)
   r=dict(item,path=str(dest),sha256=sha.hexdigest(),gzip_CRC_verified=True,gzip_layers=layers,tar_member_names=members)
   record.write_text(json.dumps(r,indent=2));print('VERIFIED',item['name'],size,flush=True);return r
  except Exception as e:
   last=str(e);print('RETRY',attempt+1,item['name'],last,flush=True)
 return dict(item,error=last)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:records=list(pool.map(download,items))
(OUT/'download_manifest.json').write_text(json.dumps(records,indent=2))
failed=[r for r in records if 'error' in r]
print('COMPLETE',len(records)-len(failed),'FAILED',len(failed),flush=True)
if failed:raise RuntimeError(failed)
print('EXAMPLE_GEX',[r for r in records if r['kind']=='GEX'][0]['tar_member_names'],flush=True)
print('EXAMPLE_VDJ',[r for r in records if r['kind']=='VDJ'][0]['tar_member_names'],flush=True)
