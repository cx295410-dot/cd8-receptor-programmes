from pathlib import Path
import urllib.request,re,json,hashlib,sys
sys.stdout.reconfigure(encoding='utf-8',errors='replace')
import pandas as pd
OUT=Path('D:/CD8_Immunotherapy_20261023/analysis_extensions/GSE176021_antigen_validation');S=OUT/'supplements';S.mkdir(exist_ok=True)
html=(OUT/'nature_original.html').read_text(encoding='utf-8')
urls=sorted(set(v for v in re.findall(r'href=[\"\']([^\"\']+)',html) if 'media.springernature' in v and any(s in v for s in ['MOESM8','MOESM4','MOESM1'])))
records=[]
for url in urls:
 p=S/url.split('/')[-1]
 if not p.exists():
  with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0'}),timeout=60) as r:b=r.read()
  p.write_bytes(b)
 records.append(dict(url=url,path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
 if p.suffix=='.xlsx':
  xf=pd.ExcelFile(p);print(p.name,'SHEETS',xf.sheet_names,flush=True)
  for sheet in xf.sheet_names:
   d=pd.read_excel(p,sheet_name=sheet,header=None)
   print('SHEET',sheet,d.shape,'PREVIEW\n',d.iloc[:9,:12].to_string(index=False,header=False),flush=True)
   d.to_csv(S/(p.stem+'_'+re.sub('[^A-Za-z0-9_-]','_',sheet)+'.csv'),index=False,header=False)
(OUT/'supplements_manifest.json').write_text(json.dumps(records,indent=2))
