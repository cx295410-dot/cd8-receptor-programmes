from pathlib import Path
import urllib.request,json,hashlib,io,zipfile
import pandas as pd
URL='https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs43018-021-00292-8/MediaObjects/43018_2021_292_MOESM3_ESM.xlsx'
OUT=Path('E:/CD8_Data/GSE179994/clinical')
META=Path('D:/CD8_Immunotherapy_20261023/analysis_extensions/independent_GSE179994')
with urllib.request.urlopen(URL,timeout=45) as response: data=response.read(20_000_001)
if len(data)>20_000_000: raise ValueError('Unexpected supplement size')
assert zipfile.ZipFile(io.BytesIO(data)).testzip() is None
OUT.mkdir(exist_ok=True)
target=OUT/'43018_2021_292_MOESM3_ESM.xlsx'; target.write_bytes(data)
manifest=dict(url=URL,path=str(target),bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),xlsx_zip_CRC_passed=True)
(META/'clinical_supplement_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
book=pd.ExcelFile(io.BytesIO(data))
print('SHEETS',book.sheet_names,flush=True)
audit=[]
for sheet in book.sheet_names:
    tab=pd.read_excel(book,sheet_name=sheet,header=None,nrows=9)
    audit.append(dict(sheet=sheet,preview=tab.fillna('').to_numpy().tolist()))
    print('SHEET',sheet,tab.to_string(index=False,header=False),flush=True)
(META/'clinical_supplement_sheet_preview.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2,default=str),encoding='utf-8')
print('CLINICAL_SOURCE_DOWNLOADED',json.dumps(manifest),flush=True)
