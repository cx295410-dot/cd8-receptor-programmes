from pathlib import Path
import sys,json,re,hashlib
sys.stdout.reconfigure(encoding='utf-8')
import pandas as pd
B=Path('D:/CD8_Immunotherapy_20261023/analysis_extensions')
O=Path('D:/Downloads/Virtual_Biotech_Study_20261025/biology_rebuild_20261026/clinical_validation_20261026')
X=B/'GSE176021_antigen_validation/supplements/41586_2021_3752_MOESM8_ESM.xlsx'
for sheet,header in [('Supplementary Table 1',1),('Supplementary Table 2',3),('Supplementary Table 3',1)]:
 d=pd.read_excel(X,sheet_name=sheet,header=header)
 print('\nSHEET',sheet,'columns',d.columns.tolist());print(d.to_string(index=False))
 d.to_csv(O/'data'/f'{sheet.replace(" ","_")}_raw.csv',index=False)
a=pd.read_csv(B/'oncology_letters_extension/new_data/GSE176021_CD8_annotations.csv.gz',index_col=0)
print('ANNOTATION_COLUMNS',a.columns.tolist());print(a.groupby('imid').size().to_string());print('CELL_EXAMPLE',a.barcode.head(6).tolist())
rows=[]
for line in (B/'oncology_letters_extension/new_data/GSE176021_filelist.txt').read_text().splitlines():
 p=line.split('\t')
 if len(p)<5 or p[0]!='File':continue
 m=re.match(r'(GSM\d+)_((?:MD01|MD043|NY016)-\d{3})_(.+)_(\d+)\.tar\.gz$',p[1])
 if m:rows.append(dict(gsm=m[1],patient=m[2],tissue=m[3],lane=int(m[4]),name=p[1],bytes=int(p[3])))
d=pd.DataFrame(rows);d.to_csv(O/'data/GEO_library_inventory.csv',index=False)
print('LIBRARIES',d.groupby(['patient','tissue']).agg(libraries=('gsm','size'),bytes=('bytes','sum')).to_string())
