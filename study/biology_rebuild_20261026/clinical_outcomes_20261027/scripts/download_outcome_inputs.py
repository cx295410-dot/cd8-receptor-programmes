"""Download original public inputs without modifying the original study directory."""
from pathlib import Path
import os,requests,hashlib,json
ROOT=Path(os.environ.get('CD8_STUDY_ROOT',r'D:\Downloads\Virtual_Biotech_Study_20261025\biology_rebuild_20261026'))
RAW=ROOT/'clinical_outcomes_20261027/raw';RAW.mkdir(parents=True,exist_ok=True)
urls={
'GSE78220_series_matrix.txt.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78220/matrix/GSE78220_series_matrix.txt.gz',
'GSE78220_PatientFPKM.xlsx':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78220/suppl/GSE78220_PatientFPKM.xlsx',
'GSE91061_series_matrix.txt.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/matrix/GSE91061_series_matrix.txt.gz',
'GSE91061_fpkm.csv.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/suppl/GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz',
'GSE135222_series_matrix.txt.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/matrix/GSE135222_series_matrix.txt.gz',
'GSE135222_expression.tsv.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/suppl/GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz',
'Homo_sapiens.gene_info.gz':'https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz',
'41467_2019_12159_MOESM1_ESM.pdf':'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-019-12159-9/MediaObjects/41467_2019_12159_MOESM1_ESM.pdf'}
records=[]
for name,url in urls.items():
    p=RAW/name
    if not p.exists():
        r=requests.get(url,timeout=180);r.raise_for_status();p.write_bytes(r.content)
    records.append({'file':name,'url':url,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(RAW.parent/'download_sources.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
print('Verified public input files:',len(records))
