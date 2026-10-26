"""Validate the code-only release without importing data-dependent pipelines."""
from pathlib import Path
import ast,hashlib,json,re,io,tokenize
ROOT=Path(__file__).resolve().parents[1]
records=json.loads((ROOT/'provenance/source_code_manifest.json').read_text())
for r in records:
 p=ROOT/r['repository_path'];assert hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'],p
scripts=list((ROOT/'pipelines').rglob('*.py'))+list((ROOT/'tools').glob('*.py'))
for p in scripts:ast.parse(p.read_text(encoding='utf-8-sig'),filename=str(p))
panel=json.loads((ROOT/'config/target_panel.json').read_text())['targets']
programmes=json.loads((ROOT/'config/programme_definitions.json').read_text())['programmes']
assert len(panel)==24 and len(set(panel))==24
markers=[g for genes in programmes.values() for g in genes]
assert len(programmes)==6 and len(markers)==38 and not set(markers)&set(panel)
for p in ROOT.rglob('*'):
 if not p.is_file() or any(part in {'.git','__pycache__','.venv','env','study_workspace','legacy_workspace'} for part in p.relative_to(ROOT).parts):continue
 assert p.suffix.lower() not in {'.docx','.h5ad','.h5','.rds','.png','.pdf','.tif','.zip'},p
 assert p.name not in {'.env','.git-credentials','id_rsa','id_ed25519'},p
 text=p.read_text(encoding='utf-8-sig')
 for pattern in [r'gh[pousr]_[A-Za-z0-9]{25,}',r'github_pat_[A-Za-z0-9_]{30,}',r'sk-(?:proj-)?[A-Za-z0-9_-]{35,}',r'-----BEGIN (?:RSA |OPENSSH )?PRIVATE KEY-----']:
  assert not re.search(pattern,text),str(p)+' contains a known secret pattern'
 email_text='\n'.join(t.string for t in tokenize.generate_tokens(io.StringIO(text).readline) if t.type in [tokenize.STRING,tokenize.COMMENT]) if p.suffix=='.py' else text
 assert not re.search(r'[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}',email_text),str(p)+' contains an email address'
print(json.dumps({'passed':True,'source_scripts':len(records),'python_files_parsed':len(scripts),'targets':len(panel),'programme_markers':len(markers)}))
