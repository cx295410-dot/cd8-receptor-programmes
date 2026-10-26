"""Copy the recorded scripts to configurable input/output roots without running them."""
from pathlib import Path
import argparse,json,ast
REPO=Path(__file__).resolve().parents[1]
def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--study-root',type=Path,required=True)
 p.add_argument('--legacy-root',type=Path,required=True)
 p.add_argument('--dry-run',action='store_true')
 args=p.parse_args();study=args.study_root.resolve();legacy=args.legacy_root.resolve()
 if study==legacy:raise SystemExit('Study and legacy roots must differ.')
 bio=study/'biology_rebuild_20261026'
 targets={'00_source_data':legacy/'scripts','01_target_panel':study/'scripts','02_transcriptome':bio/'scripts','03_expansion':bio/'reference_expansion_20261026/scripts','04_clinical':bio/'clinical_validation_20261026/scripts','05_multipanel':bio/'multipanel_rebuild_20261026/scripts'}
 replacements=[('D:/Downloads/Virtual_Biotech_Study_20261025',study.as_posix()),('D:/CD8_Immunotherapy_20261023',legacy.as_posix())]
 jobs=[]
 for stage,dest in targets.items():
  for source in sorted((REPO/'pipelines'/stage).glob('*.py')):
   text=source.read_text(encoding='utf-8-sig')
   for old,new in replacements:
    text=text.replace(old.replace('/','\\'),new).replace(old,new)
   ast.parse(text,filename=source.name)
   jobs.append((dest/source.name,text))
 jobs += [(study/'work/panel.json',(REPO/'config/target_panel.json').read_text()),(bio/'work/programme_definitions.json',(REPO/'config/programme_definitions.json').read_text()),(bio/'Software_Versions.json',(REPO/'config/software_versions_recorded.json').read_text()),(bio/'clinical_validation_20261026/analysis_protocol_frozen.json',(REPO/'config/clinical_analysis_protocol_frozen.json').read_text())]
 for dest,text in jobs:
  if dest.exists() and dest.read_text(encoding='utf-8-sig')!=text:raise SystemExit('Refusing to overwrite a different existing file: '+str(dest))
 if not args.dry_run:
  for root in [study,bio,bio/'reference_expansion_20261026',bio/'clinical_validation_20261026',bio/'multipanel_rebuild_20261026',legacy]:
   for name in ['scripts','data','raw','results','work','cache','references','figures','figure_source_data']:(root/name).mkdir(parents=True,exist_ok=True)
  for dest,text in jobs:dest.parent.mkdir(parents=True,exist_ok=True);dest.write_text(text,encoding='utf-8')
 print(json.dumps({'dry_run':args.dry_run,'planned_files':len(jobs),'study_root':str(study),'legacy_root':str(legacy)}))
if __name__=='__main__':main()
