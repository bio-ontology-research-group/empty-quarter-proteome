"""Run from the retained IBEX review package after the final sync."""
from pathlib import Path
import json,hashlib,datetime,os
r=Path(__file__).resolve().parents[1];work=r.parents[1]
areas={'H1':'h1_lineage','H2':'h2_context','H3':'h3_references','H4':'h4_structure'}
s=json.loads((r/'STATE.json').read_text());expected=json.loads((r/'root_checks/final_package_verification.json').read_text());out={'checked_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'slurm_job_id':os.environ.get('SLURM_JOB_ID'),'hypotheses':{}}
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for x in iter(lambda:f.read(1024*1024),b''):h.update(x)
 return h.hexdigest()
for name,area in areas.items():
 a=s['accepted'][name];assert a['status']=='VIABLE'
 snap=r/a['snapshot'];proof=expected['hypotheses'][name]
 assert sha(snap/'SNAPSHOT_MANIFEST.json')==proof['snapshot_manifest_sha256']
 vp=r/a['review']/'verdicts_final.json'
 if not vp.exists():vp=r/a['review']/'verdicts.json'
 assert sha(vp)==proof['review_sha256']
 m=json.loads((snap/'SNAPSHOT_MANIFEST.json').read_text());n=0
 for row in m['files']:
  p=snap/row['path'];assert p.is_file() and p.stat().st_size==row['bytes'] and sha(p)==row['sha256'],str(p)
  if row['path']!='LARGE_ARTIFACT_LOCATORS.json':
   q=work/area/row['path'];assert q.is_file() and sha(q)==row['sha256'],str(q)
  n+=1
 out['hypotheses'][name]={'snapshot_files_checked':n,'snapshot_manifest_sha256':sha(snap/'SNAPSHOT_MANIFEST.json'),'live_analysis_matches':True}
out['passed']=True
(r/'root_checks/ibex_final_package_verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
