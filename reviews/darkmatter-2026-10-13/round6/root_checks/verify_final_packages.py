"""Verify the four final approved snapshots against live artifacts and baseline."""
from pathlib import Path
import hashlib,json,datetime
r=Path(__file__).resolve().parents[1]
repo=r.parents[2]
state=json.loads((r/'STATE.json').read_text())
areas={'H1':'h1_lineage','H2':'h2_context','H3':'h3_references','H4':'h4_structure'}
h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
out={'verified_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'hypotheses':{},'baseline':{}}
for name,area in areas.items():
 a=state['accepted'][name];assert a['status']=='VIABLE'
 snap=r/a['snapshot'];review=r/a['review'];vp=review/'verdicts_final.json'
 if not vp.exists():vp=review/'verdicts.json'
 v=json.loads(vp.read_text());assert v.get('status',v.get('verdict'))=='VIABLE'
 assert not v.get('critical_issues',v.get('new_critical_issues',[]))
 assert not v.get('required_fixes',[])
 for doc,info in v.get('corrected_document_hashes',{}).items():
  expected=info.get('sha256',info.get('post_edit_sha256')) if isinstance(info,dict) else info
  assert h(repo/doc)==expected,doc
 live=repo/'reproducibility/darkmatter/round6_hypotheses'/area
 files=json.loads((snap/'SNAPSHOT_MANIFEST.json').read_text())['files'];checked=0
 for row in files:
  p=snap/row['path'];assert p.stat().st_size==row['bytes'] and h(p)==row['sha256'],str(p)
  q=live/row['path']
  if row['path']!='LARGE_ARTIFACT_LOCATORS.json':assert q.is_file() and h(q)==row['sha256'],str(q)
  checked+=1
 out['hypotheses'][name]={'status':'VIABLE','snapshot':str(snap.relative_to(repo)),'snapshot_manifest_sha256':h(snap/'SNAPSHOT_MANIFEST.json'),'review':str(vp.relative_to(repo)),'review_sha256':h(vp),'verified_files':checked}
for p,sha in json.loads((r/'baseline_hashes.json').read_text()).items():
 out['baseline'][p]=h(repo/p)==sha
 assert out['baseline'][p]
out['passed']=True
(r/'root_checks/final_package_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
