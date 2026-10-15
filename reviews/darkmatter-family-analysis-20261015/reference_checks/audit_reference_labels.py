"""Independent source-based check of F12 reference labels; no sequence searches."""
import csv,json,sys
from pathlib import Path
root=Path(sys.argv[1]) if len(sys.argv)>1 else Path('reproducibility/darkmatter/family_discovery_20261015/f12')
meta=root/'references_v2/verified_references_metadata.tsv'
rows=list(csv.DictReader(open(meta),delimiter='\t'))
issues=[];seen={};counts={}
for r in rows:
 acc=r['accession'];cat=r.get('category','');sub=r.get('subcategory','')
 if acc in seen: issues.append([acc,'duplicate_accession',seen[acc],cat])
 seen[acc]=cat
 raw=r.get('raw_record_path',r.get('json_path',''))
 p=root/raw
 if not p.exists(): issues.append([acc,'missing_raw_record',str(p)]);continue
 if p.suffix!='.json':continue
 j=json.loads(p.read_text());pd=j.get('proteinDescription',{});rec=pd.get('recommendedName',{})
 name=rec.get('fullName',{}).get('value','')
 ecs={e['value'] for e in rec.get('ecNumbers',[])}
 pf={x['id'].split('.')[0] for x in j.get('uniProtKBCrossReferences',[]) if x.get('database')=='Pfam'}
 target=None
 if 'MpgP' in cat:target='3.1.3.70'
 if 'GgpP' in cat:target='3.1.3.69'
 if 'TPP' in cat:target='3.1.3.12'
 if target and target not in ecs:issues.append([acc,'target_EC_absent',cat,name,sorted(ecs)])
 if cat.startswith('CPF') and not (pf&{'PF00875','PF03441','PF04244'}):issues.append([acc,'no_CPF_Pfam_in_source_requires_HMM_control',name,sorted(pf)])
 counts[cat]=counts.get(cat,0)+1
print(json.dumps({'input':str(meta),'rows':len(rows),'unique_accessions':len(seen),'categories':counts,'issues':issues},indent=2))
