from pathlib import Path
import csv,json,math,collections,sys
root=Path(sys.argv[1] if len(sys.argv)>1 else 'reproducibility/darkmatter/round7_analyses/h4_structure/results')
read=lambda name:list(csv.DictReader((root/name).open(),delimiter='\t'))
m=read('v7_model_summaries.tsv');r=read('v7_region_evaluations.tsv')
assert len(m)==5851 and len(r)==7773
assert len({x['model_id'] for x in m})==5851
assert sum(int(x['n_confident_regions']) for x in m)==len(r)
checks=0
by_model=collections.defaultdict(list)
for v in r: by_model[v["model_id"]].append(v)
for x in m:
 for t in range(5):
  rr=by_model[x['model_id']]
  assert int(x[f'crop_t{t}_matched_regions'])==sum(int(v[f't{t}_matched']) for v in rr)
  assert (x[f'crop_t{t}_status']=='MATCHED')==any(int(v[f't{t}_matched']) for v in rr)
  assert (x[f'crop_t{t}_status']=='NO_CONFIDENT_REGION')==(not rr)
  checks+=3
for x in r:
 for t in (0,1):
  assert int(x[f't{t}_aligned_pairs']) <= int(x[f't{t}_qspan']) <= int(x['reg_len'])
  assert math.isclose(float(x[f't{t}_aligned_pair_qcov']),int(x[f't{t}_aligned_pairs'])/int(x['reg_len']),abs_tol=0.000051)
  checks+=2
byh=collections.defaultdict(list)
for x in m:byh[x['design_stratum_h']].append(x)
out={'passed':True,'model_region_assertions':checks,'models':len(m),'regions':len(r),'census_design_strata':{h:len(v) for h,v in byh.items()},'no_region_models':sum(int(x['n_confident_regions'])==0 for x in m)}
Path('reviews/darkmatter-2026-10-13/round7/root_checks/h4_census_invariants.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out))
