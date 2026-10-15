from pathlib import Path
import csv,json,math,collections
root=Path('reproducibility/darkmatter/round6_hypotheses/h4_structure/results')
read=lambda name:list(csv.DictReader((root/name).open(),delimiter='\t'))
m=read('v3_model_summaries.tsv');r=read('v3_region_evaluations.tsv')
assert len(m)==300 and len(r)==413
assert len({x['model_id'] for x in m})==300
assert sum(int(x['n_confident_regions']) for x in m)==len(r)
assert math.isclose(sum(float(x['pilot_ht_weight']) for x in m),5851,abs_tol=1e-4)
checks=0
for x in m:
 assert math.isclose(float(x['pilot_ht_weight']),int(x['pilot_pop_n'])/int(x['pilot_samp_n']),rel_tol=1e-6)
 for t in range(5):
  rr=[v for v in r if v['model_id']==x['model_id']]
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
for h,v in byh.items():assert len(v)==int(v[0]['pilot_samp_n'])
out={'passed':True,'model_region_assertions':checks,'models':len(m),'regions':len(r),'design_strata':{h:{'sample':len(v),'population':int(v[0]['pilot_pop_n'])} for h,v in byh.items()},'no_region_models':sum(int(x['n_confident_regions'])==0 for x in m),'total_design_weights':sum(float(x['pilot_ht_weight']) for x in m)}
Path('reviews/darkmatter-2026-10-13/round6/root_checks/h4_pilot_invariants.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out))
