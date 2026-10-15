import csv,json,statistics,collections,math
from pathlib import Path
base=Path('reproducibility/darkmatter/round7_analyses/h2_context/results')
rows=list(csv.DictReader((base/'pilot_family_conservation.tsv').open(),delimiter='\t'))
summary=json.loads((base/'pilot_conservation_summary.json').read_text())
checks={}
for label,gate,nsp,opp,outcome in [('tier1_primary','eligible_primary','n_species_complete_window','ko_annotation_opportunity','excess_c_fam_primary'),('tier2_stop3_sensitivity','eligible_stop3_sens','n_species_stop3_window','ko_annotation_opportunity_stop3','excess_c_fam_stop3')]:
 r=[x for x in rows if int(x[gate])]
 median=statistics.median(float(x[opp]) for x in r)
 strata=collections.defaultdict(list)
 for x in r:
  key=('sp2' if int(x[nsp])==2 else 'sp3plus',x['phylum_group'],float(x[opp])>median)
  strata[key].append(x)
 numerator=denominator=0
 for key,v in strata.items():
  k=[float(x[outcome]) for x in v if x['arm']=='K'];q=[float(x[outcome]) for x in v if x['arm']=='QnotK']
  if not k or not q:continue
  weight=len(k)*len(q)/(len(k)+len(q))
  numerator+=weight*(statistics.mean(k)-statistics.mean(q));denominator+=weight
 contrast=numerator/denominator
 reported=summary['part_b_between_arm_covariate_stratified_permutation'][label]['mh_weighted_excess_c_diff']
 assert abs(contrast-reported)<1e-6
 checks[label]={'eligible_families':len(r),'strata':len(strata),'independent_mh_contrast':contrast,'reported':reported,'passed':True}
Path('reviews/darkmatter-2026-10-13/round7/root_checks/h2_contrasts.json').write_text(json.dumps(checks,indent=2)+'\n')
print(json.dumps(checks))
