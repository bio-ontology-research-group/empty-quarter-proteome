import csv,json,math
from pathlib import Path
root=Path('reviews/darkmatter-2026-10-13/round6/root_checks/h3_remote_tables')
read=lambda n:list(csv.DictReader((root/n).open(),delimiter='\t'))
a=read('v3_smag_tar_extraction_audit.tsv')
assert len(a)==213
assert all(x['status']=='COMPLETE_RETAINED' and x['declared_bytes']==x['observed_bytes'] and len(x['observed_sha256'])==64 for x in a[:212])
assert a[-1]['status']=='REJECTED_TRUNCATED_STREAM_END'
s=read('v3_matched_length_stratified_results.tsv');c=read('v3_matched_cohort_summary.tsv')
results=[]
for row in c:
 ss=[x for x in s if x['run_label']==row['run_label'] and x['cohort']==row['cohort']]
 assert len(ss)==5
 N=sum(int(x['N_h']) for x in ss);n=sum(int(x['n_h']) for x in ss)
 assert N==int(row['N_pop']) and n==int(row['n_panel'])
 delta=var=0
 for x in ss:
  nh=int(x['n_h']);Nh=int(x['N_h']);b=int(x['smag_only']);cc=int(x['gem_only']);d=(b-cc)/nh
  assert sum(int(x[k]) for k in ('both','smag_only','gem_only','neither'))==nh
  assert int(x['smag_hits'])==int(x['both'])+b
  assert int(x['gem_hits'])==int(x['both'])+cc
  delta+=Nh/N*d
  var+=(Nh/N)**2*(1-nh/Nh)*((b+cc)-nh*d*d)/(nh-1)/nh
 assert abs(delta*100-float(row['pop_wt_delta_pct']))<=.0051
 assert abs(math.sqrt(var)*100-float(row['pop_wt_se_delta_fpc_pct']))<=.0051
 results.append({'run':row['run_label'],'cohort':row['cohort'],'delta_pp':delta*100,'se_pp':math.sqrt(var)*100})
out={'passed':True,'complete_smag_genomes':212,'truncated_member_rejected':True,'weighted_contrasts_verified':results}
Path('reviews/darkmatter-2026-10-13/round6/root_checks/h3_table_audit.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out))
