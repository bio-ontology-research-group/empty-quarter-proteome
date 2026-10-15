import csv,json,math
from pathlib import Path
R=Path(__file__).resolve().parents[1]; rows=list(csv.DictReader(open(R/'snapshots/h3_completed/matched_full_census_summary.tsv'),delimiter='\t'))
assert len(rows)==16
out=[]
for r in rows:
 n,a,b,c,d=[int(r[k]) for k in ['total_queries','both_hits_a','smag_only_b','gem_only_c','neither_d']]
 assert n==106042 and a+b+c+d==n
 for k,x in [('smag_hit_pct',100*(a+b)/n),('gem_hit_pct',100*(a+c)/n),('paired_delta_pct',100*(b-c)/n),('discordant_ratio',b/c)]:assert abs(float(r[k])-x)<=.00501,(r['run_label'],k)
 out.append({'run':r['run_label'],'cohort':r['cohort'],'smag_hits':a+b,'gem_hits':a+c,'total':n,'difference_pp':100*(b-c)/n})
primary={x['cohort']:x for x in out if x['run']=='harmonized_seed20261014'}
report={'status':'PASS_COMPACT_ARITHMETIC_ONLY','rows':out,'primary_dark_minus_control_difference_pp':primary['Dark']['difference_pp']-primary['Control']['difference_pp'],'raw_validation':'separate Slurm audit pending','pvalues':'Underflow zeros not interpreted; paired effect sizes describe fixed queried families and references'}
(R/'root_checks/h3_summary_audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
