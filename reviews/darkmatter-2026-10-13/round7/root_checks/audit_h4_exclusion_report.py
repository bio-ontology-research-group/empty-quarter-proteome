import csv,json,collections
from pathlib import Path
R=Path(__file__).resolve().parents[1]
d=json.load(open(R/'h4_tm_anomaly/exclusion_sensitivity_diff.json'))
report={}
for name,rows in d.items():
 report[name]={'entries':len(rows),'example':rows[:1]}
 if rows and 'query' in rows[0]:
  report[name]['unique_queries']=len({r['query'] for r in rows})
  report[name]['all_control_queries']=all(r['query'].startswith('KF') for r in rows)
  report[name]['tier_counts']=dict(collections.Counter(r['tier'] for r in rows))
 if rows and 'model_id' in rows[0]:
  report[name]['unique_models']=len({r['model_id'] for r in rows})
  report[name]['all_controls']=all(r['model_id'].startswith('KF') for r in rows)
  report[name]['tier_counts']=dict(collections.Counter(r['tier'] for r in rows))
for name, rows in d.items():
 if rows and 'crop_id' in rows[0]:
  report[name]['unique_crops']=len({r['crop_id'] for r in rows})
  report[name]['all_control_crops']=all(r['crop_id'].startswith('KF') for r in rows)
  report[name]['tier_counts']=dict(collections.Counter(r['tier'] for r in rows))
for name in ['wc_status_changes','crop_status_changes','model_crop_status_changes','model_transition_changes']:
 assert all(r['tier']>=2 for r in d[name]),name
assert all(r['stratum']=='Controls' and r['tier'] in ['tier_2','tier_3','tier_4'] for r in d['stratum_census_changes'])
report['primary_control_counts']=[r for r in csv.DictReader(open(R/'root_checks/h4_outputs_preclosure/v7_stratum_census.tsv'),delimiter='\t') if r['stratum']=='Controls' and r['tier']=='tier_1']
report['scope']='Excluding alnTM>1 rows is a diagnostic, not a correction of below-one alnTM values; primary qTM classifications and historical alnTM-defined strata must be distinguished.'
(R/'root_checks/h4_exclusion_fact_audit.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
