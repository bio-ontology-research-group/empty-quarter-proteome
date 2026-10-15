from pathlib import Path
import subprocess,json,datetime
base=Path('/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean')
parts=list((base/'h3_references/downloads/smag_parts').glob('mag.tar.gz.??'))
raw=subprocess.check_output(['sacct','-n','-X','-P','-j','53738540,53745194,53738740,53741162,53743377,53744172,53742584,53738996,53738997','--format=JobID,State,Elapsed'],text=True)
jobs={}
for line in raw.splitlines():
 fields=line.strip().split('|')
 if len(fields)>=3:jobs[fields[0]]={'state':fields[1],'elapsed':fields[2]}
h4=base/'h4_structure/results'
searches={n:(h4/n).stat().st_size if (h4/n).exists() else None for n in ['v7_wc_pdb_hits.raw.tsv','v7_crop_pdb_hits.raw.tsv','v7_recovery_summary.json']}
print(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'smag_complete_parts':len(parts),'smag_total_parts':373,'smag_complete_GB':round(sum(p.stat().st_size for p in parts)/1e9,2),'jobs':jobs,'h4_result_bytes':searches}))
