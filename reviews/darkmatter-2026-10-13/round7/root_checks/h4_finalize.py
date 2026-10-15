from pathlib import Path
import sys,json,subprocess
base=Path('/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h4_structure')
sys.path.insert(0,str(base/'scripts_recovery_v3'))
from v7_crop_and_search_regions import generate_sha256_manifest,sha256_file
from v7_verify_closure import verify_closure
state=subprocess.check_output(['sacct','-j','53744678','-n','-X','-P','--format=JobID,State,ExitCode'],text=True)
assert '53744678|FAILED|1:0' in state
err=base/'logs/h4_v7_recovery_53744678.err'
assert "NameError: name 'SUB_STRATA_ORDER' is not defined" in err.read_text()
# Every numeric output is retained byte-for-byte; only verification receipt/manifest are sealed.
r=base/'results'
prior={p.name:sha256_file(str(p)) for p in r.glob('v7_*.tsv') if 'sha256' not in p.name and '.raw.' not in p.name}
subprocess.run([sys.executable,str(base/'scripts_recovery_v1/audit_h4_raw.py')],check=True)
generate_sha256_manifest(str(r),include_receipt=False)
verify_closure(str(r),mode='provisional')
generate_sha256_manifest(str(r),include_receipt=True)
verify_closure(str(r),mode='strict')
assert prior=={name:sha256_file(str(r/name)) for name in prior}
receipt={'status':'PASS','prior_summary_job':53744678,'prior_verifier_error_sha256':sha256_file(str(err)),'unchanged_numerical_table_sha256':prior,'final_manifest_sha256':sha256_file(str(r/'v7_results.sha256.tsv'))}
(base/'logs/h4_finalization_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print('FULL_H4_VERIFICATION_COMPLETE')
