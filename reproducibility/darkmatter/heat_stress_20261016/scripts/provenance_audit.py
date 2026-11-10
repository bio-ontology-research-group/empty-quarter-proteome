from common import *
import pandas as pd,sys,subprocess,platform
from concurrent.futures import ThreadPoolExecutor
R=N/'results';O=R/'audit';O.mkdir(exist_ok=True);s=pd.read_csv(R/'hs3/pilot_samples.tsv',sep='\t');jobs=[]
for r in s[s.domain=='Bacteria'].itertuples():
 for mate in ['R1','R2']:jobs.append((getattr(r,mate),getattr(r,mate+'_sha256_prior')))
ref=W/'rebuild_20260729/picrust2_current_20260803/coverm_current/prepared/reps.mmi';jobs.append((str(ref),'6bdff27a444e8594f5211634ab2cbac05c372a0ddec5f6ae2676f60d9c534bfc'))
def audit(x):
 path,prior=x;p=Path(path);before=p.stat();digest=sha(p);after=p.stat();assert(before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns);return {'path':path,'bytes':after.st_size,'sha256':digest,'prior_sha256':prior,'matches_prior':digest==prior}
with ThreadPoolExecutor(max_workers=4)as pool:rows=list(pool.map(audit,jobs))
write(O/'raw_mapping_source_hashes.tsv',rows);assert all(r['matches_prior']for r in rows),rows
versions=subprocess.check_output([sys.executable,'-m','pip','freeze'],text=True);(O/'python_packages.txt').write_text(versions);(O/'execution_environment.json').write_text(json.dumps({'python':sys.version,'platform':platform.platform(),'job':os.environ['SLURM_JOB_ID'],'host':socket.gethostname()},indent=2))
print('DONE source hashes',len(rows))
