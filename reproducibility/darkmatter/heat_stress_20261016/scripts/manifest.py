from common import *
rows=[]
for p in sorted(N.rglob('*')):
 if not p.is_file()or '__pycache__'in p.parts or '/tmp/'in str(p)or p.name=='artifact_manifest.tsv' or p.name==f"manifest-{os.environ['SLURM_JOB_ID']}.log":continue
 if p.suffix not in ['.py','.sbatch','.tsv','.json','.gz','.npz','.bam','.bai','.log','.txt','.md','.pdf','.png','.bed']:continue
 rows.append({'path_relative_to_remote_root':str(p.relative_to(N)),'bytes':p.stat().st_size,'sha256':sha(p)})
write(N/'results/audit/artifact_manifest.tsv',rows);print('Hashed',len(rows),'artifacts')
