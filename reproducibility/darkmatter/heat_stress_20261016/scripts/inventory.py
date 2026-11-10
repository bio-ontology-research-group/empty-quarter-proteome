from common import *
import gzip,importlib.util
R=N/'results';R.mkdir(exist_ok=True)
roots=[W/'06_darkmatter/darkfam_basis_clean/results',W/'06_darkmatter/darkfam_ecology_clean/results',W/'06_darkmatter/lineage_darkness_clean/results',W/'06_darkmatter/isdark',W/'06_darkmatter/orthogroups',W/'temporal_environment_20261014/inputs',W/'rebuild_20260729/picrust2_current_20260803/coverm_current',T/'functional',T/'hmm/domtbl',W/'06_pangenome']
lines=[]
for root in roots:
 lines.append('\nDIRECTORY '+str(root))
 if not root.exists():continue
 for p in sorted(root.iterdir())[:35]:
  lines.append(str(p)+' '+str(p.stat().st_size))
  if p.is_dir():
   lines.extend('  '+str(q)+' '+str(q.stat().st_size) for q in list(p.iterdir())[:5])
  elif p.suffix in {'.tsv','.csv','.json','.txt'}:
   with open(p)as f:lines.append(' HEADER '+f.readline()[:1400].rstrip())
# One representative's exact identifiers and HMM output.
g=next(r for r in tab(T/'index/genome_index.tsv')if r['universe']=='soil_representative')
for p in [Path(g['proteome']),T/'hmm/domtbl'/(g['genome_id']+'.domtbl.gz')]:
 lines.append('\nEXAMPLE '+str(p))
 with gzip.open(p,'rt')as f:
  for i,l in enumerate(f):
   lines.append(l.strip())
   if i>=5:break
(R/'inventory.txt').write_text('\n'.join(lines))
print('Inventory complete',os.environ['SLURM_JOB_ID'])
