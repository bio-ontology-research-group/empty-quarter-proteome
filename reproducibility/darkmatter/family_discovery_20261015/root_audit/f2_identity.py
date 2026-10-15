"""Resolve historical/current sequence identity directly from immutable source FASTAs."""
import gzip,hashlib,json,os,subprocess
from pathlib import Path
W=Path('/ibex/scratch/projects/c2014/rob/metagenomics');D=W/'06_darkmatter/family_discovery_20261015_clean';R=D/'root_audit/f2_identity';R.mkdir(exist_ok=True)
ids=['1PRr2__self_SemiBin_81::1PRr2_1391604_4','S52PRr1__single_SemiBin_53::S52PRr1_477345_9']
sources=[W/'06_darkmatter/orthogroups/orthogroups_rep_seq.fasta',W/'06_darkmatter/orthogroups/all_rep_proteins.faa',D/'f12/results/pilot_v3_results/f2_leads.faa']
sources += [W/f'06_darkmatter/drep_rep_proteomes/{pid.split("::")[0]}.faa.gz' for pid in ids]
records=[]
for p in sources:
 if not p.exists():
  records.append({'source':str(p),'status':'MISSING'});continue
 h=None;s=[];found=[]
 def emit():
  if h is None:return
  pid=h.split()[0]
  if pid in ids:
   seq=''.join(s).rstrip('*');rec={'source':str(p),'id':pid,'header':h,'length':len(seq),'sha256':hashlib.sha256(seq.encode()).hexdigest(),'sequence':seq};records.append(rec);found.append(pid)
 op=gzip.open if p.suffix=='.gz' else open
 with op(p,'rt') as f:
  for line in f:
   if line.startswith('>'):emit();h=line[1:].strip();s=[]
   elif h and h.split()[0] in ids:s.append(line.strip())
  emit()
 for pid in ids:
  if pid not in found:records.append({'source':str(p),'id':pid,'status':'NOT_FOUND'})
(R/'identity_records.json').write_text(json.dumps(records,indent=2))
# Query actual historical representatives; never substitute a similarly named gene.
historical=[r for r in records if r.get('sequence') and r['source'].endswith('/orthogroups_rep_seq.fasta')]
with open(R/'historical_leads.faa','w') as f:
 for r in historical:f.write('>'+r['id']+'\n'+r['sequence']+'\n')
commands=[]
if historical:
 subprocess.run([str(W/'../miniconda3/envs/eggnog/bin/diamond'),'makedb','--in',str(D/'f12/references_v2/phosphatase_references_dedup.faa'),'--db',str(R/'refs')],check=True)
 cmd=[str(W/'../miniconda3/envs/eggnog/bin/diamond'),'blastp','--query',str(R/'historical_leads.faa'),'--db',str(R/'refs.dmnd'),'--out',str(R/'diamond.tsv'),'--outfmt','6','qseqid','sseqid','pident','length','qstart','qend','sstart','send','evalue','bitscore','qlen','slen','--ultra-sensitive','--max-target-seqs','0','--evalue','1e-3','--threads',os.environ['SLURM_CPUS_PER_TASK']]
 with open(R/'diamond.log','w') as log:subprocess.run(cmd,check=True,stdout=log,stderr=log)
 commands.append(cmd)
summary={'job_id':os.environ['SLURM_JOB_ID'],'historical_leads_found':len(historical),'identity_by_id':{pid:len({r['sha256'] for r in records if r.get('id')==pid and 'sha256'in r})==1 for pid in ids},'commands':commands,'diamond_hit_rows':sum(1 for l in open(R/'diamond.tsv') if l.strip()) if historical else None,'interpretation':'no-hit cannot establish function or disprove it; source identity required before interpreting prior lead annotations'}
(R/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
