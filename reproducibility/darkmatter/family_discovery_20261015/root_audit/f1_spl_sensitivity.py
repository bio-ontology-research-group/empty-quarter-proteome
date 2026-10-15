"""Positive-control sensitivity audit: known short photolyase missed by profile coverage gate."""
import csv,json,os,subprocess,hashlib
from pathlib import Path
W=Path('/ibex/scratch/projects/c2014/rob/metagenomics');D=W/'06_darkmatter/family_discovery_20261015_clean';R=D/'root_audit/f1_spl_sensitivity';R.mkdir(exist_ok=True)
diamond=str(W/'../miniconda3/envs/eggnog/bin/diamond');cmds=[]
def run(cmd):
 cmds.append(cmd)
 with open(R/'commands.log','a')as f:subprocess.run(cmd,check=True,stdout=f,stderr=f)
ref=D/'f12/references_v2/raw_records/ncbi/QCI56913.1.fasta'
assert ref.exists()
run([diamond,'makedb','--in',str(D/'root_audit/f1_full/cohort.faa'),'--db',str(R/'eq_full'),'--threads',os.environ['SLURM_CPUS_PER_TASK']])
run([diamond,'blastp','--query',str(ref),'--db',str(R/'eq_full.dmnd'),'--out',str(R/'hits.tsv'),'--ultra-sensitive','--max-target-seqs','0','--evalue','1e-5','--threads',os.environ['SLURM_CPUS_PER_TASK'],'--outfmt','6','qseqid','sseqid','pident','qstart','qend','sstart','send','evalue','bitscore','qlen','slen'])
old={r['id']for r in csv.DictReader(open(D/'root_audit/f1_full/candidates.tsv'),delimiter='\t')};hits=[]
for l in open(R/'hits.tsv'):
 q,t,pid,qs,qe,ts,te,ev,bits,ql,tl=l.split();qc=(int(qe)-int(qs)+1)/int(ql);tc=(int(te)-int(ts)+1)/int(tl)
 if qc>=.7 and tc>=.7:hits.append({'id':t,'pident':float(pid),'query_coverage':qc,'target_coverage':tc,'evalue':float(ev),'bitscore':float(bits),'in_profile_census':int(t in old)})
with open(R/'passing.tsv','w')as f:
 w=csv.DictWriter(f,fieldnames=['id','pident','query_coverage','target_coverage','evalue','bitscore','in_profile_census'],delimiter='\t');w.writeheader();w.writerows(hits)
(R/'summary.json').write_text(json.dumps({'job_id':os.environ['SLURM_JOB_ID'],'positive_control':'QCI56913.1 / MK474076; Xu2021doi10.1093/molbev/msab191','profile_gate_blindspot':'PF03441 bestdomain coverage64/203 <0.7 despite iE1.3e-9; specificknownSPLnotrecovered by coreprofiledefinition','coverage_passing_sequence_homologs':len(hits),'outside_profile_census':sum(not r['in_profile_census']for r in hits),'interpretation':'additionalSPL-like sequenceleads, notverifiedfunctions; coreprofilecensusnotexhaustiveCPF survey','commands':cmds},indent=2))
