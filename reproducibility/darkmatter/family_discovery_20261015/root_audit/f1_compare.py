"""Reference-relative sequence analysis; all claims computed, all source memberships retained."""
import csv,hashlib,json,os,re,subprocess,statistics
from pathlib import Path
W=Path('/ibex/scratch/projects/c2014/rob/metagenomics');D=W/'06_darkmatter/family_discovery_20261015_clean';R=D/'root_audit/f1_comparison';R.mkdir(exist_ok=True)
def fasta(p):
 h=None;s=[]
 for l in open(p):
  if l.startswith('>'):
   if h:yield h,''.join(s).rstrip('*')
   h=l[1:].strip();s=[]
  else:s.append(l.strip())
 if h:yield h,''.join(s).rstrip('*')
def tab(p):return list(csv.DictReader(open(p),delimiter='\t'))
def write(p,rows,fields=None):
 with open(p,'w') as f:
  w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
commands=[]
def run(cmd):
 commands.append(cmd)
 with open(R/'commands.log','a') as f:
  f.write(' '.join(cmd)+'\n');f.flush();subprocess.run(cmd,check=True,stdout=f,stderr=f)
# Gate the curated reference collection with same CPF profiles; includes signalling proteins.
refs=D/'f12/references_v2/cpf_references_dedup.faa'
run(['hmmsearch','--cpu',os.environ['SLURM_CPUS_PER_TASK'],'--domtblout',str(R/'reference_domains.tsv'),'-o',str(R/'reference_hmm_report.txt'),str(D/'root_audit/f1_profile_pilot/cpf.hmm'),str(refs)])
refids=set()
for l in open(R/'reference_domains.tsv'):
 if l.startswith('#') or not l.strip():continue
 p=l.split()
 if float(p[12])<=1e-5 and (int(p[16])-int(p[15])+1)/int(p[5])>=.7:refids.add(p[0])
tax={}
for p in (W/'05_refine/_gtdbtk_r232').glob('gtdbtk.*.summary.tsv'):
 for r in tab(p):tax[r['user_genome'].removesuffix('.fa')]=r['classification']
metarefs={r['reference_id']:r for r in tab(D/'f12/references_v2/verified_references_metadata_dedup.tsv')}
files={'EQ':D/'root_audit/f1_full/candidates.faa','SMAG':D/'root_audit/f1_comparators/smag_cpf_candidates.faa','GEM':D/'root_audit/f1_comparators/gem_cpf_candidates.faa','curated_reference':refs}
rows=[];seqs={};groups={};excluded=[]
for source,p in files.items():
 for h,s in fasta(p):
  pid=h.split()[0]
  if source=='curated_reference' and pid not in refids:continue
  partial=re.search(r'partial=(\d+)',h)
  if source=='EQ' and (not partial or partial.group(1)!='00'):
   excluded.append({'source':source,'id':pid,'reason':'incomplete_or_unknown_CDS'});continue
  if '*' in s or not s or len(s)>2048:
   excluded.append({'source':source,'id':pid,'reason':'internal_stop_empty_or_over2048'});continue
  sha=hashlib.sha256(s.encode()).hexdigest();uid=source+'__'+pid
  lineage=tax.get(pid.split('::')[0],'') if source=='EQ' else metarefs.get(pid,{}).get('lineage','')
  phylum=next((x[3:] for x in lineage.split(';') if x.startswith('p__')),'Unknown');genus=next((x[3:] for x in lineage.split(';') if x.startswith('g__')),'Unknown')
  r={'id':uid,'source':source,'family_class':'CPF_profile_supported','length':len(s),'phylum':phylum or 'Unknown','genus':genus or 'Unknown','clean_sha256':hashlib.sha256(re.sub('[UZOB]','X',s).encode()).hexdigest(),'sequence_sha256':sha,'original_id':pid,'partial':partial.group(1) if partial else 'unknown'}
  rows.append(r);seqs[uid]=s;groups.setdefault(sha,[]).append(uid)
assert len(seqs)==len(rows)
# All source identities retained, including duplicate sequences. Embed unique per source; maps remain exact.
with open(R/'combined.faa','w') as f:
 for uid,s in seqs.items():f.write('>'+uid+'\n'+s+'\n')
write(R/'metadata.tsv',rows);write(R/'excluded.tsv',excluded,['source','id','reason'])
write(R/'sequence_memberships.tsv',[{'sha256':sha,'ids':';'.join(ids),'sources':';'.join(sorted({i.split('__')[0] for i in ids}))} for sha,ids in groups.items()])
diamond=str(W/'../miniconda3/envs/eggnog/bin/diamond')
run([diamond,'makedb','--in',str(R/'combined.faa'),'--db',str(R/'combined')])
run([diamond,'blastp','--query',str(R/'combined.faa'),'--db',str(R/'combined.dmnd'),'--out',str(R/'alignments.tsv'),'--ultra-sensitive','--max-target-seqs','0','--evalue','1e-3','--threads',os.environ['SLURM_CPUS_PER_TASK'],'--outfmt','6','qseqid','sseqid','pident','qstart','qend','sstart','send','evalue','bitscore','qlen','slen'])
meta={r['id']:r for r in rows};best={}
for l in open(R/'alignments.tsv'):
 q,t,pid,qs,qe,ts,te,ev,bits,ql,tl=l.split();qc=(int(qe)-int(qs)+1)/int(ql);tc=(int(te)-int(ts)+1)/int(tl)
 if q==t or qc<.7 or tc<.7:continue
 target=meta[t]['source'];key=(q,target)
 rec={'query':q,'target_source':target,'target':t,'identity':float(pid),'query_coverage':qc,'target_coverage':tc,'evalue':float(ev),'bitscore':float(bits)}
 if key not in best or rec['bitscore']>best[key]['bitscore']:best[key]=rec
write(R/'closest_hits.tsv',list(best.values()),['query','target_source','target','identity','query_coverage','target_coverage','evalue','bitscore'])
summary={'job_id':os.environ['SLURM_JOB_ID'],'sources':{},'gates':{'profile_ievalue':1e-5,'profile_coverage':.7,'alignment_both_coverage':.7},'limits':'sampled environmental panels; incomplete EQ excluded; reference curation is not experimental validation; no activity/endemicity/adaptation inference; tSNE illustrative; all source identities retained', 'commands':commands}
for src in files:
 rr=[r for r in rows if r['source']==src];summary['sources'][src]={'proteins':len(rr),'unique_sequences':len({r['sequence_sha256'] for r in rr})}
eq=[r for r in rows if r['source']=='EQ'];comparisons={}
for target in ['curated_reference','SMAG','GEM']:
 vals=[best[(r['id'],target)]['identity'] for r in eq if (r['id'],target)in best]
 comparisons[target]={'queries':len(eq),'with_coverage_passing_hit':len(vals),'without_coverage_passing_hit':len(eq)-len(vals),'identity_median':statistics.median(vals)if vals else None,'below_50_identity_among_hits':sum(x<50 for x in vals)}
summary['EQ_comparisons']=comparisons
(R/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
with open(R/'checksums.tsv','w') as f:
 for p in sorted(R.iterdir()):
  if p.is_file() and p.name!='checksums.tsv':
   h=hashlib.sha256()
   with open(p,'rb')as inp:
    for b in iter(lambda:inp.read(1048576),b''):h.update(b)
   f.write(h.hexdigest()+'\t'+p.name+'\n')
