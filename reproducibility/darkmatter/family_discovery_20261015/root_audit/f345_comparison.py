"""Independent complete-gene family comparisons and exact source-preserving embedding inputs."""
import csv,hashlib,json,os,random,re,statistics,subprocess
from pathlib import Path
W=Path('/ibex/scratch/projects/c2014/rob/metagenomics');D=W/'06_darkmatter/family_discovery_20261015_clean';F=D/'f345';R=D/'root_audit/f345_comparison';R.mkdir(exist_ok=True)
def tab(p):return list(csv.DictReader(open(p),delimiter='\t'))
def fasta(p):
 h=None;s=[]
 for l in open(p):
  if l.startswith('>'):
   if h:yield h,''.join(s).rstrip('*')
   h=l[1:].strip();s=[]
  else:s.append(l.strip())
 if h:yield h,''.join(s).rstrip('*')
def write(p,rows,fields):
 with open(p,'w')as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
refs={h.split('|')[0].split()[0]:s for h,s in fasta(F/'references/all_f345_references.faa')};refmeta=tab(F/'references/verified_references_metadata.tsv');allmeta=[];allseq={};summaries={};commands=[]
classes={'Type1':('f3_rhodopsin','f3','Type-1_Microbial_Rhodopsin','Type1_Microbial'),'Helio':('f3_rhodopsin','f3','Heliorhodopsin_Type-3','Heliorhodopsin'),'EctC':('f4_ectc','f4','Ectoine_Synthase_EctC','Ectoine_synthase'),'BetaCA':('f5_carbonic_anhydrase','f5','Beta-Carbonic_Anhydrase','Beta_CA'),'GammaCA':('f5_carbonic_anhydrase','f5','Gamma-Carbonic_Anhydrase','Gamma_CA')}
source_map={'EQ_MAG':'EQ','Environmental_SMAG':'SMAG','Environmental_GEM':'GEM'}
def run(cmd):
 commands.append(cmd)
 with open(R/'commands.log','a')as log:log.write(' '.join(cmd)+'\n');log.flush();subprocess.run(cmd,check=True,stdout=log,stderr=log)
for family,(directory,prefix,cls,cat)in classes.items():
 P=F/'results/full_v1'/directory;out=R/family;out.mkdir(exist_ok=True)
 candidates=tab(P/f'{prefix}_candidates_summary.tsv');seqs={h.split()[0]:s for h,s in fasta(P/f'{prefix}_all_candidates.faa')};rows=[];ss={};excluded=[]
 def add(pid,source,s,tax='',partial='reference',original=''):
  uid=family+'__'+source+'__'+pid;assert uid not in ss
  if '*'in s or not s or len(s)>2048:excluded.append({'id':uid,'reason':'internal_stop_empty_or_over2048'});return
  sha=hashlib.sha256(s.encode()).hexdigest();phylum=next((t[3:]for t in tax.split(';')if t.startswith('p__')),'Unknown');genus=next((t[3:]for t in tax.split(';')if t.startswith('g__')),'Unknown')
  rows.append({'id':uid,'source':source,'family_class':family,'length':len(s),'phylum':phylum or 'Unknown','genus':genus or 'Unknown','clean_sha256':hashlib.sha256(re.sub('[UZOB]','X',s).encode()).hexdigest(),'sequence_sha256':sha,'original_id':pid,'partial':partial});ss[uid]=s
 for r in candidates:
  if r['assigned_class']!=cls:continue
  pid=r['protein_id'];assert pid in seqs;seq=seqs[pid];assert hashlib.sha256(seq.encode()).hexdigest()==r['seq_sha256'];assert len(seq)==int(r['length_aa'])
  assert (r['prodigal_partial']=='00')==(r['is_complete_gene']=='1')
  if r['prodigal_partial']!='00':excluded.append({'id':pid,'reason':'partial_or_unknown'});continue
  add(pid,source_map[r['cohort']],seq,r['gtdbtk_taxonomy'],'00')
 for r in refmeta:
  if r['category']==cat:
   assert r['accession']in refs,r['accession'];add(r['accession'],'curated_reference',refs[r['accession']])
 assert any(r['source']=='curated_reference'for r in rows),family
 with open(out/'combined.faa','w')as f:
  for uid,s in ss.items():f.write('>'+uid+'\n'+s+'\n')
 fields=list(rows[0]);write(out/'metadata.tsv',rows,fields);write(out/'excluded.tsv',excluded,['id','reason'])
 diamond=str(W/'../miniconda3/envs/eggnog/bin/diamond');run([diamond,'makedb','--in',str(out/'combined.faa'),'--db',str(out/'db')])
 with open(out/'eq_queries.faa','w')as f:
  for r in rows:
   if r['source']=='EQ':f.write('>'+r['id']+'\n'+ss[r['id']]+'\n')
 run([diamond,'blastp','--query',str(out/'eq_queries.faa'),'--db',str(out/'db.dmnd'),'--out',str(out/'hits.tsv'),'--ultra-sensitive','--max-target-seqs','0','--evalue','1e-3','--threads',os.environ['SLURM_CPUS_PER_TASK'],'--outfmt','6','qseqid','sseqid','pident','qstart','qend','sstart','send','evalue','bitscore','qlen','slen'])
 meta={r['id']:r for r in rows};best={}
 for l in open(out/'hits.tsv'):
  q,t,pid,qs,qe,ts,te,ev,bits,ql,tl=l.split();qc=(int(qe)-int(qs)+1)/int(ql);tc=(int(te)-int(ts)+1)/int(tl)
  if qc<.7 or tc<.7 or q==t:continue
  key=(q,meta[t]['source']);rec={'query':q,'source':key[1],'target':t,'identity':float(pid),'query_coverage':qc,'target_coverage':tc,'bitscore':float(bits)}
  if key not in best or rec['bitscore']>best[key]['bitscore']:best[key]=rec
 write(out/'closest_hits.tsv',list(best.values()),['query','source','target','identity','query_coverage','target_coverage','bitscore'])
 neq=sum(r['source']=='EQ'for r in rows);summary={'source_counts':{src:sum(r['source']==src for r in rows)for src in ['EQ','SMAG','GEM','curated_reference']},'comparisons':{},'cluster_expansion':{}}
 for src in ['SMAG','GEM','curated_reference']:
  vals=[r['identity']for k,r in best.items()if k[1]==src];summary['comparisons'][src]={'EQ_queries':neq,'coverage_passing_hits':len(vals),'no_coverage_passing_hit':neq-len(vals),'median_identity':statistics.median(vals)if vals else None,'below50amonghits':sum(v<50 for v in vals)}
 mm=str(W/'../miniconda3/envs/mmseqs2/bin/mmseqs')
 for threshold in [.5,.8]:
  label=str(int(threshold*100));pre=out/('cluster'+label)
  run([mm,'easy-cluster',str(out/'combined.faa'),str(pre),str(out/('tmp'+label)),'--min-seq-id',str(threshold),'-c','0.8','--cov-mode','0','--threads',os.environ['SLURM_CPUS_PER_TASK'],'-v','1'])
  cmap={}
  for l in open(str(pre)+'_cluster.tsv'):
   rep,mem=l.split();assert mem in meta;cmap[mem]=rep
  assert set(cmap)==set(meta)
  refclus={cmap[r['id']]for r in rows if r['source']=='curated_reference'};source_ids={s:[]for s in ['EQ','SMAG','GEM']}
  for src in source_ids:
   seen=set()
   for r in sorted(rows,key=lambda r:r['id']):
    if r['source']==src and r['sequence_sha256']not in seen:seen.add(r['sequence_sha256']);source_ids[src].append(r['id'])
  n=min(100,*map(len,source_ids.values()));metrics={'equal_addition_n':n,'sources':{}}
  for src,ii in source_ids.items():
   draws=[len({cmap[i]for i in random.Random(20261015+k).sample(ii,n)}-refclus)for k in range(100)]if n else[]
   metrics['sources'][src]={'nonredundant_sequences':len(ii),'total_clusters':len({cmap[i]for i in ii}),'clusters_outside_reference':len({cmap[i]for i in ii}-refclus),'balanced_mean_clusters_outside_reference':statistics.mean(draws)if draws else None}
  summary['cluster_expansion'][label]=metrics
 summaries[family]=summary;allmeta.extend(rows);allseq.update(ss)
with open(R/'all_families.faa','w')as f:
 for uid,s in allseq.items():f.write('>'+uid+'\n'+s+'\n')
write(R/'metadata.tsv',allmeta,list(allmeta[0]))
# Deterministic balanced plot subset; retain full census, comparisons and all embeddings separately.
plot=[]
for family in classes:
 for source in ['EQ','SMAG','GEM','curated_reference']:
  rr=[r for r in allmeta if r['family_class']==family and r['source']==source];seen=set();unique=[]
  for r in sorted(rr,key=lambda r:hashlib.sha256(r['id'].encode()).hexdigest()):
   if r['sequence_sha256']not in seen:seen.add(r['sequence_sha256']);unique.append(r)
  plot.extend(unique if source=='curated_reference'else unique[:300])
write(R/'plot_metadata.tsv',plot,list(allmeta[0]))
(R/'summary.json').write_text(json.dumps({'job_id':os.environ['SLURM_JOB_ID'],'families':summaries,'commands':commands,'limits':['Profile membership not catalytic activity; source references include annotated homologs.','Named sampled environmental panels; no geographical exclusivity or universal novelty.','Clusters conditional on joint MMseqs2 sequence clustering thresholds and taxonomic composition.','Whole-protein alignments can miss short domains in fusions; do not interpret no-hit as no homolog.','Plot sampling up to300 unique sequences per source per family, all references; complete genes only.']},indent=2))
