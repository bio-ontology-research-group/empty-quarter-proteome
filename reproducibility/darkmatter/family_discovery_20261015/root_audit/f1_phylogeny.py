"""Domain-profile alignments and exploratory phylogenetic breadth with equal-size additions."""
import csv,json,os,random,subprocess,sys
from pathlib import Path
from Bio import AlignIO,Phylo,SeqIO
sys.setrecursionlimit(20000)
W=Path('/ibex/scratch/projects/c2014/rob/metagenomics');D=W/'06_darkmatter/family_discovery_20261015_clean';R=D/'root_audit/f1_phylogeny';R.mkdir(exist_ok=True);C=D/'root_audit/f1_comparison'
meta={r['id']:r for r in csv.DictReader(open(C/'metadata.tsv'),delimiter='\t')};seqs=SeqIO.to_dict(SeqIO.parse(C/'combined.faa','fasta'));summary={};commands=[]
def run(cmd,out=None):
 commands.append(cmd)
 with open(R/'commands.log','a')as log:
  log.write(' '.join(cmd)+'\n');log.flush()
  if out:
   with open(out,'w')as f:subprocess.run(cmd,check=True,stdout=f,stderr=log)
  else:subprocess.run(cmd,check=True,stdout=log,stderr=log)
run(['hmmsearch','--cpu',os.environ['SLURM_CPUS_PER_TASK'],'--domtblout',str(R/'domains.tsv'),'-o',str(R/'hmm.txt'),str(D/'root_audit/f1_profile_pilot/cpf.hmm'),str(C/'combined.faa')])
by={}
for l in open(R/'domains.tsv'):
 if l.startswith('#')or not l.strip():continue
 p=l.split()
 if float(p[12])<=1e-5 and (int(p[16])-int(p[15])+1)/int(p[5])>=.7:by.setdefault(p[4].split('.')[0],set()).add(p[0])
for profile,ids in by.items():
 # Within-source exact duplicates collapsed; cross-source membership preserved in metadata.
 seen=set();keep=[]
 for i in sorted(ids):
  key=(meta[i]['source'],meta[i]['sequence_sha256'])
  if key not in seen:keep.append(i);seen.add(key)
 SeqIO.write([seqs[i]for i in keep],R/f'{profile}.faa','fasta')
 run(['hmmalign','--trim','--outformat','Stockholm',str(D/f'f12/references_v2/{profile}.hmm'),str(R/f'{profile}.faa')],R/f'{profile}.sto')
 aln=AlignIO.read(R/f'{profile}.sto','stockholm');rf=aln.column_annotations['reference_annotation'];cols=[j for j,c in enumerate(rf)if c not in '.- ']
 aliases={i:f'tip{n:06d}' for n,i in enumerate(keep)}
 (R/f'{profile}.tip_mapping.json').write_text(json.dumps(aliases,indent=2))
 with open(R/f'{profile}.aligned.faa','w')as f:
  for rec in aln:f.write('>'+aliases[rec.id]+'\n'+''.join(str(rec.seq)[j].upper()for j in cols)+'\n')
 run([str(W/'../miniconda3/envs/selection/bin/FastTree'),'-wag','-gamma',str(R/f'{profile}.aligned.faa')],R/f'{profile}.nwk')
 tree=Phylo.read(R/f'{profile}.nwk','newick')
 reverse={v:k for k,v in aliases.items()}
 for tip in tree.get_terminals():tip.name=reverse[tip.name]
 alltips={x.name for x in tree.get_terminals()};edges=[]
 for cl in tree.find_clades():
  if cl.branch_length and cl.branch_length>0:edges.append(({t.name for t in cl.get_terminals()},cl.branch_length))
 def span(tips):return sum(length for descendants,length in edges if tips&descendants and tips-descendants)
 sources={s:{i for i in alltips if meta[i]['source']==s}for s in ['EQ','SMAG','GEM','curated_reference']};refs=sources['curated_reference'];base=span(refs)
 n=min(200,*[len(sources[s])for s in ['EQ','SMAG','GEM']]);draws=[]
 if n and len(refs)>=2:
  for seed in range(100):
   for s in ['EQ','SMAG','GEM']:
    chosen=set(random.Random(20261015+seed).sample(sorted(sources[s]),n));draws.append({'seed':20261015+seed,'source':s,'n':n,'added_branch_length':span(refs|chosen)-base})
 with open(R/f'{profile}.balanced_additions.tsv','w')as f:
  w=csv.DictWriter(f,fieldnames=['seed','source','n','added_branch_length'],delimiter='\t');w.writeheader();w.writerows(draws)
 summary[profile]={'sequences':len(keep),'alignment_match_columns':len(cols),'sources':{s:len(v)for s,v in sources.items()},'reference_scaffold_branch_length':base,'equal_addition_n':n,'mean_added_branch_length':{s:sum(r['added_branch_length']for r in draws if r['source']==s)/100 for s in ['EQ','SMAG','GEM']}if draws else{},'limits':'exploratory domain tree; equal-size reference-conditioned breadth, not novel function or taxonomic clades; no taxonomic matching; topology estimated jointly; FastTree local supports, no resampled alignment bootstrap'}
(R/'summary.json').write_text(json.dumps({'job_id':os.environ['SLURM_JOB_ID'],'profiles':summary,'commands':commands},indent=2))
