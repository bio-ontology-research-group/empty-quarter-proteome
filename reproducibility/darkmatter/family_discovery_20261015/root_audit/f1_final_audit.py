"""Independent compact-result checks and complete-gene phylogenetic sensitivity."""
import csv,hashlib,json,os,random,statistics
from collections import defaultdict
from pathlib import Path
import numpy as np
from scipy.stats import wilcoxon
from Bio import Phylo
W=Path('/ibex/scratch/projects/c2014/rob/metagenomics');D=W/'06_darkmatter/family_discovery_20261015_clean';R=D/'root_audit/f1_final_audit';R.mkdir(exist_ok=True)
def tab(p):return list(csv.DictReader(open(p),delimiter='\t'))
def digest(p):
 h=hashlib.sha256()
 with open(p,'rb')as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
# Verify pinned source/full-output integrity on worker, including the full FASTA.
hashchecks=[]
for l in open(D/'root_audit/f1_full/checksums.tsv'):
 expected,name=l.strip().split('\t');actual=digest(D/'root_audit/f1_full'/name);assert actual==expected,name;hashchecks.append(name)
rows=tab(D/'root_audit/f1_full/candidates.tsv');domains={}
for l in open(D/'root_audit/f1_full/domains.tsv'):
 if l.startswith('#')or not l.strip():continue
 p=l.split();cov=(int(p[16])-int(p[15])+1)/int(p[5])
 if float(p[12])<=1e-5 and cov>=.7:domains[p[0]]=True
assert set(domains)=={r['id']for r in rows}
summary={'job_id':os.environ['SLURM_JOB_ID'],'full_output_checksums_verified':hashchecks,'independent_domain_gate_count':len(domains),'complete_count':sum(r['partial']=='00'for r in rows),'carrier_genomes':len({r['mag']for r in rows})}
# Ecological effect and test independently reassembled from published compact sample table.
eco=tab(D/'ecology/results_root_v3/f1_community_by_sample.tsv');tests=tab(D/'ecology/results_root_v3/f1_compartment_tests.tsv');eco_audit=[]
for trait in ['cpf_complete','cpf_complete_exclusive']:
 cells=defaultdict(list)
 for r in eco:
  if r['trait']==trait:cells[(r['site'],r['trip'],r['compartment'])].append(float(r['abund_mapped_pct']))
 visits=defaultdict(list)
 for site,trip in sorted({k[:2]for k in cells}):
  if (site,trip,'Surface')in cells and (site,trip,'Deep')in cells:
   visits[site].append(statistics.median(cells[(site,trip,'Surface')])-statistics.median(cells[(site,trip,'Deep')]))
 differences=[statistics.median(v)for v in visits.values()];effect=statistics.median(differences);p=float(wilcoxon(differences).pvalue)
 expected=next(r for r in tests if r['trait']==trait and r['metric']=='abund_mapped_pct')
 assert len(visits)==int(expected['n_paired_sites']);assert abs(effect-float(expected['delta_median']))<.0003;assert abs(p-float(expected['p_wilcoxon']))<.00005
 eco_audit.append({'trait':trait,'sites':len(visits),'visits':sum(map(len,visits.values())),'paired_median_difference':effect,'p_wilcoxon_recomputed_from_rounded_sample_table':p})
summary['ecology_independent_audit']=eco_audit
ph=tab(D/'ecology/results_root_v3/f1_within_phylum_tests.tsv');order=sorted(range(len(ph)),key=lambda i:float(ph[i]['p_wilcoxon']));running=0
for rank,i in enumerate(order):
 running=max(running,(len(ph)-rank)*float(ph[i]['p_wilcoxon']));ph[i]['p_holm_all_phylum_tests']=min(1,running)
with open(R/'within_phylum_multiplicity.tsv','w')as f:
 w=csv.DictWriter(f,fieldnames=list(ph[0]),delimiter='\t');w.writeheader();w.writerows(ph)
# Equal-size breadth using complete-gene-only tips on frozen joint domain trees.
meta={r['id']:r for r in tab(D/'root_audit/f1_complete_panel/metadata_complete.tsv')};pd={}
for profile in ['PF00875','PF03441','PF04244']:
 P=D/'root_audit/f1_phylogeny';aliases=json.load(open(P/f'{profile}.tip_mapping.json'));rev={v:k for k,v in aliases.items()};tree=Phylo.read(P/f'{profile}.nwk','newick')
 for tip in tree.get_terminals():tip.name=rev[tip.name]
 tips={t.name for t in tree.get_terminals()}&set(meta);sources={s:{i for i in tips if meta[i]['source']==s}for s in ['EQ','SMAG','GEM','curated_reference']};refs=sources['curated_reference'];edges=[({t.name for t in c.get_terminals()},c.branch_length)for c in tree.find_clades()if c.branch_length and c.branch_length>0]
 def span(t):return sum(b for x,b in edges if t&x and t-x)
 n=min(200,*[len(sources[s])for s in ['EQ','SMAG','GEM']]);result={'source_counts':{s:len(v)for s,v in sources.items()},'equal_addition_n':n,'reference_tips':len(refs)}
 if len(refs)>=2 and n:
  base=span(refs);draws={s:[span(refs|set(random.Random(20261015+k).sample(sorted(sources[s]),n)))-base for k in range(100)]for s in ['EQ','SMAG','GEM']}
  result['mean_added_length']={s:statistics.mean(v)for s,v in draws.items()};result['descriptive_draw_2.5_97.5_percentiles']={s:list(np.percentile(v,[2.5,97.5]))for s,v in draws.items()}
 pd[profile]=result
summary['complete_gene_tree_sensitivity']=pd
summary['limitations']=['Tree topology fitted jointly before fragment exclusion; exploratory short-domain branch lengths, not new functional clades.','Ecology association not causal selection; within-phylum stratification does not remove species turnover.','551 metadata profiles checked; one zero-mapped sample excluded, 550 contribute ratios, 37 strictly paired sites/69 visits.','SPL positive-control failure means core-profile screen is not exhaustive CPF discovery.','Reference-only tSNE separation is not biological novelty.','Ecology provenance v3 phrase indicates ecological selection is rejected; association only.']
(R/'audit.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
