from common import *
import pandas as pd,numpy as np,re
R=N/'results';O=R/'audit';checks=[]
def check(name,passed,detail):
 assert passed,(name,detail)
 checks.append({'check':name,'passed':bool(passed),'detail':str(detail)})
samples=pd.read_csv(R/'hs3/pilot_samples.tsv',sep='\t');counts=[]
for r in samples[samples.domain=='Bacteria'].itertuples():
 p=R/'hs3/full'/r.sample;log=(p/'mapping.log').read_text();n=sum(map(int,re.findall(r'mapped (\d+) sequences',log)));check('full_reads_'+r.sample,n==r.full_reads,{'mapped_input_reads':n,'prior_reads':r.full_reads});check('single_competitive_index_'+r.sample,log.count('loaded/built the index for')==1,log.count('loaded/built the index for'));counts.append({'sample':r.sample,'input_reads':n})
 a=pd.read_csv(p/'allele_counts.tsv.gz',sep='\t');check('allele_sum_'+r.sample,np.array_equal(a[['A','C','G','T']].sum(1),a.depth),len(a));check('unique_locus_positions_'+r.sample,not a.duplicated(['protein','position']).any(),len(a))
d=pd.read_csv(R/'hs3/dosage_full_pair.tsv',sep='\t');m=d[d.role=='marker'];check('all14markers_coverage',len(m)==14 and m.coverage_eligible.all(),len(m));check('dosage_identity',np.allclose(d.dosage_log2_hot_winter,np.log2(d.normalised_depth_hot/d.normalised_depth_winter)),len(d))
e=pd.read_csv(R/'hs2_exact/exact_locus_links.tsv.gz',sep='\t',keep_default_na=False);t=pd.read_csv(R/'hs2/thermo_genes.tsv.gz',sep='\t',keep_default_na=False);j=t.merge(e[['protein','catalog_id','status']],on='protein',suffixes=('_used','_exact'),validate='one_to_one');check('all_source_links_exact',((j.family=='')|((j.catalog_id_used==j.catalog_id_exact)&(j.status=='exact_unique'))).all(),len(j));check('sequence_link_coverage',len(j)==777177,len(j))
c=pd.read_csv(R/'hs2/candidate_context.tsv',sep='\t');check('candidate_single_locus_per_genome',not c.duplicated(['family','genome']).any(),len(c));pf=pd.read_csv(R/'hs2/candidate_pfam_hits.tsv',sep='\t');check('candidate_annotation_rows_valid',set(pf.protein)<=set(pd.read_csv(R/'hs2/candidate_all_loci.tsv',sep='\t').protein),len(pf))
write(O/'final_checks.tsv',checks);write(O/'mapping_input_counts.tsv',counts)
# Hash all retained artifacts, including remote-only matrices/BAMs; exclude this manifest itself.
rows=[]
for p in sorted(N.rglob('*')):
 if not p.is_file()or '__pycache__'in p.parts or '/tmp/'in str(p)or p.name=='artifact_manifest.tsv' or p.name==f"final_audit-{os.environ['SLURM_JOB_ID']}.log":continue
 if p.suffix not in ['.py','.sbatch','.tsv','.json','.gz','.npz','.bam','.bai','.log','.txt','.md','.pdf','.png','.bed']:continue
 rows.append({'path_relative_to_remote_root':str(p.relative_to(N)),'bytes':p.stat().st_size,'sha256':sha(p)})
write(O/'artifact_manifest.tsv',rows);print('DONE final audit',len(checks),'artifacts',len(rows),counts)
