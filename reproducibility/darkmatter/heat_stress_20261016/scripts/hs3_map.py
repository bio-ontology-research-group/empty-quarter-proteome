from common import *
import pandas as pd,numpy as np,gzip,subprocess,sys,time,pysam,shutil
R=N/'results/hs3';mode=os.environ.get('HS3_MODE','pilot');task=int(os.environ['SLURM_ARRAY_TASK_ID']);row=pd.read_csv(R/'pilot_samples.tsv',sep='\t').iloc[task];assert row.R1_accessible and row.R2_accessible
O=R/mode/row['sample'];O.mkdir(parents=True,exist_ok=True);tmp=O/'tmp';tmp.mkdir(exist_ok=True)
BIN=Path('/ibex/scratch/projects/c2014/rob/miniconda3/envs/eco_probes/bin');ref=W/'rebuild_20260729/picrust2_current_20260803/coverm_current/prepared/reps.mmi';limit=5000000
provenance={'job':os.environ['SLURM_JOB_ID'],'mode':mode,'sample':row['sample'],'target_genome':row.genome,'reference':str(ref),'reference_size':ref.stat().st_size,'source_control':row.to_dict(),'selection':'first5Mpairedrecords technical pilot only'if mode=='pilot'else 'all paired reads','filter':'MAPQ>=30, baseQ>=20, NM/aligned_query_length<=0.05, flags3844excluded, perposition one observation per template'}
fqs=[]
for mate in ['R1','R2']:
 source=Path(row[mate]);provenance[mate+'_stat']={'bytes':source.stat().st_size,'mtime_ns':source.stat().st_mtime_ns}
 if mode=='pilot':
  dest=tmp/(mate+'.fq');n=0;h=hashlib.sha256()
  with gzip.open(source,'rb')as fi,open(dest,'wb')as fo:
   while n<limit:
    lines=[fi.readline()for _ in range(4)]
    if not lines[0]:break
    assert all(lines) and lines[0].startswith(b'@') and lines[2].startswith(b'+')
    data=b''.join(lines);fo.write(data);h.update(data);n+=1
  provenance[mate+'_pilot_records']=n;provenance[mate+'_pilot_uncompressed_sha256']=h.hexdigest();fqs.append(str(dest))
 else:fqs.append(str(source))
if mode=='pilot':assert provenance['R1_pilot_records']==provenance['R2_pilot_records']
unsorted=tmp/'regional.unsorted.bam';bam=O/'regional.bam'
cmd=[str(BIN/'minimap2'),'-ax','sr','-t','12',str(ref),*fqs];vcmd=[str(BIN/'samtools'),'view','-b','-q','30','-F','3844','-L',row.bed,'-o',str(unsorted),'-'];provenance['mapping_command']=cmd;provenance['regional_filter_command']=vcmd
with open(O/'mapping.log','w')as log:
 p=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=log);q=subprocess.run(vcmd,stdin=p.stdout,stderr=log);p.stdout.close();rc=p.wait();assert rc==0 and q.returncode==0,(rc,q.returncode)
subprocess.run([str(BIN/'samtools'),'sort','-@','2','-m','2G','-o',str(bam),str(unsorted)],check=True);subprocess.run([str(BIN/'samtools'),'index',str(bam)],check=True)
loci=pd.read_csv(R/'pilot_loci.tsv',sep='\t');loci=loci[loci.genome==row.genome];stats=[];alleles=[];depths=[]
with pysam.AlignmentFile(bam,'rb')as af:
 for loc in loci.itertuples():
  length=loc.end-loc.start+1;depth=np.zeros(length,int)
  for col in af.pileup(loc.contig,loc.start-1,loc.end,truncate=True,stepper='nofilter',min_base_quality=0,max_depth=100000,ignore_overlaps=False):
   templates={}
   for pr in col.pileups:
    if pr.is_del or pr.is_refskip:continue
    al=pr.alignment
    if al.flag&3844 or al.mapping_quality<30 or al.get_tag('NM')/max(al.query_alignment_length,1)>.05:continue
    qi=pr.query_position;qual=al.query_qualities[qi]
    if qual<20:continue
    if al.query_name not in templates or qual>templates[al.query_name][0]:templates[al.query_name]=(qual,al.query_sequence[qi])
   ac={b:0 for b in 'ACGT'}
   for _,b in templates.values():
    if b in ac:ac[b]+=1
   dep=sum(ac.values());depth[col.reference_pos-(loc.start-1)]=dep
   if dep>=3:alleles.append({'sample':row['sample'],'protein':loc.protein,'contig':loc.contig,'position':col.reference_pos+1,**ac,'depth':dep})
  stats.append({'sample':row['sample'],'genome':row.genome,'domain':row.domain,'season':row.season,'family':loc.family,'role':loc.role,'protein':loc.protein,'length':length,'mean_depth':depth.mean(),'median_depth':float(np.median(depth)),'breadth1':(depth>=1).mean(),'breadth5':(depth>=5).mean()})
  depths.append((loc.protein,depth))
pd.DataFrame(stats).to_csv(O/'locus_coverage.tsv',sep='\t',index=False);pd.DataFrame(alleles,columns=['sample','protein','contig','position','A','C','G','T','depth']).to_csv(O/'allele_counts.tsv.gz',sep='\t',index=False);np.savez_compressed(O/'locus_depths.npz',**dict(depths))
provenance['samtools_version']=subprocess.check_output([str(BIN/'samtools'),'--version'],text=True).splitlines()[0];provenance['pysam_version']=pysam.__version__;provenance['regional_bam_sha256']=sha(bam)
(O/'provenance.json').write_text(json.dumps(provenance,indent=2,default=str));shutil.rmtree(tmp)
print(pd.DataFrame(stats).groupby('role')[['median_depth','breadth1']].median().to_string());print('DONE',mode,row['sample'],flush=True)
