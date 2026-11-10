from common import *
import gzip,re,numpy as np,pandas as pd,platform,sys
R=N/'results';R.mkdir(exist_ok=True)
F={'PF00011':'HSP20','PF00012':'HSP70','PF00226':'DnaJ','PF01025':'GrpE','PF00118':'Cpn60','PF00166':'GroES','PF00183':'HSP90','PF00574':'ClpP'}
idx=pd.read_csv(T/'index/genome_index.tsv',sep='\t').query("universe == 'soil_representative'").set_index('genome_id').sort_index();assert len(idx)==5438 and idx.index.is_unique
m=pd.read_csv(W/'metadata/sample_metadata.tsv',sep='\t',dtype={'site':str}).set_index('sample').sort_index();assert len(m)==551 and m.index.is_unique
source=[T/'index/genome_index.tsv',W/'metadata/sample_metadata.tsv',T/'panel/panel_verified.tsv',T/'panel/mechanisms.hmm',T/'code/ibex/03_hmmsearch_array.sh',T/'results/mechanism_copy_number_by_family.tsv']
manifest=[{'role':'source','path':str(p),'bytes':p.stat().st_size,'sha256':sha(p)}for p in source]
old=pd.read_csv(T/'results/mechanism_copy_number_by_family.tsv',sep='\t').set_index('genome_id')
counts=np.zeros((len(idx),len(F)),int);strict=counts.copy();hits=[]
for i,(gid,row)in enumerate(idx.iterrows()):
 p=T/'hmm/domtbl'/(gid+'.domtbl.gz');manifest.append({'role':'hmm','path':str(p),'bytes':p.stat().st_size,'sha256':sha(p)})
 found={a:{}for a in F}
 with gzip.open(p,'rt')as f:
  for l in f:
   if not l.strip() or l.startswith('#'):continue
   v=l.split(maxsplit=22);a=v[4].split('.')[0]
   if a not in F:continue
   pid=v[0];assert pid.split('::')[0]==gid
   cov=(int(v[16])-int(v[15])+1)/int(v[5]);desc=v[22];partial=re.search(r'partial=(\d\d)',desc).group(1)
   z=desc.split('#');rec={'genome':gid,'protein':pid,'family':F[a],'accession':a,'hmm_coverage':cov,'i_evalue':float(v[12]),'partial':partial,'start':int(z[1]),'end':int(z[2]),'strand':int(z[3]),'length_aa':int(v[2])}
   if pid not in found[a]or cov>found[a][pid]['hmm_coverage']:found[a][pid]=rec
 for j,a in enumerate(F):
  rr=list(found[a].values());counts[i,j]=len(rr);strict[i,j]=sum(r['hmm_coverage']>=.7 and r['partial']=='00'for r in rr);hits+=rr
  assert counts[i,j]==old.loc[gid,a],(gid,a,counts[i,j],old.loc[gid,a])
 if i%500==0:print('HMM',i,flush=True)
pd.DataFrame(hits).to_csv(R/'marker_loci.tsv.gz',sep='\t',index=False)
pd.DataFrame(counts,index=idx.index,columns=list(F.values())).to_csv(R/'genome_counts.tsv',sep='\t')
pd.DataFrame(strict,index=idx.index,columns=list(F.values())).to_csv(R/'genome_counts_strict.tsv',sep='\t')
idx.to_csv(R/'genome_index.tsv',sep='\t');m.to_csv(R/'metadata.tsv',sep='\t')
A=[];excluded=[];prof=W/'rebuild_20260729/picrust2_current_20260803/coverm_current/profiles_full_recompute'
for s in m.index:
 p=prof/s/'abund.tsv';d=pd.read_csv(p,sep='\t',index_col=0).iloc[:,0];assert set(d.index)-{'unmapped'}==set(idx.index)
 v=d.reindex(idx.index).to_numpy(float)
 if np.isnan(v).all():
  assert float(d.loc['unmapped'])==100,(s,'NaN profile is not wholly unmapped')
  v=np.zeros_like(v)
 assert np.isfinite(v).all()and(v>=0).all(),(s,'invalid profile')
 if v.sum()==0:excluded.append(s)
 A.append(v);manifest.append({'role':'abundance','path':str(p),'bytes':p.stat().st_size,'sha256':sha(p)})
A=np.array(A);np.savez_compressed(R/'abundance.npz',raw=A,samples=m.index.to_numpy(str),genomes=idx.index.to_numpy(str));assert len(excluded)==1
sizes=pd.read_csv(W/'06_darkmatter/lineage_darkness_clean/results/reps_all.tsv',sep='\t').set_index('genome').reindex(idx.index);assert sizes.length.notna().all();sizes.to_csv(R/'genome_sizes.tsv',sep='\t')
for p in [W/'06_darkmatter/lineage_darkness_clean/results/reps_all.tsv',W/'temporal_environment_20261014/inputs/environment.tsv']:
 manifest.append({'role':'covariate','path':str(p),'bytes':p.stat().st_size,'sha256':sha(p)})
env=pd.read_csv(W/'temporal_environment_20261014/inputs/environment.tsv',sep='\t');env.to_csv(R/'environment.tsv',sep='\t',index=False)
write(R/'input_manifest.tsv',manifest)
(R/'prepare_summary.json').write_text(json.dumps({'job':os.environ['SLURM_JOB_ID'],'host':socket.gethostname(),'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__,'genomes':len(idx),'samples':len(m),'excluded_zero_mapped':excluded,'marker_loci':len(hits),'families':F,'count_matrix_matches_existing':True,'counts_total':dict(zip(F.values(),map(int,counts.sum(0)))),'strict_total':dict(zip(F.values(),map(int,strict.sum(0)))),'interpretation':'PF00118 Cpn60 combines groupI/groupII; no heat-specific label; no ClpB inference from AAA domains'},indent=2))
print('DONE',len(hits),excluded,flush=True)
