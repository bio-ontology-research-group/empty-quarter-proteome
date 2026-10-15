"""Audit full embeddings, report full native-space distances, plot declared deterministic subset."""
import csv,hashlib,json,os,subprocess
from pathlib import Path
import numpy as np
from scipy.spatial.distance import cdist
D=Path('/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean');C=D/'root_audit/f345_comparison';E=D/'root_audit/f345_embeddings';R=D/'root_audit/f345_plots';R.mkdir(exist_ok=True)
def tab(p):return list(csv.DictReader(open(p),delimiter='\t'))
def write(p,rr):
 with open(p,'w')as f:
  w=csv.DictWriter(f,fieldnames=list(rr[0]),delimiter='\t');w.writeheader();w.writerows(rr)
m=tab(C/'metadata.tsv');archive=np.load(E/'embeddings.npz',allow_pickle=True);z={k:archive[k]for k in archive.files};archive.close();ids=list(map(str,z['ids']));lookup={pid:i for i,pid in enumerate(ids)};assert len(lookup)==len(ids);assert set(lookup)=={r['id']for r in m};assert z['embeddings'].shape==(len(m),1024) and np.isfinite(z['embeddings']).all()
for r in m:
 i=lookup[r['id']];assert str(z['raw_sequence_sha256'][i])==r['sequence_sha256'];assert str(z['clean_sequence_sha256'][i])==r['clean_sha256'];assert int(z['lengths'][i])==int(r['length'])
full=[];summary={};families=sorted({r['family_class']for r in m})
for fam in families:
 rr=[r for r in m if r['family_class']==fam];refs=[r for r in rr if r['source']=='curated_reference'];x=z['embeddings'][[lookup[r['id']]for r in rr]].astype('float32');refx=z['embeddings'][[lookup[r['id']]for r in refs]].astype('float32');dd=cdist(x,refx,metric='cosine')
 for i,r in enumerate(rr):
  for j,ref in enumerate(refs):
   if r['id']==ref['id']:dd[i,j]=np.inf
  j=int(np.argmin(dd[i]));full.append({'id':r['id'],'family':fam,'source':r['source'],'nearest_curated_reference':refs[j]['id'],'cosine_distance':float(dd[i,j]),'length':r['length']})
 summary[fam]={src:{'n':len(vals),'median_nearest_curated_cosine':float(np.median(vals))}for src in sorted({r['source']for r in rr})if(vals:=[r['cosine_distance']for r in full if r['family']==fam and r['source']==src])}
write(R/'full_distances.tsv',full);(R/'full_distance_summary.json').write_text(json.dumps({'job_id':os.environ['SLURM_JOB_ID'],'families':summary,'limits':['Native-space descriptive distances; not trained-model novelty, catalytic function or geographic exclusivity.','Reference distances exclude the same ID but not distinct reference accessions with identical sequence.','No taxonomic or architectural adjustment; environment comparisons conditional on named panels.']},indent=2))
plots=tab(C/'plot_metadata.tsv')
for fam in families:
 out=R/fam;out.mkdir(exist_ok=True);rr=[r for r in plots if r['family_class']==fam];ii=[lookup[r['id']]for r in rr];np.savez_compressed(out/'plot_embeddings.npz',**{k:z[k][ii]for k in z});write(out/'plot_metadata.tsv',rr)
 subprocess.run(['/ibex/scratch/projects/c2014/rob/envs/gspa-neural/bin/python','/home/hohndor/eq-family-discovery-20261015/cpu/family_embedding_analysis.py','--embeddings',str(out/'plot_embeddings.npz'),'--metadata',str(out/'plot_metadata.tsv'),'--output-dir',str(out),'--seeds','42,123,20261015','--perplexities','5,30,50','--tsne-max-iter','1000'],check=True)
(R/'COMPLETE').write_text(os.environ['SLURM_JOB_ID']+'\n')
