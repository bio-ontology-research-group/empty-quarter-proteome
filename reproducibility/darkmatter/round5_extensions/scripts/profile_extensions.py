#!/usr/bin/env python3
"""Exploratory temporal-label and paired functional-profile nulls; run on IBEX."""
from pathlib import Path
import csv,gzip,json,hashlib
import numpy as np,pandas as pd
from scipy.stats import rankdata
W=Path('/ibex/scratch/projects/c2014/rob/metagenomics');O=W/'06_darkmatter/p4_review_20261013/results/round5_extensions';O.mkdir(exist_ok=True)
pins={}
def pin(p):
 p=Path(p);h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 pins[str(p)]=h.hexdigest();return p
meta=pd.read_csv(pin(W/'metadata/sample_metadata.tsv'),sep='\t',dtype=str).set_index('sample');assert meta.index.is_unique
# H2 global label permutations preserve co-occurrence across every family.
p=W/'06_darkmatter/darkfam_basis_clean/results/fam_presence.tsv';d=pd.read_csv(pin(p),sep='\t',index_col=0);assert d.shape==(106042,551),d.shape
assert set(np.unique(d.values))<={0,1};A=d.to_numpy(dtype=np.float32);cols=meta.loc[d.columns];trip=cols.trip_num.astype(int).to_numpy()-1;assert set(trip)==set(range(5))
blocks=[np.asarray(v) for v in cols.reset_index().groupby(['site','compartment']).indices.values()]
def allfive(labels):return ((A@np.eye(5,dtype=np.float32)[labels])>0).all(axis=1)
obs=allfive(trip);occupancy=A.sum(axis=1);bins=[('all',np.ones(len(d),dtype=bool))]+[(f'{a}-{b-1}',(occupancy>=a)&(occupancy<b)) for a,b in zip([1,10,25,50,100,200],[10,25,50,100,200,552])]
rng=np.random.default_rng(20261015);nulls={name:[] for name,mask in bins};counts=np.zeros(len(d),dtype=int)
for _ in range(499):
 labels=trip.copy()
 for ix in blocks:labels[ix]=rng.permutation(labels[ix])
 z=allfive(labels);counts+=z
 for name,mask in bins:nulls[name].append(float(z[mask].mean()) if mask.any() else None)
rows=[]
for name,mask in bins:
 if not mask.any():continue
 x=np.asarray(nulls[name]);v=float(obs[mask].mean());rows.append(dict(group=name,families=int(mask.sum()),observed=v,null_mean=float(x.mean()),null_low=float(np.quantile(x,.025)),null_high=float(np.quantile(x,.975)),observed_minus_null=v-float(x.mean()),p_upper=(1+int((x>=v).sum()))/500,p_lower=(1+int((x<=v).sum()))/500))
pd.DataFrame(rows).to_csv(O/'temporal_summary.tsv',sep='\t',index=False);pd.DataFrame(nulls).to_csv(O/'temporal_null_draws.tsv',sep='\t',index=False)
pd.DataFrame({'family_id':d.index,'occupied_samples':occupancy.astype(int),'observed_all_five':obs.astype(int),'null_all_five_count':counts}).to_csv(O/'temporal_per_family.tsv.gz',sep='\t',index=False,compression={'method':'gzip','mtime':0})
(O/'temporal_design.json').write_text(json.dumps({'permutations':499,'seed':20261015,'samples':len(cols),'blocks':len(blocks),'blocks_multiple_trips':sum(len(set(trip[ix]))>1 for ix in blocks),'samples_in_multiple_trip_blocks':sum(len(ix) for ix in blocks if len(set(trip[ix]))>1),'rows':rows},indent=2)+'\n');print('Temporal',rows[0],flush=True)
del d,A
# H3 reconstruct the original paired cohort and duplicate handling.
P=W/'rebuild_20260729/picrust2_current_20260803/inputs';F=W/'06_reference_bias/refbias_clean/results'
md=pd.read_csv(pin(P/'sample_metadata.tsv'),sep='\t',dtype=str);md=md[md.is_control!='True'];duplicates=md.groupby('sample_id').picrust2_col.apply(list).to_dict()
mag=pd.read_csv(pin(F/'mag_ko_relab_by_sample.tsv.gz'),sep='\t',index_col=0);pic=pd.read_csv(pin(P/'ko_metagenome_unstrat.tsv'),sep='\t',index_col=0);samples=sorted(set(mag.columns)&set(duplicates));assert len(samples)==448
pic=pd.DataFrame({s:pic[duplicates[s]].sum(axis=1) for s in samples});mag=mag[samples]
assert np.isfinite(mag.values).all() and np.isfinite(pic.values).all() and (mag.values>=0).all() and (pic.values>=0).all()
pic=pic/pic.sum(axis=0);mag=mag/mag.sum(axis=0)
common=sorted(set(mag.index)&set(pic.index));union=sorted(set(mag.index)|set(pic.index))
threshold=[]
for cut in [0,1e-6,1e-5,1e-4,1e-3]:
 a=mag.reindex(union,fill_value=0).values>cut;b=pic.reindex(union,fill_value=0).values>cut;inter=(a&b).sum(axis=0);un=(a|b).sum(axis=0)
 for j,s in enumerate(samples):threshold.append(dict(threshold=cut,sample=s,mag_kos=int(a[:,j].sum()),picrust_kos=int(b[:,j].sum()),intersection=int(inter[j]),jaccard=inter[j]/un[j] if un[j] else np.nan,recall_of_mag=inter[j]/a[:,j].sum() if a[:,j].sum() else np.nan))
t=pd.DataFrame(threshold);t.to_csv(O/'profile_thresholds_per_sample.tsv',sep='\t',index=False);t.groupby('threshold').median(numeric_only=True).to_csv(O/'profile_thresholds_summary.tsv',sep='\t')
old=pd.read_csv(pin(F/'per_sample_comparison.tsv'),sep='\t').set_index('sample');zero=t[t.threshold==0].set_index('sample');assert (old.loc[samples].n_pic_ko.values==zero.loc[samples].picrust_kos.values).all();assert np.allclose(old.loc[samples].recall_of_mag,zero.loc[samples].recall_of_mag,atol=1e-6)
def normal_ranks(x):
 x=rankdata(x,axis=0);x-=x.mean(axis=0);sd=np.sqrt((x*x).sum(axis=0));assert (sd>0).all();return x/sd
S=normal_ranks(mag.loc[common].values).T@normal_ranks(pic.loc[common].values);pd.DataFrame(S,index=samples,columns=samples).to_csv(O/'profile_pairwise_rank_correlations.tsv.gz',sep='\t',compression={'method':'gzip','mtime':0})
results=[];rng=np.random.default_rng(20261016);md=meta.loc[samples].reset_index()
for label,keys in [('trip_compartment',['trip','compartment']),('trip_compartment_site',['trip','compartment','site'])]:
 blocks=[np.asarray(v) for v in md.groupby(keys).indices.values() if len(v)>1];ix=np.concatenate(blocks);observed=float(np.diag(S)[ix].mean());null=[]
 for _ in range(4999):
  perm=np.arange(len(samples))
  for v in blocks:perm[v]=rng.permutation(v)
  null.append(float(S[ix,perm[ix]].mean()))
 null=np.array(null);results.append(dict(null=label,samples=len(samples),exchangeable_samples=len(ix),blocks=len(blocks),common_kos=len(common),observed_mean=observed,null_mean=float(null.mean()),advantage=observed-float(null.mean()),null_low=float(np.quantile(null,.025)),null_high=float(np.quantile(null,.975)),p_upper=(1+int((null>=observed).sum()))/5000))
 pd.DataFrame({'mean_correlation':null}).to_csv(O/f'profile_null_{label}.tsv',sep='\t',index=False)
 pd.DataFrame({'sample':[samples[i] for i in ix],'matched_rho':S[ix,ix]}).to_csv(O/f'profile_matched_{label}.tsv',sep='\t',index=False)
(O/'profile_summary.json').write_text(json.dumps({'seed':20261016,'permutations':4999,'fixed_common_universe':len(common),'results':results},indent=2)+'\n');print('Profiles',results,flush=True)
(O/'profile_input_sha256.json').write_text(json.dumps(pins,indent=2)+'\n')
