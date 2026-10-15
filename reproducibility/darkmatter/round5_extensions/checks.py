#!/usr/bin/env python3
"""Independent aggregation/invariants of the retained exploratory outputs."""
from pathlib import Path
import csv,hashlib,json
import numpy as np,pandas as pd
H=Path(__file__).resolve().parent;R=H/'results';n=0
def check(ok,label):
 global n
 assert ok,label
 n+=1
manifest=H/'outputs.sha256.tsv'
if manifest.exists():
 for r in csv.DictReader(manifest.open(),delimiter='\t'):
  check(hashlib.sha256((H/r['path']).read_bytes()).hexdigest()==r['sha256'],'hash '+r['path'])
d=pd.read_csv(R/'culture_caller/per_gene.tsv.gz',sep='\t',dtype={'partial':str})
check(not d.duplicated(['rep','protein','overlap_threshold']).any(),'unique culture-call keys')
check(len(d)==2*68537,'fixed eligible DNA-supported universe')
a=d[d.overlap_threshold==.5];b=d[d.overlap_threshold==.9]
check(a.rep.nunique()==22,'22 genome pairs')
j=a.merge(b,on=['rep','protein'],validate='one_to_one',suffixes=('_50','_90'))
for caller in ['bakta','pyrodigal']:
 check((~j[caller+'_50'].eq('no_major_CDS_overlap')|j[caller+'_90'].eq('no_major_CDS_overlap')).all(),'overlap threshold monotonicity '+caller)
summary=pd.read_csv(R/'culture_caller/summary.tsv',sep='\t')
for r in summary.itertuples():
 z=d[(d.overlap_threshold==r.threshold)&(d.dark==r.dark)]
 if r.subset=='complete_ge100aa':z=z[z.partial.eq('00')&(z.aa>=100)]
 if r.subset=='admitted_family':z=z[z.admitted_family.notna()]
 old=z.bakta.eq('no_major_CDS_overlap');new=z.pyrodigal.eq('no_major_CDS_overlap');same=z.pyrodigal.eq('same_orientation_start_phase')
 got=[len(z),z.rep.nunique(),old.sum(),new.sum(),(old&~new).sum(),(old&same).sum(),same.sum()]
 check(got==[r.n,r.pairs,r.bakta_no_cds,r.pyrodigal_no_cds,r.rescued_any,r.rescued_same_phase,r.pyrodigal_same_phase],f'culture subgroup {r.threshold} {r.subset} {r.dark}')
for r in pd.read_csv(R/'culture_caller/per_pair.tsv',sep='\t').itertuples():
 z=a[(a.rep==r.rep)&(a.dark==r.dark)]
 check(len(z)==r.n and np.isclose(z.bakta.eq('no_major_CDS_overlap').mean(),r.bakta_no_cds) and np.isclose(z.pyrodigal.eq('no_major_CDS_overlap').mean(),r.pyrodigal_no_cds),'culture pair '+r.rep+str(r.dark))
f=pd.read_csv(R/'temporal_per_family.tsv.gz',sep='\t');s=pd.read_csv(R/'temporal_summary.tsv',sep='\t');draw=pd.read_csv(R/'temporal_null_draws.tsv',sep='\t')
check(len(f)==106042 and f.family_id.is_unique,'temporal family identities')
check(len(draw)==499,'499 global permutations')
check(np.isclose(f.null_all_five_count.sum()/499/len(f),draw['all'].mean()),'temporal family/null closure')
for r in s.itertuples():
 z=f if r.group=='all' else f[f.occupied_samples.between(*map(int,r.group.split('-')))]
 x=draw[r.group].values;v=z.observed_all_five.mean()
 check(len(z)==r.families and np.isclose(v,r.observed),'temporal group observed '+r.group)
 check(np.allclose([x.mean(),np.quantile(x,.025),np.quantile(x,.975),v-x.mean()],[r.null_mean,r.null_low,r.null_high,r.observed_minus_null]),'temporal null summary '+r.group)
 check(np.isclose((1+(x>=v).sum())/500,r.p_upper) and np.isclose((1+(x<=v).sum())/500,r.p_lower),'temporal permutation tails '+r.group)
t=pd.read_csv(R/'profile_thresholds_per_sample.tsv',sep='\t')
check(len(t)==448*5 and not t.duplicated(['sample','threshold']).any(),'profile cohort keys')
for col in ['mag_kos','picrust_kos','intersection']:
 check((np.diff(t.pivot(index='sample',columns='threshold',values=col).values,axis=1)<=0).all(),'threshold monotonicity '+col)
u=t.mag_kos+t.picrust_kos-t.intersection
check(np.allclose(t.jaccard,t.intersection/u,equal_nan=True),'Jaccard identity')
check(np.allclose(t.recall_of_mag,t.intersection/t.mag_kos,equal_nan=True),'recall denominator and missing values')
S=pd.read_csv(R/'profile_pairwise_rank_correlations.tsv.gz',sep='\t',index_col=0)
check(S.shape==(448,448) and list(S.index)==list(S.columns),'pairwise matrix sample identities')
check(np.isfinite(S.values).all() and (np.abs(S.values)<=1+1e-12).all(),'finite rank correlations')
for r in json.loads((R/'profile_summary.json').read_text())['results']:
 label=r['null'];z=pd.read_csv(R/f'profile_matched_{label}.tsv',sep='\t');x=pd.read_csv(R/f'profile_null_{label}.tsv',sep='\t').mean_correlation.values
 check(len(z)==r['exchangeable_samples'] and len(x)==4999,'profile null sample/permutation counts '+label)
 check(np.allclose(z.matched_rho,[S.loc[k,k] for k in z['sample']]),'matched sample diagonals '+label)
 check(np.isclose(z.matched_rho.mean(),r['observed_mean']) and np.isclose(x.mean(),r['null_mean']),'profile mean effects '+label)
 check(np.isclose((1+(x>=r['observed_mean']).sum())/5000,r['p_upper']),'profile permutation tail '+label)
# Bind manuscript-facing claims to fixed reviewed values, independently of summary consistency.
ledger={r['claim_id']:float(r['expected']) for r in csv.DictReader((H/'claims.tsv').open(),delimiter='\t')}
def claim(key,value):check(value==ledger[key],'manuscript claim '+key)
primary=summary[summary.threshold==.5].set_index(['subset','dark'])
for key,subset,dark,col in [('culture.dark.pyrodigal_no_cds','all',1,'pyrodigal_no_cds'),('culture.dark.gained_overlap','all',1,'rescued_any'),('culture.dark.gained_start_phase','all',1,'rescued_same_phase'),('culture.annotated.pyrodigal_no_cds','all',0,'pyrodigal_no_cds'),('culture.dark.complete100.no_cds','complete_ge100aa',1,'pyrodigal_no_cds'),('culture.dark.complete100.n','complete_ge100aa',1,'n'),('culture.admitted.no_cds','admitted_family',1,'pyrodigal_no_cds'),('culture.admitted.n','admitted_family',1,'n')]:claim(key,int(primary.loc[(subset,dark),col]))
z=s[s.group=='all'].iloc[0]
claim('temporal.observed_percent',round(100*z.observed,1));claim('temporal.null_percent',round(100*z.null_mean,1));claim('temporal.delta_pp',round(100*z.observed_minus_null,2))
x=json.loads((R/'profile_summary.json').read_text());claim('profile.common_kos',x['fixed_common_universe']);claim('profile.advantage',round(x['results'][0]['advantage'],4));claim('profile.site.advantage',round(x['results'][1]['advantage'],4));claim('profile.site.exchangeable',x['results'][1]['exchangeable_samples']);claim('profile.threshold1e-3.zero_mag',int(t[t.threshold==.001].mag_kos.eq(0).sum()))
# Reviewer-driven whole-visit null and source-based corrections.
v=json.loads((R/'temporal_visits_summary.json').read_text());x=pd.read_csv(R/'temporal_visits_draws.tsv',sep='\t').all_five_frequency.values
check(len(x)==499 and v['visits']==414 and v['exchangeable_samples']==342,'whole-visit design')
check(np.isclose(x.mean(),v['null_mean']) and np.isclose((1+(x<=v['observed']).sum())/500,v['p_lower']),'whole-visit null summary')
check(round(100*v['null_mean'],1)==80.5 and v['p_lower']==.11,'whole-visit manuscript values')
review=H.parents[2]/'reviews/darkmatter-2026-10-13/round5';support=review/'reviewer_inputs'
for r in csv.DictReader((review/'reviewer_inputs.sha256.tsv').open(),delimiter='\t'):
 check(hashlib.sha256((H.parents[2]/r['local_path']).read_bytes()).hexdigest()==r['sha256'],'reviewer input '+r['local_path'])
pf=pd.read_csv(support/'04_structure/pfam_candidates_clean/results/per_protein_pfam.tsv',sep='\t')
check(int(pf[pf['set'].eq('K_members')&pf.is_family_rep.eq(1)].n_dom_resolved.gt(0).sum())==93,'Pfam low-query representative hits')
check(int(pf[pf['set'].eq('QnotK_reps')].n_dom_resolved.gt(0).sum())==101,'Pfam high-query representative hits')
for name,want in [('rep',.775),('control',.867)]:
 z=pd.read_csv(H.parent/f'inputs/04_structure/00_structure_clean/05_validate/per_{name}_structures.tsv.gz',sep='\t');check(round(z[z.confidence_qualified.eq(1)].mean_plddt.mean(),3)==want,'postgate confidence '+name)
profile=pd.read_csv(support/'06_darkmatter/contamination_clean/results/control_profiles.tsv',sep='\t');contig=pd.read_csv(support/'06_darkmatter/contamination_clean/results/control_contig_hits.tsv',sep='\t');ids=set(profile[profile.control.eq('SNTC')].genome)|set(contig[contig.control.eq('SNTC')].genome)
z=pd.read_csv(H.parent/'inputs/06_darkmatter/lineage_darkness_clean/results/v2_reps.tsv.gz',sep='\t');a=z[z.genome.isin(ids)];check(len(a)==12 and a.n_dark.sum()==1792 and a.n_proteins.sum()==42831,'flagged instance counts')
check(round(100*z.dark_fraction.median(),2)==6.31 and round(100*z[~z.genome.isin(ids)].dark_fraction.median(),2)==6.34,'flagged exclusion medians')
print(f'PASS: {n} round-5 checks of hashes, population keys, transitions, denominators and null summaries.')
