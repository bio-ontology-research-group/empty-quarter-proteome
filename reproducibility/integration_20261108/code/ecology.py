from common import *
from scipy import sparse,stats
from statsmodels.stats.multitest import multipletests
import sys
mode=sys.argv[1];G=read(R/f'results/genome_panels_{mode}.tsv').set_index('genome');mi=read(M/'mag_index.tsv').sort_values('mag_index');G=G.reindex(mi.mag_id);assert G.index.is_unique
q=read(A/'inputs/genome_quality.tsv').set_index('genome').reindex(mi.mag_id);s=read(E/'results/sample_metrics.tsv').reset_index(drop=True);AA=100*sparse.load_npz(M/'sample_by_mag_coverm_relative_abundance.npz').toarray()[s.sample_index.to_numpy()]
assert len(s)==546 and np.all(AA.sum(1)>0)
matched=set(read(E/'results/rund_20261106/matched_profiles.tsv').sample_id);s['matched_amplicon']=s.sample_id.isin(matched)
fields=read(E/'inputs/rund-20261106/site_dune_zones.tsv').set_index('site');s['edgell_dune_type']=s.site.map(fields.edgell_dune_type)
coords=read(E/'inputs/site_coordinates.tsv').set_index('site');s['x_km']=s.site.map(coords.x_km);s['y_km']=s.site.map(coords.y_km)
siteorder=s.groupby('site').transect_km.mean().sort_values().index;block={v:i for i,g in enumerate(np.array_split(siteorder,10)) for v in g}
rng=np.random.default_rng(20261106);records=[];dec=[];visitout=[];coverage=[]
def interval(v,sites):
 groups=[v[np.asarray([block[x]==b for x in sites])] for b in sorted({block[x] for x in sites})]
 means=np.array([np.concatenate([groups[i] for i in rng.integers(0,len(groups),len(groups))]).mean() for _ in range(1000)])
 return np.quantile(means,[.025,.975])
def paired(values,mask):
 df=pd.DataFrame(values[mask]);df[['site','trip','compartment']]=s.loc[mask,['site','trip','compartment']].to_numpy();v=df.groupby(['site','trip','compartment']).mean();return v
for qual in ['all','completeness70_contamination5']:
 eligible=np.ones(len(q),bool) if qual=='all' else ((q.checkm_completeness>=70)&(q.contamination<=5)).to_numpy()
 for col in G:
  t=G[col].to_numpy(float);known=np.isfinite(t)&eligible;trait=np.nan_to_num(t)*eligible;mass=AA[:,known].sum(1)
  # Unknown profile evidence is not scored as absence; denominators retained explicitly.
  whole=AA@trait;mapped=np.divide(100*whole,mass,out=np.full(len(s),np.nan),where=mass>0)
  if not known.any():continue
  for cohort in ['all','matched_amplicon']:
   mask=np.ones(len(s),bool) if cohort=='all' else s.matched_amplicon.to_numpy()
   for denom,score in [('whole_read_bound',whole),('known_catalogue_fraction',mapped)]:
    v=paired(score[:,None],mask)
    for c1,c0 in [('Surface','Deep'),('RootAdjacent','Deep')]:
     if c1 not in v.index.get_level_values(2):continue
     left=v.xs(c1,level=2);right=v.xs(c0,level=2);visits=left.index.intersection(right.index);d=(left.loc[visits,0]-right.loc[visits,0]).dropna();z=d.groupby(level='site').mean()
     if len(z)<5:continue
     lo,hi=interval(z.to_numpy(),z.index)
     p=stats.wilcoxon(z).pvalue if np.any(z!=0) else 1.
     records.append(dict(panel=col,quality=qual,cohort=cohort,denominator=denom,contrast=c1+'-minus-'+c0,omitted_trip='',sites=len(z),visits=len(d),mean_pp=z.mean(),median_pp=z.median(),ci_low=lo,ci_high=hi,p=p))
     for trip in sorted(s.trip.unique()):
      zz=d[d.index.get_level_values('trip')!=trip].groupby(level='site').mean()
      records.append(dict(panel=col,quality=qual,cohort=cohort,denominator=denom,contrast=c1+'-minus-'+c0,omitted_trip=trip,sites=len(zz),visits=len(d[d.index.get_level_values('trip')!=trip]),mean_pp=zz.mean(),median_pp=zz.median(),ci_low=np.nan,ci_high=np.nan,p=np.nan))
  if qual=='all':
   for i,r in s.iterrows():coverage.append(dict(sample_id=r.sample_id,panel=col,whole_read_bound=whole[i],known_catalogue_fraction=mapped[i],known_genome_abundance=mass[i],mapped_abundance=AA[i].sum(),unknown_genomes=int((~known).sum())))
# Exact decomposition for fully observed traits; no missing genes silently set absent.
P=AA/AA.sum(1)[:,None];v=paired(P,np.ones(len(s),bool));left=v.xs('Surface',level=2);right=v.xs('Deep',level=2);visits=left.index.intersection(right.index)
for rank,prefix in [('phylum','p__'),('genus','g__')]:
 labels=[];unclassified=[]
 for genome,tax in zip(q.index,q.gtdb_classification):
  value=next((x for x in tax.split(';') if x.startswith(prefix)),'')
  unknown=value in ['',prefix];labels.append(genome if unknown else value);unclassified.append(unknown)
 groups=[np.where(np.array(labels)==g)[0] for g in sorted(set(labels))]
 for col in G:
  trait=G[col].to_numpy(float)
  if not np.isfinite(trait).all():continue
  rows=[]
  for site,trip in visits:
   a=left.loc[(site,trip)].to_numpy(float);b=right.loc[(site,trip)].to_numpy(float);B=W=X=0.
   for ix in groups:
    ma,mb=a[ix].sum(),b[ix].sum()
    if ma>0 and mb>0:
     ta,tb=a[ix]@trait[ix]/ma,b[ix]@trait[ix]/mb
     B+=(ma-mb)*(ta+tb)/2;W+=(ta-tb)*(ma+mb)/2
    else:X+=a[ix]@trait[ix]-b[ix]@trait[ix]
   total=(a-b)@trait;assert np.isclose(B+W+X,total,atol=1e-10)
   row=dict(panel=col,rank=rank,site=site,trip=trip,total_pp=100*total,between_pp=100*B,within_pp=100*W,exclusive_pp=100*X,unclassified_mass=(a[np.array(unclassified)].sum()+b[np.array(unclassified)].sum())/2);rows.append(row)
  z=pd.DataFrame(rows);visitout+=rows
  zz=z.groupby('site')[['total_pp','between_pp','within_pp','exclusive_pp','unclassified_mass']].mean()
  for metric in zz:
   lo,hi=interval(zz[metric].to_numpy(),zz.index);dec.append(dict(panel=col,rank=rank,metric=metric,sites=len(zz),mean=zz[metric].mean(),ci_low=lo,ci_high=hi))
out=pd.DataFrame(records)
for _,ix in out[out.omitted_trip==''].groupby(['quality','cohort','denominator']).groups.items():out.loc[ix,'q_bh']=multipletests(out.loc[ix,'p'],method='fdr_bh')[1]
save(out,f'paired_guild_tests_{mode}.tsv');save(pd.DataFrame(coverage),f'sample_guild_scores_{mode}.tsv.gz');save(pd.DataFrame(dec),f'turnover_decomposition_{mode}.tsv');save(pd.DataFrame(visitout),f'turnover_paired_visits_{mode}.tsv.gz');save(s,f'model_design_{mode}.tsv')
print('Complete paired contrasts and exact turnover decomposition',mode,flush=True)
