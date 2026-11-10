from common import *
from scipy import stats
p=read(A/'results/competitive_archaeal_profiles.tsv.gz');q=read(A/'inputs/archaea.tsv').set_index('genome');thermo=set(q.index[q.gtdb_classification.str.contains('p__Thermoplasmatota;')]);external=set(p.genome)-set(q.index);groups={'EQ_Thermoplasmatota':thermo,'Thermoplasmatota_including_external':thermo|external,'other_EQ_archaea':set(q.index)-thermo};design=read(A/'inputs/paired_mapping_samples.tsv');order=design.groupby('site').transect_km.mean().sort_values().index;blocks={v:i for i,g in enumerate(np.array_split(order,10)) for v in g};rng=np.random.default_rng(20261108);out=[]
for mode,metric in [('all_primary','relative_abundance'),('mapq20','mean_coverage')]:
 for threshold in [0.1,0.2,0.5]:
  a=p[p['mode']==mode].copy();a.loc[a.breadth<threshold,metric]=0
  for group,genomes in groups.items():
   y=a[a.genome.isin(genomes)].groupby(['sample','site','trip','compartment'])[metric].sum().reset_index();v=y.groupby(['site','trip','compartment'])[metric].mean().unstack('compartment');d=(v.Deep-v.Surface).dropna();z=d.groupby('site').mean();bs=[z[z.index.map(blocks)==b].to_numpy() for b in sorted(set(z.index.map(blocks)))];boot=[np.mean(np.concatenate([bs[i] for i in rng.integers(0,len(bs),len(bs))])) for _ in range(2000)]
   out.append(dict(mode=mode,metric=metric,threshold=threshold,group=group,sites=len(z),visits=len(d),mean=z.mean(),median=z.median(),positive_sites=int((z>0).sum()),ci_low=np.quantile(boot,.025),ci_high=np.quantile(boot,.975),p=stats.wilcoxon(z).pvalue if (z!=0).any() else 1))
save(pd.DataFrame(out),'mapping_group_review.tsv')
c=read(A/'results/ani95_components.tsv');c['EQ']=c.genome.isin(q.index);summary=c.groupby('ani95_component').agg(n=('genome','size'),eq=('EQ','sum'));summary['external']=summary['n']-summary['eq'];summary.reset_index().to_csv(R/'results/ani_component_review.tsv',sep='\t',index=False)
print(pd.DataFrame(out).to_string(index=False),flush=True)
