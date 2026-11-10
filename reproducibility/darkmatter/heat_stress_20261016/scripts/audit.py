from common import *
import pandas as pd,numpy as np
from scipy.stats import false_discovery_control,wilcoxon
R=N/'results';O=R/'audit';O.mkdir(exist_ok=True);checks=[]
def check(name,condition,detail):
 assert condition,(name,detail)
 checks.append({'check':name,'passed':bool(condition),'detail':str(detail)})
p=pd.read_csv(R/'ecology/paired_tests.tsv',sep='\t');s=pd.read_csv(R/'ecology/site_differences.tsv',sep='\t')
for r in p[p.variant=='primary'].itertuples():
 v=s[(s.contrast==r.contrast)&(s.compartment==r.compartment)&(s.family==r.family)].difference.to_numpy()
 check('paired_count_'+r.contrast+'_'+r.compartment+'_'+r.family,len(v)==r.sites,len(v))
 if len(v):check('paired_effect_'+r.contrast+'_'+r.compartment+'_'+r.family,np.isclose(v.mean(),r.mean),v.mean())
 if np.isfinite(r.p):check('paired_p_'+r.contrast+'_'+r.compartment+'_'+r.family,np.isclose(wilcoxon(v).pvalue if np.any(v)else 1,r.p),r.p)
for variant,d in p.groupby('variant'):
 for label,cons in [('hot',['trip2-trip3','trip4-trip3']),('transition',['trip1-trip3','trip5-trip3'])]:
  dd=d[d.contrast.isin(cons)&d.p.notna()];check('BH_'+variant+'_'+label,np.allclose(false_discovery_control(dd.p),dd.q),len(dd))
for fn in ['HS4_environment_models.tsv','HS4_rainfall_sensitivity.tsv','within_thermoplasmatota_tests.tsv','HS4_depth_interaction.tsv']:
 d=pd.read_csv(R/'ecology'/fn,sep='\t');d=d[d.p.notna()];check('BH_'+fn,np.allclose(false_discovery_control(d.p),d.q),len(d))
x=pd.read_csv(R/'ecology/lineage_decomposition.tsv',sep='\t');check('lineage_additive_identity',np.allclose(x.total,x.between+x.within+x.unshared),len(x))
f=pd.read_csv(R/'ecology/field_observations.tsv',sep='\t');check('humidity_range',f.humidity.dropna().between(0,100).all(),len(f));check('July_temperature_retained',f.query("trip=='trip2'").temperature.notna().any(),f.query("trip=='trip2'").temperature.notna().sum())
# Independently evaluate exact tail under the same within-contig rotation null.
# Each tested family has one locus per genome, hence independent Bernoulli draws across genomes.
d=pd.read_csv(R/'hs2/tested_loci.tsv.gz',sep='\t');rows=[]
for family,g in d.groupby('family'):
 check('HS2_one_locus_per_genome_'+family,g.genome.is_unique,len(g))
 dist=np.array([1.])
 for prob in g.null_near_fraction:dist=np.convolve(dist,[1-prob,prob])
 observed=int(g.near_anchor.sum());rows.append({'family':family,'genomes':len(g),'near_genomes':observed,'expected_near_genomes':g.null_near_fraction.sum(),'p_exact':dist[observed:].sum()})
e=pd.DataFrame(rows);e['q_exact']=false_discovery_control(e.p_exact.clip(0,1));e.sort_values('q_exact').to_csv(O/'HS2_exact_tail.tsv',sep='\t',index=False)
write(O/'checks.tsv',checks)
(O/'summary.json').write_text(json.dumps({'job':os.environ['SLURM_JOB_ID'],'checks':len(checks),'all_passed':True,'exact_HS2_q05':int((e.q_exact<.05).sum()),'scope':'Independent table/statistical identities and exact null calculation; not an external review.'},indent=2));print(e.sort_values('q_exact').head(12).to_string(index=False));print('DONE audit',len(checks))
