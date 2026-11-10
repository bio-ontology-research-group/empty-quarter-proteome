from common import *
import numpy as np,pandas as pd,re
from scipy.stats import wilcoxon,t
R=N/'results';O=R/'ecology';rng=np.random.default_rng(20261017)
def bh(p):
 p=np.array(p);o=np.argsort(p);a=np.empty(len(p));a[o]=np.minimum(1,np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return a
# Within-focal-phylum associations: a separate exploratory family of tests.
d=pd.read_csv(O/'thermoplasmatota_cells.tsv',sep='\t',dtype={'site':str});families=list(pd.read_csv(R/'genome_counts.tsv',sep='\t',nrows=1).columns[1:]);tests=[]
cells=d.set_index(['site','trip','compartment'])
for hot in ['trip2','trip4']:
 for co in ['Surface','Deep','PlantRoot']:
  sites=sorted({s for s,tr,cc in cells.index if tr==hot and cc==co and(s,'trip3',co)in cells.index})
  for fam in families+['thermo_fraction']:
   vals=np.array([cells.loc[(s,hot,co),fam]-cells.loc[(s,'trip3',co),fam]for s in sites]);n=len(vals)
   if not n:continue
   boot=rng.choice(vals,(5000,n)).mean(1);pv=wilcoxon(vals).pvalue if n>=4 and np.any(vals)else(1 if n>=4 else np.nan)
   tests.append({'family':fam,'contrast':hot+'-trip3','compartment':co,'sites':n,'mean':vals.mean(),'lo':np.quantile(boot,.025),'hi':np.quantile(boot,.975),'p':pv})
tests=pd.DataFrame(tests);tests['q']=np.nan;sel=tests.p.notna();tests.loc[sel,'q']=bh(tests.loc[sel,'p']);tests.to_csv(O/'within_thermoplasmatota_tests.tsv',sep='\t',index=False)
# Rainfall sensitivity on the frozen441sampleenvironmentcohort. Missing records remain missing.
env=pd.read_csv(R/'environment.tsv',sep='\t',dtype={'site':str});rr=env.groupby(['site','trip','compartment'],as_index=False).nasa_rain_7d.mean();d=pd.read_csv(O/'field_joined_cells.tsv',sep='\t',dtype={'site':str}).merge(rr,on=['site','trip','compartment'],how='left',validate='one_to_one');d.to_csv(O/'rain_joined_cells.tsv',sep='\t',index=False)
models=[]
for usehumidity in [False,True]:
 for visitFE in [False,True]:
  predictors=['temp10','hour','lograin7']+(['humidity10']if usehumidity else [])
  dd=d.copy();dd['temp10']=dd.temperature/10;dd['humidity10']=dd.humidity/10;dd['lograin7']=np.log1p(dd.nasa_rain_7d);dd=dd.dropna(subset=predictors)
  nuisance=['site','compartment']+(['trip']if visitFE else [])
  X=pd.concat([pd.DataFrame({'intercept':1.},index=dd.index),dd[predictors],pd.get_dummies(dd[nuisance],drop_first=True,dtype=float)],axis=1);xx=X.to_numpy(float);n,k=xx.shape;rank=np.linalg.matrix_rank(xx);assert rank==k
  bread=np.linalg.inv(xx.T@xx);sites=dd.site.unique();G=len(sites)
  for fam in families:
   y=dd[fam].to_numpy();beta=np.linalg.lstsq(xx,y,rcond=None)[0];u=y-xx@beta;meat=np.zeros((k,k))
   for ss in sites:
    ii=dd.site.to_numpy()==ss;score=xx[ii].T@u[ii];meat+=np.outer(score,score)
   cov=bread@meat@bread*(G/(G-1))*((n-1)/(n-k));se=np.sqrt(np.maximum(0,np.diag(cov)))
   for term in ['temp10','lograin7']+(['humidity10']if usehumidity else []):
    j=list(X.columns).index(term);pv=2*t.sf(abs(beta[j]/se[j]),G-1)if se[j]>0 else np.nan;delta=t.ppf(.975,G-1)*se[j]
    models.append({'family':fam,'humidity_adjusted':usehumidity,'visit_fixed_effect':visitFE,'term':term,'cells':n,'sites':G,'rank':rank,'condition_number':np.linalg.cond(xx),'coefficient':beta[j],'lo':beta[j]-delta,'hi':beta[j]+delta,'p':pv})
models=pd.DataFrame(models);models['q']=bh(models.p);models.to_csv(O/'HS4_rainfall_sensitivity.tsv',sep='\t',index=False)
# Check actual field coordinates for repeated-site identity. No coordinate corrections guessed.
f=pd.read_csv(O/'field_observations.tsv',sep='\t',dtype={'site':str});coords=[]
for r in f.itertuples():
 vals=re.findall(r'[-+]?\d+(?:\.\d+)?',str(r.field_coordinates))
 if len(vals)==2:coords.append({'site':r.site,'trip':r.trip,'latitude':float(vals[0]),'longitude':float(vals[1])})
co=pd.DataFrame(coords);audit=[]
for site,dd in co.groupby('site'):
 pts=dd[['latitude','longitude']].to_numpy();maxdist=0
 for i in range(len(pts)):
  for j in range(i):
   p,q=np.deg2rad(pts[i]),np.deg2rad(pts[j]);v=np.sin((p[0]-q[0])/2)**2+np.cos(p[0])*np.cos(q[0])*np.sin((p[1]-q[1])/2)**2;maxdist=max(maxdist,6371*2*np.arcsin(np.sqrt(v)))
 audit.append({'site':site,'visits':len(dd),'maximum_coordinate_distance_km':maxdist})
pd.DataFrame(audit).to_csv(O/'site_coordinate_audit.tsv',sep='\t',index=False)
print(tests.sort_values('p').head(12).to_string(index=False));print('DONE complement')
