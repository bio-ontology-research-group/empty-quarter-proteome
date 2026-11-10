from common import *
import pandas as pd,numpy as np
from scipy.stats import t,wilcoxon,false_discovery_control
R=N/'results';O=R/'ecology';families=list(pd.read_csv(R/'genome_counts.tsv',sep='\t',nrows=1).columns[1:]);d=pd.read_csv(O/'field_joined_cells.tsv',sep='\t',dtype={'site':str});rows=[]
for exclude52 in [False,True]:
 for visitFE in [False,True]:
  dd=d.dropna(subset=['temperature','hour']).copy()
  if exclude52:dd=dd[dd.site!='52']
  dd['temp10']=dd.temperature/10;nuisance=['site','compartment']+(['trip']if visitFE else [])
  X=pd.concat([pd.DataFrame({'intercept':1.},index=dd.index),dd[['temp10','hour']],pd.get_dummies(dd[nuisance],drop_first=True,dtype=float)],axis=1);xx=X.to_numpy(float);n,k=xx.shape;assert np.linalg.matrix_rank(xx)==k;bread=np.linalg.inv(xx.T@xx);sites=dd.site.unique();G=len(sites)
  for fam in families:
   y=dd[fam].to_numpy();beta=np.linalg.lstsq(xx,y,rcond=None)[0];u=y-xx@beta;meat=np.zeros((k,k))
   for ss in sites:
    ii=dd.site.to_numpy()==ss;score=xx[ii].T@u[ii];meat+=np.outer(score,score)
   cov=bread@meat@bread*(G/(G-1))*((n-1)/(n-k));se=np.sqrt(max(0,cov[1,1]));delta=t.ppf(.975,G-1)*se
   rows.append({'family':fam,'exclude_site52':exclude52,'visit_fixed_effect':visitFE,'cells':n,'sites':G,'coefficient_per10C':beta[1],'lo':beta[1]-delta,'hi':beta[1]+delta,'p':2*t.sf(abs(beta[1]/se),G-1)})
m=pd.DataFrame(rows);m['q']=false_discovery_control(m.p);m.to_csv(O/'HS4_temperature_only_sensitivity.tsv',sep='\t',index=False)
s=pd.read_csv(O/'site_differences.tsv',sep='\t',dtype={'site':str});s=s[(s.site!='52')&s.contrast.isin(['trip2-trip3','trip4-trip3'])];rows=[]
for keys,g in s.groupby(['contrast','compartment','family']):
 v=g.difference.to_numpy();rows.append(dict(zip(['contrast','compartment','family'],keys),sites=len(v),mean=v.mean(),p=wilcoxon(v).pvalue if len(v)>=4 and np.any(v)else (1 if len(v)>=4 else np.nan)))
r=pd.DataFrame(rows);r['q']=np.nan;sel=r.p.notna();r.loc[sel,'q']=false_discovery_control(r.loc[sel,'p']);r.to_csv(O/'HS1_exclude_site52.tsv',sep='\t',index=False)
print(m.sort_values('p').head(8).to_string(index=False));print('DONE sensitivity')
