import pandas as pd, numpy as np, gzip
m=pd.read_csv('selection_covariates/matched_families.tsv',sep='\t')
d=pd.read_csv(gzip.open('inputs/06_darkmatter/darkfam_basis_clean/results/selection/dnds.tsv.gz'),sep='\t')
print(d.columns.tolist()); print(d.head(2).T)
cov=pd.read_csv('selection_covariates/family_covariates.tsv',sep='\t')
print(cov.groupby('set').size())
print('missing cov aa', cov.median_cds_aa.isna().sum())
# per-ds stratum medians, matched
for lab,g in m.groupby('ds'):
    a=g[g.set=='dark'].median_dNdS; b=g[g.set=='ctrl'].median_dNdS
    print('matched ds',lab,len(a),round(a.median(),3),round(b.median(),3))
x=d.merge(cov,on=['family_id','set'])
x['ds']=pd.cut(x.mean_dS,[0,.1,.25,.5,1,2.000001],right=False,labels=False)
for lab,g in x.groupby('ds'):
    print('all ds',lab,(g.set=='dark').sum(),(g.set=='ctrl').sum())
# share of dark retained by stratum
dk=x[x.set=='dark']; md=m[m.set=='dark']
print('dark retained frac', len(md)/len(dk))
print(md.groupby('ds').size()/dk.groupby('ds').size())
# within-stratum (full exact-strata) difference: weighted by dark (ATT) 
x['aa']=pd.cut(x.median_cds_aa,[0,60,100,200,400,np.inf],right=False,labels=False)
x['seq']=pd.cut(x.n_seq,[0,5,10,20,np.inf],right=False,labels=False)
x['gc']=np.floor(x.median_gc/.05)
res=[]
for k,g in x.groupby(['ds','aa','seq','gc']):
    a=g[g.set=='dark'].median_dNdS;b=g[g.set=='ctrl'].median_dNdS
    if len(a) and len(b): res.append((k,len(a),len(b),a.median()-b.median(), a.mean()-b.mean()))
r=pd.DataFrame(res,columns=['k','na','nb','dmed','dmean'])
print('ATT weighted mean diff', np.average(r.dmean,weights=r.na), 'ATC weighted', np.average(r.dmean,weights=r.nb), 'matched-n weighted', np.average(r.dmean,weights=np.minimum(r.na,r.nb)))
print('strata with dark>ctrl median, weighted share', np.average(r.dmed>0,weights=np.minimum(r.na,r.nb)))
# gc distribution
print(x.groupby('set').median_gc.describe())
# within matched, dS detail
print(m.groupby('set').mean_dS.describe())
# Ds fine strata in matched
m['dsf']=pd.cut(m.mean_dS,[0,.05,.1,.25,.5,.75,1,1.5,2.01])
print(m.groupby(['dsf','set']).median_dNdS.agg(['size','median']).unstack())
