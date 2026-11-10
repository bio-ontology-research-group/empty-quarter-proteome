from common import *
import numpy as np,pandas as pd
from scipy.stats import wilcoxon,t,spearmanr
R=N/'results';O=R/'ecology';O.mkdir(exist_ok=True);rng=np.random.default_rng(20261016)
def bh(p):
 p=np.array(p,float);order=np.argsort(p);v=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];out=np.empty(len(p));out[order]=np.minimum(1,v);return out
def test(v):
 v=np.array(v,float);n=len(v)
 if not n:return {'sites':0,'mean':np.nan,'median':np.nan,'lo':np.nan,'hi':np.nan,'p':np.nan,'positive':0}
 boot=rng.choice(v,(5000,n)).mean(1)
 return {'sites':n,'mean':v.mean(),'median':np.median(v),'lo':np.quantile(boot,.025),'hi':np.quantile(boot,.975),'p':wilcoxon(v).pvalue if n>=4 and np.any(v)else (1. if n>=4 else np.nan),'positive':int((v>0).sum())}
z=np.load(R/'abundance.npz');g=z['genomes'];samples=z['samples'];raw=z['raw'];keep=raw.sum(1)>0;raw=raw[keep];samples=samples[keep];a=raw/raw.sum(1)[:,None]
meta=pd.read_csv(R/'metadata.tsv',sep='\t',index_col=0,dtype={'site':str}).loc[samples];idx=pd.read_csv(R/'genome_index.tsv',sep='\t',index_col=0).loc[g];C=pd.read_csv(R/'genome_counts.tsv',sep='\t',index_col=0).loc[g];S=pd.read_csv(R/'genome_counts_strict.tsv',sep='\t',index_col=0).loc[g];size=pd.read_csv(R/'genome_sizes.tsv',sep='\t',index_col=0).loc[g]
comp=idx.completeness.to_numpy()/100
variants={'primary':C.to_numpy(float),'strict':S.to_numpy(float),'completeness_adjusted':C.to_numpy(float)/comp[:,None],'per_Mbp':C.to_numpy(float)/(size.length.to_numpy()[:,None]/1e6)}
rows=[];cellrows=[];diffs=[];decomp=[];tax=idx.gtdb_classification.fillna('').tolist();phy=np.array([v.split(';')[1]for v in tax]);genus=np.array([';'.join(v.split(';')[:6])if len(v.split(';')[5])>3 else v+';MAG='+gg for v,gg in zip(tax,g)])
cells={}
for key,d in meta.groupby(['site','trip','compartment']):
 ii=meta.index.get_indexer(d.index);cells[key]=a[ii].mean(0)
for variant,c in variants.items():
 metrics=pd.DataFrame(a@c,index=samples,columns=C.columns);metrics.index.name='sample';metrics.to_csv(O/(variant+'_sample_metrics.tsv'),sep='\t')
 cm={k:v@c for k,v in cells.items()}
 for k,v in cm.items():cellrows.append({'variant':variant,'site':k[0],'trip':k[1],'compartment':k[2],**dict(zip(C.columns,v))})
 for hot in ['trip2','trip4','trip1','trip5']:
  for compartment in ['Surface','Deep','PlantRoot']:
   sites=sorted({k[0]for k in cm if k[1:]==(hot,compartment)and(k[0],'trip3',compartment)in cm})
   for j,fam in enumerate(C.columns):
    ds=[cm[(s,hot,compartment)][j]-cm[(s,'trip3',compartment)][j]for s in sites]
    rows.append({'variant':variant,'contrast':hot+'-trip3','compartment':compartment,'family':fam,**test(ds)})
    if variant=='primary':
     for s,d in zip(sites,ds):diffs.append({'site':s,'contrast':hot+'-trip3','compartment':compartment,'family':fam,'difference':d})
# Correction is separate per prespecified primary vs sensitivity variant; primary hot comparisons48tests.
res=pd.DataFrame(rows);res['q']=np.nan
for variant in variants:
 for branch,contrasts in [('hot',['trip2-trip3','trip4-trip3']),('transitional',['trip1-trip3','trip5-trip3'])]:
  mask=(res.variant==variant)&res.contrast.isin(contrasts)&res.p.notna();res.loc[mask,'q']=bh(res.loc[mask,'p'])
res.to_csv(O/'paired_tests.tsv',sep='\t',index=False);pd.DataFrame(cellrows).to_csv(O/'cell_metrics.tsv',sep='\t',index=False);pd.DataFrame(diffs).to_csv(O/'site_differences.tsv',sep='\t',index=False)
# Exact additive decomposition of mean-copy change into between/within shared lineages and unshared mass.
for rank,labels in [('phylum',phy),('genus',genus)]:
 _,gi=np.unique(labels,return_inverse=True);ng=gi.max()+1
 for hot in ['trip2','trip4']:
  for compartment in ['Surface','Deep','PlantRoot']:
   sites=sorted({k[0]for k in cells if k[1:]==(hot,compartment)and(k[0],'trip3',compartment)in cells})
   for s in sites:
    ah=cells[(s,hot,compartment)];aw=cells[(s,'trip3',compartment)];wh=np.bincount(gi,weights=ah,minlength=ng);ww=np.bincount(gi,weights=aw,minlength=ng);both=(wh>0)&(ww>0)
    for j,fam in enumerate(C.columns):
     c=C.iloc[:,j].to_numpy();yh=np.bincount(gi,weights=ah*c,minlength=ng);yw=np.bincount(gi,weights=aw*c,minlength=ng);ph=np.divide(yh,wh,out=np.zeros(ng),where=wh>0);pw=np.divide(yw,ww,out=np.zeros(ng),where=ww>0)
     b=((wh-ww)*(ph+pw)/2)[both].sum();w=((wh+ww)*(ph-pw)/2)[both].sum();x=(yh-yw)[~both].sum();total=(ah-aw)@c;assert abs(total-b-w-x)<1e-10
     decomp.append({'site':s,'contrast':hot+'-trip3','compartment':compartment,'family':fam,'rank':rank,'total':total,'between':b,'within':w,'unshared':x})
pd.DataFrame(decomp).to_csv(O/'lineage_decomposition.tsv',sep='\t',index=False)
# Within-Thermoplasmatota conditional mean-copy composition and its community mass.
thermo=phy=='p__Thermoplasmatota';tr=[]
for k,v in cells.items():
 mass=v[thermo].sum()
 if mass>0:tr.append({'site':k[0],'trip':k[1],'compartment':k[2],'thermo_fraction':mass,**dict(zip(C.columns,(v[thermo]@C.to_numpy()[thermo])/mass))})
pd.DataFrame(tr).to_csv(O/'thermoplasmatota_cells.tsv',sep='\t',index=False)
# HS4 matched depth-by-visit difference-in-differences.
did=[]
cm={k:v@C.to_numpy()for k,v in cells.items()}
for hot in ['trip2','trip4']:
 sites=sorted({k[0]for k in cm if all((k[0],tr,co)in cm for tr in [hot,'trip3']for co in ['Surface','Deep'])})
 for j,fam in enumerate(C.columns):
  ds=[(cm[(s,hot,'Surface')][j]-cm[(s,hot,'Deep')][j])-(cm[(s,'trip3','Surface')][j]-cm[(s,'trip3','Deep')][j])for s in sites]
  did.append({'contrast':hot+'-trip3','family':fam,**test(ds)})
did=pd.DataFrame(did);did['q']=np.nan;mask=did.p.notna();did.loc[mask,'q']=bh(did.loc[mask,'p']);did.to_csv(O/'HS4_depth_interaction.tsv',sep='\t',index=False)
# Retain all field observations, fail on duplicate visit/site. They are collection-time environmental measurements, not soil-temperature time series.
import re
field=[];field_sources=[]
for num,year in [(1,2023),(2,2023),(3,2024),(4,2024),(5,2025)]:
 p=W/f'temporal_environment_20261014/inputs/trip{num}-{year}.tsv'
 with open(p)as fh: rawrows=list(csv.reader(fh,delimiter='\t'))
 header=rawrows[0];body=[r for r in rawrows[1:] if any(v.strip() for v in r)]
 if num==2:
  assert len(header)==9 and all(len(r)==10 and r[6]=='' for r in body)
  header=header[:6]+['pressure']+header[6:]
 assert all(len(r)==len(header)for r in body),(num,'ragged field sheet')
 d=pd.DataFrame(body,columns=header);field_sources.append({'path':str(p),'sha256':sha(p),'parser_note':'trip2 omitted pressure header restored as blank column; raw source unchanged'if num==2 else 'exact header'})
 for _,r in d.iterrows():
  if not str(r['site']).isdigit():continue
  def numval(x):
   try:return float(x)
   except:return np.nan
  tm=str(r.get('time',''));hour=np.nan
  if re.match(r'^\d{1,2}:\d{2}$',tm):hh,mm=tm.split(':');hour=int(hh)+int(mm)/60
  field.append({'site':str(int(r['site'])),'trip':f'trip{num}','field_date':r['date'],'field_coordinates':r.get('coordinates',''),'temperature':numval(r.get('temperature')),'humidity':numval(r.get('humidity')),'hour':hour})
field=pd.DataFrame(field);assert not field.duplicated(['site','trip']).any();field.to_csv(O/'field_observations_raw.tsv',sep='\t',index=False)
invalid=field[(~field.humidity.between(0,100)&field.humidity.notna())|(~field.temperature.between(-20,70)&field.temperature.notna())];invalid.to_csv(O/'field_invalid_values.tsv',sep='\t',index=False)
field.loc[~field.humidity.between(0,100),'humidity']=np.nan
field.loc[~field.temperature.between(-20,70),'temperature']=np.nan
field.to_csv(O/'field_observations.tsv',sep='\t',index=False);write(O/'field_sources.tsv',field_sources)
cmf=pd.DataFrame([{'site':k[0],'trip':k[1],'compartment':k[2],**dict(zip(C.columns,v))}for k,v in cm.items()]).merge(field,on=['site','trip'],how='left',validate='many_to_one');cmf.to_csv(O/'field_joined_cells.tsv',sep='\t',index=False)
# Fixed-site models with and without fixed visit. Cluster covariance across sites. Target coefficients are per10C and per10 humidity percentage points.
models=[]
for visitFE in [False,True]:
 d=cmf.dropna(subset=['temperature','humidity','hour']).copy();d['temp10']=d.temperature/10;d['humidity10']=d.humidity/10
 nuisance=['site','compartment']+(['trip']if visitFE else [])
 X=pd.concat([pd.DataFrame({'intercept':1.},index=d.index),d[['temp10','humidity10','hour']],pd.get_dummies(d[nuisance],drop_first=True,dtype=float)],axis=1)
 xx=X.to_numpy(float);n,k=xx.shape;rank=np.linalg.matrix_rank(xx);assert rank==k,(visitFE,rank,k)
 bread=np.linalg.inv(xx.T@xx);sites=d.site.unique();G=len(sites)
 for fam in C:
  yy=d[fam].to_numpy();beta=np.linalg.lstsq(xx,yy,rcond=None)[0];u=yy-xx@beta;meat=np.zeros((k,k))
  for s in sites:
   sel=d.site.to_numpy()==s;score=xx[sel].T@u[sel];meat+=np.outer(score,score)
  cov=bread@meat@bread*(G/(G-1))*((n-1)/(n-k));se=np.sqrt(np.maximum(0,np.diag(cov)))
  for term in ['temp10','humidity10']:
   j=list(X.columns).index(term);pv=2*t.sf(abs(beta[j]/se[j]),G-1)if se[j]>0 else np.nan;delta=t.ppf(.975,G-1)*se[j]
   models.append({'family':fam,'visit_fixed_effect':visitFE,'term':term,'cells':n,'sites':G,'rank':rank,'columns':k,'condition_number':np.linalg.cond(xx),'coefficient':beta[j],'lo':beta[j]-delta,'hi':beta[j]+delta,'p':pv})
models=pd.DataFrame(models);models['q']=bh(models.p);models.to_csv(O/'HS4_environment_models.tsv',sep='\t',index=False)
cor=cmf[['temperature','humidity','hour']].corr();cor.to_csv(O/'environment_correlations.tsv',sep='\t')
(R/'ecology_summary.json').write_text(json.dumps({'job':os.environ['SLURM_JOB_ID'],'samples':len(samples),'cells':len(cells),'thermoplasmatota_genomes':int(thermo.sum()),'field_complete_cells':len(d),'family_definitions':'eight domain families; Cpn60 includes bacterial and archaeal chaperonins','primary_tests':48,'inference':'paired site associations conditional on visits; not independent temporal replication or expression'},indent=2))
print('DONE ecology',res.query("variant=='primary' and contrast=='trip4-trip3'")[['family','compartment','sites','mean','p','q']].to_string(index=False))
