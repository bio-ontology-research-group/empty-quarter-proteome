"""Exploratory spatial spline sensitivity, not a replacement/reproduction of Rund's INLA model."""
from common import *
import statsmodels.api as sm
from patsy import dmatrix
from statsmodels.stats.multitest import multipletests
import sys
mode=sys.argv[1];d=read(R/f'results/model_design_{mode}.tsv');scores=read(R/f'results/sample_guild_scores_{mode}.tsv.gz');d['log_reads']=np.log1p(d.post_qc_reads);out=[];skipped=[]
for cohort in ['all','matched_amplicon']:
 sub=d if cohort=='all' else d[d.matched_amplicon]
 sub=sub.dropna(subset=['log_reads','x_km','y_km','edgell_dune_type']).copy()
 for spatial in ['design_only','spatial_spline','spatial_spline_and_dune']:
  form="1 + C(compartment, Treatment(reference='Deep')) + C(trip) + log_reads"
  if spatial!='design_only':form+=' + bs(x_km,df=4) + bs(y_km,df=4)'
  if spatial=='spatial_spline_and_dune':form+=' + C(edgell_dune_type)'
  X=dmatrix(form,sub,return_type='dataframe');assert np.linalg.matrix_rank(X)==X.shape[1],(cohort,spatial)
  for panel,g in scores.groupby('panel'):
   for endpoint in ['whole_read_bound','known_catalogue_fraction']:
    y=sub.sample_id.map(g.set_index('sample_id')[endpoint]);ok=y.notna().to_numpy()
    if sum(ok)<30:continue
    xx=X.loc[ok]
    if np.linalg.matrix_rank(xx)<xx.shape[1]:skipped.append(dict(panel=panel,cohort=cohort,model=spatial,endpoint=endpoint,reason='rank_deficient_after_unknown_exclusion'));continue
    fit=sm.OLS(y.loc[ok].to_numpy(),xx).fit(cov_type='cluster',cov_kwds={'groups':sub.loc[ok,'site']},use_t=True)
    for term in xx:
     if 'C(compartment' not in term:continue
     ci=fit.conf_int().loc[term];out.append(dict(panel=panel,cohort=cohort,model=spatial,endpoint=endpoint,term=term,n=int(ok.sum()),sites=sub.loc[ok,'site'].nunique(),coefficient_pp=fit.params[term],ci_low=ci.iloc[0],ci_high=ci.iloc[1],p=fit.pvalues[term],status='exploratory_linear_mean_spline_sensitivity_not_INLA'))
o=pd.DataFrame(out)
for _,ix in o.groupby(['cohort','model','endpoint']).groups.items():o.loc[ix,'q_bh']=multipletests(o.loc[ix,'p'].fillna(1),method='fdr_bh')[1]
save(o,f'spatial_sensitivity_{mode}.tsv');save(pd.DataFrame(skipped),f'spatial_skipped_{mode}.tsv');print('Spatial spline sensitivity complete',mode,flush=True)
