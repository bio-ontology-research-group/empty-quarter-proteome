"""Exact site-level decomposition; known F1 result is not an independent replication."""
from common import *
import importlib.util,collections,json,numpy as np
from scipy.stats import wilcoxon
R=N/'results/ecology';R.mkdir(parents=True,exist_ok=True)
spec=importlib.util.spec_from_file_location('old',D/'ecology/f1_ecology_root_v3.py');old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
meta=old.load_metadata(W/'metadata/sample_metadata.tsv');traits=old.load_rep_traits(W/'06_metabolism/uv_survival_clean/results/rep_traits.tsv');mags=sorted(traits['all_reps']);ix={g:i for i,g in enumerate(mags)}
f4=tab(N/'results/f4_integrity/loci.tsv');cpf={r['mag']for r in tab(D/'root_audit/f1_full/candidates.tsv')if r['partial']=='00'}
sets={'CPF':cpf,'K01669':traits['k01669_reps'],'EctABC_strict':{r['mag']for r in f4 if r['source']=='EQ_MAG'and r['strict_ABC']=='1'},'EctC_WEYH':{r['mag']for r in f4 if r['source']=='EQ_MAG'and r['WEYH']=='1'}}
sets['CPF_without_K01669']=cpf-traits['k01669_reps'];sets['EctABC_cov70']={r['mag']for r in f4 if r['source']=='EQ_MAG'and r['strict_ABC_cov70']=='1'}
tax={}
for p in [W/'05_refine/_gtdbtk_r232/gtdbtk.bac120.summary.tsv',W/'05_refine/_gtdbtk_r232/gtdbtk.ar53.summary.tsv']:
    for r in tab(p):tax[r['user_genome'].removesuffix('.fa')]=r['classification']
assert set(mags)<=set(tax)
lineages={}
for rank,rankidx in [('phylum',1),('genus',5)]:
    labels=[]
    for g in mags:
        t=tax[g].split(';');label=';'.join(t[:rankidx+1]);unknown=len(t[rankidx])<=3
        # Unique MAG fallback avoids treating unresolved unrelated genomes as one genus.
        if unknown:label+=';unclassified_MAG='+g
        labels.append(label)
    keys=sorted(set(labels));ii={s:i for i,s in enumerate(keys)};lineages[rank]=(keys,np.array([ii[x]for x in labels]),np.array(['unclassified_MAG='in x for x in keys]))
profiles=W/'rebuild_20260729/picrust2_current_20260803/coverm_current/profiles_full_recompute';ab={};inputs=[];skipped=[]
for s in sorted(meta):
    p=profiles/s/'abund.tsv';a,unmapped=old.read_coverm_profile(p);assert set(a)==set(mags);v=np.array([a[g]for g in mags]);inputs.append({'path':str(p),'sha256':sha(p)})
    if not np.any(v>0):
        assert unmapped==100,(s,unmapped,'undefined catalogue fractions')
        skipped.append(s);continue
    assert np.isfinite(v).all()and np.all(v>=0),(s,'invalid abundance')
    ab[s]=v/v.sum()
assert len(ab)==550,(len(ab),skipped)
cells=collections.defaultdict(list)
for s,v in ab.items():
    m=meta[s];cells[(m['site'],m['trip'],m['compartment'])].append(v)
# Arithmetic means preserve the probability simplex and exact component additivity.
cells={k:np.mean(v,axis=0)for k,v in cells.items()};pairs=[(site,trip)for site,trip,c in cells if c=='Surface'and(site,trip,'Deep')in cells]
assert len(pairs)==69,len(pairs)
visits=[];legacy=[]
for name,carriers in sets.items():
    mask=np.array([g in carriers for g in mags],dtype=float)
    for site,trip in sorted(pairs):
        vs=cells[(site,trip,'Surface')];vd=cells[(site,trip,'Deep')];total=float((vs-vd)@mask)
        for rank,(keys,idx,unknown)in lineages.items():
            ws=np.bincount(idx,weights=vs,minlength=len(keys));wd=np.bincount(idx,weights=vd,minlength=len(keys));ys=np.bincount(idx,weights=vs*mask,minlength=len(keys));yd=np.bincount(idx,weights=vd*mask,minlength=len(keys));common=(ws>0)&(wd>0);ps=np.divide(ys,ws,out=np.zeros_like(ws),where=ws>0);pd=np.divide(yd,wd,out=np.zeros_like(wd),where=wd>0)
            between=float(np.sum((ws[common]-wd[common])*(ps[common]+pd[common])/2));within=float(np.sum((ws[common]+wd[common])*(ps[common]-pd[common])/2));exclusive=float(np.sum((ys-yd)[~common]));assert abs(total-between-within-exclusive)<1e-12
            multi=np.bincount(idx,minlength=len(keys))>1
            visits.append({'trait':name,'rank':rank,'site':site,'trip':trip,'total_pp':100*total,'between_pp':100*between,'within_pp':100*within,'exclusive_pp':100*exclusive,'common_mass':float(((ws+wd)/2)[common].sum()),'multi_MAG_mass':float(((ws+wd)/2)[multi].sum()),'unclassified_mass':float(((ws+wd)/2)[unknown].sum())})
    # Original replicate-median estimand, used for the new F4 soil contrast and F1 reconstruction.
    cellmetrics=collections.defaultdict(list)
    for s,v in ab.items():
        m=meta[s];cellmetrics[(m['site'],m['trip'],m['compartment'])].append(float(v@mask)*100)
    bysite=collections.defaultdict(list)
    for site,trip in pairs:bysite[site].append(np.median(cellmetrics[(site,trip,'Surface')])-np.median(cellmetrics[(site,trip,'Deep')]))
    vals=np.array([np.median(bysite[s])for s in sorted(bysite)]);rng=np.random.default_rng(20261015);boot=np.median(rng.choice(vals,(2000,len(vals))),axis=1);p=wilcoxon(vals).pvalue if np.any(vals)else 1
    legacy.append({'trait':name,'carriers':len(carriers),'sites':len(vals),'visits':len(pairs),'median_difference_pp':float(np.median(vals)),'ci_low':float(np.quantile(boot,.025)),'ci_high':float(np.quantile(boot,.975)),'p_raw':float(p)})
    write(R/(name+'_median_site_effects.tsv'),[{'site':s,'difference_pp':float(np.median(bysite[s]))}for s in sorted(bysite)])
summary=[];siterows=[]
for name in sets:
    for rank in lineages:
        groups=collections.defaultdict(list)
        for r in visits:
            if r['trait']==name and r['rank']==rank:groups[r['site']].append(r)
        metrics=['total_pp','between_pp','within_pp','exclusive_pp','common_mass','multi_MAG_mass','unclassified_mass']
        rr=[{'trait':name,'rank':rank,'site':s,**{k:float(np.mean([r[k]for r in v]))for k in metrics}}for s,v in sorted(groups.items())];siterows+=rr
        for r in rr:assert abs(r['total_pp']-r['between_pp']-r['within_pp']-r['exclusive_pp'])<1e-9
        for k in metrics:
            v=np.array([r[k]for r in rr]);rng=np.random.default_rng(20261015);b=rng.choice(v,(2000,len(v)))
            summary.append({'trait':name,'rank':rank,'metric':k,'sites':len(v),'mean':float(v.mean()),'mean_ci_low':float(np.quantile(b.mean(axis=1),.025)),'mean_ci_high':float(np.quantile(b.mean(axis=1),.975)),'median':float(np.median(v)),'median_ci_low':float(np.quantile(np.median(b,axis=1),.025)),'median_ci_high':float(np.quantile(np.median(b,axis=1),.975))})
# Conservative family correction across all six displayed contrasts, including known F1 tests.
adj=old.holm_bonferroni([r['p_raw']for r in legacy])
for r,p in zip(legacy,adj):r['p_Holm_six_contrasts']=p
write(R/'paired_visit_decomposition.tsv',visits);write(R/'site_mean_decomposition.tsv',siterows);write(R/'decomposition_summary.tsv',summary);write(R/'median_ecology_tests.tsv',legacy);write(R/'profile_inputs.tsv',inputs)
write(R/'taxonomy.tsv',[{'mag':g,'taxonomy':tax[g]}for g in mags])
(R/'summary.json').write_text(json.dumps({'job':os.environ['SLURM_JOB_ID'],'samples':len(ab),'zero_mapped':skipped,'visits':len(pairs),'sites':len({s for s,t in pairs}),'traits':{k:len(v)for k,v in sets.items()},'decomposition':'replicate means, exact B+W+R, site means over visits, equally weighted sites','ecology_test':'separate original replicate medians, site median paired differences','unknown_genus':'unique MAG fallback; unclassified bins never pooled','limitations':['genome traits fixed at species-representative resolution','within-genus term reflects species composition, not strain adaptation','bootstrap intervals descriptive, follow-up exploratory']},indent=2))
