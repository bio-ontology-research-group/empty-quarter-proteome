from pathlib import Path
import pandas as pd,json,hashlib,csv,shutil
R=Path(__file__).resolve().parent;H=R.parents[2]/'reproducibility/darkmatter';S=H/'source_snapshot';out={}
paths=['04_structure/pfam_candidates_clean/results/per_protein_pfam.tsv','06_darkmatter/specificity_clean/results/s14_persistence_by_length.tsv','06_darkmatter/specificity_clean/results/content_vs_depth.tsv']
for name in paths:
 p=S/name;dest=R/'reviewer_inputs'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
for arm,file in [('dark','rep'),('control','control')]:
 d=pd.read_csv(H/f'inputs/04_structure/00_structure_clean/05_validate/per_{file}_structures.tsv.gz',sep='\t');a=d[d.confidence_qualified.eq(1)];out[arm+'_mean_plddt']=float(a.mean_plddt.mean())
d=pd.read_csv(R/'reviewer_inputs/04_structure/pfam_candidates_clean/results/per_protein_pfam.tsv',sep='\t');out['pfam']={}
for name,a in [('low_qtm_reps',d[d['set'].eq('K_members')&d.is_family_rep.eq(1)]),('high_qtm_reps',d[d['set'].eq('QnotK_reps')])]:out['pfam'][name]={'n':len(a),'hit':int(a.n_dom_resolved.gt(0).sum())}
c=pd.read_csv(H/'round5_extensions/results/culture_caller/per_gene.tsv.gz',sep='\t',dtype={'partial':str});c=c[c.dark.eq(1)&c.overlap_threshold.eq(.5)];gain=c[c.bakta.eq('no_major_CDS_overlap')&~c.pyrodigal.eq('no_major_CDS_overlap')];out['culture_gains']={'n':len(gain),'complete':int(gain.partial.eq('00').sum()),'median_aa':float(gain.aa.median()),'top_pair':gain.groupby('rep').size().sort_values(ascending=False).head(3).to_dict()}
pd0=pd.read_csv(R/'reviewer_inputs/06_darkmatter/contamination_clean/results/control_profiles.tsv',sep='\t');ct=pd.read_csv(R/'reviewer_inputs/06_darkmatter/contamination_clean/results/control_contig_hits.tsv',sep='\t');pd0=pd0[pd0.control.eq('SNTC')];ct=ct[ct.control.eq('SNTC')]
rep=pd.read_csv(H/'inputs/06_darkmatter/lineage_darkness_clean/results/v2_reps.tsv.gz',sep='\t');k=next(x for x in ['genome','rep','bin'] if x in ct);ids=set(pd0.genome)|set(ct[k]);a=rep[rep.genome.isin(ids)];out['flagged']={'ids':sorted(ids),'n':len(a),'proteins':int(a.n_proteins.sum()),'dark':int(a.n_dark.sum()),'total_proteins':int(rep.n_proteins.sum()),'total_dark':int(rep.n_dark.sum()),'median_dark_all':float(rep.dark_fraction.median()),'median_dark_excluding_flagged':float(rep[~rep.genome.isin(ids)].dark_fraction.median())}
(R/'root_review_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
with (R/'reviewer_inputs.sha256.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['local_path','ibex_path','sha256'])
 for p in sorted((R/'reviewer_inputs').rglob('*')):
  if p.is_file():w.writerow([p.relative_to(R.parents[2]),'/ibex/scratch/projects/c2014/rob/metagenomics/'+str(p.relative_to(R/'reviewer_inputs')),hashlib.sha256(p.read_bytes()).hexdigest()])
