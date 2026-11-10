from common import *
rng=np.random.default_rng(20261106);q=read(A/'inputs/thermoplasmatota.tsv').set_index('genome');q['order']=q.gtdb_classification.str.split(';').str[3];q['genus']=q.gtdb_classification.str.split(';').str[5]
b=read(A/'results/masked_context_benchmark.tsv').join(q[['order','genus']],on='genome');b['lift_top5']=b.top5-b.within_contig_null_top5
g=b.groupby(['genome','order','genus'])[['top1','top5','within_contig_null_top5','lift_top5']].mean().reset_index();save(g,'masked_benchmark_per_genome.tsv');rows=[]
for omitted in ['']+sorted(g.order.unique()):
 x=g if omitted=='' else g[g.order!=omitted];v=x.lift_top5.to_numpy();boot=np.mean(v[rng.integers(0,len(v),(2000,len(v)))],axis=1)
 rows.append(dict(omitted_order=omitted,genomes=len(x),mean_lift=v.mean(),ci_low=np.quantile(boot,.025),ci_high=np.quantile(boot,.975),interpretation='fixed_context_heuristic_genome_bootstrap_not_heldout_learning_or_function_accuracy'))
save(pd.DataFrame(rows),'masked_benchmark_order_sensitivity.tsv')
c=read(A/'results/gap_context_candidates.tsv.gz');genes=read(A/'inputs/archaeal_genes.tsv').set_index('gene');assert genes.index.is_unique
c=c.join(q[['order','genus','checkm_completeness','contamination']],on='genome');c=c.join(genes[['length','internal_stop','dark_strict','seq_cluster_rep','cat_cluster_rep']],on='candidate_gene')
assert c.length.notna().all();c['candidate_complete']=~c.partial.astype(bool);c['quality_pass']=(c.checkm_completeness>=70)&(c.contamination<=5)
c['review_priority']=(c.candidate_complete & (c.internal_stop==0) & c.quality_pass & (c.context_score>0))
c['evidence_tier']=np.where(c.subthreshold_sequence_hit & (c.context_score>0),'subthreshold_sequence_plus_context',np.where(c.context_score>0,'context_only','subthreshold_sequence_only'))
c['function_status']='unvalidated';c['distinct_pathway_count']=c.groupby(['genome','missing_reaction','candidate_gene']).pathway.transform('nunique')
collapsed=c.sort_values(['review_priority','subthreshold_sequence_hit','context_score'],ascending=False).drop_duplicates(['genome','missing_reaction','candidate_gene']);save(collapsed,'candidate_reaction_gene_review.tsv.gz')
# Sequence-cluster IDs identify recurrence only at this existing clustering resolution, not biochemical equivalence.
recurrent=collapsed[collapsed.review_priority].groupby(['seq_cluster_rep','missing_reaction'],dropna=False).agg(genomes=('genome','nunique'),genera=('genus','nunique'),orders=('order','nunique'),context_max=('context_score','max'),dark_strict=('dark_strict','max')).reset_index();save(recurrent,'candidate_recurrence_review.tsv')
models=read(A/'results/metabolic_model_summary.tsv');save(models,'metabolic_model_inventory.tsv')
(R/'results/candidate_audit_summary.json').write_text(json.dumps({'input_candidate_rows':len(c),'unique_genome_reaction_gene':len(collapsed),'priority_rows':int(collapsed.review_priority.sum()),'priority_dark_rows':int((collapsed.review_priority & collapsed.dark).sum()),'masked_genomes':len(g),'new_function_assignments':0,'notes':'No habitat screen. Context recovery is pathway-conditioned, not calibrated probability. Candidate recurrence is not independent of phylogeny.'},indent=2));print('Candidate review complete',flush=True)
