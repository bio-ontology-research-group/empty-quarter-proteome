"""Independent compact-output join/count audit; preserves original faulty tables."""
import csv,json,statistics
from pathlib import Path
B=Path(__file__).resolve().parent.parent;O=Path(__file__).resolve().parent
summary={}
for tr,sub in [('f3','f3_rhodopsin'),('f4','f4_ectc'),('f5','f5_carbonic_anhydrase')]:
 d=B/'f345/results/pilot_v1'/sub
 rows=list(csv.DictReader(open(d/f'{tr}_candidates_summary.tsv'),delimiter='\t'))
 ids={r['protein_id'] for r in rows};assert len(ids)==len(rows)
 hits={}
 for line in open(d/f'{tr}_diamond_reciprocal.tsv'):
  p=line.rstrip().split('\t');qid=p[0].split('|')[0];assert qid in ids
  assert qid not in hits
  hits[qid]={'raw_query_id':p[0],'reference_id':p[1].split('|')[0],'pident':float(p[2]),'alignment_length':int(p[3]),'evalue':float(p[4]),'bitscore':float(p[5])}
 out=[]
 for r in rows:
  h=hits.get(r['protein_id']);out.append({'protein_id':r['protein_id'],'class':r['assigned_class'],'raw_hit_present':h is not None,'reference_id':h['reference_id'] if h else '', 'pident':h['pident'] if h else '', 'original_reported_reference':r['top_ref_id']})
 with open(O/f'{tr}_pilot_reference_join_audit.tsv','w') as f:
  w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
 summary[tr]={'candidates':len(rows),'carrier_genomes':len({r['mag'] for r in rows}),'complete_genes':sum(r['is_complete_gene']=='1' for r in rows),'raw_reference_hits':len(hits),'false_NO_HIT_rows':sum(r['original_reported_reference']=='NO_HIT' and r['raw_hit_present'] for r in out),'identity_min_median_max':[min(h['pident'] for h in hits.values()),statistics.median(h['pident'] for h in hits.values()),max(h['pident'] for h in hits.values())] if hits else [],'note':'One-way closest curated hit, alignment coverage not retained; no novelty inference.'}
(O/'f345_pilot_reference_audit.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
