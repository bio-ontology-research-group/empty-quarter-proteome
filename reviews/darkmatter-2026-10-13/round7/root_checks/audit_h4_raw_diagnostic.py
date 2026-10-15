"""Independent full-row audit of retained Foldseek evidence; execute on Slurm."""
import csv,collections,hashlib,json,math,re,time
from pathlib import Path
base=Path('/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h4_structure');r=base/'results'
manifest=list(csv.DictReader((r/'v7_full5851_cropped_regions_manifest.tsv').open(),delimiter='\t'))
whole={x['model_id']:int(x['whole_chain_length']) for x in manifest}
crops={x['cropped_pdb'][:-4]:int(x['reg_len']) for x in manifest if x['cropped_pdb']!='NONE'}
assert len(whole)==5851 and len(crops)==7773
output={'status':'DIAGNOSTIC_NOT_ACCEPTANCE','checked_utc_unix':time.time(),'files':{}}
for name,queries in [('v7_wc_pdb_hits.raw.tsv',whole),('v7_crop_pdb_hits.raw.tsv',crops),('v7_wc_rerun_hits.raw.tsv',whole),('v7_crop_rerun_hits.raw.tsv',crops)]:
 p=r/name
 if not p.exists():
  assert 'rerun' in name
  continue
 counts=collections.Counter();h=hashlib.sha256();rows=0;anomalies=[];anomaly_count=0
 with p.open('rb') as f:
  for raw in f:
   h.update(raw);parts=raw.decode().rstrip('\r\n').split('\t');assert len(parts)==21,(name,rows)
   q=parts[0].removesuffix('.pdb');assert q in queries,q
   numeric=[float(v) for v in parts[2:15]];assert all(math.isfinite(v) for v in numeric)
   alntm,qtm,ttm=map(float,parts[2:5])
   if not all(0<=v<=1.00001 for v in [alntm,qtm,ttm]):
    anomaly_count+=1
    if len(anomalies)<20:anomalies.append({'line':rows+1,'columns':parts})
   assert 0<=float(parts[7])<=10.00001 and 0<=float(parts[9])<=1.00001
   qlen,tlen,alnlen,qs,qe,ts,te=map(int,parts[13:20]);assert qlen==queries[q],(q,qlen,queries[q])
   assert 1<=qs<=qe<=qlen and 1<=ts<=te<=tlen
   cigar=parts[20];ops=re.findall(r'(\d+)([MID])',cigar)
   assert ''.join(n+op for n,op in ops)==cigar
   sums=collections.Counter()
   for n,op in ops:sums[op]+=int(n)
   assert sum(sums.values())==alnlen
   assert sums['M']==qe-qs+1+te-ts+1-alnlen
   counts[q]+=1;rows+=1
 output['files'][name]={'rows':rows,'query_ids_with_hits':len(counts),'query_census':len(queries),'max_rows_per_query':max(counts.values(),default=0),'sha256':h.hexdigest(),'bytes':p.stat().st_size,'malformed_or_unknown_rows':0,'tm_score_range_anomalies':anomaly_count,'tm_anomaly_examples':anomalies}
(r/'v7_raw_diagnostic.json').write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps(output,indent=2))
