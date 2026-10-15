"""Harmonize gene-completeness strata across environmental comparison panels."""
import csv,hashlib,json,os,statistics
from pathlib import Path
import numpy as np
W=Path('/ibex/scratch/projects/c2014/rob/metagenomics');D=W/'06_darkmatter/family_discovery_20261015_clean';C=D/'root_audit/f1_comparison';R=D/'root_audit/f1_complete_panel';R.mkdir(exist_ok=True);H=W/'06_darkmatter/p7_analyses_clean/h3_references/results'
rows=list(csv.DictReader(open(C/'metadata.tsv'),delimiter='\t'));meta={r['id']:r for r in rows};closure={}
for src,fn in [('SMAG','smag_protein_coordinate_mapping.tsv'),('GEM','gem_harmonized_protein_coordinate_mapping_seed20261014.tsv')]:
 wanted={r['original_id']:r for r in rows if r['source']==src};found=set()
 with open(H/fn)as f:
  for r in csv.DictReader(f,delimiter='\t'):
   pid=r['canonical_seq_id']
   if pid not in wanted:continue
   assert pid not in found;found.add(pid);m=wanted[pid]
   assert r['canonical_seq_sha256']==m['sequence_sha256'],(src,pid,'sequence hash mismatch')
   assert int(r['length_aa'])==int(m['length'])
   m['partial']=r['partial']
 assert found==set(wanted),(src,'missingcoordinatejoin',len(set(wanted)-found))
 closure[src]={'candidate_rows':len(wanted),'exact_coordinate_hash_joins':len(found),'complete':sum(r['partial']=='00'for r in wanted.values())}
keep=[r for r in rows if r['source']=='curated_reference'or r['partial']=='00'];ids={r['id']for r in keep}
with open(R/'metadata_complete.tsv','w')as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(keep)
with open(R/'completeness_all.tsv','w')as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
z=np.load(D/'root_audit/f1_embeddings/embeddings.npz',allow_pickle=True);ix=np.array([i for i,uid in enumerate(z['ids'])if uid in ids]);assert len(ix)==len(ids)
np.savez_compressed(R/'embeddings_complete.npz',**{k:z[k][ix]for k in z.files})
best={}
for l in open(C/'alignments.tsv'):
 q,t,pid,qs,qe,ts,te,ev,bits,ql,tl=l.split()
 if q not in ids or t not in ids or q==t:continue
 qc=(int(qe)-int(qs)+1)/int(ql);tc=(int(te)-int(ts)+1)/int(tl)
 if qc<.7 or tc<.7:continue
 key=(q,meta[t]['source']);rec={'query':q,'target_source':key[1],'target':t,'identity':float(pid),'query_coverage':qc,'target_coverage':tc,'bitscore':float(bits)}
 if key not in best or rec['bitscore']>best[key]['bitscore']:best[key]=rec
with open(R/'closest_hits_complete.tsv','w')as f:
 w=csv.DictWriter(f,fieldnames=list(next(iter(best.values()))),delimiter='\t');w.writeheader();w.writerows(best.values())
eq=[r for r in keep if r['source']=='EQ'];comparisons={}
for source in ['curated_reference','SMAG','GEM']:
 vals=[best[(r['id'],source)]['identity']for r in eq if (r['id'],source)in best]
 comparisons[source]={'queries':len(eq),'with_hit':len(vals),'without_coverage_passing_hit':len(eq)-len(vals),'identity_median':statistics.median(vals)if vals else None,'below50amonghits':sum(v<50for v in vals)}
summary={'job_id':os.environ['SLURM_JOB_ID'],'coordinate_closure':closure,'complete_panel_sources':{s:sum(r['source']==s for r in keep)for s in ['EQ','SMAG','GEM','curated_reference']},'comparison':comparisons,'limits':'curated references exempt from caller partial flags; no newGPUrun, subsetofverifiedembeddings; profilescreenmissesknownSPLcontrol; samplesnotentireenvironmentalcatalogues'}
(R/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
