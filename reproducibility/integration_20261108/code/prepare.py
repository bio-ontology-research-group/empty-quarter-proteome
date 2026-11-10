from common import *
import hashlib
q=read(A/'inputs/genome_quality.tsv');assert len(q)==5438 and q.genome.is_unique
manifest=[];count={g:0 for g in q.genome};D=R/'proteins';D.mkdir(exist_ok=True)
assignment={g:i%50 for i,g in enumerate(sorted(q.genome))};handles=[open(D/f'shard_{i:02d}.faa','w') for i in range(50)]
source=Path('/project/k10018/empty-quarter/catalog/proteins/orthogroups/all_rep_proteins.faa');h=hashlib.sha256();active=None
with source.open() as f:
 for line in f:
  h.update(line.encode())
  if line.startswith('>'):
   g=line[1:].split('::',1)[0];assert g in assignment,g;active=assignment[g];count[g]+=1
  handles[active].write(line)
for f in handles:f.close()
assert count==q.set_index('genome').n_genes.to_dict(),'Protein count mismatch'
for g,shard in assignment.items():manifest.append(dict(genome=g,shard=shard,proteins=count[g]))
save(pd.DataFrame(manifest),'protein_shards.tsv')
(R/'results/protein_source.json').write_text(json.dumps({'path':str(source),'sha256':h.hexdigest(),'proteins':sum(count.values()),'genomes':len(count)},indent=2))
print('Prepared',len(count),'genomes',sum(count.values()),'proteins',flush=True)
