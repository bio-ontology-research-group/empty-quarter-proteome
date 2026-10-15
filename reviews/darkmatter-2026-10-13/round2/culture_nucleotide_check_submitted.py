#!/usr/bin/env python3
"""Test nucleotide support for protein calls in 22 culture/MAG comparisons.

Run on IBEX under Slurm. Excludes the confirmed E. coli contaminant. Extracts
catalogue gene intervals from original MAG assemblies and aligns to the paired
culture-derived assembly using minimap2 short-read mode with short seeds.
A single local alignment must cover >=90% of the query at >=95% identity.
This is a presence sensitivity, not a definitive absence assay.
"""
import csv,hashlib,json,subprocess,sys,importlib.util
from pathlib import Path
from collections import Counter
W=Path('/ibex/scratch/projects/c2014/rob/metagenomics');P=W/'06_darkmatter/cultured_comparison_clean';out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True)
def rows(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def fasta(p):
 d={};key=None
 for l in p.open():
  if l.startswith('>'):key=l[1:].split()[0];d[key]=[]
  else:d[key].append(l.strip())
 return {k:''.join(v) for k,v in d.items()}
inputs={}
def pin(p):inputs[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
script=P/'scripts/step3_06_pair_anatomy.py';pin(script)
spec=importlib.util.spec_from_file_location('original_anatomy',script);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
pairfile=P/'results/pair_anatomy/pair_anatomy_per_pair.tsv';metafile=P/'results/cultured_meta.tsv';gene_file=P/'tmp/pair_anatomy/mag_proteins.tsv';darkfile=P/'tmp/pair_anatomy/mag_protein_dark.tsv'
for p in [pairfile,metafile,gene_file,darkfile]:pin(p)
pairs=[r for r in rows(pairfile) if r['mag_rep']!='V42PRr2__single_SemiBin_33'];culture={r['genome_id']:Path(r['path']) for r in rows(metafile)}
genes=rows(gene_file);dark={r['protein']:r['dark']=='1' for r in rows(darkfile)}
exe='/home/hohndor/bin/minimap2';version=subprocess.check_output([exe,'--version'],text=True).strip();details=[]
for pair in pairs:
 rep=pair['mag_rep'];assembly=W/f'05_refine/_drep/out/dereplicated_genomes/{rep}.fa';target=culture[pair['cultured']]
 pin(assembly);pin(target);seq=fasta(assembly);seq.update({k.split('::')[-1]:v for k,v in list(seq.items())})
 subset=[r for r in genes if r['rep']==rep];query=out/f'{rep}.fna';lengths={}
 with query.open('w') as f:
  for r in subset:
   contig=r['contig'];nt=seq.get(contig,seq.get(contig.split('::')[-1]));assert nt is not None,contig
   a,b=int(r['start']),int(r['stop']);assert 1<=a<=b<=len(nt),(rep,a,b,len(nt))
   nt=nt[a-1:b];lengths[r['protein']]=len(nt)
   if len(nt)>=90:f.write('>'+r['protein']+'\n'+nt+'\n')
 paf=out/f'{rep}.paf';cmd=[exe,'-x','sr','-k','11','-w','5','-n','2','-m','20','-t','4','--secondary=yes','-N','20',str(target),str(query)]
 with paf.open('w') as f,(out/f'{rep}.stderr').open('w') as e:subprocess.run(cmd,stdout=f,stderr=e,check=True)
 supported=set();best={}
 for line in paf.open():
  c=line.rstrip().split('\t');coverage=(int(c[3])-int(c[2]))/int(c[1]);identity=int(c[9])/int(c[10]);pid=c[0]
  if coverage>=.9 and identity>=.95:supported.add(pid)
  if coverage*identity>best.get(pid,(0,0))[0]*best.get(pid,(0,0))[1]:best[pid]=(coverage,identity)
 fwd=P/f'tmp/pair_anatomy/dia/{rep}__m2c.m8';rev=P/f'tmp/pair_anatomy/dia/{rep}__c2m.m8';pin(fwd);pin(rev)
 shared=m.rbh(m.best_hits(str(fwd)),m.best_hits(str(rev)),50.,50.)
 for r in subset:
  pid=r['protein'];eligible=lengths[pid]>=90;cov,ident=best.get(pid,(0,0))
  details.append(dict(rep=rep,cultured=pair['cultured'],protein=pid,nt=lengths[pid],partial=r['partial'],dark=int(dark[pid]),protein_shared=int(pid in shared),eligible=int(eligible),dna_supported=int(pid in supported),best_coverage=cov,best_identity=ident))
 print(rep,len(subset),len(supported),flush=True)
with (out/'per_gene.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(details[0]),delimiter='\t');w.writeheader();w.writerows(details)
summary=[]
for d in [0,1]:
 for sh in [0,1]:
  a=[r for r in details if r['eligible'] and r['dark']==d and r['protein_shared']==sh]
  summary.append(dict(dark=d,protein_shared=sh,eligible=len(a),dna_supported=sum(r['dna_supported'] for r in a)))
(out/'summary.json').write_text(json.dumps({'pairs':len(pairs),'minimap2_version':version,'min_query_nt':90,'min_query_coverage':.9,'min_identity':.95,'groups':summary},indent=2)+'\n')
(out/'input_sha256.json').write_text(json.dumps(inputs,indent=2)+'\n')
