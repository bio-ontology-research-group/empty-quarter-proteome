#!/usr/bin/env python3
"""Summarise exact nucleotide alignments and overlapping Bakta CDS loci.
Phase classification uses the first aligned base and is not a full codon alignment.
"""
import csv,json,sys,hashlib,importlib.util
from pathlib import Path
from collections import defaultdict,Counter
W=Path('/ibex/scratch/projects/c2014/rob/metagenomics');P=W/'06_darkmatter/cultured_comparison_clean';R=W/'06_darkmatter/p4_review_20261013/results';out=R/'culture_locus';out.mkdir(exist_ok=True)
def rows(p):return list(csv.DictReader(p.open(),delimiter='\t'))
spec=importlib.util.spec_from_file_location('original',P/'scripts/step3_06_pair_anatomy.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
genes={r['protein']:r for r in rows(P/'tmp/pair_anatomy/mag_proteins.tsv')};calls=rows(R/'culture_nucleotide_cigar/per_gene.tsv');group=defaultdict(list)
for r in calls:group[r['rep']].append(r)
records=[];threshold=Counter();pins={}
for rep,rs in group.items():
 culture=rs[0]['cultured'];base=Path(m.bakta_dir(culture)[0]);gffs=list(base.glob('*.gff3'));assert len(gffs)==1,(culture,base,gffs)
 gff=gffs[0];pins[str(gff)]=hashlib.sha256(gff.read_bytes()).hexdigest();cds=defaultdict(list)
 for line in gff.open():
  if line.startswith('##FASTA'):break
  if line.startswith('#'):continue
  c=line.rstrip().split('\t')
  if len(c)==9 and c[2]=='CDS':cds[c[0]].append((int(c[3])-1,int(c[4]),c[6],c[8]))
 hits=defaultdict(list);thresholdhits=defaultdict(set)
 for line in (R/f'culture_nucleotide_cigar/{rep}.paf').open():
  c=line.rstrip().split('\t');cov=(int(c[3])-int(c[2]))/int(c[1]);ident=int(c[9])/int(c[10])
  if cov<.9:continue
  for cut in [.85,.9,.95]:
   if ident>=cut:thresholdhits[cut].add(c[0])
  if ident>=.95:hits[c[0]].append(c)
 for r in rs:
  if r['eligible']!='1':continue
  pid=r['protein']
  for cut in [.85,.9,.95]:threshold[(r['dark'],r['protein_shared'],cut,'n')]+=1;threshold[(r['dark'],r['protein_shared'],cut,'hit')]+=int(pid in thresholdhits[cut])
  if r['protein_shared']=='1' or pid not in hits:continue
  cats=set();overlaps=0
  for c in hits[pid]:
   qlen=int(c[1]);qpos=int(c[2]) if c[4]=='+' else int(c[3])-1;tpos=int(c[7]);end=int(c[8]);orientation=1 if c[4]=='+' else -1
   gs=genes[pid];strand=int(gs.get('strand',1));targetstrand='+' if strand*orientation==1 else '-'
   for a,b,s,attrs in cds[c[5]]:
    ov=max(0,min(end,b)-max(tpos,a))
    if ov<.5*(end-tpos):continue
    overlaps+=1
    phaseq=qpos if strand==1 else qlen-1-qpos;phaset=tpos-a if s=='+' else b-1-tpos
    cats.add('same_orientation_start_phase' if s==targetstrand and phaseq%3==phaset%3 else 'other_frame_or_orientation')
  cat='same_orientation_start_phase' if 'same_orientation_start_phase' in cats else 'other_frame_or_orientation' if cats else 'no_major_CDS_overlap'
  records.append(dict(rep=rep,cultured=culture,protein=pid,dark=int(r['dark']),classification=cat,overlapping_cds_hits=overlaps))
with (out/'per_gene.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter='\t');w.writeheader();w.writerows(records)
summary={'classifications':{str(d):dict(Counter(r['classification'] for r in records if r['dark']==d)) for d in [0,1]},'thresholds':[dict(dark=int(d),protein_shared=int(sh),identity=cut,n=threshold[(d,sh,cut,'n')],supported=threshold[(d,sh,cut,'hit')]) for d in ['0','1'] for sh in ['0','1'] for cut in [.85,.9,.95]],'caveat':'First-aligned-base phase, >=50% of target aligned span overlapping CDS; not full codon correspondence. Any qualifying alignment can support classification.'}
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');(out/'input_sha256.json').write_text(json.dumps(pins,indent=2)+'\n')
