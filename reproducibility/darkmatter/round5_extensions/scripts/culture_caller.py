#!/usr/bin/env python3
"""Controlled caller substitution on 22 culture assemblies; fixed nucleotide hits."""
from pathlib import Path
import csv,json,hashlib,subprocess,importlib.util,bisect
from collections import defaultdict
import numpy as np,pandas as pd
W=Path('/ibex/scratch/projects/c2014/rob/metagenomics');R=W/'06_darkmatter/p4_review_20261013/results';P=W/'06_darkmatter/cultured_comparison_clean';O=R/'round5_extensions/culture_caller';O.mkdir(parents=True,exist_ok=True)
pins={}
def pin(p):
 p=Path(p);h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 pins[str(p)]=h.hexdigest();return p
def rows(p):return list(csv.DictReader(pin(p).open(),delimiter='\t'))
script=pin(P/'scripts/step3_06_pair_anatomy.py');spec=importlib.util.spec_from_file_location('old',script);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
culture={r['genome_id']:Path(r['path']) for r in rows(P/'results/cultured_meta.tsv')};genes={r['protein']:r for r in rows(P/'tmp/pair_anatomy/mag_proteins.tsv')};calls=rows(R/'culture_nucleotide_cigar/per_gene.tsv');group=defaultdict(list)
for r in calls:
 if r['eligible']=='1' and r['dna_supported']=='1':group[r['rep']].append(r)
exe=Path('/ibex/scratch/projects/c2014/rob/miniconda3/envs/pyrodigal/bin/pyrodigal');pin(exe);version=subprocess.check_output([str(exe),'--version'],text=True,stderr=subprocess.STDOUT).strip();details=[]
def index_gff(p):
 cds=defaultdict(list)
 for line in pin(p).open():
  if line.startswith('##FASTA'):break
  if line.startswith('#'):continue
  c=line.rstrip().split('\t')
  if len(c)==9 and c[2]=='CDS':cds[c[0]].append((int(c[3])-1,int(c[4]),c[6]))
 cds={k:sorted(v) for k,v in cds.items()};starts={k:[x[0] for x in v] for k,v in cds.items()};maxlen={k:max(b-a for a,b,s in v) for k,v in cds.items()};return cds,starts,maxlen
for rep,rs in group.items():
 name=rs[0]['cultured'];assembly=pin(culture[name]);gff=O/f'{name}.gff'
 cmd=[str(exe),'-p','meta','-i',str(assembly),'-o',str(gff),'-f','gff','-j','4','--no-stop-codon']
 with (O/f'{name}.stderr').open('w') as f:subprocess.run(cmd,stdout=f,stderr=f,check=True)
 bakta=list(Path(m.bakta_dir(name)[0]).glob('*.gff3'));assert len(bakta)==1
 indexes={'bakta':index_gff(bakta[0]),'pyrodigal':index_gff(gff)};hits=defaultdict(list)
 for line in pin(R/f'culture_nucleotide_cigar/{rep}.paf').open():
  c=line.rstrip().split('\t')
  if (int(c[3])-int(c[2]))/int(c[1])>=.9 and int(c[9])/int(c[10])>=.95:hits[c[0]].append(c)
 for r in rs:
  pid=r['protein'];assert hits[pid]
  for threshold in [.5,.9]:
   item=dict(rep=rep,cultured=name,protein=pid,dark=int(r['dark']),partial=r['partial'],nt=int(r['nt']),aa=int(genes[pid]['aa']),overlap_threshold=threshold)
   for caller,(cds,starts,maxlen) in indexes.items():
    cats=set()
    for c in hits[pid]:
     qlen=int(c[1]);qpos=int(c[2]) if c[4]=='+' else int(c[3])-1;tpos=int(c[7]);end=int(c[8]);orientation=1 if c[4]=='+' else -1;strand=int(genes[pid].get('strand',1));targetstrand='+' if strand*orientation==1 else '-';target=c[5]
     left=bisect.bisect_left(starts.get(target,[]),tpos-maxlen.get(target,0));right=bisect.bisect_left(starts.get(target,[]),end)
     for a,b,s in cds.get(target,[])[left:right]:
      if max(0,min(end,b)-max(tpos,a))<threshold*(end-tpos):continue
      phaseq=qpos if strand==1 else qlen-1-qpos;phaset=tpos-a if s=='+' else b-1-tpos
      cats.add('same_orientation_start_phase' if s==targetstrand and phaseq%3==phaset%3 else 'other_frame_or_orientation')
    item[caller]='same_orientation_start_phase' if 'same_orientation_start_phase' in cats else 'other_frame_or_orientation' if cats else 'no_major_CDS_overlap'
   details.append(item)
 print(rep,len(rs),'DNA-supported calls',flush=True)
d=pd.DataFrame(details);old=pd.read_csv(pin(R/'culture_locus_all/per_gene.tsv'),sep='\t');j=d[d.overlap_threshold==.5].merge(old,on=['rep','protein'],validate='one_to_one');assert len(j)==68537 and (j.bakta==j.classification).all()
f=pd.read_csv(pin(R/'culture_family_membership/per_gene.tsv'),sep='\t')[['rep','protein','admitted_family']];d=d.merge(f,on=['rep','protein'],how='left',validate='many_to_one');d.to_csv(O/'per_gene.tsv.gz',sep='\t',index=False,compression={'method':'gzip','mtime':0})
rows_out=[]
for threshold in [.5,.9]:
 for label,mask in [('all',np.ones(len(d),bool)),('complete_ge100aa',(d.partial.astype(str).eq('00')&(d.aa>=100))),('admitted_family',d.admitted_family.notna())]:
  for dark in [0,1]:
   a=d[(d.overlap_threshold==threshold)&mask&(d.dark==dark)]
   if not len(a):continue
   no=a.bakta.eq('no_major_CDS_overlap');rescued=no&~a.pyrodigal.eq('no_major_CDS_overlap');same=a.pyrodigal.eq('same_orientation_start_phase')
   rows_out.append(dict(threshold=threshold,subset=label,dark=dark,n=len(a),pairs=a.rep.nunique(),bakta_no_cds=int(no.sum()),pyrodigal_no_cds=int(a.pyrodigal.eq('no_major_CDS_overlap').sum()),rescued_any=int(rescued.sum()),rescued_same_phase=int((no&same).sum()),pyrodigal_same_phase=int(same.sum())))
pd.DataFrame(rows_out).to_csv(O/'summary.tsv',sep='\t',index=False)
pair=[]
for (rep,dark),a in d[d.overlap_threshold==.5].groupby(['rep','dark']):pair.append(dict(rep=rep,dark=int(dark),n=len(a),bakta_no_cds=float(a.bakta.eq('no_major_CDS_overlap').mean()),pyrodigal_no_cds=float(a.pyrodigal.eq('no_major_CDS_overlap').mean())))
pd.DataFrame(pair).to_csv(O/'per_pair.tsv',sep='\t',index=False)
(O/'run.json').write_text(json.dumps({'pyrodigal_version':version,'pairs':len(group),'configuration':'-p meta -f gff -j 4 --no-stop-codon','conditioning':'fixed >=90%coverage >=95%identity DNA alignments; any qualifying alignment; phase at first aligned base only'},indent=2)+'\n');(O/'input_sha256.json').write_text(json.dumps(pins,indent=2)+'\n');print(rows_out,flush=True)
