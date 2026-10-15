import sys,glob,os
P='/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/cultured_comparison_clean/tmp/pair_anatomy/mag/'
Q='/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p4_review_20261013/results/culture_nucleotide/'
import csv
strand={};part={}
for r in csv.DictReader(open('/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/cultured_comparison_clean/tmp/pair_anatomy/mag_proteins.tsv'),delimiter='\t'):
    strand[r['protein']]=r['strand'];part[r['protein']]=r['partial']
code={}
b='TCAG';aa='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
i=0
for x in b:
 for y in b:
  for z in b: code[x+y+z]=aa[i];i+=1
def rc(s): return s[::-1].translate(str.maketrans('ACGTN','TGCAN'))
def fa(p):
 d={};k=None
 for l in open(p):
  if l[0]=='>':k=l[1:].split()[0];d[k]=[]
  else:d[k].append(l.strip())
 return {k:''.join(v) for k,v in d.items()}
print(os.listdir(P)[:5])
tot=ok=0;bad=[]
from collections import Counter
c=Counter()
for f in sorted(glob.glob(Q+'*.fna')):
  rep=os.path.basename(f)[:-4]
  cand=[p for p in glob.glob(P+rep+'*') if p.endswith(('.faa','.fasta','.fa'))]
  if not cand: c['noprot']+=1; continue
  prot=fa(cand[0]); nt=fa(f)
  for k,s in nt.items():
    if k not in prot: c['missing']+=1; continue
    s2=rc(s) if strand[k] in ('-1','-') else s
    # determine frame: try offsets 0,1,2
    p=prot[k].rstrip('*')
    best=0
    for off in range(3):
      t=''.join(code.get(s2[j:j+3],'X') for j in range(off,len(s2)-2,3))
      if p[1:30] in t: best=1;break
    tot+=1; ok+=best; c[(part[k],strand[k],best)]+=1
print(tot,ok); print(c)
