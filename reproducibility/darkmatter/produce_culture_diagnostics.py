#!/usr/bin/env python3
"""Produce per-pair diagnostic outputs separately from verification."""
import json
from pathlib import Path
import pandas as pd
H=Path(__file__).resolve().parent
c=pd.read_csv(H/'culture_nucleotide_cigar/per_gene.tsv',sep='\t');c=c[c.eligible==1]
per=c[c.protein_shared==1].groupby('rep').dna_supported.agg(['sum','count']);per['recovery']=per['sum']/per['count'];per.to_csv(H/'culture_nucleotide_cigar/per_pair_positive_control.tsv',sep='\t');keep=per.index[per.recovery>=.9];post=[]
for (d,sh),g in c[c.rep.isin(keep)].groupby(['dark','protein_shared']):post.append(dict(dark=int(d),protein_shared=int(sh),n=len(g),supported=int(g.dna_supported.sum())))
(H/'culture_nucleotide_cigar/posthoc_recovery_subset.json').write_text(json.dumps({'pairs':len(keep),'criterion':'post hoc >=90% recovery of protein-shared calls','per_pair_min':per.recovery.min(),'per_pair_max':per.recovery.max(),'groups':post},indent=2)+'\n')
