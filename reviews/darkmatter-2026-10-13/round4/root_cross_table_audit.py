"""Independent key-level cross-table audit, without rewriting evidence."""
from pathlib import Path
import json
import pandas as pd
H=Path(__file__).resolve().parents[3]/'reproducibility/darkmatter'
c=pd.read_csv(H/'culture_nucleotide_cigar/per_gene.tsv',sep='\t')
l=pd.read_csv(H/'culture_locus_all/per_gene.tsv',sep='\t')
f=pd.read_csv(H/'culture_family_membership/per_gene.tsv',sep='\t')
s=pd.read_csv(H/'structural_raw_audit/per_model.tsv',sep='\t')
key=['rep','protein'];checks=[]
for name,d in [('culture',c),('locus',l),('membership',f)]:
 assert not d.duplicated(key).any(),name
 checks.append(name+': unique genome/protein keys')
expected=c[(c.eligible==1)&(c.dna_supported==1)]
joined=expected.merge(l,on=key,how='outer',indicator=True,validate='one_to_one',suffixes=('_call','_locus'))
assert joined._merge.eq('both').all()
for col in ['cultured','dark','protein_shared']:assert joined[col+'_call'].eq(joined[col+'_locus']).all(),col
checks.append('Every eligible DNA-supported call has exactly one locus classification with matching metadata')
expected=c[(c.eligible==1)&(c.dark==1)]
joined=expected.merge(f,on=key,how='outer',indicator=True,validate='one_to_one',suffixes=('_call','_family'))
assert joined._merge.eq('both').all()
for col in c.columns.difference(key):assert joined[col+'_call'].fillna('NA').eq(joined[col+'_family'].fillna('NA')).all(),col
checks.append('Family tracing covers exactly the eligible dark culture population and preserves every call field')
assert (s.raw_rows.eq(0)==((s.qtm==0)&(s.alntm==0))).all()
checks.append('Zero scores used by Figure 4 exactly identify zero raw retrieved rows in these inputs')
assert not s.duplicated(['arm','id']).any()
checks.append('Structural model identifiers unique within arms')
print(json.dumps({'passed':checks,'culture_calls':len(c),'dna_supported_loci':len(l),'family_traced_calls':len(f),'structural_models':len(s)},indent=2))
