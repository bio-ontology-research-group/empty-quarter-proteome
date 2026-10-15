#!/usr/bin/env python3
import json
from pathlib import Path
import pandas as pd
from round2_checks import read
H=Path(__file__).resolve().parent
f=pd.read_csv(H/'shadow_overlap/family_flags.tsv',sep='\t');overlap=set(f.loc[f.opposite_nondark_overlap_members>0,'family_id']);anti=set((H/'antifam/flagged_families.txt').read_text().split());d=pd.DataFrame(read('06_darkmatter/darkfam_basis_clean/results/selection/dnds.tsv'));out=[]
for label,flags in [('antifam',anti),('opposite_overlap',overlap),('union',anti|overlap)]:
 a=d[(d['set']=='dark')&~d.family_id.isin(flags)].median_dNdS.astype(float)
 b=d[(d['set']=='dark')&d.family_id.isin(flags)].median_dNdS.astype(float)
 out.append(dict(exclusion=label,flagged_families=len(flags),remaining_evaluable=len(a),median_dnds=a.median(),fraction_below_one=(a<1).mean(),flagged_evaluable=len(b),flagged_median=b.median(),flagged_q1=b.quantile(.25),flagged_q3=b.quantile(.75),flagged_fraction_below_one=(b<1).mean()))
(H/'shadow_overlap/selection_sensitivity.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
