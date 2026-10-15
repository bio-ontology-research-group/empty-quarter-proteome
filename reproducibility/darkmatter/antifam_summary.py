#!/usr/bin/env python3
import csv,json,sys
from pathlib import Path
W=Path('/ibex/scratch/projects/c2014/rob/metagenomics');R=W/'06_darkmatter/darkfam_basis_clean/results';out=W/'06_darkmatter/p4_review_20261013/results/antifam'
def hits(name):return {l.split()[0] for l in (out/name).open() if not l.startswith('#') and l.strip()}
dark=hits('dark.tbl');ctrl=hits('ctrl.tbl');fams=set((R/'FAMILIES_FINAL.txt').read_text().split());flags=set();members=set()
for line in (R/'fam_cluster.tsv').open():
 f,m=line.rstrip().split('\t')
 if f in fams and m in dark:flags.add(f);members.add(m)
dnds=list(csv.DictReader((R/'selection/dnds.tsv').open(),delimiter='\t'));print(dnds[0].keys())
import numpy as np
key='family_id' if 'family_id' in dnds[0] else 'family'
remaining=[float(r['median_dNdS']) for r in dnds if r['set']=='dark' and r[key] not in flags]
summary={'dark_proteins_with_hit':len(dark),'control_proteins_with_hit':len(ctrl),'admitted_members_with_hit':len(members),'admitted_families_with_hit':len(flags),'evaluable_unflagged_dark_n':len(remaining),'evaluable_unflagged_dark_median_dnds':float(np.median(remaining)),'evaluable_unflagged_dark_fraction_below_one':float(np.mean(np.array(remaining)<1))}
(out/'flagged_families.txt').write_text('\n'.join(sorted(flags))+'\n');(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary)
