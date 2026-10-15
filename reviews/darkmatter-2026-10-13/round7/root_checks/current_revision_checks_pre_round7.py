#!/usr/bin/env python3
"""Recompute review-revision claims from retained inputs; no upstream search claim."""
import csv,json,hashlib
from pathlib import Path
import pandas as pd
from round2_checks import read
H=Path(__file__).resolve().parent;R=H.parents[1];checks=[]
def record(key,value,expected,source):
 assert value==expected,(key,value,expected)
 checks.append(dict(claim_id=key,value=value,expected=expected,source=source,validation='recomputed from retained table'))
m=pd.DataFrame(read('mimag_rrna_trna_clean/results/mimag_tiers.tsv'))
record('quality.mimag_high',int((m.mimag_tier=='high').sum()),126,'inputs/mimag_rrna_trna_clean/results/mimag_tiers.tsv.gz')
record('quality.mimag_medium',int((m.mimag_tier=='medium').sum()),5312,'inputs/mimag_rrna_trna_clean/results/mimag_tiers.tsv.gz')
record('quality.all_rrnas',int(m[['has_5S','has_16S','has_23S']].astype(int).all(axis=1).sum()),157,'inputs/mimag_rrna_trna_clean/results/mimag_tiers.tsv.gz')
record('quality.trna18',int((m.n_tRNA_unique_aa.astype(int)>=18).sum()),3629,'inputs/mimag_rrna_trna_clean/results/mimag_tiers.tsv.gz')
c=pd.read_csv(H/'culture_nucleotide_cigar/per_gene.tsv',sep='\t');c=c[c.eligible==1]
for d,sh,n,k in [(0,0,16056,4877),(0,1,65477,62925),(1,0,1105,516),(1,1,236,219)]:
 g=c[(c.dark==d)&(c.protein_shared==sh)];record(f'culture.dna.{d}.{sh}.n',len(g),n,'culture_nucleotide_cigar/per_gene.tsv');record(f'culture.dna.{d}.{sh}.supported',int(g.dna_supported.sum()),k,'culture_nucleotide_cigar/per_gene.tsv')
lin=pd.DataFrame(read('06_darkmatter/lineage_darkness_clean/results/v2_reps.tsv'))
for k,v in {1:10,2:560,3:272,4:1006,5:1613,6:1855,7:122}.items():record(f'taxonomy.first_digit_or_empty.{k}',int((lin.first_unnamed_rank.astype(int)==k).sum()),v,'inputs/06_darkmatter/lineage_darkness_clean/results/v2_reps.tsv.gz')
st=pd.read_csv(H/'structural_raw_audit/per_model.tsv',sep='\t')
for arm,expected in [('dark',(2838,22,1031,95,1009,73,97,797)),('control',(3013,66,101,68,35,2,69,26))]:
 a=st[st.arm==arm];retrieved=a.raw_rows>0
 values=[len(a),int((~retrieved).sum()),int((a.qtm<.5).sum()),int((a.alntm<.5).sum()),int(((a.qtm<.5)&retrieved).sum()),int(((a.alntm<.5)&retrieved).sum()),int((a.alntm_min50<.5).sum()),int(((a.qtm<.5)&(a.alntm_min100_cov50>=.5)).sum())]
 for key,value,want in zip(['qualified','no_rows','qtm_below','alntm_below','retrieved_qtm_below','retrieved_alntm_below','min50_alntm_below','qtm_below_but_min100_cov50_alntm_hit'],values,expected):record(f'structure.raw.{arm}.{key}',value,want,'structural_raw_audit/per_model.tsv')
composition=pd.DataFrame(read('06_darkmatter/lineage_darkness_clean/results/v2_dark_composition_by_phylum.tsv'))
t=composition[composition.phylum=='Thermoplasmatota'].iloc[0]
record('lineage.thermoplasmatota_family_enrichment_1dp',round(float(t.enrichment_families),1),4.6,'inputs/06_darkmatter/lineage_darkness_clean/results/v2_dark_composition_by_phylum.tsv.gz')
flags=set((H/'antifam/flagged_families.txt').read_text().split());ov=pd.read_csv(H/'shadow_overlap/family_flags.tsv',sep='\t');flags.update(ov.loc[ov.opposite_nondark_overlap_members>0,'family_id'])
d=pd.DataFrame(read('06_darkmatter/darkfam_basis_clean/results/selection/dnds.tsv'));v=d[(d['set']=='dark')&d.family_id.isin(flags)].median_dNdS.astype(float)
record('selection.flagged_evaluable',len(v),91,'inputs/06_darkmatter/darkfam_basis_clean/results/selection/dnds.tsv.gz + antifam/shadow flags')
record('selection.flagged_median',round(float(v.median()),4),.7622,'inputs/06_darkmatter/darkfam_basis_clean/results/selection/dnds.tsv.gz + antifam/shadow flags')
record('selection.flagged_below_one',int((v<1).sum()),70,'inputs/06_darkmatter/darkfam_basis_clean/results/selection/dnds.tsv.gz + antifam/shadow flags')
for dark,total,union in [(0,81533,70354),(1,1341,752)]:
 a=c[c.dark==dark];record(f'culture.union.{dark}.eligible',len(a),total,'culture_nucleotide_cigar/per_gene.tsv');record(f'culture.union.{dark}.supported',int(((a.protein_shared==1)|(a.dna_supported==1)).sum()),union,'culture_nucleotide_cigar/per_gene.tsv')
locus=pd.read_csv(H/'culture_locus_all/per_gene.tsv',sep='\t')
for dark,n,absent in [(0,67802,607),(1,735,344)]:
 a=locus[locus.dark==dark];record(f'culture.locus_all.{dark}.dna_supported',len(a),n,'culture_locus_all/per_gene.tsv');record(f'culture.locus_all.{dark}.no_major_cds',int((a.classification=='no_major_CDS_overlap').sum()),absent,'culture_locus_all/per_gene.tsv')
member=pd.read_csv(H/'culture_family_membership/per_gene.tsv',sep='\t');record('culture.admitted_family_members',int(member.admitted_family.notna().sum()),65,'culture_family_membership/per_gene.tsv');record('culture.distinct_admitted_families',int(member.admitted_family.nunique()),47,'culture_family_membership/per_gene.tsv')
import sys
if '--write-ledger' in sys.argv:
 with (H/'current_claims.tsv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(checks[0]),delimiter='\t');w.writeheader();w.writerows(checks)
else:
 ledger=list(csv.DictReader((H/'current_claims.tsv').open(),delimiter='\t'))
 assert len(ledger)==len(checks)
 for actual,saved in zip(checks,ledger):
  assert all(str(v)==saved[k] for k,v in actual.items()),(actual,saved)
text=(R/'manuscript/main.tex').read_text()
for invalid in ['Ribosomal and transfer RNA content was not','The largest single\nowner','1.38 with an interval','Fig.~\\ref{fig:design}d','Fig.~\\ref{fig:design}e','No pH measurement exists','cannot be rebuilt from its sources','fig_specificity_controls.pdf','fig_environment.pdf']:
 assert invalid not in text,invalid
for required in ['126 as high quality','516 of 1{,}105','62{,}925 of 65{,}477','797 have an alignment-normalised','These means and the pooled medians are different estimands','Supplementary Note S1']:
 assert required in text,required
print('PASS:',len(checks),'revision claims and manuscript regression guards; historical whole-paper gates remain separate.')
