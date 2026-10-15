from pathlib import Path
import csv,json,collections
root=Path('reproducibility/darkmatter/round6_hypotheses/h2_context/results')
read=lambda n:list(csv.DictReader((root/n).open(),delimiter='\t'))
f=read('pilot_family_conservation.tsv');n=read('pilot_neighbourhoods.tsv');b=read('bin_contig_verification.tsv')
assert len(f)==300 and len(n)==1518
bc=collections.Counter(x['verification_status'] for x in b)
for x in b:
 if x['verification_status'].startswith('verified'):
  assert x['n_raw_contigs']==x['n_bakta_contigs']==x['n_sha256_matched']
  assert int(x['n_rc_mismatches'])==int(x['n_hash_mismatches'])==0
eligible=[x for x in n if int(x['is_primary_locus']) and int(x['eligible_locus_primary'])]
assert len({(x['family'],x['species']) for x in eligible})==len(eligible)
for x in eligible:
 assert int(x['contig_len_bp'])>0
 assert all(x[p+'_audit_category'] in ('mapped_KO','mapped_noKO') for p in ('m2','m1','p1','p2'))
for x in f:
 ns={v['species'] for v in eligible if v['family']==x['family']}
 assert len(ns)==int(x['n_species_complete_window'])
 if int(x['eligible_primary']):
  assert len(ns)>=2
  assert 0<=float(x['r_fam_primary'])<=1
  assert int(x['best_ko_species_count'])<=len(ns)
 else:assert int(x['c_fam_primary'])==-1
out={'passed':True,'families':len(f),'loci':len(n),'selected_eligible_loci':len(eligible),'primary_eligible_families':sum(int(x['eligible_primary']) for x in f),'bin_verification_counts':dict(bc)}
Path('reviews/darkmatter-2026-10-13/round6/root_checks/h2_table_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
