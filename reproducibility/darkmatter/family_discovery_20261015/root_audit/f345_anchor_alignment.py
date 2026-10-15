"""Alignment-anchored reference-site plausibility checks; not enzyme-activity validation."""
import csv,json,hashlib,os,subprocess
from pathlib import Path
from Bio import AlignIO,SeqIO
W=Path('/ibex/scratch/projects/c2014/rob/metagenomics');D=W/'06_darkmatter/family_discovery_20261015_clean';C=D/'root_audit/f345_comparison';R=D/'root_audit/f345_anchor_alignment';R.mkdir(exist_ok=True);F=D/'f345/references'
spec={'Type1':('P02945','PF01036_Bac_rhodopsin.hmm','retinal'),'BetaCA':('P61517','PF00484_Pro_CA.hmm','Zn(2+)'),'EctC':('A0A385TEI2','PF06339_Ectoine_synth.hmm','EctC_Czech2019'),'GammaCA':('P40881','NF050158_gamma_CA.hmm','Zn(2+)')};summary={}
for family,(acc,profile,kind)in spec.items():
 meta={r['id']:r for r in csv.DictReader(open(C/family/'metadata.tsv'),delimiter='\t')};seqs=SeqIO.to_dict(SeqIO.parse(C/family/'combined.faa','fasta'));rid=family+'__curated_reference__'+acc;assert rid in seqs
 record=json.load(open(F/'raw_uniprot_json'/f'{acc}.json'));assert str(seqs[rid].seq)==record['sequence']['value']
 features=[]
 for feat in record.get('features',[]):
  keep=(kind=='retinal'and'retinylidene'in feat.get('description',''))or(kind=='Zn(2+)'and feat.get('ligand',{}).get('name')=='Zn(2+)')
  if keep:
   a=feat['location']['start']['value'];b=feat['location']['end']['value'];assert a==b;features.append({'position':a,'residue':record['sequence']['value'][a-1],'type':feat['type'],'description':feat.get('description',''),'ligand':feat.get('ligand'), 'evidences':feat.get('evidences',[])})
 if kind=='EctC_Czech2019':
  features=[{'position':p,'residue':aa,'type':'primary-paper catalytic-site plausibility anchor','description':'Czech et al. 2019 doi:10.1038/s41598-018-36247-w, Fig4 and Trp21 section; native 130aa A0A385TEI2 / 5ONN','source_doi':'10.1038/s41598-018-36247-w'}for p,aa in [(21,'W'),(57,'E'),(84,'Y'),(92,'H')]]
  assert all(record['sequence']['value'][f['position']-1]==f['residue']for f in features)
 assert features
 with open(R/f'{family}.sto','w')as out,open(R/f'{family}.log','w')as log:subprocess.run(['hmmalign','--outformat','Stockholm',str(F/'profiles'/profile),str(C/family/'combined.faa')],check=True,stdout=out,stderr=log)
 aln=AlignIO.read(R/f'{family}.sto','stockholm');rf=aln.column_annotations['reference_annotation'];byid={r.id:str(r.seq)for r in aln};ref=byid[rid];pos=0;poscol={}
 for col,char in enumerate(ref):
  if char not in '.-':pos+=1;poscol[pos]=col
 assert pos==record['sequence']['length']
 for feat in features:
  feat['alignment_column']=poscol[feat['position']];feat['profile_match_column']=rf[feat['alignment_column']]not in '.- '
 rows=[]
 for pid,s in byid.items():
  if meta[pid]['source']!='EQ':continue
  residues=''.join(s[x['alignment_column']].upper()for x in features);expected=''.join(x['residue']for x in features);assessable=all(x['profile_match_column']for x in features)
  rows.append({'id':pid,'family':family,'reference':acc,'reference_positions':','.join(str(x['position'])for x in features),'reference_residues':expected,'aligned_candidate_residues':residues,'all_sites_profile_aligned':int(assessable),'exact_reference_site_pattern':int(assessable and residues==expected)})
 with open(R/f'{family}.site_checks.tsv','w')as f:
  wr=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');wr.writeheader();wr.writerows(rows)
 summary[family]={'reference':acc,'profile':profile,'reference_json_sha256':hashlib.sha256((F/'raw_uniprot_json'/f'{acc}.json').read_bytes()).hexdigest(),'features':features,'EQ_complete_candidates':len(rows),'exact_reference_site_pattern':sum(r['exact_reference_site_pattern']for r in rows),'all_sites_profile_aligned':all(x['profile_match_column']for x in features),'limits':'homologous alignment and reference-site identity only; alternatives/compensatory substitutions possible; no ligand occupancy, activity, substrate specificity or thermal properties demonstrated'}
(R/'summary.json').write_text(json.dumps({'job_id':os.environ['SLURM_JOB_ID'],'families':summary},indent=2))
