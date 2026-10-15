# Independent small local reference-record audit (no biological bulk computation).
import csv,json,re,hashlib
from pathlib import Path
R=Path('reproducibility/darkmatter/family_discovery_20261015/f12')
rows=list(csv.DictReader(open(R/'references_v2/verified_references_metadata.tsv'),delimiter='\t'))
issues=[]
for r in rows:
 p=R/r['raw_fasta_path']
 lines=p.read_text().splitlines()
 s=''.join(x.strip() for x in lines if not x.startswith('>'))
 if not lines[0].startswith('>') or len(s)!=int(r['length_aa']):issues.append([r['accession'],'invalid fasta or length'])
 if r['source_repo'].startswith('UniProt'):
  j=json.loads((R/r['raw_record_path']).read_text())
  if s!=j['sequence']['value']:issues.append([r['accession'],'json/fasta sequence mismatch'])
  if j['primaryAccession']!=r['accession']:issues.append([r['accession'],'accession mismatch'])
gb=(R/'references_v2/raw_records/ncbi/MK474076.gb').read_text()
pid=re.search('/protein_id="([^"]+)"',gb).group(1)
trans=re.sub(r'\s','',re.search('/translation="([^"]+)"',gb,re.S).group(1))
seq=''.join(x.strip() for x in (R/'references_v2/raw_records/ncbi/QCI56913.1.fasta').read_text().splitlines() if not x.startswith('>'))
assert pid=='QCI56913.1' and trans==seq
print(json.dumps({'n_reference_rows':len(rows),'issues':issues,'spl_nucleotide_protein_translation_match':True,'spl_length':len(seq),'note':'Integrity/source-record identity only; does not validate class_assignment or evidence tiers.'},indent=2))
