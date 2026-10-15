from pathlib import Path
import csv,hashlib,json,random,sys,zipfile,xml.etree.ElementTree as ET
b=Path('reproducibility/darkmatter/round7_analyses/h3_references/results')
source=Path('reviews/darkmatter-2026-10-13/round7/smag_source_audit/evidence/Supplementary_Data_2.xlsx')
with zipfile.ZipFile(source) as z:
 strings=[''.join(x.itertext()) for x in ET.fromstring(z.read('xl/sharedStrings.xml'))]
 tree=ET.fromstring(z.read('xl/worksheets/sheet1.xml'))
 vals=[];id_column=None
 for row in tree.findall('.//{*}row'):
  decoded={}
  for c in row:
   v=c.find('{*}v')
   if v is None:continue
   decoded[c.get('r','').rstrip('0123456789')]=strings[int(v.text)] if c.get('t')=='s' else v.text
  if id_column is None:
   matches=[k for k,v in decoded.items() if v=='user_genome']
   if matches:
    assert len(matches)==1;id_column=matches[0]
  elif id_column in decoded:vals.append(decoded[id_column])
assert id_column is not None;xls_ids=vals
frame=(b/'smag_catalogue_40039_ids.txt').read_text().splitlines()
assert len(xls_ids)==len(frame)==len(set(frame))==40039 and set(frame)==set(xls_ids)
observed=(b/'sampled_smag_1500_ids_seed20261014.txt').read_text().splitlines()
expected=sorted(random.Random(20261014).sample(sorted(frame),1500))
assert observed==expected and len(set(observed))==1500
out={'status':'PASS','python_version':sys.version,'smag_source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'smag_frame_matches_all_40039_xlsx_ids':True,'smag_exact_replay_rule':'sorted(random.Random(20261014).sample(sorted(frame),1500))','gem':{}}
gem=sorted((b/'gem_catalogue_52515_ids.txt').read_text().splitlines());assert len(gem)==len(set(gem))==52515
for seed in [20261014,20261015,20261016]:
 shuffled=list(gem);random.Random(seed).shuffle(shuffled);expected=set(shuffled[:2055]);members={}
 for caller in ['harmonized','published']:
  rows=list(csv.DictReader((b/f'gem_{caller}_manifest_seed{seed}.tsv').open(),delimiter='\t'))
  assert len(rows)==2055 and {r['genome_id'] for r in rows}==expected
  members[caller]=len(rows)
 out['gem'][str(seed)]={'same_exact_genome_set_both_callers':True,'counts':members,'selection_rule':'random.Random(seed).shuffle(sorted(frame)); first2055'}
Path('reviews/darkmatter-2026-10-13/round7/root_checks/h3_sampling_audit.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out))
