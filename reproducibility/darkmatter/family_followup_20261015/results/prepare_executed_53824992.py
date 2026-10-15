from common import *
import json,subprocess,urllib.request,zipfile,xml.etree.ElementTree as ET
from Bio import SeqIO
R=N/'results';I=N/'inputs';I.mkdir(exist_ok=True);R.mkdir(exist_ok=True)
urls={'gem_metadata.tsv':'https://portal.nersc.gov/GEM/genomes/genome_metadata.tsv','6SU3.pdb':'https://files.rcsb.org/download/6SU3.pdb','6SU3.json':'https://data.rcsb.org/rest/v1/core/entry/6SU3'}
for name,url in urls.items():
    if not (I/name).exists():
        subprocess.run(['curl','-fL','--retry','2','--max-time','180',url,'-o',str(I/name)],check=True)
write(R/'sources.tsv',[{'file':str(I/name),'url':url,'sha256':sha(I/name)}for name,url in urls.items()])
# Convert publisher-supplied workbook using standard-library XML, preserving original cells.
ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
with zipfile.ZipFile(I/'smag_extended_table2.xlsx') as z:
    strings=[]
    if 'xl/sharedStrings.xml' in z.namelist():
        strings=[''.join(e.itertext())for e in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',ns)]
    sheets=[]
    for fn in sorted(x for x in z.namelist()if x.startswith('xl/worksheets/sheet')and x.endswith('.xml')):
        rr=[]
        for row in ET.fromstring(z.read(fn)).findall('.//s:row',ns):
            vals={}
            for c in row.findall('s:c',ns):
                letters=''.join(x for x in c.attrib['r']if x.isalpha());idx=0
                for ch in letters:idx=idx*26+ord(ch)-64
                v=c.find('s:v',ns);value=v.text if v is not None else ''
                if c.attrib.get('t')=='s':value=strings[int(value)]
                elif c.attrib.get('t')=='inlineStr':value=''.join(c.find('s:is',ns).itertext())
                vals[idx-1]=value
            if vals:rr.append([vals.get(i,'')for i in range(max(vals)+1)])
        dest=I/(Path(fn).stem+'.tsv')
        with open(dest,'w')as f:csv.writer(f,delimiter='\t').writerows(rr)
        sheets.append({'file':str(dest),'rows':len(rr),'first_rows':rr[:4]})
(R/'metadata_headers.json').write_text(json.dumps({'smag':sheets,'gem':tab(I/'gem_metadata.tsv')[:2]},indent=2))
# All quality-filtered heliorhodopsins, exact-sequence deduplicated, retain source mapping.
meta=tab(D/'root_audit/f3_site_conserving/Helio/metadata.tsv')
seqs=SeqIO.to_dict(SeqIO.parse(D/'root_audit/f345_comparison/Helio/combined.faa','fasta'))
bysha={};records=[];mapping=[]
for row in sorted(meta,key=lambda x:x['id']):
    seq=str(seqs[row['id']].seq);digest=hashlib.sha256(seq.encode()).hexdigest();assert digest==row['sequence_sha256']
    if digest not in bysha:
        pid='HER_'+digest[:16];bysha[digest]=pid;records.append((pid,seq))
    mapping.append({**row,'structure_id':bysha[digest]})
with open(I/'helio_structures.faa','w')as f:
    for pid,seq in records:f.write('>'+pid+'\n'+seq+'\n')
write(I/'helio_structure_mapping.tsv',mapping)
anchor=next(r['structure_id']for r in mapping if r['original_id']=='AVZ43932.1')
pilot=[r for r in records if r[0]==anchor]+[next(r for r in records if r[0]==m['structure_id'])for src in ['EQ','SMAG','GEM']for m in [next(m for m in mapping if m['source']==src)]]
with open(I/'helio_pilot.faa','w')as f:
    for pid,seq in dict(pilot).items():f.write('>'+pid+'\n'+seq+'\n')
(R/'prepare_summary.json').write_text(json.dumps({'job':os.environ['SLURM_JOB_ID'],'helio_records':len(meta),'unique_sequences':len(records),'pilot_anchor':anchor},indent=2))
