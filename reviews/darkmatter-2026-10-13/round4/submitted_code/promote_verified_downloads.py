#!/usr/bin/env python3
"""Promote retained XML with exact embedded DOI identity into the canonical manifest."""
import csv,hashlib,json,shutil,xml.etree.ElementTree as ET
from pathlib import Path
base=Path('resources/literature/darkmatter-2026-10-13');bydoi={}
for p in (base/'papers').rglob('*.xml'):
 try:r=ET.parse(p)
 except ET.ParseError:continue
 for x in r.findall('.//article-meta/article-id'):
  if x.get('pub-id-type')=='doi' and x.text:bydoi.setdefault(x.text.lower().strip(),p)
manifest=base/'verified_sources.tsv';rows=list(csv.DictReader(manifest.open(),delimiter='\t'));added=[]
for row in rows:
 if row['status']=='fulltext_doi_verified':continue
 p=bydoi.get(row['doi'].lower())
 if p is None:continue
 key=hashlib.sha256(row['doi'].lower().encode()).hexdigest()[:16];out=base/'verified'/(key+'.xml');shutil.copyfile(p,out)
 row.update(status='fulltext_doi_verified',path=str(out),sha256=hashlib.sha256(out.read_bytes()).hexdigest(),error='',source_url='Retained download: '+str(p))
 meta=base/'verified'/(key+'.json');info=json.loads(meta.read_text()) if meta.exists() else {};info.update(row);info['retained_source_path']=str(p);meta.write_text(json.dumps(info,indent=2)+'\n');added.append(row['doi'])
with manifest.open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
print('Promoted',added)
