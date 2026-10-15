from pathlib import Path
import re,json,xml.etree.ElementTree as ET
p=Path('manuscript/references.bib');bib=p.read_text();audit=[]
for key in ['chan2021trnascan','eberhardt2012antifam','mugal2014time','nelson2020olgenie','orr2026spurious','shaw2023skani']:
 m=re.search(r'@\w+\{'+key+r',([\s\S]*?)(?=\n@|\Z)',bib);old=m[0];doi=re.search(r'doi\s*=\s*\{([^}]+)',old).group(1)
 matches=[]
 for f in Path('resources/literature/darkmatter-2026-10-13/verified').glob('*.json'):
  d=json.loads(f.read_text())
  if d.get('doi','').lower()==doi.lower():matches.append(d)
 assert len(matches)==1
 source=Path(matches[0]['path']);a=ET.parse(source).getroot().find('.//article-meta')
 assert any(x.text.lower()==doi.lower() for x in a.findall("article-id[@pub-id-type='doi']"))
 pages=a.findtext('page-range') or a.findtext('elocation-id') or a.findtext('fpage')
 if a.findtext('lpage'):pages=a.findtext('fpage')+'--'+a.findtext('lpage')
 pages=pages.replace('–','--')
 fields={k:v for k,v in [('volume',a.findtext('volume')),('number',a.findtext('issue')),('pages',pages)] if v}
 added={k:v for k,v in fields.items() if not re.search(r'\b'+k+r'\s*=',old)}
 body=old.rstrip();assert body.endswith('}')
 body=body[:-1].rstrip().rstrip(',')+',\n'+''.join(' '+k+'={'+v+'},\n' for k,v in added.items())+'}\n'
 bib=bib.replace(old,body);audit.append({'key':key,'doi':doi,'source':str(source),'added':added})
p.write_text(bib);Path('paper/references.bib').write_text(bib)
Path('reviews/darkmatter-2026-10-13/round5/bibliography_additions.json').write_text(json.dumps(audit,indent=2)+'\n')
