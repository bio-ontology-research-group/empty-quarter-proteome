#!/usr/bin/env python3
"""Retain free full-text endpoints for DOI leads lacking canonical OA XML.
Access failures/challenges are recorded, never called downloaded papers.
"""
import csv,json,hashlib,urllib.request,concurrent.futures,re
from pathlib import Path
B=Path('resources/literature/darkmatter-2026-10-13');O=B/'papers/additional_free';O.mkdir(exist_ok=True)
def one(r):
 key=hashlib.sha256(r['doi'].lower().encode()).hexdigest()[:16];m=json.loads((B/'verified'/(key+'.json')).read_text()).get('metadata',{});urls=[]
 if m.get('pmcid'):urls.append('https://pmc.ncbi.nlm.nih.gov/articles/'+m['pmcid']+'/')
 for x in m.get('fullTextUrlList',{}).get('fullTextUrl',[]):
  if x.get('availabilityCode') in ['F','OA'] and x.get('documentStyle') in ['pdf','html']:urls.append(x['url'])
 results=[]
 for i,u in enumerate(dict.fromkeys(urls)):
  if i>=3:break
  rec={'doi':r['doi'],'title':r['title'],'url':u,'path':'','status':'failed','error':''}
  try:
   req=urllib.request.Request(u,headers={'User-Agent':'Mozilla/5.0 (research literature archive)'})
   with urllib.request.urlopen(req,timeout=25) as f:data=f.read()
   ext='.pdf' if data.startswith(b'%PDF') else '.html';p=O/(key+'_'+str(i)+ext);p.write_bytes(data);rec['path']=str(p)
   if ext=='.pdf':rec['status']='pdf_needs_identity_check'
   else:
    t=data.decode(errors='replace');doi=re.findall(r'<meta[^>]+name=["\']citation_doi["\'][^>]+content=["\']([^"\']+)',t,re.I)
    rec['status']='fulltext_doi_verified' if any(x.lower()==r['doi'].lower() for x in doi) and len(data)>20000 else 'retained_unverified_html'
   rec['sha256']=hashlib.sha256(data).hexdigest()
  except Exception as e:rec['error']=str(e)
  results.append(rec)
  if rec['status']=='fulltext_doi_verified':break
 return results
rows=[r for r in csv.DictReader((B/'verified_sources.tsv').open(),delimiter='\t') if r['status']=='doi_metadata_verified']
with concurrent.futures.ThreadPoolExecutor(max_workers=5) as ex:results=[r for group in ex.map(one,rows) for r in group]
(O/'retrievals.json').write_text(json.dumps(results,indent=2)+'\n');print([(r['doi'],r['status']) for r in results])
