#!/usr/bin/env python3
"""Auxiliary SMAG transfer, reverse order, checked HTTP ranges and final public hashes."""
import concurrent.futures,csv,hashlib,json,os,re,time,urllib.request,urllib.error
from pathlib import Path
BASE=Path('downloads/smag_parts');TMP=BASE/'legacy_resumed_tmp';TMP.mkdir(parents=True,exist_ok=True)
def hashes(p):
 m=hashlib.md5();s=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(2**20):m.update(b);s.update(b)
 return m.hexdigest(),s.hexdigest()
def fetch(x):
 key=x['key'];size=x['size'];expected=x['checksum'].split(':')[-1];dest=BASE/key;tmp=TMP/(key+'.tmp');started=time.time()
 for attempt in range(1,31):
  # Other downloader may have completed this part while this task was waiting/retrying.
  if dest.exists() and dest.stat().st_size==size:
   m,s=hashes(dest)
   if m==expected:return {'part_key':key,'byte_size':size,'md5_checksum':m,'sha256_checksum':s,'status':'VERIFIED','seconds':round(time.time()-started,2)}
  offset=tmp.stat().st_size if tmp.exists() else 0
  if offset==size:
   m,s=hashes(tmp)
   if m==expected:
    os.replace(tmp,dest);print(key,'verified',round(time.time()-started,1),'s',flush=True)
    return {'part_key':key,'byte_size':size,'md5_checksum':m,'sha256_checksum':s,'status':'RESUMED_VERIFIED','seconds':round(time.time()-started,2)}
   tmp.rename(tmp.with_name(tmp.name+f'.bad_md5_{time.time_ns()}'));offset=0
  elif offset>size:
   raise RuntimeError(f'Oversized partial file: {tmp}')
  url='https://zenodo.org/records/8223844/files/'+key+'?download=1'
  headers={'User-Agent':'KAUST-research-reproducibility/1.0','Accept-Encoding':'identity'}
  if offset:headers['Range']=f'bytes={offset}-'
  try:
   with urllib.request.urlopen(urllib.request.Request(url,headers=headers),timeout=120) as response:
    status=response.status
    if status==206:
     cr=response.headers.get('Content-Range','');match=re.fullmatch(r'bytes (\d+)-(\d+)/(\d+)',cr)
     if not match or int(match[1])!=offset or int(match[3])!=size:raise RuntimeError(f'Invalid range: {cr}, expected offset{offset},total{size}')
     mode='ab' if offset else 'wb'
    elif status==200:
     mode='wb' # Server ignored Range: replace the partial, never append a complete response.
    else:raise RuntimeError(f'Unexpected HTTP status{status}')
    with tmp.open(mode) as out:
     while b:=response.read(2**20):out.write(b)
   print(key,'received',tmp.stat().st_size,'/',size,'attempt',attempt,flush=True)
   time.sleep(1)
  except (OSError,urllib.error.HTTPError) as exc:
   delay=min(120,5*attempt)
   if isinstance(exc,urllib.error.HTTPError) and exc.code==429:
    value=exc.headers.get('Retry-After','')
    if value.isdigit():delay=max(delay,int(value))
   print(key,'retry',attempt,type(exc).__name__,str(exc),'partial_bytes',tmp.stat().st_size if tmp.exists() else 0,'wait',delay,flush=True);time.sleep(delay)
 raise RuntimeError(f'{key}: failed after30 attempts; partial retained')
meta=json.load(open('results/zenodo_8223844.json'));parts=sorted([x for x in meta['files'] if x['key'].startswith('mag.tar.gz.')],key=lambda x:x['key'],reverse=True);assert len(parts)==373
parts=sorted(parts,key=lambda x:x['key'])[200:300]
rows=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
 for result in pool.map(fetch,parts):rows.append(result)
with open('results/smag_legacy_resumed_verified_parts_manifest.tsv','w') as out:
 w=csv.DictWriter(out,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(sorted(rows,key=lambda x:x['part_key']))
print('LEGACY_SUBSET_VERIFIED',sum(x['byte_size'] for x in rows),flush=True)
