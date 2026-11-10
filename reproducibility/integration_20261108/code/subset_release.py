from pathlib import Path
import json,gzip,hashlib,re
I=Path(__file__).resolve().parents[1]/'inputs';rules=json.loads((I/'panels.json').read_text())['rules'];wanted=set(re.findall(r"'([^']+)'",' '.join(rules.values())));data=gzip.open(I/'microtrait.hmmdb.gz','rt').read();found={};aliases={}
for block in re.split(r'^//\s*$',data,flags=re.M):
 if not block.strip():continue
 m=re.search(r'^NAME\s+(\S+)',block,re.M);assert m
 name=m[1]
 if name in wanted:
  assert name not in found;assert re.search(r'^TC\s',block,re.M),name
  found[name]=(block.strip()+'\n//\n').encode()
print('Release matches',len(found),'of',len(wanted),'missing',sorted(wanted-set(found)),flush=True)
# Keep an immutable release subset distinct from the earlier source-tree pilot.
(I/'panels-release.hmm').write_bytes(b''.join(found[n] for n in sorted(found)))
(I/'hmm-release-manifest.json').write_text(json.dumps([dict(leaf=n,name=n,sha256=hashlib.sha256(found[n]).hexdigest(),source='microtrait.hmmdb.gz release latest',release_asset_sha256=hashlib.sha256((I/'microtrait.hmmdb.gz').read_bytes()).hexdigest()) for n in sorted(found)],indent=2))
(I/'missing_release_hmms.json').write_text(json.dumps(sorted(wanted-set(found))))
