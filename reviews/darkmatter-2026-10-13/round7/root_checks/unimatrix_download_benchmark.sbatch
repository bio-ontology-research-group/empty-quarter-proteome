#!/bin/bash
#SBATCH --job-name=smag_net_bench
#SBATCH --partition=debug
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:30:00
#SBATCH --output=/data/hohndor/darkmatter_round7_download/benchmark_%j.log
set -euo pipefail
cd /data/hohndor/darkmatter_round7_download
python3 -u - <<'PY'
import urllib.request,hashlib,json,time,os
x=json.load(open('unimatrix_download_benchmark.json'));t=time.time();dest=x['key'];tmp=dest+'.tmp'
with urllib.request.urlopen(x['links']['self'],timeout=120) as r,open(tmp,'wb') as out:
 while b:=r.read(1024*1024):out.write(b)
b=open(tmp,'rb').read();assert len(b)==x['size'];assert hashlib.md5(b).hexdigest()==x['checksum'].split(':')[-1];os.replace(tmp,dest)
result={'part':dest,'bytes':len(b),'seconds':time.time()-t,'sha256':hashlib.sha256(b).hexdigest()};json.dump(result,open('benchmark_result.json','w'),indent=2);print(json.dumps(result))
PY
