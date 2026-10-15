from common import *
import json,subprocess,zipfile
R=N/'inputs/bulzu2021';R.mkdir(parents=True,exist_ok=True)
api=R/'figshare_v3.json'
subprocess.run(['curl','-fL','--retry','2','--max-time','120','https://api.figshare.com/v2/articles/13286486/versions/3','-o',str(api)],check=True)
j=json.load(open(api));assert j['version']==3
files=[]
for entry in j['files']:
    p=R/entry['name']
    if not p.exists():subprocess.run(['curl','-fL','--retry','2','--max-time','180',entry['download_url'],'-o',str(p)],check=True)
    assert p.stat().st_size==entry['size'];assert hashlib.md5(p.read_bytes()).hexdigest()==entry['computed_md5'];files.append({'path':str(p),'url':entry['download_url'],'sha256':sha(p),'md5':entry['computed_md5']})
write(R/'downloads.tsv',files)
with zipfile.ZipFile(R/'Sequence_files.zip')as z:
    names=z.namelist();(R/'sequence_members.txt').write_text('\n'.join(names))
    for name in names:
        p=(R/'sequences'/name).resolve();assert (R/'sequences').resolve()in p.parents
    z.extractall(R/'sequences')
print((R/'Readme.txt').read_text(),flush=True);print('\n'.join(names),flush=True)
