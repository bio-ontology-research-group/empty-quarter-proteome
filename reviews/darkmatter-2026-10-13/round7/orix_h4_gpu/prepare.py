import csv, hashlib, json, os, pathlib, shutil, subprocess, tarfile, urllib.request
P=pathlib.Path
base=P('/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean')
h4=base/'h4_structure'; root=h4/'orix_gpu_20261014'; out=root/'payload'
assert os.environ.get('SLURM_JOB_ID')
assert not out.exists(), 'Immutable preparation directory already exists'
out.mkdir(parents=True)
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''): h.update(b)
 return h.hexdigest()
def run(*a): subprocess.run([str(x) for x in a],check=True)
# The release matches the original CPU run, while the GPU prefilter remains a different algorithm.
u='https://github.com/steineggerlab/foldseek/releases/download/10-941cd33/foldseek-linux-gpu.tar.gz'
a=root/'foldseek-linux-gpu.tar.gz'; urllib.request.urlretrieve(u,a)
with tarfile.open(a) as t: t.extractall(out/'software',filter='data')
f=out/'software/foldseek/bin/foldseek'
v=subprocess.check_output([str(f),'version'],text=True).strip()
assert v=='941cd33ff0771cd2e3f144e3293e22a2b87e9fda',v
pin=json.load(open(root/'pdb_reference_verification.json'))
reference=[]
for e in pin['files']:
 p=P(e['path']); assert sha(p)==e['expected'],str(p)
 reference.append({'path':str(p),'sha256':e['expected']})
manifest=h4/'results/v7_full5851_cropped_regions_manifest.tsv'
rows=list(csv.DictReader(open(manifest),delimiter='\t'))
wc={r['model_id']+'.pdb':r['source_pdb_sha256'] for r in rows}
cr={r['cropped_pdb']:r['crop_pdb_sha256'] for r in rows if r['cropped_pdb']!='NONE'}
assert len(wc)==5851 and len(cr)==7773,(len(wc),len(cr))
for name,folder,expected in [('wc','v7_wc_queries',wc),('crop','v7_crops',cr)]:
 d=h4/folder
 assert {p.name for p in d.glob('*.pdb')}==set(expected)
 for n,h in expected.items(): assert sha(d/n)==h,n
 (out/name).mkdir()
 run(f,'createdb',d,out/name/'query','--threads','8')
 shutil.copy2(manifest,out/'query_manifest.tsv')
(out/'pdb').mkdir()
run(f,'makepaddedseqdb',P(reference[0]['path']),out/'pdb/target','--threads','8')
# Materialize any source symlinks before transfer: all references must be self-contained.
for p in list(out.rglob('*')):
 if p.is_symlink():
  target=p.resolve(); p.unlink(); shutil.copy2(target,p)
meta={'job':os.environ['SLURM_JOB_ID'],'host':os.uname().nodename,'foldseek_version':v,'download_url':u,'archive_sha256':sha(a),'reference_components':reference,'whole_models':5851,'crop_regions':7773,'query_manifest_sha256':sha(manifest),'comparison':'GPU ungapped prefilter sensitivity; not an exact CPU replay'}
(out/'provenance.json').write_text(json.dumps(meta,indent=2)+'\n')
with open(out/'payload.sha256','w') as w:
 for p in sorted(out.rglob('*')):
  if p.is_file() and p.name!='payload.sha256': w.write(sha(p)+'  '+str(p.relative_to(out))+'\n')
(root/'PREPARED').write_text(sha(out/'payload.sha256')+'\n')
print('PREPARED',json.dumps(meta))
