"""Independent profile-first F1 pilot, distinct from Flash's DIAMOND-first pipeline."""
import csv,gzip,hashlib,json,os,subprocess
from pathlib import Path
W=Path('/ibex/scratch/projects/c2014/rob/metagenomics')
D=W/'06_darkmatter/family_discovery_20261015_clean'
R=D/'root_audit/f1_profile_pilot';R.mkdir(parents=True,exist_ok=True)
rows=list(csv.DictReader(open(W/'05_refine/_drep/out/data_tables/Wdb.csv')))
assert len(rows)==5438 and len({r['cluster'] for r in rows})==5438
rows=sorted(rows,key=lambda r:hashlib.sha256(r['genome'].encode()).hexdigest())[:256]
manifest=[];n=0
with open(R/'cohort.faa','w') as out:
 for r in rows:
  mag=r['genome'].removesuffix('.fa');p=W/f'06_darkmatter/drep_rep_proteomes/{mag}.faa.gz'
  h=hashlib.sha256(p.read_bytes()).hexdigest();ng=0
  with gzip.open(p,'rt') as f:
   for line in f:
    if line.startswith('>'):ng+=1;assert line[1:].startswith(mag+'::')
    out.write(line)
  manifest.append({'mag':mag,'cluster':r['cluster'],'path':str(p),'sha256':h,'proteins':ng});n+=ng
with open(R/'cohort.tsv','w') as f:
 w=csv.DictWriter(f,fieldnames=list(manifest[0]),delimiter='\t');w.writeheader();w.writerows(manifest)
# Original primary-source HMMs, independently copied and checksum pinned.
profiles=['PF00875','PF03441','PF04244']
with open(R/'cpf.hmm','wb') as f:
 for acc in profiles:
  p=D/f'f12/references_v2/{acc}.hmm';b=p.read_bytes();assert ('ACC   '+acc).encode() in b;f.write(b)
cmd=['hmmsearch','--cpu',os.environ['SLURM_CPUS_PER_TASK'],'--domE','1e-3','-E','1e-3','--domtblout',str(R/'domains.tsv'),'-o',str(R/'hmm_report.txt'),str(R/'cpf.hmm'),str(R/'cohort.faa')]
subprocess.run(cmd,check=True)
hits={}
for line in open(R/'domains.tsv'):
 if line.startswith('#') or not line.strip():continue
 p=line.split();assert len(p)>=22
 acc=p[4].split('.')[0];assert acc in profiles
 cov=(int(p[16])-int(p[15])+1)/int(p[5])
 if float(p[12])<=1e-5 and cov>=.70:
  rec={'id':p[0],'profile':acc,'domain_ievalue':float(p[12]),'domain_bitscore':float(p[13]),'profile_coverage':cov}
  if p[0] not in hits or rec['domain_bitscore']>hits[p[0]]['domain_bitscore']:hits[p[0]]=rec
records=[]
def emit(h,s,out):
 pid=h.split()[0]
 if pid not in hits:return
 r=hits[pid].copy();r.update(mag=pid.split('::')[0],length=len(s),sha256=hashlib.sha256(s.encode()).hexdigest(),partial='unknown',header=h)
 for item in h.split(';'):
  if 'partial=' in item:r['partial']=item.split('partial=')[-1].strip()
 records.append(r);out.write('>'+h+'\n'+s+'\n')
with open(R/'candidates.faa','w') as out:
 h=None;s=[]
 for line in open(R/'cohort.faa'):
  if line.startswith('>'):
   if h:emit(h,''.join(s),out)
   h=line[1:].strip();s=[]
  else:s.append(line.strip())
 if h:emit(h,''.join(s),out)
assert len(records)==len(hits)
with open(R/'candidates.tsv','w') as f:
 w=csv.DictWriter(f,fieldnames=list(records[0]) if records else ['id'],delimiter='\t');w.writeheader();w.writerows(records)
summary={'job_id':os.environ['SLURM_JOB_ID'],'genomes':256,'proteins':n,'profile_supported_candidates':len(records),'complete_partial00':sum(r['partial']=='00' for r in records),'carrier_genomes':len({r['mag'] for r in records}),'command':cmd,'interpretation':'profile-supported CPF candidates; not verified functions, subfamilies or novelty'}
(R/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
with open(R/'checksums.tsv','w') as f:
 for p in sorted(R.iterdir()):
  if p.is_file() and p.name!='checksums.tsv':f.write(hashlib.sha256(p.read_bytes()).hexdigest()+'\t'+p.name+'\n')
