#!/usr/bin/env python3
"""Round-4 corrections on IBEX. Preserve original upstream files and outputs."""
from pathlib import Path
import csv,hashlib,importlib.util,json
import numpy as np
from scipy.stats import t as tdist
W=Path('/ibex/scratch/projects/c2014/rob/metagenomics')
O=W/'06_darkmatter/p4_review_20261013/results/round4_corrections';O.mkdir(exist_ok=True)
pins={}
def pin(p):
 p=Path(p);h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 pins[str(p)]=h.hexdigest();return p
def module(name,p):
 spec=importlib.util.spec_from_file_location(name,pin(p));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
m=module('eco1',W/'06_darkmatter/darkfam_ecology_clean/scripts/eco_1_indicators.py')
env=m.load_env(pin(m.ENVF));samples,fams,M=m.load_matrix(pin(m.ABUND))
keep=[i for i,s in enumerate(samples) if env.get(s,{}).get('compartment') in ['Surface','Deep','PlantRoot']]
excluded=[s for i,s in enumerate(samples) if i not in keep];sub=M[:,keep];y=np.array([float(env[samples[i]]['compartment']=='PlantRoot') for i in keep]);present=(sub>0).sum(axis=1);mask=present>=m.MIN_PRESENT
r,p=m.spearman_vec(sub[mask],y);q=m.bh_fdr(p)
(O/'indicators').mkdir(exist_ok=True)
with (O/'indicators/plantroot_vs_other_all.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['family_id','spearman_r','p','q','n_present'])
 for row in sorted(zip(np.array(fams)[mask],r,p,q,present[mask]),key=lambda a:(a[3],-abs(a[1]))):w.writerow([row[0],f'{row[1]:.6f}',f'{row[2]:.6e}',f'{row[3]:.6e}',row[4]])
summary={'n_samples':len(keep),'root_samples':int(y.sum()),'excluded':excluded,'minimum_present':m.MIN_PRESENT,'n_tested':int(mask.sum()),'n_significant':int((q<.05).sum()),'positive':int(((q<.05)&(r>0)).sum()),'negative':int(((q<.05)&(r<0)).sum())}
(O/'compartment_cohort.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary,flush=True)
del M,sub
c=module('eco4',W/'06_darkmatter/darkfam_ecology_clean/scripts/eco_4_compartment_families.py');c.RES=str(O);c.AXIS=str(O/'indicators/plantroot_vs_other_all.tsv')
for p in [c.PHYLUM,c.PRESENCE,c.SITETRIP]:pin(p)
c.main()
# Correct residual degrees of freedom; preserve trait matrices and original table.
src=pin(W/'06_metabolism/osmolyte_clean/scripts/osmo_5_gradients.py').read_text()
old='return pearsonr(rx - C @ bx, ry - C @ by)'
new='''rho = pearsonr(rx - C @ bx, ry - C @ by)[0]
    df = len(x) - np.linalg.matrix_rank(C) - 1
    p = 2 * tdist.sf(abs(rho) * np.sqrt(df / max(1 - rho*rho, np.finfo(float).tiny)), df)
    return rho, p'''
assert src.count(old)==1
src=src.replace(old,new).replace('from scipy.stats import rankdata, spearmanr, pearsonr','from scipy.stats import rankdata, spearmanr, pearsonr, t as tdist').replace('with open(f"{OSMO}/gradient_tests.tsv", "w")',f'with open({str(O/"gradient_tests.tsv")!r}, "w")')
(O/'osmo_5_gradients_corrected.py').write_text(src)
for p in [W/'06_metabolism/osmolyte_clean/results/community_by_sample.tsv',W/'06_metabolism/uv_survival_clean/results/community_by_sample.tsv']:pin(p)
exec(compile(src,str(O/'osmo_5_gradients_corrected.py'),'exec'),{'__name__':'__main__'})
# Determine whether preferring qualifying complete rRNA hits changes any retained call.
a=module('mimag',W/'mimag_rrna_trna_clean/scripts/02_summarise.py');changes=[];nfiles=0
for p in sorted(a.RRNA_DIR.glob('*.gff')):
 pin(p);nfiles+=1;old=a.parse_barrnap_gff(p);good={}
 for line in p.open():
  if line.startswith('#'):continue
  f=line.rstrip().split('\t')
  if len(f)<9:continue
  match=a.re.search(r'Name=([^;]+)',f[8])
  if not match:continue
  name=match.group(1);length=abs(int(f[4])-int(f[3]))+1;partial='(partial)' in f[8] or 'partial=true' in f[8].lower()
  good[name]=good.get(name,False) or (not partial and length>=a.RRNA_LEN_FLOOR.get(name,0))
 for name in a.RRNA_LEN_FLOOR:
  before=a.has_full_rrna(old,name);after=good.get(name,False)
  if before!=after:changes.append({'rep':p.stem,'type':name,'old':before,'corrected':after})
(O/'rrna_parser_audit.json').write_text(json.dumps({'files':nfiles,'changed_calls':changes},indent=2)+'\n')
(O/'input_sha256.json').write_text(json.dumps(pins,indent=2)+'\n')
print('rRNA parser audit:',nfiles,'files;',len(changes),'changed calls',flush=True)
