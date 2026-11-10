from common import *
import numpy as np,pandas as pd,gzip,re
R=N/'results';O=R/'hs3';O.mkdir(exist_ok=True)
z=np.load(R/'abundance.npz');raw=z['raw'];samples=z['samples'];g=z['genomes'];ix=pd.read_csv(R/'genome_index.tsv',sep='\t',index_col=0).loc[g];m=pd.read_csv(R/'metadata.tsv',sep='\t',index_col=0,dtype={'site':str}).loc[samples];c=pd.read_csv(R/'genome_counts_strict.tsv',sep='\t',index_col=0).loc[g]
candidates=[]
for site in sorted(m.site.unique()):
 for co in ['Surface','Deep','PlantRoot']:
  winter=np.flatnonzero(((m.site==site)&(m.compartment==co)&(m.trip=='trip3')).to_numpy());hot=np.flatnonzero(((m.site==site)&(m.compartment==co)&(m.trip=='trip4')).to_numpy())
  if not len(winter)or not len(hot):continue
  mw=raw[winter].max(0);mh=raw[hot].max(0)
  for j in np.flatnonzero((mw>0)&(mh>0)&(c.HSP20.to_numpy()>0)):
   ws=winter[raw[winter,j].argmax()];hs=hot[raw[hot,j].argmax()];domain='Thermoplasmatota'if'p__Thermoplasmatota;'in ix.iloc[j].gtdb_classification else 'Bacteria'if'd__Bacteria;'in ix.iloc[j].gtdb_classification else 'OtherArchaea'
   candidates.append({'genome':g[j],'site':site,'compartment':co,'domain':domain,'winter_sample':samples[ws],'hot_sample':samples[hs],'winter_raw_pct':mw[j],'hot_raw_pct':mh[j],'minimum_raw_pct':min(mw[j],mh[j])})
d=pd.DataFrame(candidates).sort_values('minimum_raw_pct',ascending=False);d.to_csv(O/'pair_eligibility.tsv.gz',sep='\t',index=False)
# Technical pilot: highest shared raw abundance per focal lineage, not selected by change magnitude or p value.
selected=pd.concat([d[d.domain==dom].head(1)for dom in ['Thermoplasmatota','Bacteria']],ignore_index=True);assert len(selected)==2
markers=pd.read_csv(R/'marker_loci.tsv.gz',sep='\t');targets=[]
for row in selected.itertuples():
 r=ix.loc[row.genome];pp=Path(r.proteome);allgenes=[]
 with gzip.open(pp,'rt')as f:
  for l in f:
   if not l.startswith('>'):continue
   x=l[1:].strip().split(' # ');pid=x[0];a,b=int(x[1]),int(x[2]);partial=re.search(r'partial=(\d\d)',x[4]).group(1)
   allgenes.append({'genome':row.genome,'protein':pid,'contig':pid.split('::')[1].rsplit('_',1)[0],'start':a,'end':b,'partial':partial})
 allgenes=pd.DataFrame(allgenes);mk=markers[(markers.genome==row.genome)&(markers.hmm_coverage>=.7)&(markers.partial.astype(str).str.zfill(2)=='00')]
 loci=allgenes.merge(mk[['protein','family']],on='protein',how='inner');loci['role']='marker'
 ctrl=allgenes[(allgenes.partial=='00')&(~allgenes.protein.isin(markers.protein))&((allgenes.end-allgenes.start).between(500,2000))].sample(n=min(60,len(allgenes)),random_state=20261016);ctrl['family']='background_CDS';ctrl['role']='background'
 target=pd.concat([loci,ctrl]);target['domain']=row.domain;targets.append(target)
targets=pd.concat(targets,ignore_index=True);targets.to_csv(O/'pilot_loci.tsv',sep='\t',index=False)
sr=[]
for row in selected.itertuples():
 bed=O/(row.domain+'.bed');dd=targets[targets.genome==row.genome]
 bed.write_text(''.join(f'{r.contig}\t{r.start-1}\t{r.end}\t{r.protein}\n'for r in dd.itertuples()))
 for season,s in [('winter',row.winter_sample),('hot',row.hot_sample)]:
  p=W/f'rebuild_20260729/picrust2_current_20260803/coverm_current/profiles_full_recompute/{s}/run_provenance.json';prov=json.loads(p.read_text());rr={'sample':s,'genome':row.genome,'domain':row.domain,'season':season,'bed':str(bed),'full_reads':prov['total_reads'],'R1':prov['R1_path'],'R2':prov['R2_path'],'R1_sha256_prior':prov['R1_sha256'],'R2_sha256_prior':prov['R2_sha256']}
  for mate in ['R1','R2']:
   try:
    path=Path(rr[mate]);rr[mate+'_bytes']=path.stat().st_size
    with open(path,'rb')as f:f.read(1)
    rr[mate+'_accessible']=True
   except Exception as e:rr[mate+'_accessible']=False;rr[mate+'_error']=str(e)
  sr.append(rr)
pd.DataFrame(sr).to_csv(O/'pilot_samples.tsv',sep='\t',index=False);selected.to_csv(O/'selected_pairs.tsv',sep='\t',index=False)
(O/'selection_summary.json').write_text(json.dumps({'job':os.environ['SLURM_JOB_ID'],'eligible_species_site_compartment_pairs':len(d),'pilot_pairs':selected.to_dict('records'),'selection':'highest minimum raw relative abundance across a matched August/February pair, separately Thermoplasmatota and Bacteria; HSP20 strict carrier','pilot':'first5000000pairedrecords, technical coverage gate only, no seasonal inference from read prefix','full_gate':'at least5xprojected full depth and at least0.5 observed pilot breadth for marker and background; real full-depth gate >=5xmedian and>=.9breadth; expand site replication only if feasible'},indent=2))
print(selected.to_string(index=False));print(pd.DataFrame(sr)[['sample','R1_accessible','R2_accessible']].to_string(index=False))
