from common import *
import pandas as pd,numpy as np,gzip,re,pickle,collections,time
R=N/'results/hs2';R.mkdir(exist_ok=True)
idx=pd.read_csv(T/'index/genome_index.tsv',sep='\t').query("universe=='soil_representative'")
idx=idx[idx.gtdb_classification.str.contains('p__Thermoplasmatota;')];assert len(idx)==237
markers=pd.read_csv(N/'results/marker_loci.tsv.gz',sep='\t');markers=markers[(markers.hmm_coverage>=.7)&(markers.partial.astype(str).str.zfill(2)=='00')];anchors=set(markers.protein)
genes=[];sources=[]
for r in idx.itertuples():
 p=Path(r.proteome);sources.append({'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size})
 with gzip.open(p,'rt')as f:
  for l in f:
   if not l.startswith('>'):continue
   z=l[1:].strip().split(' # ');pid=z[0];raw=pid.split('::')[1];key=r.source_sample+'::'+raw
   genes.append({'genome':r.genome_id,'protein':pid,'catalog_id':key,'contig':raw.rsplit('_',1)[0],'start':int(z[1]),'end':int(z[2]),'strand':int(z[3]),'partial':re.search(r'partial=(\d\d)',z[4]).group(1),'anchor':int(pid in anchors),'taxonomy':r.gtdb_classification})
d=pd.DataFrame(genes);assert d.protein.is_unique
exact_path=N/'results/hs2_exact/exact_locus_links.tsv.gz';exact=pd.read_csv(exact_path,sep='\t',keep_default_na=False).set_index('protein');assert set(d.protein)==set(exact.index);d['name_based_catalog_id']=d.catalog_id;d['catalog_id']=d.protein.map(exact.catalog_id);d['sequence_link_status']=d.protein.map(exact.status);sources.append({'path':str(exact_path),'sha256':sha(exact_path),'bytes':exact_path.stat().st_size});targets=set(x.encode()for x in d.catalog_id if x);B=W/'06_darkmatter/darkfam_basis_clean/results';p=B/'member_fam_final.tsv';mf={}
with open(p,'rb')as f:
 for l in f:
  member,fam=l.rstrip(b'\r\n').split(b'\t');assert member not in mf;mf[member]=fam
sources.append({'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size})
# Both columns verified against original family builder; map original loci through90%membership.
p=W/'04_protein_catalog/all_linclust90/all_cluster90.tsv';links={};h=hashlib.sha256();n=0
with open(p,'rb',buffering=8*1024*1024)as f:
 for l in f:
  h.update(l);rep,member=l.rstrip(b'\r\n').split(b'\t');n+=1
  if member in targets:
   fam=mf.get(rep)
   if fam:
    value=(rep.decode(),fam.decode())
    if member in links:assert links[member]==value
    links[member]=value
  if n%50000000==0:print('cluster_rows',n,'mapped_dark_genes',len(links),flush=True)
sources.append({'path':str(p),'sha256':h.hexdigest(),'bytes':p.stat().st_size})
d['cluster90']=d.catalog_id.map(lambda x:links.get(x.encode(),('',''))[0]);d['family']=d.catalog_id.map(lambda x:links.get(x.encode(),('',''))[1]);d.to_csv(R/'thermo_genes.tsv.gz',sep='\t',index=False);write(R/'mapping_sources.tsv',sources)
(R/'mapping_summary.json').write_text(json.dumps({'job':os.environ['SLURM_JOB_ID'],'genomes':len(idx),'genes':len(d),'unknown_family_loci':int((d.family!='').sum()),'unknown_families':d.loc[d.family!='','family'].nunique(),'anchor_loci':int(d.anchor.sum()),'known_anchor_overlap':int(((d.family!='')&(d.anchor==1)).sum()),'cluster90_rows':n,'linkage':'exact amino-acid sequence on same contig -> unique original source locus ->90%representative -> admitted darkfamily'},indent=2))
print('DONE HS2 mapping',len(links),flush=True)
