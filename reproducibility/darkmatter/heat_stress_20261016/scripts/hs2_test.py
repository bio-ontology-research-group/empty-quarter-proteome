from common import *
import pandas as pd,numpy as np
from scipy.stats import wilcoxon
R=N/'results';O=R/'hs2';rng=np.random.default_rng(20261016);B=9999
def bh(p):
 p=np.array(p);o=np.argsort(p);a=np.empty(len(p));a[o]=np.minimum(1,np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return a
d=pd.read_csv(O/'thermo_genes.tsv.gz',sep='\t',keep_default_na=False);d['partial']=d.partial.astype(str).str.zfill(2);d=d.sort_values(['genome','contig','start','protein']).reset_index(drop=True)
d['pos']=d.groupby(['genome','contig']).cumcount();d['n_contig']=d.groupby(['genome','contig']).protein.transform('size');d['eligible']=(d.family!='')&(d.partial=='00')&(d.anchor==0)&(d.pos>=2)&(d.pos<d.n_contig-2)
# One prespecified first eligible locus per family per species representative avoids paralog pseudoreplication.
e=d[d.eligible].drop_duplicates(['family','genome'],keep='first').copy();ng=e.groupby('family').genome.nunique();fams=sorted(ng[ng>=5].index);e=e[e.family.isin(fams)];fi={f:i for i,f in enumerate(fams)}
null=np.zeros((len(fams),B),dtype=np.uint16);obs=np.zeros(len(fams),int);expected=np.zeros(len(fams));detail=[]
blocks=d.groupby(['genome','contig'],sort=False)
for (g,c),selected in e.groupby(['genome','contig'],sort=False):
 block=blocks.get_group((g,c))
 flags=block.anchor.to_numpy(bool);n=len(flags);near=np.zeros(n,bool)
 for delta in [-2,-1,1,2]:near|=np.roll(flags,delta)
 shifts=rng.integers(0,n,B)
 for r in selected.itertuples():
  k=fi[r.family];actual=bool(near[r.pos]);obs[k]+=actual;expected[k]+=near.mean();null[k]+=near[(r.pos-shifts)%n]
  anchors=block.iloc[max(0,r.pos-2):r.pos+3];anchors=anchors[anchors.anchor==1]
  detail.append({'family':r.family,'genome':g,'protein':r.protein,'contig':c,'pos':r.pos,'n_contig':n,'near_anchor':int(actual),'nearby_anchors':';'.join(anchors.protein),'null_near_fraction':near.mean(),'taxonomy':r.taxonomy})
p=((null>=obs[:,None]).sum(1)+1)/(B+1);q=bh(p);out=pd.DataFrame({'family':fams,'genomes':[ng[f]for f in fams],'near_genomes':obs,'expected_near_genomes':expected,'p':p,'q':q});out['enrichment']=out.near_genomes/out.expected_near_genomes.replace(0,np.nan);out=out.sort_values(['q','p','near_genomes'],ascending=[True,True,False]);out.to_csv(O/'neighbourhood_tests.tsv',sep='\t',index=False);pd.DataFrame(detail).to_csv(O/'tested_loci.tsv.gz',sep='\t',index=False);np.savez_compressed(O/'null_counts.npz',families=np.array(fams),counts=null)
# Candidate ecological follow-up is an exploratory second stage, selection solely on neighbourhood evidence.
selected=out[out.q<=.05].head(20)
if selected.empty:selected=out.head(10)
chosen=set(selected.family);z=np.load(R/'abundance.npz');g=z['genomes'];samples=z['samples'];raw=z['raw'];ok=raw.sum(1)>0;raw=raw[ok];samples=samples[ok];a=raw/raw.sum(1)[:,None]
meta=pd.read_csv(R/'metadata.tsv',sep='\t',index_col=0,dtype={'site':str}).loc[samples];idx=pd.read_csv(R/'genome_index.tsv',sep='\t',index_col=0).loc[g];thermo=idx.gtdb_classification.str.contains('p__Thermoplasmatota;').to_numpy();mass=a[:,thermo].sum(1)
rr=[];sampleout=[]
for fam in sorted(chosen):
 hosts=set(d.loc[d.family==fam,'genome']);v=a[:,np.isin(g,list(hosts))].sum(1);values={'community_fraction':v,'within_Thermoplasmatota_fraction':np.divide(v,mass,out=np.full(len(v),np.nan),where=mass>0)}
 for norm,y in values.items():
  df=meta[['site','trip','compartment']].copy();df['value']=y;df['sample']=df.index;df['family']=fam;df['normalisation']=norm;sampleout.append(df)
  cells=df.groupby(['site','trip','compartment']).value.mean()
  for hot in ['trip2','trip4']:
   for co in ['Surface','Deep','PlantRoot']:
    sites=sorted({s for s,tr,cc in cells.index if tr==hot and cc==co and(s,'trip3',co)in cells.index});vals=np.array([cells.loc[(s,hot,co)]-cells.loc[(s,'trip3',co)]for s in sites]);vals=vals[np.isfinite(vals)]
    if not len(vals):continue
    boot=rng.choice(vals,(2000,len(vals))).mean(1);pv=wilcoxon(vals).pvalue if len(vals)>=4 and np.any(vals)else(1 if len(vals)>=4 else np.nan)
    rr.append({'family':fam,'normalisation':norm,'contrast':hot+'-trip3','compartment':co,'sites':len(vals),'mean_difference_pp':100*vals.mean(),'lo_pp':100*np.quantile(boot,.025),'hi_pp':100*np.quantile(boot,.975),'p':pv})
rr=pd.DataFrame(rr);rr['q']=np.nan;mask=rr.p.notna();rr.loc[mask,'q']=bh(rr.loc[mask,'p']);rr.to_csv(O/'candidate_ecology.tsv',sep='\t',index=False);pd.concat(sampleout).to_csv(O/'candidate_sample_metrics.tsv.gz',sep='\t',index=False)
(O/'summary.json').write_text(json.dumps({'job':os.environ['SLURM_JOB_ID'],'families_tested':len(fams),'eligible_family_genome_loci':len(e),'q05':int((out.q<=.05).sum()),'q10':int((out.q<=.1).sum()),'null':'9999 independent circular rotations of anchor positions within each contig; preserves within-contig spacing and genome/lineage membership; excludes target anchors, complete target genes >=2positions from boundary; firsteligiblelocus per family/genome','limitations':['association is genomic context, not function','shared ancestry can conserve neighbourhoods, no independent-origin claim','rotations do not match gene length','ecological follow-up uses same observations and is exploratory','carrier abundance is static gene-content projection, not host-independent gene coverage']},indent=2))
print(out.head(15).to_string(index=False));print('DONE HS2')
