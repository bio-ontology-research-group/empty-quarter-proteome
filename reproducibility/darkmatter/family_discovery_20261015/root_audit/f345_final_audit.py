"""Independent audit of raw family gates, original sequences, coordinates and EctC context."""
import csv,hashlib,json,os,re,collections
from pathlib import Path
W=Path('/ibex/scratch/projects/c2014/rob/metagenomics');D=W/'06_darkmatter/family_discovery_20261015_clean';F=D/'f345';R=D/'root_audit/f345_final_audit';R.mkdir(exist_ok=True)
def tab(p):return list(csv.DictReader(open(p),delimiter='\t'))
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb')as f:
  while b:=f.read(4*1024*1024):h.update(b)
 return h.hexdigest()
def fasta(p):
 h=None;s=[]
 for l in open(p):
  if l.startswith('>'):
   if h:yield h,''.join(s).rstrip('*')
   h=l[1:].strip();s=[]
  else:s.append(l.strip())
 if h:yield h,''.join(s).rstrip('*')
def dom(p):
 for l in open(p):
  if l.startswith('#'):continue
  x=l.split()
  if len(x)<22:continue
  yield x[0],x[4].split('.')[0],float(x[13]),float(x[12]),(int(x[16])-int(x[15])+1)/int(x[5]),int(x[17]),int(x[18])
tracks={'f3':('f3_rhodopsin',{'PF01036','PF18761','NF038020'}),'f4':('f4_ectc',{'PF06339','NF009806'}),'f5':('f5_carbonic_anhydrase',{'NF050158','NF040597','PF00484','TIGR02701'})};sources=['EQ_MAG','Environmental_SMAG','Environmental_GEM'];rows=[];summary={};inputs={};markers={s:collections.defaultdict(dict)for s in sources}
for pre,(directory,targets)in tracks.items():
 P=F/'results/full_v1'/directory;rr=tab(P/f'{pre}_candidates_summary.tsv');rows+=rr;result={}
 for src in sources:
  raw=P/f'{pre}_domtbl_{src}.tsv';inputs[str(raw)]=sha(raw);targetids={v[0]for v in dom(raw)if v[1]in targets};best=collections.defaultdict(dict)
  for v in dom(raw):
   pid,acc,score,ev,cov,a,b=v
   if pre=='f4'and acc in {'TIGR02406','TIGR02407','TIGR02408'}and ev<=1e-5 and cov>=.5:
    role={'TIGR02406':'EctA','TIGR02407':'EctB','TIGR02408':'EctD'}[acc];markers[src][pid][role]=v
   if pid in targetids and(acc not in best[pid]or score>best[pid][acc][2]):best[pid][acc]=v
  def top(by,accs):return max((v for a,v in by.items()if a in accs),key=lambda v:v[2],default=None)
  def score(v):return v[2]if v else 0
  accepted={}
  for pid,by in best.items():
   if pre=='f3':
    t=top(by,{'PF01036'});h=top(by,{'PF18761','NF038020'});g=top(by,{'PF00001','PF00002'});b=t if score(t)>=score(h)else h
    cls='Type-1_Microbial_Rhodopsin'if score(t)>=score(h)else'Heliorhodopsin_Type-3'
    ok=b and b[3]<=1e-5 and b[4]>=.7 and b[2]-score(g)>=20
   elif pre=='f4':
    b=top(by,targets);c=top(by,{'PF00190','PF07883','PF02678'});cls='Ectoine_Synthase_EctC';ok=b and b[3]<=1e-10 and b[4]>=.7 and b[2]-score(c)>=25
   else:
    g=top(by,{'NF050158','NF040597'});bta=top(by,{'PF00484'});c=top(by,{'TIGR02701'});mx=max(score(g),score(bta),score(c));ctrl=top(by,{'PF00132','PF00583'})
    if c and mx==score(c):b=c;cls='Carboxysome_Shell_CA_CsoSCA';extra=True
    elif bta and mx==score(bta):b=bta;cls='Beta-Carbonic_Anhydrase';extra=b[2]>=30
    else:b=g;cls='Gamma-Carbonic_Anhydrase';extra=b and b[2]-score(ctrl)>=20
    ok=b and b[3]<=1e-7 and b[4]>=.7 and extra
   if ok:accepted[pid]=cls
  reported={r['protein_id']:r['assigned_class']for r in rr if r['cohort']==src}
  assert accepted==reported,(pre,src,len(accepted),len(reported),list(set(accepted)^set(reported))[:10]);result[src]={'raw_gate_recount':len(accepted),'complete':sum(r['cohort']==src and r['prodigal_partial']=='00'for r in rr),'classes':dict(collections.Counter(accepted.values()))}
 summary[pre]=result
print('Raw gates closed',flush=True)
# Independently stream original cohort sequences and coordinates, retaining only candidate/marker records.
H=W/'06_darkmatter/p7_analyses_clean/h3_references/results';paths={'EQ_MAG':D/'root_audit/f1_full/cohort.faa','Environmental_SMAG':H/'smag_pool_canonical.faa','Environmental_GEM':H/'gem_harmonized_pool_canonical_seed20261014.faa'};coordpaths={'Environmental_SMAG':H/'smag_protein_coordinate_mapping.tsv','Environmental_GEM':H/'gem_harmonized_protein_coordinate_mapping_seed20261014.tsv'};coords={};sourceclosure={}
for src,p in paths.items():
 by={r['protein_id']:r for r in rows if r['cohort']==src};assert len(by)==sum(r['cohort']==src for r in rows),'cross track repeated ID';wanted=set(by)|set(markers[src]);found={};count=0
 for h,s in fasta(p):
  count+=1;pid=h.split()[0]
  if pid not in wanted:continue
  assert pid not in found;digest=hashlib.sha256(s.encode()).hexdigest();found[pid]=(digest,len(s))
  if pid in by:
   r=by[pid];assert digest==r['seq_sha256'] and len(s)==int(r['length_aa']),(src,pid,'source sequence mismatch')
  if src=='EQ_MAG':
   mag,rawid=pid.split('::');fields=h.split('#');start=int(fields[1]);end=int(fields[2]);strand=fields[3].strip();partial=re.search(r'partial=(\d\d)',h).group(1);ctg=rawid.rsplit('_',1)[0];coords[(src,pid)]=(mag,ctg,start,end,strand,partial)
 assert set(found)==wanted,(src,'missing original records',list(wanted-set(found))[:5])
 if src!='EQ_MAG':
  cp=coordpaths[src];inputs[str(cp)]=sha(cp)
  for r in csv.DictReader(open(cp),delimiter='\t'):
   pid=r['canonical_seq_id']
   if pid not in wanted:continue
   assert found[pid]==(r['canonical_seq_sha256'],int(r['length_aa'])),(src,pid,'mapping SHA mismatch');key=(src,pid);assert key not in coords
   coords[key]=(r.get('mag_id',r.get('genome_id')),r['contig_id'],int(r['start_coord']),int(r['end_coord']),r['strand'],r['partial'])
 for pid,r in by.items():
  co=coords[(src,pid)];assert co[0]==r['mag_or_genome'] and co[1]==r['contig'] and co[2]==int(r['start']) and co[3]==int(r['stop']) and co[4]==r['strand'] and co[5]==r['prodigal_partial'],(src,pid,co,r['contig'],r['strand']);assert (co[5]=='00')==(r['is_complete_gene']=='1')
 sourceclosure[src]={'source_fasta':str(p),'source_fasta_sha256':sha(p),'proteins_streamed':count,'candidates_verified':len(by),'candidate_and_marker_coordinates':len(wanted)};print('Source closure',src,len(by),flush=True)
# Context independently recalculate nearest non-self role hit with correct intergenic distance.
context=[]
def dist(a,b):return max(0,max(a[2],b[2])-min(a[3],b[3])-1)
for r in rows:
 if r['assigned_class']!='Ectoine_Synthase_EctC':continue
 src=r['cohort'];pid=r['protein_id'];co=coords[(src,pid)];flags={};same={}
 for role in ['EctA','EctB','EctD']:
  pool=[(mp,coords[(src,mp)])for mp,v in markers[src].items()if role in v and mp!=pid and coords[(src,mp)][:2]==co[:2]];nearest=min((dist(co,mc)for mp,mc in pool),default=None);expected=nearest is not None and nearest<=5000;assert int(r['has_nearby_ect'+role[-1]])==int(expected),(pid,role,nearest,r['nearby_ect'+role[-1]+'_dist_bp'])
  flags[role]=expected;same[role]=any(dist(co,mc)<=5000 and mc[4]==co[4]for mp,mc in pool)
  chosen=r['nearby_ect'+role[-1]+'_id']
  if chosen!='None':assert chosen!=pid and chosen in markers[src] and role in markers[src][chosen] and dist(co,coords[(src,chosen)])==nearest
 fusion=r['is_same_orf_ectcd_fusion']=='1'
 if fusion:
  a,b=map(int,re.split(r'\.\.|-',r['same_orf_ectc_domain_coords']));c,d=map(int,re.split(r'\.\.|-',r['same_orf_ectd_domain_coords']));assert b<c or d<a;assert 'EctD'in markers[src][pid]
 context.append({'id':pid,'source':src,'complete':int(r['prodigal_partial']=='00'),'ABC_neighborhood':int(flags['EctA']and flags['EctB']),'ABC_same_strand':int(same['EctA']and same['EctB']),'ABCD_neighborhood':int(flags['EctA']and flags['EctB']and(flags['EctD']or fusion)),'same_orf_EctCD':int(fusion),'reported_context':r['pathway_organization']})
with open(R/'ectc_context_independent.tsv','w')as f:
 w=csv.DictWriter(f,fieldnames=list(context[0]),delimiter='\t');w.writeheader();w.writerows(context)
counts={}
for src in sources:
 rr=[r for r in context if r['source']==src and r['complete']];counts[src]={'complete_EctC':len(rr),**{key:sum(r[key]for r in rr)for key in ['ABC_neighborhood','ABC_same_strand','ABCD_neighborhood','same_orf_EctCD']}}
for p in (F/'results/full_v1/frozen_code').glob('*.py'):inputs[str(p)]=sha(p)
(R/'audit.json').write_text(json.dumps({'job_id':os.environ['SLURM_JOB_ID'],'raw_gates':summary,'source_closure':sourceclosure,'ectc_complete_context':counts,'input_checksums':inputs,'limits':['Complete means Prodigal partial=00, not complete genome or experimentally validated ORF.','Neighborhoods within 5kb are not demonstrated operons; opposite strands separately checked.','Profile thresholds are operational and class selection can miss fusions/multidomain proteins.','No catalytic activity, sensory function or exceptional EQ novelty follows from census alone.']},indent=2))
print('Independent audit complete',flush=True)
