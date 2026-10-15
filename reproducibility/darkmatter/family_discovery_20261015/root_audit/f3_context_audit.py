"""Independent raw-marker, source-sequence and window audit; prespecified matched sensitivity."""
import csv,collections,gzip,hashlib,json,os
from pathlib import Path
from scipy.stats import binomtest
W=Path('/ibex/scratch/projects/c2014/rob/metagenomics');D=W/'06_darkmatter/family_discovery_20261015_clean';F=D/'f3_context';P=F/'results';R=D/'root_audit/f3_context_audit';R.mkdir(exist_ok=True)
def tab(p):return list(csv.DictReader(open(p),delimiter='\t'))
def fasta(handle):
 h=None;s=[]
 for l in handle:
  if l.startswith('>'):
   if h:yield h,''.join(s)
   h=l[1:].strip();s=[]
  else:s.append(l.strip())
 if h:yield h,''.join(s)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
rr=tab(P/'full_loci_manifest.tsv');reported={r['rhodopsin_id']:r for r in tab(P/'full_locus_marker_summary.tsv')};profile=F/'profiles/f3_sensory_context.hmm';ga={};name=None
for l in open(profile):
 if l.startswith('NAME'):name=l.split()[1]
 if l.startswith('GA'):ga[name]=tuple(float(x.strip(';'))for x in l.split()[1:3])
markers=collections.defaultdict(set);nraw=0
for l in open(P/'full_domtbl_sensory_context.tsv'):
 if l.startswith('#'):continue
 x=l.split()
 if not x:continue
 assert float(x[7])>=ga[x[3]][0]-.01 and float(x[13])>=ga[x[3]][1]-.01,(x[0],x[3],'GA failure');nraw+=1
 if x[3]not in {'Bac_rhodopsin','Heliorhodopsin'}:markers[x[0]].add(x[3])
seqs={h.split()[0]:s for h,s in fasta(open(P/'full_combined_search_proteins.faa'))};original={};coords={};bycontig=collections.defaultdict(list);source_hashes={}
for mag in sorted({r['mag']for r in rr}):
 p=W/'06_darkmatter/drep_rep_proteomes'/(mag+'.faa.gz');source_hashes[str(p)]=sha(p)
 for h,s in fasta(gzip.open(p,'rt')):
  fields=h.split('#');pid=fields[0].strip();ctg=pid.split('::')[1].rsplit('_',1)[0];co=(mag,ctg,int(fields[1]),int(fields[2]),fields[3].strip());coords[pid]=co;bycontig[co[:2]].append(pid)
  if pid in seqs:assert seqs[pid]==s,(pid,'source mismatch');original[pid]=s
assert set(original)==set(seqs)
results=[]
for r in rr:
 pid=r['rhodopsin_id'];cp=r['control_gene_id'];rho=coords[pid];assert rho[2]==int(r['rho_start'])and rho[3]==int(r['rho_stop']);a=max(1,rho[2]-5000);b=rho[3]+5000;neighbors={p for p in bycontig[rho[:2]]if p!=pid and coords[p][2]<=b and coords[p][3]>=a};assert neighbors==set(filter(None,r['neighbor_ids'].split(';')))
 ctrl=set()
 if cp:
  co=coords[cp];a=max(1,co[2]-5000);b=co[3]+5000;ctrl={p for p in bycontig[co[:2]]if p!=cp and coords[p][2]<=b and coords[p][3]>=a};assert ctrl==set(filter(None,r['control_neighbor_ids'].split(';')))
  if r['control_scope']=='same_contig':assert rho[:2]==co[:2];assert co[2]-5000>rho[3]+5000 or co[3]+5000<rho[2]-5000
 assert len(neighbors)==int(r['n_neighbors'])and len(ctrl)==int(r['control_n_neighbors'])
 rv=int(any(markers[p]for p in neighbors));cv=int(any(markers[p]for p in ctrl));tr=int(any(markers[p]&{'MCPsignal','HAMP','CheW'}for p in neighbors));assert rv==int(reported[pid]['rho_any_sensory'])and cv==int(reported[pid]['ctrl_any_sensory'])
 results.append({'id':pid,'family':r['assigned_class'],'mag':r['mag'],'rho_any':rv,'ctrl_any':cv,'rho_transducer':tr,'same_contig':int(r['control_scope']=='same_contig'),'rho_neighbors':len(neighbors),'ctrl_neighbors':len(ctrl),'eligible':int(r['control_scope']=='same_contig'and abs(len(neighbors)-len(ctrl))<=2),'taxonomy':r['gtdbtk_taxonomy']})
summary={};retained=[]
for fam in sorted({r['family']for r in results}):
 allr=[r for r in results if r['family']==fam];eligible=[r for r in allr if r['eligible']];seen=set();selected=[]
 for r in sorted(eligible,key=lambda r:hashlib.sha256(r['id'].encode()).hexdigest()):
  if r['mag']not in seen:seen.add(r['mag']);selected.append(r)
 retained+=selected;b=sum(r['rho_any']==1 and r['ctrl_any']==0 for r in selected);c=sum(r['rho_any']==0 and r['ctrl_any']==1 for r in selected);n=len(selected);p=binomtest(b,b+c,.5).pvalue if b+c else 1.0
 summary[fam]={'all_loci':len(allr),'raw_any_signaling':sum(r['rho_any']for r in allr),'raw_control_any':sum(r['ctrl_any']for r in allr),'raw_transducer':sum(r['rho_transducer']for r in allr),'same_contig_pairs':sum(r['same_contig']for r in allr),'neighbor_matched_pairs_before_MAG_dedup':len(eligible),'primary_one_per_MAG_pairs':n,'primary_rho_positive':sum(r['rho_any']for r in selected),'primary_control_positive':sum(r['ctrl_any']for r in selected),'discordant_rho_only':b,'discordant_control_only':c,'paired_difference_percentage_points':100*(b-c)/n if n else None,'p_exact_paired':p,'transducer_phyla':dict(collections.Counter(next((x[3:]for x in r['taxonomy'].split(';')if x.startswith('p__')),'Unknown')for r in allr if r['rho_transducer']))}
last=0
for k,(fam,v)in enumerate(sorted(summary.items(),key=lambda kv:kv[1]['p_exact_paired'])):last=max(last,min(1,v['p_exact_paired']*(2-k)));v['p_Holm_two_classes']=last
for name,rows in [('all_loci.tsv',results),('primary_matched_loci.tsv',retained)]:
 with open(R/name,'w')as f:
  wr=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');wr.writeheader();wr.writerows(rows)
(R/'audit.json').write_text(json.dumps({'job_id':os.environ['SLURM_JOB_ID'],'raw_GA_passing_domains':nraw,'original_sequences_verified':len(original),'original_MAGs':len(source_hashes),'source_hashes':source_hashes,'profile_sha256':sha(profile),'families':summary,'limits':['Same-contig +/-2-neighbor-count matching and one SHA-ranked locus per MAG/class fixed before inspecting outcomes.','Generic signaling-domain proximity does not prove rhodopsin sensory function.','Gene-coordinate extents are lower-bound assembled-contig extents, not actual contig lengths.','Absence of nearby markers does not identify proteins as ion pumps or refute sensing.','No external-cohort context comparison or functional experiment; Type1 sensory architecture may be established haloarchaeal biology.']},indent=2))
print(json.dumps(summary,indent=2))
