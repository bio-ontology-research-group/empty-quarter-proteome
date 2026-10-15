"""Heliorhodopsin anchor check and independent best-reference retinal-site projection."""
import csv,json,hashlib,os,subprocess,collections
from pathlib import Path
from Bio import AlignIO,SeqIO
W=Path('/ibex/scratch/projects/c2014/rob/metagenomics');D=W/'06_darkmatter/family_discovery_20261015_clean';C=D/'root_audit/f345_comparison';F=D/'f345/references';R=D/'root_audit/f3_reference_sensitivity';R.mkdir(exist_ok=True)
# Primary source Pushkarev 2018 Extended Data Fig1: 48C12 retinal lysine 241.
seqs=SeqIO.to_dict(SeqIO.parse(C/'Helio/combined.faa','fasta'));rid='Helio__curated_reference__AVZ43932.1';gp=SeqIO.read(F/'raw_ncbi/AVZ43932.1.gp','genbank');assert str(seqs[rid].seq)==str(gp.seq);assert str(gp.seq)[240]=='K'
with open(R/'helio.sto','w')as f:subprocess.run(['hmmalign','--outformat','Stockholm',str(F/'profiles/PF18761_Heliorhodopsin.hmm'),str(C/'Helio/combined.faa')],check=True,stdout=f)
aln=AlignIO.read(R/'helio.sto','stockholm');seqaln={r.id:str(r.seq)for r in aln};s=seqaln[rid];pos=0;col=None
for i,c in enumerate(s):
 if c not in '.-':pos+=1
 if pos==241:col=i;break
assert aln.column_annotations['reference_annotation'][col]not in '.- '
hr=[{'id':pid,'aligned_K241_site':s[col].upper(),'exact_K':int(s[col].upper()=='K')}for pid,s in seqaln.items()if'__EQ__'in pid]
# Reference profile coverage diagnostic (not an activity test). Missing positive controls exposes sensitivity limits.
refrows=[]
for fam,prof in [('Type1','PF01036_Bac_rhodopsin.hmm'),('Helio','PF18761_Heliorhodopsin.hmm'),('EctC','PF06339_Ectoine_synth.hmm'),('BetaCA','PF00484_Pro_CA.hmm'),('GammaCA','NF050158_gamma_CA.hmm')]:
 records=[r for r in SeqIO.parse(C/fam/'combined.faa','fasta')if'__curated_reference__'in r.id];SeqIO.write(records,R/(fam+'_refs.faa'),'fasta')
 with open(R/(fam+'_hmmsearch.log'),'w')as f:subprocess.run(['hmmsearch','--cpu','2','--domtblout',str(R/(fam+'_refs.domtbl')),str(F/'profiles'/prof),str(R/(fam+'_refs.faa'))],check=True,stdout=f)
 best={}
 for l in open(R/(fam+'_refs.domtbl')):
  if l.startswith('#'):continue
  x=l.split();sc=float(x[13]);cov=(int(x[16])-int(x[15])+1)/int(x[5]);rec={'id':x[0],'family':fam,'length':int(x[2]),'best_profile_coverage':cov,'domain_ievalue':float(x[12]),'score':sc}
  if x[0]not in best or sc>best[x[0]]['score']:best[x[0]]=rec
 for r in records:refrows.append(best.get(r.id,{'id':r.id,'family':fam,'length':len(r.seq),'best_profile_coverage':0,'domain_ievalue':None,'score':None}))
# Pairwise projection uses each target reference's actual annotated retinylidene position.
di=str(W/'../miniconda3/envs/eggnog/bin/diamond')
with open(R/'diamond.log','w')as log:
 subprocess.run([di,'makedb','--in',str(R/'Type1_refs.faa'),'--db',str(R/'Type1_refdb')],check=True,stdout=log,stderr=log)
 subprocess.run([di,'blastp','--query',str(C/'Type1/eq_queries.faa'),'--db',str(R/'Type1_refdb.dmnd'),'--out',str(R/'Type1_reference_alignments.tsv'),'--ultra-sensitive','--max-target-seqs','0','--evalue','1e-3','--threads','4','--outfmt','6','qseqid','sseqid','pident','qstart','qend','sstart','send','evalue','bitscore','qlen','slen','qseq','sseq'],check=True,stdout=log,stderr=log)
features={};calibration=[]
a=AlignIO.read(D/'root_audit/f345_anchor_alignment/Type1.sto','stockholm');aa={r.id:str(r.seq)for r in a};anchor=json.load(open(D/'root_audit/f345_anchor_alignment/summary.json'))['families']['Type1']['features'][0]['alignment_column']
for r in SeqIO.parse(R/'Type1_refs.faa','fasta'):
 acc=r.id.rsplit('__',1)[1];x=json.load(open(F/'raw_uniprot_json'/(acc+'.json')));ff=[f['location']['start']['value']for f in x.get('features',[])if'retinylidene'in f.get('description','')];assert len(ff)==1;features[r.id]=ff[0];native_position=sum(c not in '.-'for c in aa[r.id][:anchor+1]);calibration.append({'reference':acc,'annotated_retinal_position':ff[0],'HMM_projected_position':native_position,'match':native_position==ff[0] and aa[r.id][anchor]=='K'})
best={}
for l in open(R/'Type1_reference_alignments.tsv'):
 x=l.rstrip().split('\t');q,t,pid,qs,qe,ts,te,ev,bit,ql,tl,qseq,tseq=x
 if (int(qe)-int(qs)+1)/int(ql)<.7 or (int(te)-int(ts)+1)/int(tl)<.7:continue
 if q in best and float(bit)<=best[q]['bitscore']:continue
 spos=int(ts)-1;qpos=int(qs)-1;res='UNALIGNED';querypos=None
 # DIAMOND qseq/sseq are full sequences, not gapped traceback! Exact projection requires qseq_translated/sseq? Do not misuse them.
 best[q]={'id':q,'reference':t,'identity':float(pid),'bitscore':float(bit),'reference_retinal_position':features[t]}
# Use source pairwise aligner against the highest-bitscore qualifying curated hit; retain disagreement, no rescue by selecting another hit.
from Bio.Align import PairwiseAligner,substitution_matrices
aligner=PairwiseAligner();aligner.mode='global';aligner.substitution_matrix=substitution_matrices.load('BLOSUM62');aligner.open_gap_score=-10;aligner.extend_gap_score=-.5
qs=SeqIO.to_dict(SeqIO.parse(C/'Type1/eq_queries.faa','fasta'));ts=SeqIO.to_dict(SeqIO.parse(R/'Type1_refs.faa','fasta'));tr=[]
for q,rec in best.items():
 target=str(ts[rec['reference']].seq);query=str(qs[q].seq);al=aligner.align(target,query)[0];pos=rec['reference_retinal_position']-1;res='GAP';qp=None
 for (ta,tb),(qa,qb)in zip(*al.aligned):
  if ta<=pos<tb:qp=int(qa+pos-ta);res=query[qp];break
 rec.update({'pairwise_projected_residue':res,'query_position':qp+1 if qp is not None else None,'exact_K':int(res=='K'),'HMM_projected_residue':aa[q][anchor].upper(),'agreement_with_HMM':int(res==aa[q][anchor].upper())});tr.append(rec)
for name,rr in [('helio_sites.tsv',hr),('reference_profile_coverage.tsv',refrows),('type1_pairwise_sites.tsv',tr),('type1_reference_calibration.tsv',calibration)]:
 with open(R/name,'w')as f:w=csv.DictWriter(f,fieldnames=list(rr[0]),delimiter='\t');w.writeheader();w.writerows(rr)
(R/'summary.json').write_text(json.dumps({'job_id':os.environ['SLURM_JOB_ID'],'Helio':{'reference':'AVZ43932.1','reference_length':len(gp.seq),'reference_source_sha256':hashlib.sha256((F/'raw_ncbi/AVZ43932.1.gp').read_bytes()).hexdigest(),'site':241,'source':'Pushkarev2018 doi:10.1038/s41586-018-0225-9 Extended Data Fig1','EQ_complete':len(hr),'exact_K':sum(r['exact_K']for r in hr)},'Type1':{'curated_reference_site_calibration':all(r['match']for r in calibration),'references_calibrated':len(calibration),'qualified_best_reference_queries':len(tr),'exact_K_pairwise':sum(r['exact_K']for r in tr),'same_as_HMM':sum(r['agreement_with_HMM']for r in tr)},'reference_profile_under70':[r for r in refrows if r['best_profile_coverage']<.7],'limits':['Site identity is plausibility, not retinal binding or biochemical function.','Pairwise global alignment is an alternative projection, not independent functional validation.','Profile-limited references retained as sensitivity diagnostics, not automatically assigned no function.']},indent=2))
