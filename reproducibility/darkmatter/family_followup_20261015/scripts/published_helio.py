from common import *
import zipfile,json,subprocess,collections,numpy as np
from Bio import SeqIO,AlignIO
R=N/'results/published_helio';R.mkdir(parents=True,exist_ok=True);I=N/'inputs/bulzu2021'
with zipfile.ZipFile(I/'Alignment_files.zip')as z:
    names=[s for s in z.namelist()if s.endswith('HeR_phylogenetic_tree_rep_seq_alignment.fas')];assert len(names)==1,names
    (R/'published_alignment.fas').write_bytes(z.read(names[0]))
published=[];mapping=[]
for i,r in enumerate(SeqIO.parse(R/'published_alignment.fas','fasta')):
    s=str(r.seq).replace('-','').replace('.','').upper();assert set(s)<=set('ACDEFGHIKLMNPQRSTVWYBXZ'),(r.id,set(s));pid='BULZU_'+str(i);published.append((pid,s));mapping.append({'id':pid,'original_header':r.description,'length':len(s),'sha256':hashlib.sha256(s.encode()).hexdigest()})
write(R/'reference_mapping.tsv',mapping)
with open(R/'published_ungapped.faa','w')as f:
    for pid,s in published:f.write('>'+pid+'\n'+s+'\n')
profile=D/'f345/references/profiles/PF18761_Heliorhodopsin.hmm'
with open(R/'hmmsearch.log','w')as f:subprocess.run(['hmmsearch','--cpu','4','--domtblout',str(R/'domains.tsv'),str(profile),str(R/'published_ungapped.faa')],check=True,stdout=f)
keep=set()
for l in open(R/'domains.tsv'):
    if l.startswith('#'):continue
    x=l.split()
    if float(x[12])<=1e-5 and(int(x[16])-int(x[15])+1)/int(x[5])>=.7:keep.add(x[0])
ref='Helio__curated_reference__AVZ43932.1';anchor=next(r for r in SeqIO.parse(D/'root_audit/f345_comparison/Helio/combined.faa','fasta')if r.id==ref)
with open(R/'profile_pass.faa','w')as f:
    for pid,s in published:
        if pid in keep:f.write('>'+pid+'\n'+s+'\n')
    f.write('>'+ref+'\n'+str(anchor.seq)+'\n')
with open(R/'profile.sto','w')as f:subprocess.run(['hmmalign','--outformat','Stockholm',str(profile),str(R/'profile_pass.faa')],check=True,stdout=f)
al=AlignIO.read(R/'profile.sto','stockholm');aa={r.id:str(r.seq)for r in al};pos=0
for col,c in enumerate(aa[ref]):
    if c not in '.-':pos+=1
    if pos==241:break
assert aa[ref][col]=='K';kkeep={pid for pid in keep if aa[pid][col].upper()=='K'}
with open(R/'site_conserving_references.faa','w')as f:
    for pid,s in published:
        if pid in kkeep:f.write('>'+pid+'\n'+s+'\n')
eqmeta=[r for r in tab(D/'root_audit/f3_site_conserving/Helio/metadata.tsv')if r['source']=='EQ'];eqids={r['id']for r in eqmeta};eqseq={r.id:str(r.seq)for r in SeqIO.parse(D/'root_audit/f345_comparison/Helio/combined.faa','fasta')if r.id in eqids}
with open(R/'EQ.faa','w')as f:
    for pid,s in eqseq.items():f.write('>'+pid+'\n'+s+'\n')
diamond='/ibex/scratch/projects/c2014/rob/miniconda3/envs/eggnog/bin/diamond'
with open(R/'diamond.log','w')as log:
    subprocess.run([diamond,'makedb','--in',str(R/'site_conserving_references.faa'),'--db',str(R/'reference')],check=True,stdout=log,stderr=log)
    subprocess.run([diamond,'blastp','--query',str(R/'EQ.faa'),'--db',str(R/'reference'),'--out',str(R/'hits.tsv'),'--ultra-sensitive','--max-target-seqs','0','--evalue','1e-3','--threads','4','--outfmt','6','qseqid','sseqid','pident','qstart','qend','sstart','send','evalue','bitscore','qlen','slen'],check=True,stdout=log,stderr=log)
best={}
for l in open(R/'hits.tsv'):
    q,t,pid,qs,qe,ts,te,ev,bs,ql,tl=l.split();qcov=(int(qe)-int(qs)+1)/int(ql);tcov=(int(te)-int(ts)+1)/int(tl)
    if qcov<.7 or tcov<.7:continue
    r={'id':q,'target':t,'identity':float(pid),'bitscore':float(bs),'query_coverage':qcov,'target_coverage':tcov}
    if q not in best or r['bitscore']>best[q]['bitscore']:best[q]=r
rows=[]
for r in eqmeta:rows.append({**r,**best.get(r['id'],{'target':'','identity':'','bitscore':'','query_coverage':'','target_coverage':''})})
write(R/'closest_published_hits.tsv',rows)
(R/'summary.json').write_text(json.dumps({'job':os.environ['SLURM_JOB_ID'],'dataset':'10.6084/m9.figshare.13286486.v3','input_alignment_sequences':len(published),'profile_pass':len(keep),'site_conserving_pass':len(kkeep),'EQ_queries':len(eqmeta),'EQ_coverage_passing_hits':len(best),'EQ_median_identity':float(np.median([r['identity']for r in best.values()]))if best else None,'EQ_below50_identity':sum(r['identity']<50 for r in best.values()),'by_phylum':{p:{'queries':len(rr),'hits':len(vals),'median_identity':float(np.median(vals))if vals else None}for p in sorted({r['phylum']for r in eqmeta})for rr in [[r for r in eqmeta if r['phylum']==p]]for vals in [[best[r['id']]['identity']for r in rr if r['id']in best]]},'limits':['Published phylogeny representatives, not an exhaustive current database','Alignment-derived reference sequences can have trimmed ends','Gene completeness metadata is unavailable for this literature reference panel','No model training-data novelty or activity inference']},indent=2))
