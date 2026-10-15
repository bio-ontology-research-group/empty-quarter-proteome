from common import *
import collections,re,json,gzip
from Bio import SeqIO,AlignIO
from Bio.Seq import Seq
R=N/'results/f4_integrity';R.mkdir(parents=True,exist_ok=True)
F=D/'f345/results/full_v1/f4_ectc';H=W/'06_darkmatter/p7_analyses_clean/h3_references'
candidates=tab(F/'f4_candidates_summary.tsv')
anchor=json.load(open(D/'root_audit/f345_anchor_alignment/summary.json'))['families']['EctC']['features']
assert all(x['profile_match_column']for x in anchor)
sites={r.id.split('__',2)[2]:''.join(str(r.seq)[x['alignment_column']].upper()for x in anchor)=='WEYH'for r in AlignIO.read(D/'root_audit/f345_anchor_alignment/EctC.sto','stockholm')}
proteins={r.id:str(r.seq).rstrip('*')for r in SeqIO.parse(F/'f4_all_candidates.faa','fasta')}
rows=[];coords={};markers={};inputfiles=[]
for src in ['EQ_MAG','Environmental_SMAG','Environmental_GEM']:
    marks=collections.defaultdict(dict);raw=F/f'f4_domtbl_{src}.tsv';inputfiles.append(raw)
    for l in open(raw):
        if l.startswith('#'):continue
        x=l.split();acc=x[4].split('.')[0]
        if acc in {'TIGR02406','TIGR02407','TIGR02408'}and float(x[12])<=1e-5 and (int(x[16])-int(x[15])+1)/int(x[5])>=.5:
            role={'TIGR02406':'A','TIGR02407':'B','TIGR02408':'D'}[acc];rec=(float(x[13]),(int(x[16])-int(x[15])+1)/int(x[5]))
            if role not in marks[x[0]]or rec[0]>marks[x[0]][role][0]:marks[x[0]][role]=rec
    cr=[r for r in candidates if r['cohort']==src];wanted=set(marks)|{r['protein_id']for r in cr};co={}
    if src=='EQ_MAG':
        path=D/'root_audit/f1_full/cohort.faa';inputfiles.append(path)
        # Header-only streaming avoids storing the 9GB cohort.
        for l in open(path):
            if not l.startswith('>'):continue
            pid=l[1:].split()[0]
            if pid not in wanted:continue
            fields=l.split('#');mag,rawid=pid.split('::');co[pid]={'mag':mag,'contig':rawid.rsplit('_',1)[0],'start':int(fields[1]),'end':int(fields[2]),'strand':int(fields[3]),'partial':re.search(r'partial=(\d\d)',l).group(1)}
    else:
        cp=H/'results'/('smag_protein_coordinate_mapping.tsv'if src=='Environmental_SMAG'else'gem_harmonized_protein_coordinate_mapping_seed20261014.tsv');inputfiles.append(cp)
        for r in csv.DictReader(open(cp),delimiter='\t'):
            if r['canonical_seq_id']in wanted:co[r['canonical_seq_id']]={'mag':r.get('mag_id',r.get('genome_id')),'contig':r['contig_id'],'start':int(r['start_coord']),'end':int(r['end_coord']),'strand':{'-':-1,'+':1,'-1':-1,'1':1}[r['strand']],'partial':r['partial']}
    assert wanted==set(co),(src,len(wanted),len(co));coords[src]=co;markers[src]=marks
    byctg=collections.defaultdict(list)
    for pid in marks:byctg[(co[pid]['mag'],co[pid]['contig'])].append(pid)
    for r in cr:
        if r['prodigal_partial']!='00':continue
        pid=r['protein_id'];c=co[pid];rec={'id':pid,'source':src,'mag':c['mag'],'contig':c['contig'],'start':c['start'],'end':c['end'],'strand':c['strand'],'WEYH':int(sites[pid]),'fusion_CD':int(r['is_same_orf_ectcd_fusion'])}
        chosen={};partial_near=False
        for role in ['A','B','D']:
            eligible=[]
            for q in byctg[(c['mag'],c['contig'])]:
                if q==pid or role not in marks[q]:continue
                a=co[q];dist=max(0,max(a['start'],c['start'])-min(a['end'],c['end'])-1)
                if dist<=5000:
                    partial_near|=a['partial']!='00'
                    if a['strand']==c['strand']and a['partial']=='00':eligible.append((dist,q))
            chosen[role]=min(eligible)[1]if eligible else ''
            rec['complete_'+role]=chosen[role]
        rec['strict_ABC']=int(rec['WEYH']and bool(chosen['A'])and bool(chosen['B']))
        rec['strict_ABC_cov70']=int(rec['strict_ABC']and all(marks[chosen[k]][k][1]>=.7 for k in ['A','B']))
        rec['partial_marker_near']=int(partial_near)
        ordered=[(c['start'],'C')]+[(co[q]['start'],k)for k,q in chosen.items()if q]
        rec['gene_order']='-'.join(k for _,k in sorted(ordered,reverse=c['strand']==-1));rows.append(rec)
    print('Coordinates',src,len(co),flush=True)
# Actual genome FASTA lengths, no use of last predicted gene as contig length.
smagfiles={p.stem:p for p in (H/'genomes/smag_extracted').rglob('*')if p.suffix in {'.fa','.fna','.fasta'}}
for src in coords:
    rr=[r for r in rows if r['source']==src];bygen=collections.defaultdict(list)
    for r in rr:bygen[r['mag']].append(r)
    for mag,gr in bygen.items():
        if src=='EQ_MAG':p=W/'05_refine/_drep/out/dereplicated_genomes'/(mag+'.fa')
        elif src=='Environmental_GEM':p=H/'genomes/gem_nucleotide'/(mag+'.fna.gz')
        else:p=smagfiles.get(mag)
        if p is None or not p.exists():raise FileNotFoundError((src,mag,p))
        inputfiles.append(p);op=gzip.open if str(p).endswith('.gz')else open
        with op(p,'rt')as f:seqs={x.id:str(x.seq)for x in SeqIO.parse(f,'fasta')}
        for r in gr:
            ctg=r['contig'];seq=seqs.get(ctg)
            if seq is None:seq=seqs.get(ctg.split('_',1)[-1])
            if seq is None:raise KeyError((p,ctg,list(seqs)[:2]))
            r['contig_length']=len(seq);r['nearest_edge_bp']=min(r['start']-1,len(seq)-r['end']);assert r['nearest_edge_bp']>=0
            dna=Seq(seq[r['start']-1:r['end']]);dna=dna.reverse_complement()if r['strand']==-1 else dna
            assert len(dna)%3==0
            translated=str(dna.translate(table=11)).rstrip('*');expected=proteins[r['id']]
            r['translation_matches']=int(len(translated)==len(expected)and translated[1:]==expected[1:]and expected[0]=='M')
            assert r['translation_matches'],('Translation mismatch',r['id'],translated[:20],expected[:20])
            r['context_status']='strict_ABC'if r['strict_ABC']else('edge_censored'if r['nearest_edge_bp']<=5000 or r['partial_marker_near']else'internal_without_complete_ABC')
write(R/'loci.tsv',rows)
summary={s:{'complete_EctC':len(rr),'WEYH':sum(r['WEYH']for r in rr),'strict_ABC':sum(r['strict_ABC']for r in rr),'strict_carrier_MAGs':len({r['mag']for r in rr if r['strict_ABC']}),'cov70_strict_ABC':sum(r['strict_ABC_cov70']for r in rr),'WEYH_context':dict(collections.Counter(r['context_status']for r in rr if r['WEYH']))}for s in coords for rr in [[r for r in rows if r['source']==s]]}
write(R/'inputs.tsv',[{'path':str(p),'sha256':sha(p)}for p in sorted(set(inputfiles))])
(R/'summary.json').write_text(json.dumps({'job':os.environ['SLURM_JOB_ID'],'sources':summary},indent=2))
