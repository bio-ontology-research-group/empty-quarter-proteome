from common import *
import collections,json,subprocess,random,statistics,numpy as np
from Bio import SeqIO,AlignIO,Phylo
from Bio.Align import PairwiseAligner,substitution_matrices
from scipy.spatial.distance import cdist
R=N/'results/f3_comparison';R.mkdir(parents=True,exist_ok=True);C=D/'root_audit/f345_comparison/Helio'
meta=tab(D/'root_audit/f3_site_conserving/Helio/metadata.tsv');candidates={r['protein_id']:r for r in tab(D/'f345/results/full_v1/f3_rhodopsin/f3_candidates_summary.tsv')}
smag={r['user_genome']:r for r in tab(N/'inputs/sheet1.tsv')};gem={r['genome_id']:r for r in tab(N/'inputs/gem_metadata.tsv')}
alias={'Actinobacteriota':'Actinomycetota','Actinobacteria':'Actinomycetota','Proteobacteria':'Pseudomonadota','Candidatus Patescibacteria':'Patescibacteriota','Patescibacteria':'Patescibacteriota','Chloroflexi':'Chloroflexota','Deinococcus-Thermus':'Deinococcota'}
def canonical(s):return alias.get(s,s)
for r in meta:
    src=r['source'];orig=r['original_id']
    if src=='EQ':raw=candidates[orig]['gtdbtk_taxonomy'];r['genome']=candidates[orig]['mag_or_genome'];ph=next(x[3:]for x in raw.split(';')if x.startswith('p__'))
    elif src in {'SMAG','GEM'}:
        g=candidates[orig]['mag_or_genome'];r['genome']=g
        if src=='SMAG':d=smag[g];raw=';'.join(d[k]for k in ['classification','phylum','class','order','family','genus','species']);ph=d['phylum'].removeprefix('p__')
        else:
            raw=gem[g]['ecosystem'];assert raw.startswith('d__')and';p__'in raw;ph=next(x[3:]for x in raw.split(';')if x.startswith('p__'))
    else:
        r['genome']='reference_'+orig
        if orig=='AVZ43932.1':raw='Actinomycetota;48C12 experimental metagenomic clone';ph='Actinomycetota'
        else:
            j=json.load(open(D/'f345/references/raw_uniprot_json'/(orig+'.json')));lin=j['organism'].get('lineage',[]);raw=';'.join(lin);known={'Actinomycetota','Patescibacteriota','Chloroflexota','Halobacteriota','Deinococcota','Pseudomonadota','Bacteroidota','Bacillota','Cyanobacteriota','Acidobacteriota','Planctomycetota','Verrucomicrobiota'};ph=next((canonical(x)for x in lin if canonical(x)in known),'Unresolved')
    r['published_taxonomy']=raw;r['phylum']=canonical(ph)
write(R/'taxonomy.tsv',meta)
seqs=SeqIO.to_dict(SeqIO.parse(C/'combined.faa','fasta'));al=AlignIO.read(D/'root_audit/f3_reference_sensitivity/helio.sto','stockholm');als={x.id:str(x.seq)for x in al};rf=al.column_annotations['reference_annotation'];match=[i for i,c in enumerate(rf)if c not in '.- '];core={};arch=[]
for r in meta:
    a=als[r['id']];ix=[i for i in match if a[i]not in '.-'];lo,hi=min(ix),max(ix);s=''.join(c for c in a[lo:hi+1]if c not in '.-').upper();core[r['id']]=s;full=str(seqs[r['id']].seq)
    assert s in full;arch.append({'id':r['id'],'source':r['source'],'phylum':r['phylum'],'full_length':len(full),'core_length':len(s),'N_extension':full.index(s),'C_extension':len(full)-full.index(s)-len(s),'match_occupancy':len(ix)/len(match)})
write(R/'architecture.tsv',arch)
with open(R/'core.faa','w')as f:
    for pid,s in core.items():f.write('>'+pid+'\n'+s+'\n')
with open(R/'full_same_set.faa','w')as f:
    for pid in core:f.write('>'+pid+'\n'+str(seqs[pid].seq)+'\n')
# Joint fixed core partition, with every source treated identically.
mm='/ibex/scratch/projects/c2014/rob/miniconda3/envs/mmseqs2/bin/mmseqs'
maps={}
for kind in ['frozen_full','full','core']:
    for lev in ['50','80']:
        if kind=='frozen_full':p=C/('cluster'+lev+'_cluster.tsv')
        else:
            prefix=R/(kind+lev);p=Path(str(prefix)+'_cluster.tsv')
            with open(R/('mmseqs_'+kind+lev+'.log'),'w')as log:subprocess.run([mm,'easy-cluster',str(R/('core.faa'if kind=='core'else'full_same_set.faa')),str(prefix),str(R/('tmp_'+kind+lev)),'--min-seq-id',str(int(lev)/100),'-c','0.8','--cov-mode','0','--threads','4'],check=True,stdout=log,stderr=log)
        maps[(kind,lev)]={mem:rep for rep,mem in (l.split()for l in open(p))}
# Sequence-SHA dedup within source, with explicit conflicting-taxonomy assertion.
unique=[];seen={}
for r in sorted(meta,key=lambda r:r['id']):
    k=(r['source'],r['sequence_sha256'])
    if k in seen:
        assert seen[k]['phylum']==r['phylum'],('identical sequence in different phyla',k)
        continue
    seen[k]=r;unique.append(r)
sources=['EQ','SMAG','GEM'];counts=[{'source':s,'phylum':ph,'n_unique':sum(r['source']==s and r['phylum']==ph for r in unique)}for ph in sorted({r['phylum']for r in unique})for s in sources+['curated_reference']];write(R/'source_phylum_counts.tsv',counts)
raref=[];unsupported=[]
for kind,lev in maps:
    cm=maps[(kind,lev)];refclus={cm[r['id']]for r in unique if r['source']=='curated_reference'}
    for ph in ['ALL']+sorted({r['phylum']for r in unique if r['source']=='EQ'}):
        groups={s:[r['id']for r in unique if r['source']==s and(ph=='ALL'or r['phylum']==ph)]for s in sources}
        for comp in [('EQ','SMAG','GEM'),('EQ','SMAG'),('EQ','GEM')]:
            n=min(len(groups[s])for s in comp)
            if n<5:
                unsupported.append({'kind':kind,'identity':lev,'phylum':ph,'comparison':'+'.join(comp),'n':n});continue
            n=min(n,100)
            for s in comp:
                draws=np.array([len({cm[x]for x in random.Random(20261015+i).sample(groups[s],n)}-refclus)for i in range(1000)])
                raref.append({'kind':kind,'identity':lev,'phylum':ph,'comparison':'+'.join(comp),'source':s,'equal_n':n,'mean_external_clusters':float(draws.mean()),'sampling_q025':float(np.quantile(draws,.025)),'sampling_q975':float(np.quantile(draws,.975))})
write(R/'rarefaction.tsv',raref);write(R/'unsupported_strata.tsv',unsupported)
# Native-space proximity: environmental comparator union, leave-genome-out same-source baseline.
z=np.load(D/'root_audit/f345_embeddings/embeddings.npz',allow_pickle=True);ids=list(map(str,z['ids']));emb=z['embeddings'];z.close();zi={x:i for i,x in enumerate(ids)};dist=cdist(emb[[zi[r['id']]for r in unique]],emb[[zi[r['id']]for r in unique]],'cosine')
pa=PairwiseAligner();pa.mode='global';pa.substitution_matrix=substitution_matrices.load('BLOSUM62');pa.open_gap_score=-10;pa.extend_gap_score=-.5
near=[]
for i,r in enumerate(unique):
    if r['source']=='curated_reference':continue
    for mode in ['external_environment','within_source_other_genome','curated']:
        jj=[j for j,t in enumerate(unique)if t['phylum']==r['phylum']and((mode=='external_environment'and t['source']in sources and t['source']!=r['source'])or(mode=='within_source_other_genome'and t['source']==r['source']and t['genome']!=r['genome'])or(mode=='curated'and t['source']=='curated_reference'))]
        if not jj:continue
        j=min(jj,key=lambda j:dist[i,j]);t=unique[j];a=pa.align(core[r['id']],core[t['id']])[0];counts_a=a.counts();pid=100*counts_a.identities/a.shape[1]
        near.append({'id':r['id'],'source':r['source'],'phylum':r['phylum'],'mode':mode,'n_comparators':len(jj),'nearest_embedding_id':t['id'],'cosine':float(dist[i,j]),'core_global_identity_to_embedding_neighbor':float(pid)})
write(R/'native_neighbors.tsv',near)
# HMM match-column phylogeny. Safe IDs and map retained; support is FastTree SH-like local support.
with open(R/'match_alignment.faa','w')as f:
    for i,r in enumerate(unique):f.write('>T'+str(i)+'\n'+''.join(als[r['id']][j].upper().replace('.','-')for j in match)+'\n')
write(R/'tree_tips.tsv',[{'tip':'T'+str(i),**r}for i,r in enumerate(unique)])
with open(R/'helio.tree','w')as f,open(R/'FastTree.log','w')as log:subprocess.run(['/ibex/scratch/projects/c2014/rob/miniconda3/envs/selection/bin/FastTree','-wag','-gamma',str(R/'match_alignment.faa')],check=True,stdout=f,stderr=log)
(R/'summary.json').write_text(json.dumps({'job':os.environ['SLURM_JOB_ID'],'records':len(meta),'unique_source_sequences':len(unique),'match_columns':len(match),'metadata_sources':{'SMAG':sha(N/'inputs/smag_extended_table2.xlsx'),'GEM':sha(N/'inputs/gem_metadata.tsv')},'phylum_aliases':alias,'reference_taxonomy_unresolved':sum(r['phylum']=='Unresolved'for r in meta if r['source']=='curated_reference'),'limits':['Published taxonomy versions differ; phylum name harmonization is not full GTDB reclassification','Rarefaction quantiles are conditional subsampling ranges, not independent biological confidence intervals','One sampled panel per source; no desert endemism inference','Core-global identity is to the embedding-nearest neighbor, not a maximum identity search','ESMFold and ProstT5 predictions share sequence-derived evidence; neither validates activity']},indent=2))
