from common import *
import json,re,subprocess,collections,numpy as np
R=N/'results/read_support';R.mkdir(parents=True,exist_ok=True)
sam='/ibex/scratch/projects/c2014/rob/miniconda3/envs/instrain/bin/samtools'
f3={r['protein_id']:r for r in tab(D/'f345/results/full_v1/f3_rhodopsin/f3_candidates_summary.tsv')};targets=[]
for r in tab(D/'root_audit/f3_site_conserving/Helio/metadata.tsv'):
    if r['source']=='EQ':targets.append({**f3[r['original_id']],'track':'F3_Helio','junction':None})
for r in tab(D/'f345/results/full_v1/f4_ectc/f4_candidates_summary.tsv'):
    if r['cohort']=='EQ_MAG'and r['is_same_orf_ectcd_fusion']=='1'and r['prodigal_partial']=='00':
        domains=sorted([tuple(map(int,re.split(r'\.\.|-',r[k])))for k in ['same_orf_ectc_domain_coords','same_orf_ectd_domain_coords']]);mid=(domains[0][1]+domains[1][0])//2;pos=int(r['start'])+3*(mid-1)if r['strand']in {'+','1'}else int(r['stop'])-3*(mid-1);targets.append({**r,'track':'F4_EctCD','junction':pos})
rows=[];headers={};sources={}
for r in targets:
    sample=r['mag_or_genome'].split('__')[0];bam=Path('/ibex/scratch/projects/c2014/EmptyQuarter_Data/soil/wgs_assemblies/binning')/sample/(sample+'_metaG.reads.sorted.bam');assert bam.exists(),bam
    idx=Path(str(bam)+'.bai')
    if not idx.exists():
        idx=bam.with_suffix('.bai')
    if not idx.exists():
        idx=R/(sample+'.bai');subprocess.run([sam,'index','-@','4',str(bam),str(idx)],check=True)
    if sample not in headers:
        h=subprocess.check_output([sam,'view','-H',str(bam)],text=True);(R/(sample+'.header.sam')).write_text(h);ctgs={}
        for l in h.splitlines():
            if l.startswith('@SQ'):
                fields=dict(x.split(':',1)for x in l.split('\t')[1:]);ctgs[fields['SN']]=int(fields['LN'])
        headers[sample]=ctgs;sources[sample]={'sample':sample,'bam':str(bam),'bytes':bam.stat().st_size,'mtime_ns':bam.stat().st_mtime_ns,'index':str(idx),'index_sha256':sha(idx),'header_sha256':sha(R/(sample+'.header.sam'))}
    ctg=r['contig'];ctgs=headers[sample]
    if ctg not in ctgs:ctg=ctg.split('_',1)[-1]
    assert ctg in ctgs,(sample,ctg)
    a,b=int(r['start']),int(r['stop']);region=f'{ctg}:{a}-{b}';pid=r['protein_id'];stem=hashlib.sha256(pid.encode()).hexdigest()[:16]
    # BAM##idx##path permits an index stored separately from an immutable BAM.
    bamarg=str(bam)+'##idx##'+str(idx)
    out=subprocess.check_output([sam,'depth','-aa','-s','-q','20','-Q','20','-G','3844','-r',region,bamarg],text=True);(R/(stem+'.depth.tsv')).write_text(out);depth={int(x[1]):int(x[2])for x in(l.split()for l in out.splitlines())};assert all(a<=p<=b for p in depth);v=np.array([depth.get(p,0)for p in range(a,b+1)])
    subset=R/(stem+'.bam');subprocess.run([sam,'view','-b','-o',str(subset),bamarg,region],check=True);subprocess.run([sam,'index',str(subset)],check=True)
    replay=subprocess.check_output([sam,'depth','-aa','-s','-q','20','-Q','20','-G','3844','-r',region,str(subset)],text=True);assert out==replay,('regional BAM replay differs',pid)
    rec={'id':pid,'track':r['track'],'sample':sample,'region':region,'length_nt':b-a+1,'median_depth':float(np.median(v)),'mean_depth':float(v.mean()),'breadth_ge1':float((v>=1).mean()),'breadth_ge5':float((v>=5).mean()),'coverage_support_gate':int(np.median(v)>=5 and(v>=1).mean()>=.95),'junction_position':r['junction']or'','junction_spanning_templates':''}
    if r['junction']:
        p=r['junction'];window=f'{ctg}:{max(1,p-100)}-{min(ctgs[ctg],p+100)}';samtext=subprocess.check_output([sam,'view','-q','20','-F','3844',bamarg,window],text=True);names=set()
        for line in samtext.splitlines():
            x=line.split('\t');loc=int(x[3]);blocks=[]
            for n,op in re.findall(r'(\d+)([MIDNSHP=X])',x[5]):
                n=int(n)
                if op in 'M=X':blocks.append((loc,loc+n-1));loc+=n
                elif op in 'DN':loc+=n
            if any(lo<=p-20 and hi>=p+20 for lo,hi in blocks):names.add(x[0])
        rec['junction_spanning_templates']=len(names)
        # Retain the exact qualifying alignments for audit; read names remain on project storage.
        (R/(stem+'.junction.sam')).write_text(samtext)
    rec['regional_bam']=str(subset);rec['regional_bam_sha256']=sha(subset);rec['depth_replay_exact']=1
    rows.append(rec);print(json.dumps(rec),flush=True)
write(R/'loci.tsv',rows);write(R/'alignment_sources.tsv',list(sources.values()))
(R/'summary.json').write_text(json.dumps({'job':os.environ['SLURM_JOB_ID'],'samtools_version':subprocess.check_output([sam,'--version'],text=True).splitlines()[0],'mapping_quality':20,'base_quality':20,'excluded_flags':3844,'overlapping_mates':'samtools depth -s','gate':'median depth >=5 and breadth at depth>=1 >=0.95','junction':'at least20reference bases on each side of domain-gap midpoint in one contiguous M/= /X CIGAR block; distinct query names','tracks':{t:{'loci':len(rr),'coverage_gate_pass':sum(r['coverage_support_gate']for r in rr),'median_locus_depth':float(np.median([r['median_depth']for r in rr]))}for t in sorted({r['track']for r in rows})for rr in [[r for r in rows if r['track']==t]]},'limits':['Support in original assembly mappings, not an independent assembly or functional assay','BAMs identified by exact path, size, mtime, header and index hashes; full BAM hashes are not recomputed here','One source sample per locus; mapping specificity conditional on original assembly and MAPQ']},indent=2))
