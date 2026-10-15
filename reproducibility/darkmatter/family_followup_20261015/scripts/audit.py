"""Independent algebra/count checks and checksum closure for accepted follow-up outputs."""
from common import *
import json,math,collections,statistics,re,subprocess
from scipy.stats import wilcoxon
R=N/'results/audit';R.mkdir(parents=True,exist_ok=True)
checks=[]
def passed(name,n):checks.append({'check':name,'n':n,'status':'PASS'})
f4=tab(N/'results/f4_integrity/loci.tsv');s4=json.load(open(N/'results/f4_integrity/summary.json'))
for src,reported in s4['sources'].items():
    rr=[r for r in f4 if r['source']==src];assert len(rr)==reported['complete_EctC'];strict=[r for r in rr if r['WEYH']=='1'and r['complete_A']and r['complete_B']];assert len(strict)==reported['strict_ABC'];assert len({r['mag']for r in strict})==reported['strict_carrier_MAGs'];assert all(r['translation_matches']=='1'for r in rr);assert all(r['context_status']!='internal_without_complete_ABC'or(int(r['nearest_edge_bp'])>5000 and r['partial_marker_near']=='0')for r in rr)
    passed('F4 '+src+' counts and translation closure',len(rr))
ec=tab(N/'results/ecology/paired_visit_decomposition.tsv');site=tab(N/'results/ecology/site_mean_decomposition.tsv')
for r in ec+site:assert abs(float(r['total_pp'])-sum(float(r[k])for k in ['between_pp','within_pp','exclusive_pp']))<1e-9
passed('Exact decomposition per visit and site',len(ec)+len(site))
for r in tab(N/'results/ecology/median_ecology_tests.tsv'):
    values=[float(x['difference_pp'])for x in tab(N/'results/ecology'/(r['trait']+'_median_site_effects.tsv'))];assert len(values)==int(r['sites']);assert abs(statistics.median(values)-float(r['median_difference_pp']))<1e-10;assert abs(wilcoxon(values).pvalue-float(r['p_raw']))<1e-10
passed('Independent median and Wilcoxon reconstruction',6)
meta=tab(N/'results/f3_finer/taxonomy_unique.tsv');expected=[]
for row in tab(N/'results/f3_finer/rarefaction.tsv'):
    cm={m:r for r,m in(l.split()for l in open(N/'results/f3_comparison'/(row['kind']+row['identity']+'_cluster.tsv')))};refs={cm[r['id']]for r in meta if r['source']=='curated_reference'};ids=[r['id']for r in meta if r['source']==row['source']and r[row['rank']]==row['taxon']];counts=collections.Counter(cm[i]for i in ids);n=int(row['equal_n']);Nseq=len(ids);den=math.comb(Nseq,n);ex=sum(1-(math.comb(Nseq-k,n)/den if Nseq-k>=n else 0)for c,k in counts.items()if c not in refs);reported=float(row['mean_external_clusters']);assert abs(ex-reported)<.15,(row,ex)
    expected.append({**row,'analytical_expected_clusters':ex,'MC_minus_expectation':reported-ex})
write(R/'analytical_rarefaction.tsv',expected);passed('Analytical hypergeometric check of rarefaction',len(expected))
models=tab(N/'results/structure_audit/structures.tsv');assert len(models)==178;assert len({r['structure_id']for r in models})==175
for r in models:
    raw=(N/'results/structure_audit'/(r['structure_id']+'.tmalign.txt')).read_text();assert 'TM-score='in raw;assert 0<=float(r['query_TM'])<=1 and 0<=float(r['reference_TM'])<=1;assert 0<=float(r['mean_CA_pLDDT'])<=100
passed('Structure sequence closure and score bounds',len(models))
readrows=tab(N/'results/read_support/loci.tsv')
for r in readrows:
    bam=Path(r['regional_bam']);stem=bam.stem;assert sha(bam)==r['regional_bam_sha256']
    depthfile=bam.with_suffix('.depth.tsv');lines=depthfile.read_text();depth=[int(l.split()[2])for l in lines.splitlines()];assert len(depth)==int(r['length_nt'])
    med=statistics.median(depth);breadth=sum(d>=1 for d in depth)/len(depth)
    assert med==float(r['median_depth']) and abs(breadth-float(r['breadth_ge1']))<1e-12
    assert int(med>=5 and breadth>=.95)==int(r['coverage_support_gate'])
    replay=subprocess.check_output(['/ibex/scratch/projects/c2014/rob/miniconda3/envs/instrain/bin/samtools','depth','-aa','-s','-q','20','-Q','20','-G','3844','-r',r['region'],str(bam)],text=True);assert replay==lines
    if r['junction_position']:
        pos=int(r['junction_position']);names=set()
        for line in bam.with_suffix('.junction.sam').read_text().splitlines():
            fields=line.split('\t');assert int(fields[4])>=20 and not(int(fields[1])&3844)
            cursor=int(fields[3])
            for length,op in re.findall(r'(\d+)([MIDNSHP=X])',fields[5]):
                length=int(length)
                if op in 'M=X'and cursor<=pos-20 and cursor+length-1>=pos+20:names.add(fields[0])
                if op in 'MDN=X':cursor+=length
        assert len(names)==int(r['junction_spanning_templates'])
passed('Regional BAM depth replay, coverage gates and junction recount',len(readrows))
bundle=N/'expanded_embeddings'
assert json.load(open(bundle/'verification_report.json'))['status']=='PASSED'
for line in (bundle/'final_bundle.sha256').read_text().splitlines():
    digest,name=line.split(maxsplit=1);assert sha(bundle/name)==digest,name
passed('Expanded ProstT5 independently verified bundle checksums',1744)
expanded=[]
pubids={l[1:].strip()for l in open(N/'results/published_helio/site_conserving_references.faa')if l.startswith('>')}
for row in tab(N/'results/published_clusters/rarefaction.tsv'):
    cm={m:r for r,m in(l.split()for l in open(N/'results/published_clusters'/(row['kind']+row['identity']+'_cluster.tsv')))}
    refs={cm[i]for i in pubids}|{cm[r['id']]for r in meta if r['source']=='curated_reference'}
    ids=[r['id']for r in meta if r['source']==row['source']and r[row['rank']]==row['taxon']];n=int(row['equal_n']);cnt=collections.Counter(cm[i]for i in ids)
    exact=sum(1-(math.comb(len(ids)-k,n)/math.comb(len(ids),n)if len(ids)-k>=n else 0)for c,k in cnt.items()if c not in refs)
    assert abs(exact-float(row['analytical_expected_clusters']))<1e-12
    assert abs(exact-float(row['mean_external_clusters']))<.15
    expanded.append({**row,'independent_expectation':exact})
write(R/'expanded_analytical_rarefaction.tsv',expanded);passed('Expanded reference rarefaction exact expectation',len(expanded))
for row in tab(N/'results/review_checks/reference_free_rarefaction.tsv'):
    cm={m:r for r,m in(l.split()for l in open(N/'results/review_checks'/(row['kind']+row['identity']+'_cluster.tsv')))}
    ids=[r['id']for r in meta if r['source']==row['source']and r['class']=='Saccharimonadia'];n=int(row['equal_n']);cnt=collections.Counter(cm[i]for i in ids)
    exact=sum(1-(math.comb(len(ids)-k,n)/math.comb(len(ids),n)if len(ids)-k>=n else 0)for k in cnt.values())
    assert abs(exact-float(row['mean_total_clusters']))<.15
passed('Reference-free rarefaction analytical expectations',8)
dual=tab(N/'results/review_checks/dual_AB_domains.tsv')
assert len(dual)==16 and all(int(r['domain_overlap_aa'])==0 and float(r['profile_coverage'])>=.7 for r in dual)
passed('Distinct A/B domain positions in dual-profile ORFs',8)
# Reference-file source metadata and original GPU bundle stay accessible on IBEX.
required=[N/'inputs/smag_extended_table2.xlsx',N/'inputs/gem_metadata.tsv',N/'inputs/6SU3.pdb',N/'structures/provenance.json',N/'results/published_helio/summary.json']
assert all(p.is_file()for p in required)
files=[]
for name in ['README.md','REPRODUCIBILITY.txt','DESIGN.txt','verify_family.py','plot_family.py','IBEX_JOB_RECEIPTS.tsv','ORIX_JOB_RECEIPTS.tsv']:
    p=N/name
    if p.is_file():files.append({'path':name,'bytes':p.stat().st_size,'sha256':sha(p)})
for root in [N/'scripts',N/'inputs',N/'structures',N/'expanded_embeddings',N/'results',N/'literature',N/'manuscript_snapshot',N/'reproducibility_snapshot']:
    for p in sorted(root.rglob('*')):
        if not p.is_file()or R in p.parents or any(part.startswith('tmp')or part=='__pycache__'for part in p.parts):continue
        files.append({'path':str(p.relative_to(N)),'bytes':p.stat().st_size,'sha256':sha(p)})
write(R/'artifact_manifest.tsv',files)
(R/'summary.json').write_text(json.dumps({'job':os.environ['SLURM_JOB_ID'],'checks':checks,'retained_files':len(files),'bytes':sum(r['bytes']for r in files),'scope':'Reconstruction from accepted compact outputs plus retained source identities; upstream family-search raw audits remain in the prior frozen bundle.','limitations':['No wet-lab activity measurements','No claim of novel function, desert endemism or absence from model training data','Raw-read support only if separately documented in accepted read-support outputs']},indent=2))
