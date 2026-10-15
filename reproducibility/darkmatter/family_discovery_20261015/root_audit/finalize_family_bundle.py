"""Close retained family artifacts and diagnose implemented reference-gate sensitivity."""
import csv,hashlib,importlib.util,json,os,platform,subprocess,sys
from pathlib import Path
D=Path('/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean');A=D/'root_audit';F=D/'f345';R=A/'final_bundle';R.mkdir(exist_ok=True)
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb')as f:
  while b:=f.read(4*1024*1024):h.update(b)
 return h.hexdigest()
def tab(p):return list(csv.DictReader(open(p),delimiter='\t'))
assert (A/'f345_plots/COMPLETE').exists();assert(A/'f3_site_conserving/COMPLETE').exists()
# Probe exact implemented F5 selector against known references. This is a sensitivity diagnostic.
spec=importlib.util.spec_from_file_location('f5_gate',F/'scripts/10_f5_carbonic_anhydrase_full.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
with open(R/'f5_reference_hmm.log','w')as log:subprocess.run(['hmmsearch','--cpu','2','--domtblout',str(R/'f5_references.domtbl'),str(F/'references/profiles/f5_carbonic_anhydrase.hmm'),str(F/'references/f5_carbonic_anhydrase_references.faa')],check=True,stdout=log)
by={}
for h in mod.parse_domtblout(R/'f5_references.domtbl'):by.setdefault(h['target_name'],[]).append(h)
accepted,rejected=mod.evaluate_f5_gates(by);positive={'Beta_CA':'Beta-Carbonic_Anhydrase','Gamma_CA':'Gamma-Carbonic_Anhydrase','CsoSCA_Outgroup':'Carboxysome_Shell_CA_CsoSCA'};records=[]
for r in tab(F/'references/verified_references_metadata.tsv'):
 if r['category']not in positive:continue
 candidates=[pid for pid in by if pid.split('|')[0]==r['accession']];pid=candidates[0]if candidates else None;found=accepted.get(pid,{});records.append({'accession':r['accession'],'reference_category':r['category'],'expected_class':positive[r['category']],'selected_class':found.get('assigned_class','NOT_ACCEPTED'),'passes_implemented_gate':int(found.get('assigned_class')==positive[r['category']]),'experimental_function_flag':r['has_function_experimental_evidence'],'best_profile':found.get('best_hit',{}).get('profile_acc','NA'),'best_profile_coverage':found.get('best_hit',{}).get('hmm_coverage','NA')})
with open(R/'f5_reference_gate_diagnostic.tsv','w')as f:w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter='\t');w.writeheader();w.writerows(records)
# Close ORIX-generated embedding checksum bundle after complete replication to IBEX.
verifications={}
for path in [A/'f345_embeddings/final_bundle.sha256']:
 cp=subprocess.run(['sha256sum','-c',str(path)],cwd=path.parent,capture_output=True,text=True);(R/'embedding_checksum_verification.txt').write_text(cp.stdout+cp.stderr);assert cp.returncode==0;verifications[str(path.relative_to(D))]='PASS'
# EctC evidence intersection avoids conflating family membership, site pattern and context.
site={r['id'].split('__EQ__',1)[1]:r for r in tab(A/'f345_anchor_alignment/EctC.site_checks.tsv')};ect=[r for r in tab(A/'f345_final_audit/ectc_context_independent.tsv')if r['source']=='EQ_MAG'and r['complete']=='1'];joint={k:sum(r[k]=='1'and site[r['id']]['exact_reference_site_pattern']=='1'for r in ect)for k in ['ABC_neighborhood','ABC_same_strand','same_orf_EctCD']}
# Exhaustive manifest of the retained analysis branches; skip rebuildable search indices, tmp directories and the already independently hash-verified 9GB EQ concatenation.
roots=[A,D/'ecology/results_root_v3',D/'f345/references',D/'f345/results/full_v1',D/'f345/scripts',D/'f3_context',D/'f12/references_v2',D/'f12/results/pilot_v3_results',D/'f12/scripts']
manifest=[]
for root in roots:
 for p in sorted(root.rglob('*')):
  if not p.is_file()or R in p.parents or any(part.startswith('tmp')or part=='__pycache__'for part in p.relative_to(root).parts)or p.name=='cohort.faa'or p.suffix in {'.dmnd','.pyc'}:continue
  manifest.append({'path':str(p.relative_to(D)),'bytes':p.stat().st_size,'sha256':sha(p)})
with open(R/'artifact_manifest.tsv','w')as f:w=csv.DictWriter(f,fieldnames=['path','bytes','sha256'],delimiter='\t');w.writeheader();w.writerows(manifest)
versions={}
for name in ['numpy','scipy','sklearn','matplotlib','Bio']:
 mod0=__import__(name);versions[name]=getattr(mod0,'__version__','unknown')
(R/'closure.json').write_text(json.dumps({'job_id':os.environ['SLURM_JOB_ID'],'host':platform.node(),'python':sys.version,'libraries':versions,'manifest_files':len(manifest),'manifest_bytes':sum(r['bytes']for r in manifest),'checksum_checks':verifications,'EctC_reference_site_and_context':joint,'F5_reference_gate_failures':[r for r in records if not r['passes_implemented_gate']],'external_source_inventory':['root_audit/f1_full/cohort.tsv:5438originalFAAhashes','root_audit/f345_final_audit/audit.json:fullsourceFASTAandcoordinatehashes','ecology/results_root_v3/profile_inputs.tsv:551CoverMprofilehashes','f345/references/references_manifest.tsv:retainedreferenceprovenance','f3_context/metadata:originalMAGmetadata','root_audit/f345_embeddings/model_cache.sha256:ORIXcheckpointfiles'],'limits':['Manifest closes the analysis state; raw sequencing-read remapping and wet-lab validation were not performed for this family track.','Known reference gate failures limit census sensitivity; no exhaustive enzyme-class counts.','CPUjob53821447 and site-sensitivity53822479 complete before closure; failed/superseded jobs retained separately.']},indent=2))
print('Closure complete',len(manifest),flush=True)
