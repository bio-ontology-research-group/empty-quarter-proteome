from common import *
from scipy import sparse
import sys
mode=sys.argv[1];mi=read(M/'mag_index.tsv').sort_values('mag_index');genomes=mi.mag_id.tolist();n=len(genomes)
rules=json.loads((R/'inputs/panels.json').read_text())['rules'];out=pd.DataFrame({'genome':genomes});audit=[]
if mode=='ko':
 K=sparse.load_npz(M/'mag_by_ko_copy.npz').tocsc();ki=read(M/'ko_index.tsv').sort_values('ko_index');lookup={k:i for i,k in enumerate(ki.ko)};values={}
 for row in read(R/'inputs/microtrait_hmm.txt').fillna('').to_dict('records'):
  leaf=row['microtrait_hmm-name'];kos=re.findall(r'K\d{5}',row['microtrait_hmm-dbxref_kegg'])
  # KO projections are a sensitivity analysis, not equivalent to profile-specific microTrait calls.
  if kos:values[leaf]=np.maximum.reduce([np.asarray(K[:,lookup[k]].toarray()).ravel()>0 if k in lookup else np.zeros(n,bool) for k in kos]).astype(float)
 for name,expr in rules.items():
  v=evaluate(expr,values,n);out['KO_projection:'+name]=v;audit.append(dict(panel=name,source='published_rules_KO_projection',unknown=int(np.isnan(v).sum()),positive=int(np.nansum(v))))
 # Existing curated repair flags preserve legacy definitions for comparison only.
 uv=read(E/'inputs/ibex/06_metabolism/uv_survival_clean/results/rep_traits.tsv').set_index('rep').reindex(genomes)
 assert not uv.annotated.isna().any()
 for c in ['photolyase','nhej','carotenoid','SOS_translesion','base_excision_repair','dna_protection_oxidative','mismatch_repair','nucleotide_excision_repair']:out['legacy:'+c]=uv[c].to_numpy()
 for name,expr in {'trehalose_strict':"'K16055' | ('K00697' & 'K01087')",'ectoine_marker':"'K06720'",'glucosylglycerol_strict':"'K03692' & 'K05978'",'mannosylglycerate_strict':"'K18818' | ('K05947' & 'K07026')"}.items():
  kv={k:(np.asarray(K[:,lookup[k]].toarray()).ravel()>0).astype(float) if k in lookup else np.zeros(n) for k in re.findall(r'K\d{5}',expr)};out['curated:'+name]=evaluate(expr,kv,n)
else:
 mf=json.loads((R/'inputs/hmm-manifest.json').read_text());byname={x['name']:x['leaf'] for x in mf};assert len(byname)==len(mf)
 values={x['leaf']:np.zeros(n) for x in mf};gi={g:i for i,g in enumerate(genomes)};hits=[]
 for i in range(50):
  assert (R/f'results/hmm/{i:02d}.DONE').exists(),i
  with (R/f'results/hmm/{i:02d}.tbl').open() as f:
   for line in f:
    if line.startswith('#'):continue
    x=line.split();gene=x[0];leaf=byname[x[2]];g=gene.split('::',1)[0];values[leaf][gi[g]]=1
    hits.append(dict(genome=g,gene=gene,profile=leaf,evalue=float(x[4]),score=float(x[5])))
 save(pd.DataFrame(hits),'profile_gene_hits.tsv.gz')
 for name,expr in rules.items():
  v=evaluate(expr,values,n);out['HMM:'+name]=v;audit.append(dict(panel=name,source='microtrait_trusted_cutoff',unknown=int(np.isnan(v).sum()),positive=int(np.nansum(v))))
save(out,f'genome_panels_{mode}.tsv');save(pd.DataFrame(audit),f'panel_audit_{mode}.tsv');print(mode,out.shape,flush=True)
