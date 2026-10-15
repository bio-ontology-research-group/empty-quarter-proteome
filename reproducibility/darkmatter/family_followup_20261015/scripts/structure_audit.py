from common import *
import json,re,subprocess,collections,numpy as np
from Bio import SeqIO,AlignIO
from Bio.PDB import PDBParser,PDBIO,Select
from Bio.SeqUtils import seq1
from Bio.Seq import Seq
R=N/'results/structure_audit';R.mkdir(parents=True,exist_ok=True);S=N/'structures';I=N/'inputs'
assert (S/'COMPLETE').exists()
parser=PDBParser(QUIET=True);structure=parser.get_structure('reference',I/'6SU3.pdb')
chain=next(c for c in structure[0]if len([r for r in c if r.id[0]==' '])>200)
class ChainSelect(Select):
    def accept_chain(self,c):return c.id==chain.id
    def accept_residue(self,r):return r.id[0]==' '
io=PDBIO();io.set_structure(structure);io.save(str(R/'6SU3_chain.pdb'),ChainSelect())
tool=N/'tools/USalign/TMalign';assert tool.exists()
mapping=tab(I/'helio_structure_mapping.tsv');seqs={r.id:str(r.seq)for r in SeqIO.parse(I/'helio_structures.faa','fasta')};summary={};rows=[]
for pid,seq in seqs.items():
    path=S/(pid+'.pdb');p=parser.get_structure(pid,path);res=[r for r in p.get_residues()if r.id[0]==' '];assert len(res)==len(seq);assert ''.join(seq1(r.resname)for r in res)==seq
    ca=np.load(S/(pid+'.ca_plddt.npy'));pdbca=np.array([r['CA'].bfactor for r in res]);assert np.max(np.abs(pdbca*100-ca))<=.501
    output=subprocess.check_output([str(tool),str(path),str(R/'6SU3_chain.pdb')],text=True);(R/(pid+'.tmalign.txt')).write_text(output)
    scores=re.findall(r'TM-score=\s*([0-9.]+).*(?:Chain|Structure)_([12])',output);assert len(scores)==2,output
    lengths=re.findall(r'Length of (?:Chain|Structure)_([12]):\s*(\d+)',output);assert dict(lengths)['1']==str(len(seq))
    alignment=re.search(r'Aligned length=\s*(\d+), RMSD=\s*([0-9.]+)',output);assert alignment
    summary[pid]={'structure_id':pid,'length':len(seq),'mean_CA_pLDDT':float(ca.mean()),'fraction_CA_ge70':float((ca>=70).mean()),'query_TM':float(dict((n,v)for v,n in scores)['1']),'reference_TM':float(dict((n,v)for v,n in scores)['2']),'aligned_length':int(alignment[1]),'RMSD':float(alignment[2])}
for r in mapping:rows.append({**r,**summary[r['structure_id']]})
write(R/'structures.tsv',rows)
# Verify EQ full-length translations from recovered nucleotide contigs, independently of source protein hash.
candidates={r['protein_id']:r for r in tab(D/'f345/results/full_v1/f3_rhodopsin/f3_candidates_summary.tsv')};genes=[]
for r in mapping:
    if r['source']!='EQ':continue
    c=candidates[r['original_id']];mag=c['mag_or_genome'];path=W/'05_refine/_drep/out/dereplicated_genomes'/(mag+'.fa');records={x.id:str(x.seq)for x in SeqIO.parse(path,'fasta')};contig=c['contig'];nt=records.get(contig,records.get(contig.split('_',1)[-1]));assert nt is not None,(path,contig,list(records)[:2]);a,b=int(c['start']),int(c['stop']);dna=Seq(nt[a-1:b]);dna=dna.reverse_complement()if c['strand']in {'-','-1'}else dna;tr=str(dna.translate(table=11)).rstrip('*');expected=seqs[r['structure_id']];match=len(tr)==len(expected)and tr[1:]==expected[1:]and expected[0]=='M';assert match,(r['id'],'translation')
    genes.append({'id':r['id'],'mag':mag,'contig':contig,'start':a,'end':b,'contig_length':len(nt),'nearest_edge_bp':min(a-1,len(nt)-b),'translation_match':int(match),'assembly_path':str(path),'assembly_sha256':sha(path)})
write(R/'EQ_assembly_verification.tsv',genes)
ss={s:{'records':len(rr),'median_CA_pLDDT':float(np.median([r['mean_CA_pLDDT']for r in rr])),'median_query_TM':float(np.median([r['query_TM']for r in rr])),'query_TM_ge05':sum(r['query_TM']>=.5 for r in rr),'mean_CA_ge70':sum(r['mean_CA_pLDDT']>=70 for r in rr)}for s in sorted({r['source']for r in rows})for rr in [[r for r in rows if r['source']==s]]}
(R/'summary.json').write_text(json.dumps({'job':os.environ['SLURM_JOB_ID'],'reference':'6SU3','chain':chain.id,'sources':ss,'EQ_assembly_translations_verified':len(genes),'TMalign_binary_sha256':sha(tool),'ESMFold_provenance':json.load(open(S/'provenance.json')),'limits':['Predicted structural compatibility, not experimental validation','Full-query TM normalization includes terminal extensions','Assembly translation agreement does not establish raw-read support or exclude assembly errors']},indent=2))
