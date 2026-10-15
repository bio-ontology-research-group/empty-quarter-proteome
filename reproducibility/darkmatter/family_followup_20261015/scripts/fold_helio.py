"""Pinned ESMFold architecture check; no claim of biochemical validation."""
import os,sys,json,csv,hashlib,time,platform
from pathlib import Path
import numpy as np,torch,transformers
from transformers import AutoTokenizer,EsmForProteinFolding
assert os.environ.get('SLURM_JOB_ID') and torch.cuda.is_available()
inp,out=map(Path,sys.argv[1:3]);out.mkdir(parents=True,exist_ok=True)
modelpath=Path('/mnt/data/u/hohndor/folding/hf/hub/models--facebook--esmfold_v1/snapshots/75a3841ee059df2bf4d56688166c8fb459ddd97a')
tok=AutoTokenizer.from_pretrained(modelpath,local_files_only=True)
model=EsmForProteinFolding.from_pretrained(modelpath,low_cpu_mem_usage=True,local_files_only=True).cuda().eval();model.esm=model.esm.half();model.trunk.set_chunk_size(64)
records=[]
for line in open(inp):
    if line.startswith('>'):records.append([line[1:].split()[0],''])
    else:records[-1][1]+=line.strip()
rows=[]
for pid,seq in records:
    assert set(seq)<=set('ACDEFGHIKLMNPQRSTVWY'),(pid,'noncanonical residue')
    p=out/(pid+'.json')
    if p.exists()and(out/(pid+'.pdb')).exists():rows.append(json.load(open(p)));continue
    token=tok([seq],return_tensors='pt',add_special_tokens=False)['input_ids'].cuda();assert token.shape[1]==len(seq)
    start=time.time()
    with torch.inference_mode():pred=model(token,num_recycles=4)
    pdb=model.output_to_pdb(pred)[0];(out/(pid+'.pdb')).write_text(pdb)
    ca=pred['plddt'][0,:,1].float().cpu().numpy();assert ca.shape==(len(seq),);assert np.isfinite(ca).all()
    # This pinned transformers implementation uses categorical_lddt in [0,1].
    assert ca.max()<=1.00001 and ca.min()>=0,('Unexpected confidence scale',ca.min(),ca.max())
    np.save(out/(pid+'.ca_plddt.npy'),ca*100)
    r={'id':pid,'length':len(seq),'sequence_sha256':hashlib.sha256(seq.encode()).hexdigest(),'mean_CA_pLDDT':float(ca.mean()*100),'fraction_CA_pLDDT_ge70':float((ca>=.7).mean()),'seconds':time.time()-start,'job':os.environ['SLURM_JOB_ID']};p.write_text(json.dumps(r,indent=2));rows.append(r);print(json.dumps(r),flush=True)
    del pred;torch.cuda.empty_cache()
with open(out/'summary.tsv','w')as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
(out/'provenance.json').write_text(json.dumps({'job':os.environ['SLURM_JOB_ID'],'model':str(modelpath),'torch':torch.__version__,'transformers':transformers.__version__,'python':platform.python_version(),'gpu':torch.cuda.get_device_name(),'recycles':4,'chunk_size':64,'input_sha256':hashlib.sha256(inp.read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'plddt':'CA atom index 1; native [0,1] scaled to [0,100]','n':len(rows)},indent=2))
(out/'COMPLETE').write_text(os.environ['SLURM_JOB_ID']+'\n')
