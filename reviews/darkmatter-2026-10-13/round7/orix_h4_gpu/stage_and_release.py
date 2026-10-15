#!/usr/bin/env python3
"""Wait for CPU preparation; transfer completed payload; release held ORIX jobs."""
import datetime,json,pathlib,shlex,subprocess,sys,time
HERE=pathlib.Path(__file__).resolve().parent
IBEX='/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h4_structure/orix_gpu_20261014'
ORIX='/mnt/data/u/hohndor/eq-current-structures/h4-round7-gpu-20261014'
CACHE=pathlib.Path('/home/leechuck/.cache/darkmatter-h4-orix-20261014/payload')
PREP='53760480'; JOB=sys.argv[1]
assert JOB.isdigit()
def ssh(host,cmd):
 return subprocess.check_output(['ssh','-o','BatchMode=yes','-o','LogLevel=ERROR','-o','ConnectTimeout=15',host,cmd],text=True,timeout=60).strip()
def status(stage,**kw):
 d={'stage':stage,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'preparation_job':PREP,'orix_job':JOB,**kw}
 p=HERE/'transfer_status.json'; tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(d,indent=2)+'\n');tmp.replace(p)
 print(json.dumps(d),flush=True)
def transfer(src,dst):
 for attempt in range(1,4):
  p=subprocess.run(['rsync','-aL','--partial','--timeout=120','-e','ssh -o BatchMode=yes -o LogLevel=ERROR -o ConnectTimeout=15',src,dst])
  if p.returncode==0:return
  status('TRANSFER_RETRY',attempt=attempt);time.sleep(15)
 raise RuntimeError('Transfer failed three times; ORIX job remains held')
try:
 status('WAITING_FOR_CPU_PREPARATION')
 deadline=time.monotonic()+7200
 while time.monotonic()<deadline:
  state=ssh('ibex',f'sacct -nX -j {PREP} -o State -P').strip()
  if state=='COMPLETED':break
  if state.startswith(('FAILED','CANCELLED','TIMEOUT','OUT_OF_MEMORY','NODE_FAIL')):raise RuntimeError('Preparation '+state)
  time.sleep(45)
 else:raise TimeoutError('Preparation not complete within two hours; ORIX job remains held')
 expected=ssh('ibex',f'cat {IBEX}/PREPARED')
 assert len(expected)==64 and all(c in '0123456789abcdef' for c in expected)
 CACHE.mkdir(parents=True,exist_ok=True)
 status('TRANSFERRING_IBEX_TO_LOCAL')
 transfer(f'ibex:{IBEX}/payload/',str(CACHE)+'/')
 import hashlib
 assert hashlib.sha256((CACHE/'payload.sha256').read_bytes()).hexdigest()==expected
 status('TRANSFERRING_LOCAL_TO_ORIX')
 transfer(str(CACHE)+'/',f'orix:{ORIX}/payload/')
 # Full file hashes are checked by the guarded Slurm job before any search.
 got=ssh('orix',f'sha256sum {ORIX}/payload/payload.sha256').split()[0]
 assert got==expected
 ssh('orix',f'printf "%s\\n" {shlex.quote(expected)} > {ORIX}/TRANSFER_COMPLETE.tmp && mv {ORIX}/TRANSFER_COMPLETE.tmp {ORIX}/TRANSFER_COMPLETE && scontrol release {JOB}')
 status('RELEASED_TO_ORIX_FREECYCLE',manifest_sha256=expected)
except Exception as e:
 status('ERROR_JOB_REMAINS_HELD',error=str(e));raise
