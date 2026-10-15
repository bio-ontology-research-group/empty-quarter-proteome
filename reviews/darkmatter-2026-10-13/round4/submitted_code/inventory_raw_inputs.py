#!/usr/bin/env python3
"""Check the frozen run table's raw paths over SSH; metadata only, no writes on IBEX."""
import csv
import json
from pathlib import Path
import shlex
import subprocess

here = Path(__file__).resolve().parent
repo = here.parents[1]
with (repo / 'code/ena_deposit/wgs_runs.tsv').open() as source:
    rows = list(csv.DictReader(source, delimiter='\t'))
REMOTE = "import json,os,csv,sys\nrows=json.load(sys.stdin)\nw=csv.writer(sys.stdout,delimiter='\\t');w.writerow(['sample','run_alias','read','path','resolved_path','expected_bytes','observed_bytes','status'])\nfor r in rows:\n for mate in ('r1','r2'):\n  p=r[mate+'_path'];expected=r[mate+'_bytes']\n  try:\n   s=os.stat(p);status='size_match' if s.st_size==int(expected) else 'size_mismatch'\n   w.writerow([r['sample_id'],r['run_alias'],mate,p,os.path.realpath(p),expected,s.st_size,status])\n  except OSError as e:w.writerow([r['sample_id'],r['run_alias'],mate,p,'',expected,'',type(e).__name__+':'+str(e)])\n"
with (here / 'raw_input_manifest.tsv').open('w') as out:
    subprocess.run(['ssh', 'ibex', 'python3 -c ' + shlex.quote(REMOTE)], input=json.dumps(rows), text=True, stdout=out, check=True)
