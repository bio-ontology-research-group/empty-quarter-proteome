#!/usr/bin/env python3
"""run_h1_v2_pipeline.py - Master runner for Hypothesis 1 v2 ancestry-aware pipeline."""

import os
import sys
import subprocess
import time

def run_step(step_name, cmd):
    print(f"\n{'='*75}\nSTARTING: {step_name}\nCMD: {' '.join(cmd)}\n{'='*75}")
    t0 = time.time()
    res = subprocess.run(cmd)
    dt = time.time() - t0
    if res.returncode != 0:
        print(f"FAILED: {step_name} exited with code {res.returncode} in {dt:.2f}s")
        sys.exit(res.returncode)
    print(f"COMPLETED: {step_name} in {dt:.2f}s")

def main():
    here = os.path.dirname(os.path.abspath(__file__))
    
    # Environment paths on IBEX
    py_drep = "/ibex/scratch/projects/c2014/rob/miniconda3/envs/drep/bin/python"
    py_meta = "/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python"
    
    if not os.path.exists(py_drep):
        py_drep = sys.executable
    if not os.path.exists(py_meta):
        py_meta = sys.executable

    steps = [
        ("Step 1: Extract Phylogenetic Covariance Matrix (DendroPy)", 
         [py_drep, os.path.join(here, "01_extract_phylogenetic_vcv.py")]),
        ("Step 2: Fit PGLS Ancestry-Aware Models & Lambda Profiling", 
         [py_meta, os.path.join(here, "02_pgls_ancestry_models.py")]),
        ("Step 3: Reconciled Coupling & Descriptive Accessory Audit", 
         [py_meta, os.path.join(here, "03_reconciled_coupling_and_accessory.py")])
    ]

    t_start = time.time()
    for name, cmd in steps:
        run_step(name, cmd)

    total_time = time.time() - t_start
    print(f"\n{'='*75}\nH1 V2 PIPELINE COMPLETED SUCCESSFULLY in {total_time:.2f}s\n{'='*75}")

if __name__ == "__main__":
    main()
