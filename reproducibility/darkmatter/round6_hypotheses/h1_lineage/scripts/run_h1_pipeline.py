#!/usr/bin/env python3
"""run_h1_pipeline.py - Master execution runner for Hypothesis 1 feasibility pilot."""

import os
import sys
import subprocess
import time

def run_step(step_name, cmd):
    print(f"\n{'='*70}\nSTARTING: {step_name}\nCMD: {' '.join(cmd)}\n{'='*70}")
    t0 = time.time()
    res = subprocess.run(cmd)
    dt = time.time() - t0
    if res.returncode != 0:
        print(f"FAILED: {step_name} exited with code {res.returncode} in {dt:.2f}s")
        sys.exit(res.returncode)
    print(f"COMPLETED: {step_name} in {dt:.2f}s")

def main():
    here = os.path.dirname(os.path.abspath(__file__))
    py_exec = sys.executable

    steps = [
        ("Step 1: Denominator & Gene-Caller Audit", [py_exec, os.path.join(here, "01_denominator_caller_audit.py")]),
        ("Step 2: Decoupling & Statistical Models", [py_exec, os.path.join(here, "02_decoupling_and_confounding_models.py")]),
        ("Step 3: Phylogenetic Ancestry Audit", [py_exec, os.path.join(here, "03_phylogenetic_ancestry_audit.py")]),
        ("Step 4: Accessory Orthogroup Prototype", [py_exec, os.path.join(here, "04_accessory_orthogroup_prototype.py")]),
    ]

    t_start = time.time()
    for name, cmd in steps:
        run_step(name, cmd)

    total_time = time.time() - t_start
    print(f"\n{'='*70}\nALL H1 PIPELINE STEPS COMPLETED SUCCESSFULLY in {total_time:.2f}s\n{'='*70}")

if __name__ == "__main__":
    main()
