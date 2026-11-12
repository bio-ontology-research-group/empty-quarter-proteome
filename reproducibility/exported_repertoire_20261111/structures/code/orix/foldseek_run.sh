#!/bin/bash
# Foldseek GPU search of monomer models (whole chains and confident-region crops) against PDB and AFDB v6.
# Runs as a side task inside a folding worker job (fold_worker_multi.sbatch, EXTRA_CMD) so that the GPU is
# kept busy by structure prediction while Foldseek loads the 80-GB AFDB target (I/O) and aligns (CPU).
# PDB: --alignment-type 1 -e 10 -s 9.5 --max-seqs 1000 --gpu 1 --prefilter-mode 1; match = qTM >= 0.5.
# AFDB v6 (3Di+AA only on ORIX): --alignment-type 2 -s 9.5 -e 1e-3 --max-seqs 50 --gpu 1; match = E < 1e-3.
set -euo pipefail
D=/mnt/data/u/hohndor/eq-exported-dark-20261111
FS=/mnt/data/u/hohndor/eq-prostt5/bin/foldseek
PDB=/mnt/data/u/hohndor/eq-current-references/pdb-current-20260802/padded/pdb_pad
AFDB=/mnt/data/u/hohndor/eq-current-references/afdb-v6-uniprot-2025-03/gpu_portable_8dc75c74/padded/afdb_pad
W=$D/foldseek
T=${FS_THREADS:-6}
case "$(hostname -s)" in *login*) echo "refusing to run on a login node"; exit 64;; esac
echo "[$(date -Is)] job=$SLURM_JOB_ID node=$(hostname)"; nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader
"$FS" version | tee $W/foldseek_version.txt
FMT1='query,target,alntmscore,qtmscore,ttmscore,prob,lddt,evalue,bits,fident,qcov,tcov,qlen,tlen,alnlen,qstart,qend,tstart,tend,theader'
FMT2='query,target,evalue,bits,fident,qcov,tcov,qlen,tlen,alnlen,qstart,qend,tstart,tend,theader'
for ARM in whole crops; do
  Q=$W/db/monomer_$ARM
  mkdir -p $W/db $W/tmp_$ARM $W/hits
  [ -s $Q.dbtype ] || "$FS" createdb $W/queries/monomer_$ARM $Q --threads $T -v 2
  if [ ! -s $W/hits/pdb_$ARM.tsv ]; then
    rm -f $W/db/pdb_${ARM}_res*
    "$FS" search $Q $PDB $W/db/pdb_${ARM}_res $W/tmp_$ARM/pdb --alignment-type 1 --gpu 1 --prefilter-mode 1 \
        -e 10 -s 9.5 --max-seqs 1000 --threads $T -a -v 2
    "$FS" convertalis $Q $PDB $W/db/pdb_${ARM}_res $W/hits/pdb_$ARM.tsv.tmp --format-output "$FMT1" --threads $T -v 2
    mv $W/hits/pdb_$ARM.tsv.tmp $W/hits/pdb_$ARM.tsv
  fi
  echo "[$(date -Is)] PDB $ARM rows=$(wc -l < $W/hits/pdb_$ARM.tsv)"
  if [ ! -s $W/hits/afdb_$ARM.tsv ]; then
    rm -f $W/db/afdb_${ARM}_res*
    "$FS" search $Q $AFDB $W/db/afdb_${ARM}_res $W/tmp_$ARM/afdb --alignment-type 2 --sort-by-structure-bits 0 \
        --gpu 1 -e 1e-3 -s 9.5 --max-seqs 50 --threads $T --db-load-mode 2 -v 2
    "$FS" convertalis $Q $AFDB $W/db/afdb_${ARM}_res $W/hits/afdb_$ARM.tsv.tmp --format-output "$FMT2" --threads $T -v 2
    mv $W/hits/afdb_$ARM.tsv.tmp $W/hits/afdb_$ARM.tsv
  fi
  echo "[$(date -Is)] AFDB $ARM rows=$(wc -l < $W/hits/afdb_$ARM.tsv)"
  rm -rf $W/tmp_$ARM
done
echo "[$(date -Is)] FOLDSEEK_DONE"
