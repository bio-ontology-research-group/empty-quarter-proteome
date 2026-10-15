#!/bin/bash
#SBATCH --job-name=eq-family-read-inventory
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch,debug
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
set -euo pipefail
cd /ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_followup_20261015_clean
find /ibex/scratch/projects/c2014/EmptyQuarter_Data/soil/wgs_assemblies -maxdepth 5 \( -name '*.bam' -o -name '*.cram' \) > results/source_alignment_paths.txt
mkdir -p literature
curl -fL --retry 2 --max-time 120 https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7049168/fullTextXML -o literature/helio_structure_PMC7049168.xml
sha256sum literature/helio_structure_PMC7049168.xml > literature/checksums.txt
