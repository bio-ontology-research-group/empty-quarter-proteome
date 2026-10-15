#!/bin/bash
# §2.5 gene-environment rerun on the 106,042 universe: eco_0 sequencing depth -> eco_1 indicators
# -> eco_2 space-partial -> eco_3 community PERMANOVA/varpart/PCoA -> eco_4 plant-root
# compartment families.
# Methods = June phase_e/phase_b verbatim (see README).
#SBATCH --job-name=ale_eco25
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_ecology_clean/logs/eco_%j.out
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --partition=batch
set -euo pipefail
S=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_ecology_clean/scripts
cd "$S"
echo "[eco] library selftest"; python3 eco_stats.py
echo "[eco] 0/5 sequencing depth"; python3 eco_0_depth.py
echo "[eco] 1/5 indicators";  python3 eco_1_indicators.py
echo "[eco] 2/5 space-partial"; python3 eco_2_partial.py
echo "[eco] 3/5 community";   python3 eco_3_community.py
echo "[eco] 4/5 compartment families"; python3 eco_4_compartment_families.py
echo "[eco] all done"
