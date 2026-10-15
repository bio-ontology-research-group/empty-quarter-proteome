#!/bin/bash
# darkclean_7 — (A) length/completeness profile of ALL 449M representatives (the denominator,
# so the dark fraction can be stated on a >=100 aa / complete basis); (B) the same profile for
# the MAG-mapped dark subset, reusing darkclean_6 outputs (cheap intersection).
#SBATCH --job-name=dark_denom
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/logs/dark_denom_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/logs/dark_denom_%j.err
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --partition=batch
set -euo pipefail

PC=/ibex/scratch/projects/c2014/rob/metagenomics/04_protein_catalog
FAA=$PC/all_linclust90/all_reps90.faa.gz
MAGDARK=$PC/darkmatter_clean/05_mag_dark/mag_dark_rep_ids.txt.gz
LEN6=$PC/darkmatter_clean/06_dark_lengths/dark_rep_len.tsv.gz
SIZE6=$PC/darkmatter_clean/06_dark_lengths/dark_rep_size.tsv.gz
OUT=$PC/darkmatter_clean/06_dark_lengths
mkdir -p "$OUT"

echo "[$(date -Iseconds)] A: length/completeness of ALL representatives (denominator, counters only)"
zcat "$FAA" | awk '
    /^>/{
        n=split($0,a," # ");
        if(n>=5){ aa=int((a[3]-a[2]+1)/3);
                  pf="NA"; if(match(a[5],/partial=[0-9]+/)) pf=substr(a[5],RSTART+8,RLENGTH-8) }
        else    { aa=-1; pf="NA" }
        tot++;
        if(aa<0){ bad++; next }
        if(aa<60)L["lt60"]++; else if(aa<100)L["60_99"]++; else L["ge100"]++;
        if(pf=="00"){ comp++; if(aa>=100) ge100comp++ } else inc++;
    }
    END{
        printf "all_reps\t%d\n", tot;
        printf "all_lt60\t%d\n", L["lt60"];
        printf "all_60_99\t%d\n", L["60_99"];
        printf "all_ge100\t%d\n", L["ge100"];
        printf "all_complete\t%d\n", comp;
        printf "all_incomplete\t%d\n", inc;
        printf "all_ge100_and_complete\t%d\n", ge100comp;
        printf "all_bad\t%d\n", bad;
    }' > "$OUT/all_reps_length_summary.txt"
cat "$OUT/all_reps_length_summary.txt"

echo "[$(date -Iseconds)] B: MAG-mapped dark subset profile (intersection)"
zcat "$MAGDARK" > "$OUT/mag_dark_ids.txt"; wc -l "$OUT/mag_dark_ids.txt"
zcat "$LEN6" | awk -v IDS="$OUT/mag_dark_ids.txt" '
    BEGIN{ while((getline i<IDS)>0) m[i]=1 }
    ($1 in m){ aa=$2; pf=$3; nrep++;
        if(aa<0){ bad++; next }
        if(aa<60)L["lt60"]++; else if(aa<100)L["60_99"]++; else L["ge100"]++;
        if(pf=="00")comp++; else inc++ }
    END{
        printf "mag_dark_reps\t%d\n", nrep;
        printf "mag_lt60\t%d\n", L["lt60"];
        printf "mag_60_99\t%d\n", L["60_99"];
        printf "mag_ge100\t%d\n", L["ge100"];
        printf "mag_complete\t%d\n", comp;
        printf "mag_incomplete\t%d\n", inc;
    }' > "$OUT/mag_dark_length_summary.txt"
zcat "$SIZE6" | awk -v IDS="$OUT/mag_dark_ids.txt" '
    BEGIN{ while((getline i<IDS)>0) m[i]=1 }
    ($1 in m){ tot+=$2; n++ } END{ printf "mag_dark_reps_check\t%d\nmag_dark_proteins\t%d\n", n, tot }' \
    >> "$OUT/mag_dark_length_summary.txt"
cat "$OUT/mag_dark_length_summary.txt"
rm -f "$OUT/mag_dark_ids.txt"
echo "[$(date -Iseconds)] done -> all_reps_length_summary.txt , mag_dark_length_summary.txt"
