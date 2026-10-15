#!/bin/bash
# darkclean_6 — exact dark PROTEIN count + length/completeness distribution for the
# complete 69,938,628-rep dark set. Answers: (a) exact proteins with no UniRef50 hit,
# (b) how many dark reps/proteins are short (<100 aa) or incomplete (prodigal partial != 00).
#SBATCH --job-name=dark_len
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/logs/dark_len_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/logs/dark_len_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --partition=batch
set -euo pipefail

PC=/ibex/scratch/projects/c2014/rob/metagenomics/04_protein_catalog
DARK=$PC/darkmatter_clean/04_catalog_dark/dark_rep_ids.txt.gz   # 69,938,628 dark rep IDs
CLU=$PC/all_linclust90/all_cluster90.tsv                        # rep<TAB>member (753,305,533 lines)
FAA=$PC/all_linclust90/all_reps90.faa.gz                        # rep sequences (449,265,306), prodigal headers
OUT=$PC/darkmatter_clean/06_dark_lengths
mkdir -p "$OUT"

echo "[$(date -Iseconds)] step0: materialise dark id list"
zcat "$DARK" > "$OUT/dark_ids.txt"
wc -l "$OUT/dark_ids.txt"

echo "[$(date -Iseconds)] step1: dark protein count = cluster sizes for dark reps (col1)"
awk 'NR==FNR{d[$1]=1; next} ($1 in d){c[$1]++} END{for(r in c) print r"\t"c[r]}' \
    "$OUT/dark_ids.txt" "$CLU" | gzip > "$OUT/dark_rep_size.tsv.gz"

echo "[$(date -Iseconds)] step2: per dark rep, aa length (from coords) + partial flag (from header)"
zcat "$FAA" | awk -v IDS="$OUT/dark_ids.txt" '
    BEGIN{ while((getline l < IDS)>0) d[l]=1 }
    /^>/{
        n=split($0, a, " # ");
        id=substr(a[1],2);
        if(!(id in d)) next;
        if(n>=5){ aa=int((a[3]-a[2]+1)/3);
                  pf="NA"; if(match(a[5],/partial=[0-9]+/)) pf=substr(a[5],RSTART+8,RLENGTH-8); }
        else    { aa=-1; pf="NA" }
        print id"\t"aa"\t"pf;
    }' | gzip > "$OUT/dark_rep_len.tsv.gz"

echo "[$(date -Iseconds)] step3: join size+length, summarise"
zcat "$OUT/dark_rep_len.tsv.gz" > "$OUT/_len.tsv"
zcat "$OUT/dark_rep_size.tsv.gz" | awk -v LENF="$OUT/_len.tsv" '
    BEGIN{ while((getline l < LENF)>0){ split(l,a,"\t"); AA[a[1]]=a[2]; PF[a[1]]=a[3] } }
    { rep=$1; sz=$2; nrep++; tot+=sz; aa=AA[rep]; pf=PF[rep];
      if(aa=="" || aa<0){ nolen++; next }
      if(aa<60) b="lt60"; else if(aa<100) b="60_99"; else b="ge100";
      R[b]++; P[b]+=sz;
      if(pf=="00"){ compR++; compP+=sz } else { partR++; partP+=sz }
    }
    END{
      printf "dark_reps\t%d\n", nrep;
      printf "dark_proteins\t%d\n", tot;
      printf "reps_lt60\t%d\tproteins_lt60\t%d\n",   R["lt60"],  P["lt60"];
      printf "reps_60_99\t%d\tproteins_60_99\t%d\n", R["60_99"], P["60_99"];
      printf "reps_ge100\t%d\tproteins_ge100\t%d\n", R["ge100"], P["ge100"];
      printf "reps_complete(partial00)\t%d\tproteins_complete\t%d\n", compR, compP;
      printf "reps_incomplete(partial!=00)\t%d\tproteins_incomplete\t%d\n", partR, partP;
      printf "reps_missing_length\t%d\n", nolen;
    }' | tee "$OUT/SUMMARY.txt"

rm -f "$OUT/_len.tsv"
echo "[$(date -Iseconds)] done -> $OUT/SUMMARY.txt"
