#!/usr/bin/env python3
"""RQ-C4 — structure-based annotation rescue of sequence-dark proteins.

Documents the IBEX commands and processes the small summaries pulled back to this dir.
No heavy data lives here; the 7.8 GB hit table stays on IBEX.

IBEX (M=/ibex/scratch/projects/c2014/rob/metagenomics, dir 04_structure/foldseek_results):

  # Step 1 — counts
  wc -l all_query_ids.txt afdb_hit_ids.txt struct_nohit_ids.txt

  # Step 2 — SwissProt (=characterised) accession set, then streaming join (sbatch c4_characterised.sh)
  awk -F'\\t' '{a=$2;sub(/^AF-/,"",a);sub(/-F1-model_v[0-9]+$/,"",a);print a}' \\
      ../fsdb_swissprot.lookup | sort -u > sp_acc.txt          # 590,183
  awk -F'\\t' 'FNR==NR{S[$1]=1;next}
    {q=$1;t=$2;sub(/^AF-/,"",t);sub(/-F1-model_v[0-9]+$/,"",t);
     seen[q]=1;rows++; if(t in S){qc[q]=1;sprows++}}
    END{for(k in seen)nq++;for(k in qc)nc++;
        print "rescued_queries_unique\\t"nq; ...}' sp_acc.txt afdb_hits.tsv

Pulled here: structure_rescue_summary.tsv  (key/value summary; the file this script reads).
"""
import csv, os

HERE = os.path.dirname(os.path.abspath(__file__))
SUMMARY = os.path.join(HERE, "structure_rescue_summary.tsv")


def load(path):
    d = {}
    with open(path) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) >= 2 and not row[0].startswith("#"):
                d[row[0]] = row[1]
    return d


def main():
    d = load(SUMMARY)
    nq = int(d["N_query"])
    nhit = int(d["N_afdb_hit"])
    nnohit = int(d["N_nohit"])
    assert nhit + nnohit == nq, "counts must partition the query set"
    rescue = 100.0 * nhit / nq
    print(f"Sequence-dark proteins queried structurally : {nq:,}")
    print(f"  gain an AFDB STRUCTURAL hit (rescued)      : {nhit:,} ({rescue:.2f}%)")
    print(f"  remain structurally novel (no AFDB hit)    : {nnohit:,} ({100-rescue:.2f}%)")
    if "rescued_queries_with_swissprot_target" in d:
        ncfg = int(d["rescued_queries_with_swissprot_target"])
        frac = 100.0 * ncfg / nhit
        print(f"  rescued whose match is a reviewed SwissProt (characterised) target:")
        print(f"    {ncfg:,} of {nhit:,} rescued ({frac:.1f}%)  "
              f"= {100.0*ncfg/nq:.2f}% of all sequence-dark proteins")
        print("  (SwissProt membership = curated function; a lower bound on 'interpretable'.)")


if __name__ == "__main__":
    main()
