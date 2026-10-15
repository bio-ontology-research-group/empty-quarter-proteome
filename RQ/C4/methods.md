# RQ-C4 — Methods

**Question.** Can structure-based annotation rescue function that sequence/KO annotation misses?

## Inputs (IBEX)
`M=/ibex/scratch/projects/c2014/rob/metagenomics`, dir `M/04_structure/foldseek_results/`:

- `all_query_ids.txt` — all structurally-queried dark proteins (the ~2.84M UniRef50-no-hit /
  sequence-dark mmseqs reps that were submitted to ProstT5+FoldSeek vs AFDB).
- `afdb_hit_ids.txt` — the subset that obtained an AFDB structural hit.
- `struct_nohit_ids.txt` — the subset with no AFDB structural hit (structurally novel).
- `afdb_hits.tsv` → symlink to `M/04_structure/prostt5_foldseek_orix/merged_afdb_hits.tsv`
  (7.8 GB, 80,231,849 rows). FoldSeek m8, 12 cols:
  `query, target, fident, alnlen, mismatch, gapopen, qstart, qend, tstart, tend, evalue, bits`.
  `target` = AFDB accession `AF-<UniProtAcc>-F1-model_v4`.

## Step 1 — structure-rescue rate
Counts (login node, `wc -l`):

```
wc -l all_query_ids.txt afdb_hit_ids.txt struct_nohit_ids.txt
```

`N_query` = sequence-dark proteins structurally queried; `N_afdb_hit` = those that gain an AFDB
structural hit despite NO UniRef50 sequence hit. Structure-rescue rate = `N_afdb_hit / N_query`.

## Step 2 — characterised vs hypothetical target
Whether an AFDB structural match carries an interpretable function cannot be read off the
accession alone. We use SwissProt membership as a conservative proxy for "characterised":
SwissProt entries are manually reviewed and carry a curated protein name/function, whereas the
bulk of AFDB/TrEMBL targets are themselves predicted ("Uncharacterized protein").

- SwissProt accession set from the FoldSeek SwissProt DB lookup
  `M/04_structure/fsdb_swissprot.lookup` (590,183 entries), stripping `AF-`/`-F1-model_vN`:

```
awk -F'\t' '{a=$2; sub(/^AF-/,"",a); sub(/-F1-model_v[0-9]+$/,"",a); print a}' \
    fsdb_swissprot.lookup | sort -u > sp_acc.txt          # 590,183
```

- Single streaming pass over `afdb_hits.tsv` (sbatch `c4_characterised.sh`, batch partition,
  4 cpu / 16 GB) joining each hit's target accession to `sp_acc.txt`; a rescued query is counted
  "characterised" if ANY of its AFDB targets is a reviewed SwissProt entry:

```
awk -F'\t' 'FNR==NR{S[$1]=1;next}
  {q=$1;t=$2;sub(/^AF-/,"",t);sub(/-F1-model_v[0-9]+$/,"",t);
   seen[q]=1;rows++; if(t in S){qc[q]=1;sprows++}}
  END{for(k in seen)nq++; for(k in qc)nc++;
      print "rescued_queries_unique\t"nq;
      print "rescued_queries_with_swissprot_target\t"nc;
      print "frac_rescued_characterised\t"nc/nq}' sp_acc.txt afdb_hits.tsv
```

Version note: the SwissProt FoldSeek DB is `model_v6`, `afdb_hits.tsv` targets are `model_v4`;
matching is on the version-stripped UniProt core accession (version differs, accession does not).

## Caveat
Being a SwissProt target is a *lower bound* on "interpretable function" (some TrEMBL targets are
also informative) but the structural match still must be inspected to assign function; a FoldSeek
hit to a SwissProt entry of moderate TM/lddt is a *hint*, not an annotation. An AFDB hit whose
target is itself a hypothetical/uncharacterized protein conveys no function.
