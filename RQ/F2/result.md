# RQ-F2 — How novel are the Empty Quarter BGCs vs known biosynthesis (MIBiG / BiG-FAM)?

**Status:** historical exploratory upper-bound analysis; not promoted to the active manuscript

## Answer
The archived BiG-SLICE analysis (IBEX job 47874840, completed 2026-06-29) co-clustered **119,090** antiSMASH
regions with **2,502 MIBiG 3.1** reference BGCs (BiG-SLICE v2.0.2, BIRCH, T=0.4). It produced **42,365**
gene-cluster families (GCFs) containing an EQ region, of which:

| Metric | Value |
|--------|------:|
| EQ GCFs that are **novel** (zero MIBiG members) | **42,064 / 42,365 = 99.29%** |
| EQ GCFs shared with a MIBiG reference | 301 (0.71%) |
| EQ BGCs falling in a novel GCF | 116,382 / 119,090 = **97.73%** |
| Total centroids (EQ+MIBiG) | 43,858 |

In that historical output, the no-MIBiG GCF fraction by product label was: terpene-precursor 99.95%, T3PKS 99.55%,
RiPP-like 99.54%, NRPS-like 99.28%, T1PKS 99.16%, ectoine 98.76%, lanthipeptide 98.22%, siderophore 98.77%
(`novelty_by_product.tsv`, figure `fig_bgc_novelty.{pdf,png}`). These values describe the clustering output only;
they are not a conservative estimate of biological BGC novelty.

## ⚠️ Critical caveat — fragmentation inflates apparent novelty (this is an upper bound)
**107,316 / 119,090 (90.1%) of the EQ BGCs are on contig edges (fragmented).** A partial BGC carries fewer
biosynthetic domains, so its BiG-SLICE feature vector lands far from any *complete* MIBiG reference even when
the true cluster is a known type — i.e. fragmentation **systematically inflates** the novel-GCF count. The
99.3% is therefore a fragmentation-sensitive **upper bound**, not a conservative novelty estimate. It is not
promoted to the active manuscript. Required follow-ups (planned, in methods.md):
1. **`--complete`-only rerun** (the ~11,774 non-edge BGCs) → the conservative novelty rate.
2. **BiG-FAM `--query` (Mode 2)** against ~1.2 M predicted BGCs (not just 2.5 k characterized MIBiG) — the
   stricter "novel vs anything ever predicted" test.
3. **A↔B link:** do novel BGCs carry UniRef50-dark core biosynthetic genes? (Needs an antiSMASH-CDS ↔ Bakta
   ID crosswalk; not yet wired.)

## Evidence / reproduce
- IBEX pipeline: `87_bgc_gcf_novelty.sbatch` + `87_bgc_novelty_analyze.py` (copies here and in `$M/scripts/ibex/`);
  job **47874840** (32 cpu / 300 G). BiG-SLICE DB `$M/06_metabolism/bgc_novelty/bigslice_out/result/data.db`
  (2.2 GB, stays on IBEX). Curated tables pulled to `RQ/F2/results/` (`result_summary.json`,
  `novelty_by_product.tsv`, `novelty_by_chemclass.tsv`, `gcf_composition.tsv`).
- Figure: `python3 RQ/F2/make_fig_bgc_novelty.py` → `fig_bgc_novelty.{pdf,png}`.

**References:** kautsar2020bigslice (BiG-SLICE); kautsar2021bigfam (BiG-FAM); blin2023antismash (antiSMASH);
terlouw2023mibig (MIBiG 3.0). Keys in `manuscript/references.bib`.
