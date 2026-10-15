#!/usr/bin/env python3
"""Turn the control runs into the three tables the paper needs.

  results/control_read_stats.tsv   reads before and after fastp, reads mapped to the catalogue
  results/control_profiles.tsv     every representative with relative abundance > 0 in a control,
                                   with its R232 taxonomy and its behaviour in the 550 soil samples
  results/control_contig_hits.tsv  which representatives the control CONTIGS hit, from the minimap2
                                   PAF against the same concatenated reference CoverM mapped to
  results/sntc_ungated_counts.tsv  the blank's alignment records per representative with CoverM's
                                   default 10 % covered-fraction gate switched off (step 08)

Read-only except for those three files. Run after jobs 04 and 05 have finished.
"""
import json
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats

W = "/ibex/scratch/projects/c2014/rob/metagenomics"
P = "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/contamination_clean"
READS = "/ibex/scratch/projects/c2014/EmptyQuarter_Data/soil/wgs_assemblies/processed_reads"
PREP = os.path.join(W, "rebuild_20260729/picrust2_current_20260803/coverm_current/prepared")
RES = os.path.join(P, "results")
CONTROLS = ["SNTC", "SPositive"]

# ---------------------------------------------------------------- taxonomy of the 5,438
tax = {}
for f in ("gtdbtk.bac120.summary.tsv", "gtdbtk.ar53.summary.tsv"):
    d = pd.read_csv(os.path.join(W, "05_refine/_gtdbtk_r232", f), sep="\t")
    for _, r in d.iterrows():
        tax[r["user_genome"]] = r["classification"]
assert len(tax) == 5438, f"expected 5438 representatives, got {len(tax)}"


def short(c):
    """'d__Bacteria;...;s__Escherichia coli' -> 'p__Pseudomonadota | s__Escherichia coli'."""
    parts = {t[:1]: t[3:] for t in c.split(";")}
    for rank in ("s", "g", "f", "o", "c"):
        if parts.get(rank):
            return f"p__{parts.get('p', '')} | {rank}__{parts[rank]}"
    return f"p__{parts.get('p', '')} | unresolved below phylum"


# ---------------------------------------------------------------- soil cohort context
ab = pd.read_csv(os.path.join(W, "10_thermal/abundance/rep_relative_abundance.tsv"),
                 sep="\t", index_col=0).fillna(0.0)
des = pd.read_csv(os.path.join(W, "10_thermal/abundance/sample_design.tsv"),
                  sep="\t", index_col="sample_id").reindex(ab.columns)
dna = pd.read_csv(os.path.join(W, "metadata/dna_conc.tsv"), sep="\t",
                  index_col="sample")["dna_ng_ul"].reindex(ab.columns)
mapped_pct = 100.0 - des["unmapped_percent"]
usable = mapped_pct > 0
has_dna = usable & dna.notna()

# ---------------------------------------------------------------- read stats
rows = []
for s in CONTROLS:
    j = json.load(open(os.path.join(READS, f"{s}_report.json")))
    rec = {"control": s,
           "reads_before_fastp": j["summary"]["before_filtering"]["total_reads"],
           "bases_before_fastp": j["summary"]["before_filtering"]["total_bases"],
           "reads_after_fastp_paired": j["summary"]["after_filtering"]["total_reads"],
           "bases_after_fastp_paired": j["summary"]["after_filtering"]["total_bases"],
           "q30_rate_after": j["summary"]["after_filtering"]["q30_rate"],
           "duplication_rate": j.get("duplication", {}).get("rate", np.nan)}
    log = os.path.join(RES, "coverm", s, "coverm.log")
    m = None
    if os.path.exists(log):
        for line in open(log):
            g = re.search(r"found (\d+) reads mapped out of (\d+) total \(([\d.]+)%\)", line)
            if g:
                m = g
    if m:
        rec.update(coverm_reads_total=int(m.group(2)),
                   coverm_reads_mapped=int(m.group(1)),
                   coverm_mapped_pct=float(m.group(3)))
    rows.append(rec)
read_stats = pd.DataFrame(rows)
read_stats.to_csv(os.path.join(RES, "control_read_stats.tsv"), sep="\t", index=False)
print("=== read stats ===")
print(read_stats.to_string(index=False))

# ---------------------------------------------------------------- profiles
prof_rows = []
for s in CONTROLS:
    f = os.path.join(RES, "coverm", s, "abund.tsv")
    if not os.path.exists(f):
        print(f"WARNING: {f} missing, skipped", file=sys.stderr)
        continue
    d = pd.read_csv(f, sep="\t")
    d.columns = ["genome", "rel_abund_pct"]
    unmapped = float(d.loc[d.genome == "unmapped", "rel_abund_pct"].iloc[0])
    d = d[(d.genome != "unmapped") & (d.rel_abund_pct > 0)].sort_values("rel_abund_pct",
                                                                       ascending=False)
    for _, r in d.iterrows():
        g = r.genome
        v = ab.loc[g] if g in ab.index else None
        rec = {"control": s, "genome": g,
               "rel_abund_pct": r.rel_abund_pct,
               "share_of_mapped_pct": 100.0 * r.rel_abund_pct / (100.0 - unmapped),
               "classification": tax.get(g, "NOT IN R232 SUMMARY"),
               "taxon": short(tax[g]) if g in tax else "NA",
               "is_ecoli_representative": g == "V42PRr2__single_SemiBin_33"}
        if v is not None:
            det = v > 0
            rec.update(soil_n_detected=int(det.sum()),
                       soil_prevalence_pct=100.0 * det.mean(),
                       soil_median_rel_abund_detected=float(v[det].median()) if det.any() else 0.0,
                       soil_rho_dna=stats.spearmanr(v[has_dna], dna[has_dna])[0],
                       soil_rho_mapped=stats.spearmanr(v[usable], mapped_pct[usable])[0])
        prof_rows.append(rec)
prof = pd.DataFrame(prof_rows)
prof.to_csv(os.path.join(RES, "control_profiles.tsv"), sep="\t", index=False,
            float_format="%.6g")
print("\n=== control profiles (representatives with relative abundance > 0) ===")
for s in CONTROLS:
    sub = prof[prof.control == s]
    print(f"\n-- {s}: {len(sub)} representatives > 0")
    print(sub.head(20)[["genome", "rel_abund_pct", "share_of_mapped_pct", "taxon",
                        "soil_n_detected", "soil_rho_dna"]].to_string(index=False))

# ---------------------------------------------------------------- contig hits
gdef = pd.read_csv(os.path.join(PREP, "genome_def.tsv"), sep="\t",
                   header=None, names=["genome", "contig"])
c2g = dict(zip(gdef.contig, gdef.genome))
hit_rows = []
for s in CONTROLS:
    paf = os.path.join(RES, "identify", f"{s}_contigs_vs_reps.paf")
    if not os.path.exists(paf):
        print(f"WARNING: {paf} missing, skipped", file=sys.stderr)
        continue
    cols = ["qname", "qlen", "qstart", "qend", "strand", "tname", "tlen", "tstart", "tend",
            "nmatch", "alen", "mapq"]
    d = pd.read_csv(paf, sep="\t", header=None, usecols=range(12), names=cols)
    if d.empty:
        print(f"{s}: PAF empty, no contig aligned to the catalogue")
        continue
    d["genome"] = d.tname.map(c2g)
    d["identity"] = d.nmatch / d.alen
    agg = (d.groupby("genome")
             .agg(n_alignments=("qname", "size"),
                  n_query_contigs=("qname", "nunique"),
                  aligned_bp=("alen", "sum"),
                  matching_bp=("nmatch", "sum"),
                  max_identity=("identity", "max"))
             .reset_index())
    agg["mean_identity"] = agg.matching_bp / agg.aligned_bp
    agg["control"] = s
    agg["taxon"] = agg.genome.map(lambda g: short(tax[g]) if g in tax else "NA")
    hit_rows.append(agg.sort_values("matching_bp", ascending=False))
if hit_rows:
    hits = pd.concat(hit_rows)
    hits = hits[["control", "genome", "taxon", "n_query_contigs", "n_alignments",
                 "aligned_bp", "matching_bp", "mean_identity", "max_identity"]]
    hits.to_csv(os.path.join(RES, "control_contig_hits.tsv"), sep="\t", index=False,
                float_format="%.4f")
    print("\n=== control contigs vs the catalogue (top 15 per control by matching bp) ===")
    for s in CONTROLS:
        sub = hits[hits.control == s]
        print(f"\n-- {s}: {len(sub)} representatives hit by contigs")
        print(sub.head(15).to_string(index=False))

# ---------------------------------------------------------------- blank, breadth gate off
# CoverM's count method sums to more than the library (130.4 % here), so these are alignment
# records, not reads: a read chimeric across genomes contributes to each. They are a diagnostic of
# where blank sequence piles up, not an abundance.
cf = os.path.join(RES, "coverm", "SNTC_counts", "counts.tsv")
if os.path.exists(cf):
    d = pd.read_csv(cf, sep="\t")
    d.columns = ["genome", "alignment_records"]
    d = d[d.alignment_records > 0].sort_values("alignment_records", ascending=False)
    d["pct_of_library_reads"] = 100.0 * d.alignment_records / 335190
    d["taxon"] = d.genome.map(lambda g: short(tax[g]) if g in tax else "NA")
    d["passes_default_breadth_gate"] = d.genome.isin(
        prof.loc[prof.control == "SNTC", "genome"])
    d.to_csv(os.path.join(RES, "sntc_ungated_counts.tsv"), sep="\t", index=False,
             float_format="%.4f")
    print(f"\n=== SNTC with the 10 % breadth gate off: {len(d)} representatives hold at least one "
          f"alignment record, {d.alignment_records.sum()} records = "
          f"{100 * d.alignment_records.sum() / 335190:.1f} % of the library's reads ===")
    print(d.head(12).to_string(index=False))
print("\nwrote", RES)
