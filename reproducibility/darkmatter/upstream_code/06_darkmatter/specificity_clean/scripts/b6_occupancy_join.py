#!/usr/bin/env python3
"""
§3.1 / B6 — occupancy joined to selection and to reference-absence, on the 106,042 universe.

THE DEFINITION TRAP THIS SCRIPT EXISTS TO AVOID
  There are two incompatible occupancy definitions on this basis and they differ by ~7x in the
  "restricted" class:
    * READ-MAPPING PRESENCE  = fam_presence.tsv (CoverM covf>=0.5) x sample_site_trip.tsv.
      A family is present at a site if reads from a sample at that site map to it. This is June's
      analogue and THE ONE WE USE — it is what "occurs at this site" means to a reader.
    * ASSEMBLY DETECTION = fam_recurrence.tsv col5 n_sites. A family counts at a site only if it
      ASSEMBLED there. Much stricter, and an artefact of assembly depth as much as of biology.
  The figure panel and the prose must use the SAME definition or they disagree by that factor. This
  script reports both so the discrepancy is documented rather than rediscovered.

FRAMING: DESCRIPTIVE ONLY. RQ/B6/result.md already walked back the constraint claim ("does not show
that one class is under stronger evolutionary constraint"); that walk-back stands. Report the
association, do not causally interpret it.

USAGE
  b6_occupancy_join.py <basis_results_dir> <specificity_results_dir> <out_prefix> [TAG]
  TAG selects which no-hit set defines the reference-absent class (added 2026-09-27, same pattern
  as spec_5's TAG). Unset = the original run: input dark_nohit_union_all.ids, outputs without
  suffix. TAG=with_smag = dark_nohit_union_with_smag.ids, every output suffixed _with_smag. The
  "has a hit" class is always the complement within the full 106,042-family universe.
    b6_occupancy_join.py <B> <S> results/b6 with_smag
"""
import sys, os
from statistics import median

def main():
    B, S, outp = sys.argv[1:4]
    TAG = sys.argv[4] if len(sys.argv) > 4 else ""
    SUF = f"_{TAG}" if TAG else ""
    IDS = TAG if TAG else "all"

    # ---- sample -> site, restricted to the 551 matrix columns ----
    with open(os.path.join(B, "fam_presence.tsv")) as fh:
        cols = fh.readline().rstrip("\n").split("\t")[1:]
    site_of = {}
    with open(os.path.join(B, "sample_site_trip.tsv")) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                site_of[f[0]] = f[1]
    col_site = [site_of.get(c) for c in cols]
    missing = sum(1 for s in col_site if s is None)
    sites_all = sorted({s for s in col_site if s})
    sys.stderr.write(f"columns={len(cols)} unmapped={missing} sites={len(sites_all)}\n")
    assert missing == 0, "unmapped sample column — fix the map before trusting any occupancy number"
    site_idx = {s: i for i, s in enumerate(sites_all)}
    col_site_i = [site_idx[s] for s in col_site]
    NS = len(sites_all)

    # ---- occupancy from the read-mapping matrix ----
    occ = {}          # family -> (n_samples, n_sites)
    with open(os.path.join(B, "fam_presence.tsv")) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            fam = f[0]
            seen = bytearray(NS)
            ns = 0
            for j, v in enumerate(f[1:]):
                if v == "1":
                    ns += 1
                    seen[col_site_i[j]] = 1
            occ[fam] = (ns, sum(seen))
    sys.stderr.write(f"families with occupancy: {len(occ)}\n")

    # ---- assembly-detection occupancy, for the same families only ----
    rec = {}
    with open(os.path.join(B, "fam_recurrence.tsv")) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] in occ:
                rec[f[0]] = int(f[4])          # n_sites, assembly detection
    sys.stderr.write(f"families with fam_recurrence rows: {len(rec)}\n")

    # ---- selection, dark set only ----
    dnds = {}
    with open(os.path.join(B, "selection", "dnds.tsv")) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ci = {k: i for i, k in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[ci["set"]] != "dark":
                continue
            try:
                dnds[f[ci["family_id"]]] = float(f[ci["mean_dNdS"]])
            except ValueError:
                pass
    sys.stderr.write(f"dark families with dN/dS: {len(dnds)}\n")

    # ---- rep protein lengths: needed because read-mapping occupancy is LENGTH-BIASED (see below) ----
    L, cur, n = {}, None, 0
    with open(os.path.join(B, "fam_rep_seq.fasta")) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None: L[cur] = n
                cur = line[1:].split()[0]; n = 0
            else:
                n += len(line.strip())
    if cur is not None: L[cur] = n

    # ---- reference absence (union across all named catalogues searched so far) ----
    nohit = set()
    p = os.path.join(S, f"dark_nohit_union_{IDS}.ids")
    if os.path.exists(p):
        with open(p) as fh:
            nohit = {l.split()[0] for l in fh if l.strip()}
    sys.stderr.write(f"dark families absent from every named reference: {len(nohit)}\n")

    def bin_of(n):
        if n <= 3:  return "1-3"
        if n <= 10: return "4-10"
        if n <= 30: return "11-30"
        if n <= 55: return "31-55"
        return "56-61"
    ORDER = ["1-3", "4-10", "11-30", "31-55", "56-61"]

    rows = []
    for fam, (nsamp, nsite) in occ.items():
        rows.append((fam, nsamp, nsite, dnds.get(fam), fam in nohit, rec.get(fam)))

    # ---- table 1: occupancy x selection x absence ----
    with open(outp + "_occupancy_bins" + SUF + ".tsv", "w") as out:
        out.write("site_bin\tn_families\tpct_universe\tmedian_dNdS\tn_with_dNdS\tpct_dNdS_lt1\t"
                  "pct_absent_from_all_refs\tmedian_n_samples\n")
        for b in ORDER:
            sub = [r for r in rows if bin_of(r[2]) == b]
            if not sub: continue
            d = [r[3] for r in sub if r[3] is not None]
            out.write("%s\t%d\t%.2f\t%s\t%d\t%s\t%.2f\t%d\n" % (
                b, len(sub), 100.0 * len(sub) / len(rows),
                ("%.4f" % median(d)) if d else "NA", len(d),
                ("%.1f" % (100.0 * sum(1 for x in d if x < 1) / len(d))) if d else "NA",
                100.0 * sum(1 for r in sub if r[4]) / len(sub),
                int(median([r[1] for r in sub]))))

    # ---- table 2: the definition discrepancy, documented ----
    with open(outp + "_definition_check" + SUF + ".tsv", "w") as out:
        out.write("definition\tn_families_le3_sites\tpct_universe\tmedian_n_sites\n")
        for name, vals in (("read_mapping_presence", [r[2] for r in rows]),
                           ("assembly_detection",    [r[5] for r in rows if r[5] is not None])):
            out.write("%s\t%d\t%.2f\t%d\n" % (
                name, sum(1 for v in vals if v <= 3),
                100.0 * sum(1 for v in vals if v <= 3) / len(vals), int(median(vals))))

    # ---- table 3: the §3.1 join — are reference-absent families more restricted? ----
    with open(outp + "_absence_vs_occupancy" + SUF + ".tsv", "w") as out:
        out.write("class\tn_families\tmedian_n_sites\tmedian_n_samples\tpct_le3_sites\tmedian_dNdS\n")
        for name, sel in (("absent_from_all_refs", [r for r in rows if r[4]]),
                          ("has_hit_somewhere",    [r for r in rows if not r[4]])):
            if not sel: continue
            d = [r[3] for r in sel if r[3] is not None]
            out.write("%s\t%d\t%d\t%d\t%.2f\t%s\n" % (
                name, len(sel), int(median([r[2] for r in sel])), int(median([r[1] for r in sel])),
                100.0 * sum(1 for r in sel if r[2] <= 3) / len(sel),
                ("%.4f" % median(d)) if d else "NA"))

    # ---- table 4: THE LENGTH CONFOUND. Read this before quoting anything from table 1. ----
    # CoverM presence at covf>=0.5 is easier to satisfy for a SHORT gene, so "widespread" on the
    # read-mapping definition partly measures shortness rather than biogeography. And the dN/dS
    # pipeline (fam_6c_dnds.sh) drops sequences with len(nt)<90 or protein<25aa and then drops any
    # family left with <3 usable sequences — i.e. it drops short families. The two biases compose, so
    # the per-bin dN/dS medians in table 1 are computed on differentially filtered subsets and the
    # apparent occupancy->dN/dS trend is NOT interpretable. This table is the evidence for that.
    with open(outp + "_length_confound" + SUF + ".tsv", "w") as out:
        out.write("site_bin\tn_families\tmedian_rep_len\tpct_lt100aa\tmedian_len_with_dNdS\t"
                  "median_len_dropped\tpct_with_dNdS\n")
        for b in ORDER:
            sub = [r for r in rows if bin_of(r[2]) == b and r[0] in L]
            if not sub: continue
            lh = [L[r[0]] for r in sub if r[3] is not None]
            ld = [L[r[0]] for r in sub if r[3] is None]
            out.write("%s\t%d\t%d\t%.1f\t%s\t%s\t%.1f\n" % (
                b, len(sub), int(median([L[r[0]] for r in sub])),
                100.0 * sum(1 for r in sub if L[r[0]] < 100) / len(sub),
                int(median(lh)) if lh else "NA", int(median(ld)) if ld else "NA",
                100.0 * len(lh) / len(sub)))

    # ---- table 5: the §3.1 claim, LENGTH-STRATIFIED so the confound above cannot drive it ----
    STRATA = [(0, 60, "lt60aa"), (60, 100, "60-99"), (100, 150, "100-149"),
              (150, 250, "150-249"), (250, 10**9, "ge250aa")]
    with open(outp + "_absence_occupancy_by_length" + SUF + ".tsv", "w") as out:
        out.write("len_stratum\tn_absent\tmedian_sites_absent\tn_has_hit\tmedian_sites_has_hit\t"
                  "median_samples_absent\tmedian_samples_has_hit\n")
        for lo, hi, lab in STRATA:
            a = [r for r in rows if r[4] and r[0] in L and lo <= L[r[0]] < hi]
            h = [r for r in rows if not r[4] and r[0] in L and lo <= L[r[0]] < hi]
            if not a or not h: continue
            out.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\n" % (
                lab, len(a), int(median([r[2] for r in a])), len(h), int(median([r[2] for r in h])),
                int(median([r[1] for r in a])), int(median([r[1] for r in h]))))

    sys.stderr.write("wrote %s_{occupancy_bins,definition_check,absence_vs_occupancy,"
                     "length_confound,absence_occupancy_by_length}%s.tsv\n" % (outp, SUF))

if __name__ == "__main__":
    main()
