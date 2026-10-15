#!/usr/bin/env python
"""
mimag_rrna_trna_clean/scripts/02_summarise.py

Run AFTER the rRNA/tRNA array job (01_rrna_trna_array.sbatch) has finished (or partially
finished — genomes with no output pair are simply skipped and reported as missing).

Reads, per representative genome:
  - results/rrna/<rep>.gff            barrnap GFF3 (5S/16S/23S rRNA hits)
  - results/trna/<rep>.trnascan.txt   tRNAscan-SE tabular output (tRNA hits by isotype)
  - results/domain_map.tsv            rep -> Bacteria/Archaea (from GTDB-Tk r232)
  - 05_refine/_checkm2_refined/quality_report.tsv   completeness/contamination (read-only)

Writes:
  - results/rrna_trna_per_rep.tsv   rep, domain, has_5S, has_16S, has_23S,
                                     n_tRNA_unique_aa, n_tRNA_total
  - results/mimag_tiers.tsv         rep, domain, completeness, contamination,
                                     has_5S, has_16S, has_23S, n_tRNA_unique_aa,
                                     old_tier (completeness/contamination only, matches the
                                     2,259 HQ / 3,179 MQ catalogue split), mimag_tier (full
                                     MIMAG rule incl. rRNA + >=18 tRNA)

Prints a tier count summary and the HQ -> MQ drop count to stdout.

Idempotent to rerun: always recomputed from the current state of results/rrna and results/trna,
so partial array progress can be summarised at any time and rerun once the array finishes.
"""
import csv
import re
import sys
from pathlib import Path

BASE = Path("/ibex/scratch/projects/c2014/rob/metagenomics/mimag_rrna_trna_clean")
REFINE = Path("/ibex/scratch/projects/c2014/rob/metagenomics/05_refine")

GENOME_LIST = BASE / "results" / "genome_list.txt"
DOMAIN_MAP = BASE / "results" / "domain_map.tsv"
RRNA_DIR = BASE / "results" / "rrna"
TRNA_DIR = BASE / "results" / "trna"
QUALITY_REPORT = REFINE / "_checkm2_refined" / "quality_report.tsv"

REVISION_OUT = BASE.parent / "06_darkmatter/p4_review_20261013/results/round4_corrections/mimag"
OUT_PER_REP = REVISION_OUT / "rrna_trna_per_rep.tsv"
OUT_TIERS = REVISION_OUT / "mimag_tiers.tsv"

RRNA_LEN_FLOOR = {"5S_rRNA": 80, "16S_rRNA": 900, "23S_rRNA": 1900}  # barrnap partial-hit sanity


def load_genome_list():
    return [l.strip() for l in GENOME_LIST.read_text().splitlines() if l.strip()]


def load_domain_map():
    dom = {}
    with open(DOMAIN_MAP) as fh:
        for line in fh:
            rep, domain = line.rstrip("\n").split("\t")
            dom[rep] = domain
    return dom


def load_quality(reps):
    q = {}
    with open(QUALITY_REPORT) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            name = row["Name"]
            if name in reps:
                q[name] = (float(row["Completeness"]), float(row["Contamination"]))
    return q


def parse_barrnap_gff(path):
    """Return dict rRNA_type -> best (longest, non-partial-preferred) length, or {} if no file."""
    hits = {}
    if not path.exists() or path.stat().st_size == 0:
        return hits
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            start, end, attrs = f[3], f[4], f[8]
            m = re.search(r"Name=([^;]+)", attrs)
            if not m:
                continue
            rtype = m.group(1)
            length = abs(int(end) - int(start)) + 1
            partial = "(partial)" in attrs or "partial=true" in attrs.lower()
            prev = hits.get(rtype)
            if prev is None or ((not partial, length) > (not prev[1], prev[0])):
                hits[rtype] = (length, partial)
    return hits


def has_full_rrna(hits, rtype):
    if rtype not in hits:
        return False
    length, partial = hits[rtype]
    floor = RRNA_LEN_FLOOR.get(rtype, 0)
    return (not partial) and length >= floor


def parse_trnascan(path):
    """Return (n_unique_isotype_aa, n_total_tRNA) from tRNAscan-SE tabular output.
    Undetermined ('Undet', 'Sup', pseudo-flagged) calls are excluded from both counts —
    MIMAG counts tRNAs covering the standard amino acids."""
    if not path.exists() or path.stat().st_size == 0:
        return 0, 0
    aa_seen = set()
    total = 0
    with open(path) as fh:
        lines = fh.readlines()
    for line in lines[3:]:  # first 3 lines are header/dashes
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 10:
            continue
        isotype = parts[4].strip()
        note = parts[9].strip().lower() if len(parts) > 9 else ""
        if "pseudo" in note:
            continue
        if isotype in ("Undet", "Sup", "") or isotype.startswith("Undet"):
            continue
        # collapse iMet/fMet/Ile2 etc. to their standard amino acid bucket
        aa = re.sub(r"[0-9]", "", isotype)
        aa = "Met" if aa in ("fMet", "iMet") else aa
        aa_seen.add(aa)
        total += 1
    return len(aa_seen), total


def mimag_tier(completeness, contamination, has5, has16, has23, n_aa):
    if completeness >= 90 and contamination <= 5 and has5 and has16 and has23 and n_aa >= 18:
        return "high"
    if completeness >= 50 and contamination <= 10:
        return "medium"
    return "low"


def old_tier(completeness, contamination):
    if completeness >= 90 and contamination <= 5:
        return "high"
    if completeness >= 50 and contamination <= 10:
        return "medium"
    return "low"


def main():
    reps = load_genome_list()
    domain = load_domain_map()
    quality = load_quality(set(reps))

    missing_quality = [r for r in reps if r not in quality]
    if missing_quality:
        print(f"WARNING: {len(missing_quality)} reps missing from CheckM2 quality_report.tsv "
              f"(first 5: {missing_quality[:5]})", file=sys.stderr)

    per_rep_rows = []
    tier_rows = []
    n_missing_rrna_out = 0
    n_missing_trna_out = 0

    for rep in reps:
        dom = domain.get(rep, "NA")
        gff_path = RRNA_DIR / f"{rep}.gff"
        trna_path = TRNA_DIR / f"{rep}.trnascan.txt"

        if not gff_path.exists():
            n_missing_rrna_out += 1
        if not trna_path.exists():
            n_missing_trna_out += 1

        rrna_hits = parse_barrnap_gff(gff_path)
        has5 = has_full_rrna(rrna_hits, "5S_rRNA")
        has16 = has_full_rrna(rrna_hits, "16S_rRNA")
        has23 = has_full_rrna(rrna_hits, "23S_rRNA")
        n_aa, n_total = parse_trnascan(trna_path)

        per_rep_rows.append([rep, dom, int(has5), int(has16), int(has23), n_aa, n_total])

        if rep in quality:
            comp, cont = quality[rep]
            old = old_tier(comp, cont)
            new = mimag_tier(comp, cont, has5, has16, has23, n_aa)
            tier_rows.append([rep, dom, comp, cont, int(has5), int(has16), int(has23),
                               n_aa, old, new])

    REVISION_OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT_PER_REP, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["rep", "domain", "has_5S", "has_16S", "has_23S",
                    "n_tRNA_unique_aa", "n_tRNA_total"])
        w.writerows(per_rep_rows)

    with open(OUT_TIERS, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["rep", "domain", "completeness", "contamination", "has_5S", "has_16S",
                    "has_23S", "n_tRNA_unique_aa", "old_tier", "mimag_tier"])
        w.writerows(tier_rows)

    # --- summary ---
    print(f"Representatives: {len(reps)}")
    print(f"Missing rRNA output (array not yet run/failed): {n_missing_rrna_out}")
    print(f"Missing tRNA output (array not yet run/failed): {n_missing_trna_out}")

    from collections import Counter
    old_counts = Counter(r[8] for r in tier_rows)
    new_counts = Counter(r[9] for r in tier_rows)
    print(f"\nOld tier (completeness/contamination only), n={len(tier_rows)}:")
    for k in ("high", "medium", "low"):
        print(f"  {k}: {old_counts.get(k, 0)}")
    print(f"\nFull MIMAG tier (+ rRNA + >=18 tRNA), n={len(tier_rows)}:")
    for k in ("high", "medium", "low"):
        print(f"  {k}: {new_counts.get(k, 0)}")

    old_high = [r for r in tier_rows if r[8] == "high"]
    dropped = [r for r in old_high if r[9] != "high"]
    dropped_to_medium = [r for r in dropped if r[9] == "medium"]
    dropped_to_low = [r for r in dropped if r[9] == "low"]
    print(f"\nOf the {len(old_high)} completeness/contamination-only 'high' representatives:")
    print(f"  drop to medium under full MIMAG rule: {len(dropped_to_medium)}")
    print(f"  drop to low under full MIMAG rule: {len(dropped_to_low)}")
    print(f"  remain high: {len(old_high) - len(dropped)}")

    print(f"\nWrote {OUT_PER_REP}")
    print(f"Wrote {OUT_TIERS}")


if __name__ == "__main__":
    main()
