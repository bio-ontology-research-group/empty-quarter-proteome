#!/usr/bin/env python
"""Turn the Pfam-A 38.2 domtblout into per-protein architecture metrics and set summaries.

Per protein
  n_dom_raw        reported domains passing the family gathering threshold
  n_dom_resolved   after greedy overlap resolution (best bit score wins; a domain
                   is dropped when its envelope overlaps a kept envelope by more
                   than 50 % of the shorter of the two)
  n_families       distinct Pfam accessions among the resolved domains
  covered_aa       residues inside the union of the reported envelopes
  covered_frac     covered_aa / protein length
  max_gap          longest contiguous run of residues outside every envelope
  has_gap80        max_gap >= 80

Per set: hit share, hit-and-gap share, covered-fraction distribution, domain-count
distribution, length distribution. Plus the top 20 Pfam families among K members and
the share of the 1,031 K families with at least one Pfam-carrying member.
"""
import collections
import gzip
import json
import os

W = "/ibex/scratch/projects/c2014/rob/metagenomics"
BASE = os.path.join(W, "04_structure/pfam_candidates_clean")
RES = os.path.join(BASE, "results")
DOMTBL = os.path.join(RES, "pfam_domtbl.txt")
PFAM_DAT = os.path.join(W, "rebuild_20260729/references/pfam_38_2/source/Pfam-A.hmm.dat.gz")
MEMBERS = os.path.join(RES, "k_family_members.tsv")
SETS = os.path.join(RES, "set_membership.tsv")

GAP = 80
OVERLAP_FRAC = 0.5


def load_pfam_dat():
    """accession (no version) -> (id, description, clan)."""
    info = {}
    acc = idn = desc = clan = None
    with gzip.open(PFAM_DAT, "rt") as fh:
        for line in fh:
            if line.startswith("#=GF ID"):
                idn = line.split(None, 2)[2].strip()
            elif line.startswith("#=GF AC"):
                acc = line.split(None, 2)[2].strip().split(".")[0]
            elif line.startswith("#=GF DE"):
                desc = line.split(None, 2)[2].strip()
            elif line.startswith("#=GF CL"):
                clan = line.split(None, 2)[2].strip()
            elif line.startswith("//"):
                if acc:
                    info[acc] = (idn, desc, clan or "")
                acc = idn = desc = clan = None
    return info


def load_domtbl():
    """protein -> (tlen, [ (env_from, env_to, score, pfam_acc, pfam_id) ])."""
    hits = collections.defaultdict(list)
    tlen = {}
    with open(DOMTBL) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 22:
                continue
            target, _, tl, qname, qacc = f[0], f[1], int(f[2]), f[3], f[4]
            dscore = float(f[13])
            env_from, env_to = int(f[19]), int(f[20])
            tlen[target] = tl
            acc = qacc.split(".")[0] if qacc != "-" else qname
            hits[target].append((env_from, env_to, dscore, acc, qname))
    return tlen, hits


def resolve(doms):
    """Greedy: keep the best-scoring domain, drop later ones overlapping it by >50 %."""
    kept = []
    for d in sorted(doms, key=lambda x: -x[2]):
        s, e = d[0], d[1]
        clash = False
        for k in kept:
            ov = min(e, k[1]) - max(s, k[0]) + 1
            if ov > 0 and ov > OVERLAP_FRAC * min(e - s + 1, k[1] - k[0] + 1):
                clash = True
                break
        if not clash:
            kept.append(d)
    return sorted(kept)


def merge(intervals):
    out = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def quantile(xs, q):
    if not xs:
        return float("nan")
    xs = sorted(xs)
    i = q * (len(xs) - 1)
    lo = int(i)
    hi = min(lo + 1, len(xs) - 1)
    return xs[lo] + (xs[hi] - xs[lo]) * (i - lo)


def wilson(k, n):
    if n == 0:
        return (float("nan"), float("nan"))
    z = 1.959963984540054
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * ((p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5)
    return ((c - h) / d * 100, (c + h) / d * 100)


def main():
    pfam = load_pfam_dat()
    tlen_hit, hits = load_domtbl()

    # every protein in the three sets, with its set label and length
    prot_set, prot_len, is_rep = {}, {}, {}
    with open(SETS) as fh:
        next(fh)
        for line in fh:
            p, s, r, L = line.rstrip("\n").split("\t")
            prot_set[p] = s
            prot_len[p] = int(L)
            is_rep[p] = r == "1"

    for p, L in tlen_hit.items():
        assert prot_len.get(p) == L, (p, L, prot_len.get(p))

    rows = {}
    for p, s in prot_set.items():
        L = prot_len[p]
        doms = hits.get(p, [])
        res = resolve(doms)
        cov = merge([(d[0], d[1]) for d in doms])
        covered = sum(e - s2 + 1 for s2, e in cov)
        gaps, prev = [], 0
        for s2, e in cov:
            if s2 - prev - 1 > 0:
                gaps.append(s2 - prev - 1)
            prev = e
        if L - prev > 0:
            gaps.append(L - prev)
        rows[p] = dict(
            protein=p,
            set=s,
            is_family_rep=int(is_rep[p]),
            length=L,
            n_dom_raw=len(doms),
            n_dom_resolved=len(res),
            n_families=len({d[3] for d in res}),
            covered_aa=covered,
            covered_frac=covered / L if L else 0.0,
            max_gap=max(gaps) if gaps else 0,
            has_gap80=int((max(gaps) if gaps else 0) >= GAP),
            pfam_accs=",".join(sorted({d[3] for d in res})) or "-",
            pfam_ids=",".join(sorted({d[4] for d in res})) or "-",
        )

    cols = [
        "protein", "set", "is_family_rep", "length", "n_dom_raw", "n_dom_resolved",
        "n_families", "covered_aa", "covered_frac", "max_gap", "has_gap80",
        "pfam_accs", "pfam_ids",
    ]
    with open(os.path.join(RES, "per_protein_pfam.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for p in sorted(rows):
            r = rows[p]
            fh.write(
                "\t".join(
                    ("%.4f" % r[c]) if c == "covered_frac" else str(r[c]) for c in cols
                )
                + "\n"
            )

    # ---------------- per-set summary ----------------
    groups = collections.OrderedDict()
    for label in ("K_members", "QnotK_reps", "controls"):
        groups[label] = [r for r in rows.values() if r["set"] == label]
    groups["K_reps_only"] = [
        r for r in rows.values() if r["set"] == "K_members" and r["is_family_rep"]
    ]

    summ = []
    for label, rs in groups.items():
        n = len(rs)
        hit = [r for r in rs if r["n_dom_raw"] > 0]
        hitgap = [r for r in hit if r["has_gap80"]]
        nohit = n - len(hit)
        cf = [r["covered_frac"] for r in hit]
        lo, hi = wilson(len(hit), n)
        lo2, hi2 = wilson(len(hitgap), n)
        summ.append(
            dict(
                set=label,
                n=n,
                median_length=quantile([r["length"] for r in rs], 0.5),
                mean_length=sum(r["length"] for r in rs) / n,
                n_with_hit=len(hit),
                pct_with_hit=100.0 * len(hit) / n,
                hit_wilson95=[round(lo, 2), round(hi, 2)],
                n_no_hit=nohit,
                n_hit_and_gap80=len(hitgap),
                pct_hit_and_gap80=100.0 * len(hitgap) / n,
                hitgap_wilson95=[round(lo2, 2), round(hi2, 2)],
                pct_of_hits_with_gap80=100.0 * len(hitgap) / len(hit) if hit else float("nan"),
                covfrac_hits_median=quantile(cf, 0.5),
                covfrac_hits_q1=quantile(cf, 0.25),
                covfrac_hits_q3=quantile(cf, 0.75),
                covfrac_hits_mean=sum(cf) / len(cf) if cf else float("nan"),
                pct_hits_covfrac_ge_0_8=100.0 * sum(1 for x in cf if x >= 0.8) / len(cf)
                if cf
                else float("nan"),
                pct_hits_covfrac_lt_0_5=100.0 * sum(1 for x in cf if x < 0.5) / len(cf)
                if cf
                else float("nan"),
                median_max_gap_hits=quantile([r["max_gap"] for r in hit], 0.5),
                dom0=sum(1 for r in rs if r["n_dom_resolved"] == 0),
                dom1=sum(1 for r in rs if r["n_dom_resolved"] == 1),
                dom2=sum(1 for r in rs if r["n_dom_resolved"] == 2),
                dom3plus=sum(1 for r in rs if r["n_dom_resolved"] >= 3),
                pct_multidomain_of_hits=100.0
                * sum(1 for r in hit if r["n_dom_resolved"] >= 2)
                / len(hit)
                if hit
                else float("nan"),
            )
        )

    keys = list(summ[0].keys())
    with open(os.path.join(RES, "set_summary.tsv"), "w") as fh:
        fh.write("\t".join(keys) + "\n")
        for s in summ:
            fh.write(
                "\t".join(
                    ("%.4f" % s[k]) if isinstance(s[k], float) else
                    (json.dumps(s[k]) if isinstance(s[k], list) else str(s[k]))
                    for k in keys
                )
                + "\n"
            )

    # ---------------- covered-fraction histogram ----------------
    edges = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0001]
    with open(os.path.join(RES, "covered_fraction_hist.tsv"), "w") as fh:
        fh.write("set\tbin_low\tbin_high\tn_hit_proteins\tpct_of_hit_proteins\n")
        for label, rs in groups.items():
            hit = [r["covered_frac"] for r in rs if r["n_dom_raw"] > 0]
            for i in range(len(edges) - 1):
                c = sum(1 for x in hit if edges[i] <= x < edges[i + 1])
                fh.write(
                    "%s\t%.1f\t%.1f\t%d\t%.2f\n"
                    % (label, edges[i], min(edges[i + 1], 1.0), c,
                       100.0 * c / len(hit) if hit else float("nan"))
                )

    # ---------------- length-stratified hit rate ----------------
    # K members are longer than the Q-not-K representatives, and a longer protein
    # has more chances of carrying a Pfam domain, so the raw hit rates are not
    # directly comparable. Stratify.
    bins = [(100, 149), (150, 199), (200, 299), (300, 399), (400, 599), (600, 10 ** 9)]
    with open(os.path.join(RES, "hit_rate_by_length.tsv"), "w") as fh:
        fh.write("length_bin\tset\tn\tn_with_hit\tpct_with_hit\tn_hit_and_gap80\tpct_hit_and_gap80\n")
        for lo, hi in bins:
            for label in ("K_members", "QnotK_reps", "controls"):
                rs = [r for r in groups[label] if lo <= r["length"] <= hi]
                if not rs:
                    continue
                h = [r for r in rs if r["n_dom_raw"] > 0]
                g = [r for r in h if r["has_gap80"]]
                fh.write(
                    "%d-%s\t%s\t%d\t%d\t%.2f\t%d\t%.2f\n"
                    % (lo, "inf" if hi > 10 ** 8 else hi, label, len(rs), len(h),
                       100.0 * len(h) / len(rs), len(g), 100.0 * len(g) / len(rs))
                )

    # ---------------- how many hits are themselves unknown-function families ----
    def is_duf(acc):
        idn, desc, _ = pfam.get(acc, ("", "", ""))
        return idn.startswith(("DUF", "UPF")) or "unknown function" in desc.lower()

    with open(os.path.join(RES, "duf_share.tsv"), "w") as fh:
        fh.write("set\tn_hit\tn_hit_all_duf\tpct_of_hits_all_duf\tn_hit_any_named\tpct_of_set_any_named\n")
        for label in ("K_members", "QnotK_reps", "controls"):
            hit = [r for r in groups[label] if r["n_dom_raw"] > 0]
            alld = [r for r in hit if all(is_duf(a) for a in r["pfam_accs"].split(","))]
            named = [r for r in hit if any(not is_duf(a) for a in r["pfam_accs"].split(","))]
            fh.write(
                "%s\t%d\t%d\t%.2f\t%d\t%.2f\n"
                % (label, len(hit), len(alld),
                   100.0 * len(alld) / len(hit) if hit else float("nan"),
                   len(named), 100.0 * len(named) / len(groups[label]))
            )

    # ---------------- top 20 Pfam families among K members ----------------
    for label in ("K_members", "QnotK_reps", "controls"):
        cnt = collections.Counter()
        for r in groups[label]:
            if r["pfam_accs"] != "-":
                for a in r["pfam_accs"].split(","):
                    cnt[a] += 1
        n_hit = sum(1 for r in groups[label] if r["n_dom_raw"] > 0)
        with open(os.path.join(RES, "top_pfam_%s.tsv" % label), "w") as fh:
            fh.write("rank\tpfam_acc\tpfam_id\tclan\tn_proteins\tpct_of_hit_proteins\tdescription\n")
            for i, (a, c) in enumerate(cnt.most_common(20), 1):
                idn, desc, clan = pfam.get(a, (a, "", ""))
                fh.write(
                    "%d\t%s\t%s\t%s\t%d\t%.2f\t%s\n"
                    % (i, a, idn, clan or "-", c,
                       100.0 * c / n_hit if n_hit else float("nan"), desc)
                )

    # ---------------- family-level view of K ----------------
    fam_members = collections.defaultdict(list)
    with open(MEMBERS) as fh:
        next(fh)
        for line in fh:
            rep, _fam, mem = line.rstrip("\n").split("\t")
            fam_members[rep].append(mem)
    fam_rows = []
    for rep, mems in fam_members.items():
        rr = [rows[m] for m in mems if m in rows]
        any_hit = sum(1 for r in rr if r["n_dom_raw"] > 0)
        any_gap = sum(1 for r in rr if r["n_dom_raw"] > 0 and r["has_gap80"])
        fam_rows.append(
            dict(
                modelled_rep=rep,
                n_members=len(mems),
                n_members_scanned=len(rr),
                n_members_with_hit=any_hit,
                family_has_hit=int(any_hit > 0),
                family_has_hit_and_gap80=int(any_gap > 0),
                rep_has_hit=int(rows[rep]["n_dom_raw"] > 0) if rep in rows else -1,
                rep_covered_frac=round(rows[rep]["covered_frac"], 4) if rep in rows else -1,
            )
        )
    fk = list(fam_rows[0].keys())
    with open(os.path.join(RES, "k_family_pfam.tsv"), "w") as fh:
        fh.write("\t".join(fk) + "\n")
        for r in sorted(fam_rows, key=lambda x: x["modelled_rep"]):
            fh.write("\t".join(str(r[k]) for k in fk) + "\n")

    nfam = len(fam_rows)
    fam_hit = sum(r["family_has_hit"] for r in fam_rows)
    fam_gap = sum(r["family_has_hit_and_gap80"] for r in fam_rows)
    rep_hit = sum(1 for r in fam_rows if r["rep_has_hit"] == 1)
    flo, fhi = wilson(fam_hit, nfam)
    fam_summary = dict(
        n_k_families=nfam,
        n_families_with_any_pfam_member=fam_hit,
        pct_families_with_any_pfam_member=round(100.0 * fam_hit / nfam, 2),
        families_wilson95=[round(flo, 2), round(fhi, 2)],
        n_families_with_pfam_member_and_gap80=fam_gap,
        pct_families_with_pfam_member_and_gap80=round(100.0 * fam_gap / nfam, 2),
        n_modelled_reps_with_pfam=rep_hit,
        pct_modelled_reps_with_pfam=round(100.0 * rep_hit / nfam, 2),
    )
    with open(os.path.join(RES, "k_family_summary.json"), "w") as fh:
        json.dump(fam_summary, fh, indent=1)

    # ---------------- console report ----------------
    print("=== per-set summary ===")
    hdr = ["set", "n", "med_len", "%hit", "%hit+gap80", "%of hits w gap80",
           "covfrac med", "%hits cov>=0.8", "%hits multi-dom"]
    print("\t".join(hdr))
    for s in summ:
        print(
            "%s\t%d\t%.0f\t%.2f\t%.2f\t%.2f\t%.3f\t%.2f\t%.2f"
            % (s["set"], s["n"], s["median_length"], s["pct_with_hit"],
               s["pct_hit_and_gap80"], s["pct_of_hits_with_gap80"],
               s["covfrac_hits_median"], s["pct_hits_covfrac_ge_0_8"],
               s["pct_multidomain_of_hits"])
        )
    print()
    print("=== K families ===")
    print(json.dumps(fam_summary, indent=1))


if __name__ == "__main__":
    main()
