#!/usr/bin/env python3
"""16S abundance weighting of the atmospheric trace-gas guild.

Maps the genome-resolved guild (MAGs with high-affinity hydB 1h/1l or form-I coxL) onto the
amplicon community to estimate what fraction of the community can "live on air".

Linkage: amplicon taxonomy is GTDB-style, so we link MAG GTDB taxonomy <-> ASV taxonomy at
genus (precise) and phylum (full-coverage fallback). Per taxon we compute guild PREVALENCE
(guild MAGs / total MAGs); community guild fraction per sample = sum over taxa of
relative_abundance(taxon) * prevalence(taxon). We also report a presence-based upper bound
(any taxon with >=1 guild MAG counts fully) and the covered fraction (abundance in taxa that
have >=1 MAG at all). Run for the full table and the 'alive' (non-relic) table.

Run from the amplicon repo env (pandas+pyarrow):
  cd ~/Documents/papers/empty-quarter-amplicon/repository && \
  uv run python ~/Documents/papers/empty-quarter-mags/code/make_16s_weighting.py
"""
import os, sys
import pandas as pd

MAGS = os.path.expanduser("~/Documents/papers/empty-quarter-mags")
AMP  = os.path.expanduser("~/Documents/papers/empty-quarter-amplicon/repository/cache")
GTDB = "/tmp/eq_gtdb"
OUT  = os.path.join(MAGS, "results", "atmospheric", "abundance")
os.makedirs(OUT, exist_ok=True)

# ---- 1. MAG GTDB lineage: mag_id -> {rank: taxon}
RANKS = ["phylum", "class", "order", "family", "genus"]
PFX = {"phylum": "p__", "class": "c__", "order": "o__", "family": "f__", "genus": "g__"}
mag_tax = {}
for s in ("bac120", "ar53"):
    p = os.path.join(GTDB, f"gtdbtk.{s}.summary.tsv")
    if not os.path.exists(p): continue
    df = pd.read_csv(p, sep="\t", usecols=["user_genome", "classification"])
    for mag, lin in zip(df.user_genome, df.classification):
        d = {}
        for tok in str(lin).split(";"):
            for r, pf in PFX.items():
                if tok.startswith(pf):
                    v = tok[len(pf):].strip()
                    d[r] = v if v else None
        mag_tax[mag] = d
print(f"[mags] GTDB lineage for {len(mag_tax)} MAGs")

# ---- 2. guild membership from atmo_classification.tsv
cls = pd.read_csv(os.path.join(MAGS, "results", "atmospheric", "atmo_classification.tsv"), sep="\t")
ha = cls[cls.is_high_affinity == 1]
h2_mags = set(ha[ha.marker == "hydB"].mag_id)
co_mags = set(ha[ha.marker == "coxL"].mag_id)
any_mags = h2_mags | co_mags
print(f"[guild] H2 {len(h2_mags)} | CO {len(co_mags)} | any {len(any_mags)} MAGs")

# ---- 3. per-taxon prevalence (n_total, n_H2, n_CO, n_any) at each rank
def prevalence(rank):
    tot, nh, nc, na = {}, {}, {}, {}
    for mag, d in mag_tax.items():
        t = d.get(rank)
        if not t: continue
        tot[t] = tot.get(t, 0) + 1
        if mag in h2_mags: nh[t] = nh.get(t, 0) + 1
        if mag in co_mags: nc[t] = nc.get(t, 0) + 1
        if mag in any_mags: na[t] = na.get(t, 0) + 1
    rows = []
    for t, n in tot.items():
        rows.append((t, n, nh.get(t, 0), nc.get(t, 0), na.get(t, 0),
                     nh.get(t, 0)/n, nc.get(t, 0)/n, na.get(t, 0)/n))
    pv = pd.DataFrame(rows, columns=[rank, "n_total", "n_H2", "n_CO", "n_any",
                                     "prev_H2", "prev_CO", "prev_any"]).set_index(rank)
    return pv

prev = {r: prevalence(r) for r in ("genus", "phylum")}
prev["genus"].sort_values("n_any", ascending=False).to_csv(os.path.join(OUT, "per_genus_prevalence.tsv"), sep="\t")
prev["phylum"].sort_values("n_any", ascending=False).to_csv(os.path.join(OUT, "per_phylum_prevalence.tsv"), sep="\t")

# ---- 4. amplicon taxonomy (ASV -> rank)
tax = pd.read_parquet(os.path.join(AMP, "taxonomy.parquet"))
tax = tax.rename(columns={c: c for c in tax.columns})  # ASV is the index
def asv_rank(rank):
    s = tax[rank].copy()
    return s.where(s.notna() & (s != "NA"), other=None)

def compartment(col):
    samp = col.split("_", 1)[1] if "_" in col else col
    return "PlantRoot" if "PRr" in samp else ("Deep" if "Dr" in samp else ("Surface" if "Sr" in samp else "?"))

# ---- 5. weighting per table x rank
def weight(table_file, label):
    ft = pd.read_parquet(os.path.join(AMP, table_file))   # ASV x sample
    samples = list(ft.columns)
    recs = []
    for rank in ("genus", "phylum"):
        amap = asv_rank(rank).reindex(ft.index)
        g = ft.groupby(amap).sum()                         # taxon x sample (NaN taxa dropped)
        colsum = g.sum(axis=0).replace(0, pd.NA)
        rel = g.div(colsum, axis=1)                        # relative abundance within assigned taxa
        pv = prev[rank]
        common = rel.index.intersection(pv.index)
        covered = rel.loc[common].sum(axis=0)              # frac of (assigned) abundance in taxa with >=1 MAG
        for metric, col in (("any", "prev_any"), ("H2", "prev_H2"), ("CO", "prev_CO")):
            w = rel.loc[common].mul(pv.loc[common, col], axis=0).sum(axis=0)   # prevalence-weighted
            for samp in samples:
                recs.append((label, rank, samp, compartment(samp), metric,
                             float(w.get(samp, 0)), float(covered.get(samp, 0))))
        # presence-based upper bound (any taxon with >=1 guild MAG)
        cap = pv.index[pv["n_any"] > 0]
        capcommon = rel.index.intersection(cap)
        pres = rel.loc[capcommon].sum(axis=0)
        for samp in samples:
            recs.append((label, rank, samp, compartment(samp), "any_presence",
                         float(pres.get(samp, 0)), float(covered.get(samp, 0))))
    return pd.DataFrame(recs, columns=["table", "rank", "sample", "compartment", "metric",
                                       "guild_fraction", "covered_fraction"])

allw = pd.concat([weight("feature_table.parquet", "full"),
                  weight("feature_table_alive.parquet", "alive")], ignore_index=True)
allw.to_csv(os.path.join(OUT, "guild_abundance_by_sample.tsv"), sep="\t", index=False)

# ---- 6. summary (mean across samples; overall + by compartment)
def summarise(df):
    rows = []
    for (table, rank, metric), grp in df.groupby(["table", "rank", "metric"]):
        rows.append((table, rank, metric, "overall", grp.guild_fraction.mean(),
                     grp.guild_fraction.std(), grp.covered_fraction.mean()))
        for comp, g2 in grp.groupby("compartment"):
            rows.append((table, rank, metric, comp, g2.guild_fraction.mean(),
                         g2.guild_fraction.std(), g2.covered_fraction.mean()))
    return pd.DataFrame(rows, columns=["table", "rank", "metric", "scope",
                                       "mean_guild_fraction", "sd", "mean_covered_fraction"])

summ = summarise(allw)
summ.to_csv(os.path.join(OUT, "guild_abundance_summary.tsv"), sep="\t", index=False)

# ---- headline print
ov = summ[summ.scope == "overall"]
def lookup(table, rank, metric, col="mean_guild_fraction"):
    r = ov[(ov.table == table) & (ov["rank"] == rank) & (ov.metric == metric)]
    return float(r[col].iloc[0]) if len(r) else float("nan")
print("\n=== HEADLINE: prevalence-weighted community fraction able to use atmospheric trace gases ===")
for table in ("full", "alive"):
    for rank in ("genus", "phylum"):
        print(f"  {table:5s}/{rank:6s}: any={lookup(table,rank,'any'):.1%}  "
              f"H2={lookup(table,rank,'H2'):.1%}  CO={lookup(table,rank,'CO'):.1%}  "
              f"(presence upper-bound any={lookup(table,rank,'any_presence'):.1%}; "
              f"covered={lookup(table,rank,'any','mean_covered_fraction'):.1%})")
print("\nwrote:", OUT)
