#!/usr/bin/env python3
"""Visit-matched depth contrasts, with sites as independent sampling units."""
import hashlib
import json
import os
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.stats import wilcoxon

W = Path("/ibex/scratch/projects/c2014/rob/metagenomics")
HERE = Path(__file__).resolve().parents[1]
OUT = HERE / "results"
HEAT = W / "06_darkmatter/heat_stress_20261016/results"
assert os.environ.get("SLURM_JOB_ID"), "Run on an allocated Slurm worker"
OUT.mkdir(exist_ok=True)
sources = [HEAT / x for x in ("abundance.npz", "metadata.tsv", "genome_index.tsv")]
data = np.load(sources[0])
meta = pd.read_csv(sources[1], sep="\t", dtype={"site": str}).set_index("sample")
idx = pd.read_csv(sources[2], sep="\t").set_index("genome_id")
samples, genomes, raw = data["samples"], data["genomes"], data["raw"]
meta = meta.loc[samples].copy()
idx = idx.loc[genomes]
assert raw.shape == (551, 5438)
assert idx.index.is_unique and meta.index.is_unique
mask = idx.gtdb_classification.str.contains("p__Thermoplasmatota;").to_numpy()
assert mask.sum() == 237
total = raw.sum(axis=1)
meta["usable"] = total > 0
meta["thermo_read_percent"] = raw[:, mask].sum(axis=1)
meta["thermo_mapped_percent"] = np.divide(
    meta.thermo_read_percent.to_numpy() * 100, total,
    out=np.full(len(total), np.nan), where=total > 0)
meta["thermo_species_detected"] = (raw[:, mask] > 0).sum(axis=1)
meta["total_mapped_percent"] = total
assert list(meta.index[~meta.usable]) == ["F5Dr1"]
columns = ["site", "trip", "compartment", "usable", "thermo_read_percent",
           "thermo_mapped_percent", "thermo_species_detected", "total_mapped_percent"]
meta[columns].to_csv(OUT / "ecology_samples.tsv", sep="\t")
use = meta[meta.usable]
metrics = ["thermo_read_percent", "thermo_mapped_percent"]
# Average replicate soil samples within a site, visit and compartment.
cells = use.groupby(["site", "trip", "compartment"])[metrics].mean().reset_index()
cells.to_csv(OUT / "ecology_cells.tsv", sep="\t", index=False)
summary, pairs_all, sites_all = [], [], []
for metric in metrics:
    wide = cells.pivot(index=["site", "trip"], columns="compartment", values=metric)
    for contrast, a, b in [("Deep-Surface", "Deep", "Surface"),
                           ("Root-Deep", "PlantRoot", "Deep")]:
        # The source metadata uses PlantRoot for root-adjacent samples.
        if a not in wide:
            raise ValueError(f"Missing compartment {a}: {list(wide.columns)}")
        paired = wide[[a, b]].dropna().copy()
        paired["difference"] = paired[a] - paired[b]
        paired = paired.reset_index()
        pairs_all.append(paired.assign(metric=metric, contrast=contrast))
        for omitted in ["none"] + sorted(paired.trip.unique().tolist()):
            p = paired if omitted == "none" else paired[paired.trip != omitted]
            site = p.groupby("site")[[a, b, "difference"]].mean()
            if omitted == "none":
                sites_all.append(site.reset_index().assign(metric=metric, contrast=contrast))
            differences = site.difference.to_numpy()
            rng = np.random.default_rng(20261020)
            bootstrap = np.median(rng.choice(differences, (10000, len(site)), replace=True), axis=1)
            ci = np.quantile(bootstrap, [.025, .975])
            summary.append(dict(metric=metric, contrast=contrast, omitted_visit=omitted,
                                n_site_visits=len(p), n_sites=len(site),
                                n_positive=int((differences > 0).sum()),
                                median_a=float(site[a].median()), median_b=float(site[b].median()),
                                median_paired_difference=float(np.median(differences)),
                                ci_low=float(ci[0]), ci_high=float(ci[1]),
                                wilcoxon_two_sided_p=float(wilcoxon(differences).pvalue)))
summary = pd.DataFrame(summary)
# One declared block: two compartment contrasts x two abundance denominators.
primary = summary.index[summary.omitted_visit == "none"]
ordered = summary.loc[primary].sort_values("wilcoxon_two_sided_p").index
adjusted = np.maximum.accumulate(
    summary.loc[ordered, "wilcoxon_two_sided_p"].to_numpy() * np.arange(4, 0, -1))
summary.loc[ordered, "holm_p_four_primary_tests"] = np.minimum(adjusted, 1)
summary.to_csv(OUT / "ecology_summary.tsv", sep="\t", index=False)
pd.concat(pairs_all).to_csv(OUT / "ecology_site_visit_pairs.tsv", sep="\t", index=False)
pd.concat(sites_all).to_csv(OUT / "ecology_site_contrasts.tsv", sep="\t", index=False)
manifest = [{"path": str(p), "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
             "bytes": p.stat().st_size} for p in sources]
pd.DataFrame(manifest).to_csv(OUT / "ecology_input_manifest.tsv", sep="\t", index=False)
report = dict(job=os.environ["SLURM_JOB_ID"], host=platform.node(),
              numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__,
              usable_samples=len(use), detected_samples=int((use.thermo_read_percent > 0).sum()),
              median_read_percent=float(use.thermo_read_percent.median()),
              median_detected_species=float(use.thermo_species_detected.median()),
              replicate_reduction="arithmetic mean within site x visit x compartment",
              visit_reduction="equal-weight mean of paired visit differences within site",
              inference="two-sided Wilcoxon on site differences; 10000 site bootstrap samples",
              seed=20261020)
(OUT / "ecology_run.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
