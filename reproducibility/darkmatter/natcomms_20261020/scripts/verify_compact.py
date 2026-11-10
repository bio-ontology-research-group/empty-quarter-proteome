#!/usr/bin/env python3
"""Verify retained ecological aggregates; this does not remap raw reads."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

HERE = Path(__file__).resolve().parents[1]
R = HERE / "results"
samples = pd.read_csv(R / "ecology_samples.tsv", sep="\t")
cells = pd.read_csv(R / "ecology_cells.tsv", sep="\t")
sites = pd.read_csv(R / "ecology_site_contrasts.tsv", sep="\t")
summary = pd.read_csv(R / "ecology_summary.tsv", sep="\t")
assert samples["sample"].is_unique and len(samples) == 551
assert samples.usable.sum() == 550
assert samples.loc[~samples.usable, "sample"].tolist() == ["F5Dr1"]
assert (samples.loc[samples.usable, "thermo_read_percent"] > 0).sum() == 482
assert not cells.duplicated(["site", "trip", "compartment"]).any()
verified = []
for metric in ["thermo_read_percent", "thermo_mapped_percent"]:
    for contrast, first in [("Deep-Surface", "Deep"), ("Root-Deep", "PlantRoot")]:
        second = "Surface" if first == "Deep" else "Deep"
        selected = cells.pivot(index=["site", "trip"], columns="compartment", values=metric)
        selected = selected[[first, second]].dropna()
        # Reconstruct site contrasts directly from compartment-cell values.
        difference = (selected[first] - selected[second]).groupby(level="site").mean()
        frozen = sites[(sites.metric == metric) & (sites.contrast == contrast)].set_index("site")
        assert np.allclose(difference.sort_index(), frozen.difference.sort_index())
        row = summary[(summary.metric == metric) & (summary.contrast == contrast)
                      & (summary.omitted_visit == "none")].iloc[0]
        assert row.n_sites == len(difference)
        assert row.n_site_visits == len(selected)
        assert row.n_positive == (difference > 0).sum()
        assert np.isclose(row.median_paired_difference, difference.median())
        assert np.isclose(row.wilcoxon_two_sided_p, wilcoxon(difference).pvalue, rtol=1e-8, atol=0)
        verified.append(row.to_dict())
primary = pd.DataFrame(verified).sort_values("wilcoxon_two_sided_p")
primary["holm_p_four_primary_tests"] = np.minimum(
    1, np.maximum.accumulate(primary.wilcoxon_two_sided_p.to_numpy() * np.arange(4, 0, -1)))
primary.to_csv(R / "ecology_primary_verified.tsv", sep="\t", index=False)
depth = primary[(primary.metric == "thermo_read_percent") & (primary.contrast == "Deep-Surface")].iloc[0]
assert (depth.n_sites, depth.n_site_visits, depth.n_positive) == (37, 69, 36)
assert round(depth.median_paired_difference, 2) == .46
assert (round(depth.ci_low, 2), round(depth.ci_high, 2)) == (.33, .86)
assert np.isclose(depth.holm_p_four_primary_tests, 4.0745362639427185e-10, rtol=1e-8, atol=0)
sensitivity = summary[(summary.contrast == "Deep-Surface") & (summary.omitted_visit != "none")]
assert len(sensitivity) == 10 and (sensitivity.median_paired_difference > 0).all()
manifest = HERE / "outputs.sha256.tsv"
if manifest.exists():
    for row in pd.read_csv(manifest, sep="\t").itertuples():
        assert hashlib.sha256((HERE / row.path).read_bytes()).hexdigest() == row.sha256, row.path
print("PASS: sample cohort, cell keys, four paired contrasts, Holm block, reported depth effect, and ten visit-omission directions.")
print("Scope: retained abundance summaries and derived statistics; no raw-read rerun.")
