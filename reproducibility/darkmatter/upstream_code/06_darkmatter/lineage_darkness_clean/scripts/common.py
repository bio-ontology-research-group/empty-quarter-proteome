"""Shared loaders for the lineage-darkness and size/function analyses.

Workbench W is READ ONLY. Two MAG naming bases exist:
  Wdb / genomeInformation / GTDB / rep_ko :  <sample>__self_SemiBin_<n>[_sub][.fa]
                                             <sample>__single_SemiBin_<n>[_sub][.fa]
  per_bin_dark                            :  selfsup__<sample>__SemiBin_<n>
                                             ss__<sample>__SemiBin_<n>
So per_bin_dark drops the ".fa" and the "_sub" suffix and moves the assembly
type ("self" -> "selfsup", "single" -> "ss") to the front as a prefix.
"""
import re
import pandas as pd
import numpy as np

W = "/ibex/scratch/projects/c2014/rob/metagenomics"
F_DARK = W + "/04_protein_catalog/darkmatter_clean/05_mag_dark/per_bin_dark.tsv.gz"
F_GI = W + "/05_refine/_drep/out/data_tables/genomeInformation.csv"
F_WDB = W + "/05_refine/_drep/out/data_tables/Wdb.csv"
F_BAC = W + "/05_refine/_gtdbtk_r232/gtdbtk.bac120.summary.tsv"
F_AR = W + "/05_refine/_gtdbtk_r232/gtdbtk.ar53.summary.tsv"
F_KO = W + "/04_functional/repko_clean/results/rep_ko_v2.tsv"


def strip_ext(s):
    return re.sub(r"\.(fa|fasta)$", "", str(s))


def to_dark_id(g):
    """Wdb-basis id (extension already stripped) -> per_bin_dark id, or None."""
    base = re.sub(r"_sub$", "", g)
    if "__self_SemiBin_" in base:
        return "selfsup__" + base.replace("__self_SemiBin_", "__SemiBin_")
    if "__single_SemiBin_" in base:
        return "ss__" + base.replace("__single_SemiBin_", "__SemiBin_")
    return None


def load_reps():
    """The 5,438 species representatives with quality and GTDB lineage."""
    wdb = pd.read_csv(F_WDB)
    wdb["genome"] = wdb["genome"].map(strip_ext)

    gi = pd.read_csv(F_GI)
    gi["genome"] = gi["genome"].map(strip_ext)
    gi = gi[["genome", "completeness", "contamination", "length", "N50"]]

    tax = pd.concat([
        pd.read_csv(F_BAC, sep="\t", usecols=["user_genome", "classification"]),
        pd.read_csv(F_AR, sep="\t", usecols=["user_genome", "classification"]),
    ])
    tax = tax.rename(columns={"user_genome": "genome"})
    tax["phylum"] = tax["classification"].map(lambda s: _rank(s, "p"))
    tax["genus"] = tax["classification"].map(lambda s: _rank(s, "g"))

    df = wdb[["genome"]].merge(gi, on="genome", how="left")
    df = df.merge(tax[["genome", "phylum", "genus"]], on="genome", how="left")
    df["phylum"] = df["phylum"].fillna("unclassified")
    df["genus"] = df["genus"].fillna("unclassified")
    df["log10_length"] = np.log10(df["length"])
    df["log10_N50"] = np.log10(df["N50"])
    df["corrected_size"] = df["length"] / (df["completeness"] / 100.0)
    df["log10_corrected_size"] = np.log10(df["corrected_size"])
    return df


def _rank(s, letter):
    """Pull one rank out of a GTDB classification string; blank -> unclassified."""
    if not isinstance(s, str):
        return "unclassified"
    for part in s.split(";"):
        part = part.strip()
        if part.startswith(letter + "__"):
            v = part[3:]
            return v if v else "unclassified"
    return "unclassified"


def load_dark():
    """per_bin_dark, all 37,949 bins, on its own naming basis."""
    return pd.read_csv(F_DARK, sep="\t")


def load_ko():
    """rep_ko_v2: no header, genome_id TAB comma-joined KO list."""
    ko = pd.read_csv(F_KO, sep="\t", header=None, names=["genome", "ko_list"])
    ko["ko_richness"] = ko["ko_list"].map(lambda s: len(set(str(s).split(","))))
    return ko[["genome", "ko_richness"]]


def top_phyla(df, n=8):
    return list(df["phylum"].value_counts().head(n).index)
