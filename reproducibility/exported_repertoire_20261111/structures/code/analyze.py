#!/usr/bin/env python3
"""Laptop-side summary of analysis4_structures (no cluster compute).
Reads results/ (step1 tables, alignment depths, ORIX metrics and Foldseek hits) and writes the result tables.
Confidence values are on a 0-1 scale. Qualification gate (paper): >= 100 residues and mean pLDDT >= 0.70.
PDB match: query-normalised TM-score >= 0.5 for the whole chain or any confident-region crop.
AFDB match: any hit with E < 1e-3 (3Di+AA alignment; the ORIX AFDB v6 copy has no coordinates)."""
import json, re, glob, sys
from pathlib import Path
import numpy as np, pandas as pd
R = Path(__file__).parent / "results"
O = R / "final"; O.mkdir(exist_ok=True)
S = {}

def cat(pattern):
    fs = sorted(glob.glob(str(R / "orix" / pattern)))
    out = []
    for f in fs:
        d = pd.read_csv(f, sep="\t"); d["prediction_set"] = Path(f).name.split("_metrics.")[0]
        out.append(d)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()

# ---------------- monomers ----------------
reps = pd.read_csv(R / "step1/representatives.tsv", sep="\t")
dep = pd.read_csv(R / "monomer_alignment_depth_own.tsv", sep="\t")
mm = cat("monomer*_metrics.shard_*.tsv")
mm = mm[mm.status == "ok"].drop_duplicates("id")
mono = reps.merge(dep.drop(columns=["in_981", "rep_gene", "tier", "signal_peptide_cut", "mature_aa"]),
                  on="family_rep", how="left").merge(mm.rename(columns={"id": "label"}), on="label", how="left")
mono["modelled"] = mono.status.eq("ok")
mono["qualified"] = mono.modelled & (mono.n_res >= 100) & (mono.mean_plddt >= 0.70)
mono["qualified_model_1_only"] = mono.modelled & (mono.n_res >= 100) & (mono.mean_plddt_model_1 >= 0.70)
mono["not_modelled_reason"] = np.where(mono.modelled, "", np.where(~mono.usable, "no_usable_representative_" + mono.tier,
                                       np.where(mono.mature_aa > 2500, "longer_than_2500_residues", "missing_model")))
f = mono[mono.in_981]
S["monomer"] = {
    "n_families": int(len(f)), "n_usable_representative": int(f.usable.sum()), "n_modelled": int(f.modelled.sum()),
    "not_modelled": f[~f.modelled].not_modelled_reason.value_counts().to_dict(),
    "n_modelled_ge100_residues": int((f.modelled & (f.n_res >= 100)).sum()),
    "n_modelled_lt100_residues": int((f.modelled & (f.n_res < 100)).sum()),
    "n_qualified": int(f.qualified.sum()),
    "n_qualified_model_1_only": int(f.qualified_model_1_only.sum()),
    "n_mean_plddt_ge070_any_length": int((f.modelled & (f.mean_plddt >= 0.70)).sum()),
    "median_mean_plddt": float(f.mean_plddt.median()), "median_ptm": float(f.ptm.median()),
    "median_frac_res_ge70": float(f.frac_res_plddt_ge70.median()),
    "n_ptm_ge05_among_ge100": int((f.modelled & (f.n_res >= 100) & (f.ptm >= 0.5)).sum()),
}
bins = [0, 10, 30, 100, 300, 1e9]; labs = ["<10", "10-30", "30-100", "100-300", ">=300"]
g = f[f.modelled & (f.n_res >= 100)].copy()
g["neff_bin"] = pd.cut(g.neff80_own, bins, labels=labs, right=False)
S["monomer"]["qualified_by_neff_own_ge100res"] = {
    k: {"n": int(len(v)), "qualified": int(v.qualified.sum()), "fraction": round(float(v.qualified.mean()), 3),
        "median_plddt": round(float(v.mean_plddt.median()), 3)} for k, v in g.groupby("neff_bin", observed=True)}
g["len_bin"] = pd.cut(g.n_res, [100, 150, 250, 400, 700, 1e9], labels=["100-149", "150-249", "250-399", "400-699", ">=700"], right=False)
S["monomer"]["qualified_by_length"] = {k: {"n": int(len(v)), "qualified": int(v.qualified.sum())}
                                       for k, v in g.groupby("len_bin", observed=True)}
g["sp"] = g.modal_sp_class
S["monomer"]["qualified_by_signal_peptide_class"] = {k: {"n": int(len(v)), "qualified": int(v.qualified.sum())}
                                                     for k, v in g.groupby("sp")}
from scipy.stats import spearmanr
rho = spearmanr(g.neff80_own, g.mean_plddt)
S["monomer"]["spearman_neff_own_vs_plddt_ge100res"] = [round(float(rho[0]), 3), float(rho[1])]
dd = f[f.usable]
S["alignment_depth_981"] = {
    "n": int(len(dd)),
    "median_nr95_sequences_eq_only": float(dd.n_nr95_eq_only.median()), "median_nr95_sequences_all_1957": float(dd.n_nr95_own.median()),
    "median_neff80_eq_only": float(dd.neff80_eq_only.median()), "median_neff80_all_1957": float(dd.neff80_own.median()),
    "iqr_neff80_eq_only": [float(x) for x in dd.neff80_eq_only.quantile([.25, .75])],
    "iqr_neff80_all_1957": [float(x) for x in dd.neff80_own.quantile([.25, .75])],
    "median_ratio_neff_all_over_eq": float((dd.neff80_own / dd.neff80_eq_only).median()),
    "n_ratio_ge_1.5": int(((dd.neff80_own / dd.neff80_eq_only) >= 1.5).sum()),
    "n_ratio_ge_2": int(((dd.neff80_own / dd.neff80_eq_only) >= 2).sum()),
    "n_neff_lt30_eq_only": int((dd.neff80_eq_only < 30).sum()), "n_neff_lt30_all_1957": int((dd.neff80_own < 30).sum()),
    "n_neff_ge100_eq_only": int((dd.neff80_eq_only >= 100).sum()), "n_neff_ge100_all_1957": int((dd.neff80_own >= 100).sum()),
    "median_genomes_in_alignment": float(dd.n_genomes_own.median()),
    "external_database": "none in the folded alignments (own 3,617,011 proteins only)"}

# ---------------- Foldseek ----------------
C1 = "query,target,alntmscore,qtmscore,ttmscore,prob,lddt,evalue,bits,fident,qcov,tcov,qlen,tlen,alnlen,qstart,qend,tstart,tend,theader".split(",")
C2 = "query,target,evalue,bits,fident,qcov,tcov,qlen,tlen,alnlen,qstart,qend,tstart,tend,theader".split(",")
def rd(name, cols):
    p = R / "orix_foldseek_top" / name
    if not p.exists() or p.stat().st_size == 0: return pd.DataFrame(columns=cols)
    return pd.read_csv(p, sep="\t", header=None, names=cols, quoting=3)
have_fs = (R / "orix_foldseek_top/afdb_crops.tsv").exists()
if have_fs:
    pw, pc = rd("pdb_whole.tsv", C1), rd("pdb_crops.tsv", C1)
    aw, ac = rd("afdb_whole.tsv", C2), rd("afdb_crops.tsv", C2)
    for d_ in (pw, pc, aw, ac):
        d_["query"] = d_["query"].astype(str).str.replace(r"\.pdb$", "", regex=True)
    pc["label"] = pc["query"].str.split("__crop").str[0]; ac["label"] = ac["query"].str.split("__crop").str[0]
    pw["label"] = pw["query"]; aw["label"] = aw["query"]
    def desc(h, t):
        h = str(h); return h[len(str(t)):].strip() if h.startswith(str(t)) else h
    for d_ in (pw, pc, aw, ac): d_["description"] = [desc(h, t) for h, t in zip(d_.theader, d_.target)]
    bw = pw.sort_values("qtmscore", ascending=False).drop_duplicates("label")
    bc = pc.sort_values("qtmscore", ascending=False).drop_duplicates("label")
    ba = aw.sort_values("evalue").drop_duplicates("label"); bac = ac.sort_values("evalue").drop_duplicates("label")
    mono = mono.merge(bw[["label", "target", "qtmscore", "ttmscore", "qcov", "evalue", "description"]].rename(
        columns=lambda c: c if c == "label" else "pdb_whole_" + c), on="label", how="left")
    mono = mono.merge(bc[["label", "query", "target", "qtmscore", "qcov", "description"]].rename(
        columns=lambda c: c if c == "label" else "pdb_crop_" + c), on="label", how="left")
    mono = mono.merge(ba[["label", "target", "evalue", "fident", "qcov", "description"]].rename(
        columns=lambda c: c if c == "label" else "afdb_whole_" + c), on="label", how="left")
    mono = mono.merge(bac[["label", "query", "target", "evalue", "description"]].rename(
        columns=lambda c: c if c == "label" else "afdb_crop_" + c), on="label", how="left")
    mono["pdb_match_whole"] = mono.pdb_whole_qtmscore.fillna(0) >= 0.5
    mono["pdb_match_crop"] = mono.pdb_crop_qtmscore.fillna(0) >= 0.5
    mono["pdb_match"] = mono.pdb_match_whole | mono.pdb_match_crop
    mono["afdb_match"] = mono.afdb_whole_evalue.fillna(1) < 1e-3
    mono["afdb_match_crop"] = mono.afdb_crop_evalue.fillna(1) < 1e-3
    # keyword fold/function classes from the best PDB and AFDB descriptions (heuristic; for orientation only)
    KW = [("pilin_or_flagellin_like", r"pilin|pilus|fimbri|archaell|flagell|type iv|type 4 pil"),
          ("S-layer_or_cell_surface", r"s-layer|surface layer|cell surface|cell wall|surface glycoprotein|surface protein"),
          ("adhesin_Ig_or_beta_sandwich", r"adhesin|immunoglobulin|ig-like|fibronectin|cadherin|invasin|intimin|big[_ -]|pkd|cohesin|dockerin|carbohydrate[- ]binding|cbm|lectin|thrombospondin"),
          ("solute_binding_or_transport", r"binding protein|periplasmic|abc transporter|substrate-binding|solute|transporter|permease|porin|channel"),
          ("protease_or_peptidase", r"proteas|peptidas|subtilis|thermolysin|proteinase"),
          ("glycosyl_hydrolase_or_transferase", r"glycos|glucan|chitin|cellul|xylan|amylase|galactos|mannos|sialid|lyase|pectat"),
          ("other_enzyme", r"ase\b|enzyme|reductase|oxidase|kinase|hydrolase|transferase|synthase|dehydrogenase|esterase|lipase|phosphatase|nuclease|cytochrome"),
          ("DUF_or_uncharacterised", r"uncharacteri|hypothetical|duf\d|unknown function|domain of unknown|putative")]
    def kw(s):
        s = str(s).lower()
        if s in ("nan", ""): return "no_hit"
        for k, pat in KW:
            if re.search(pat, s): return k
        return "other_named"
    mono["pdb_best_description"] = np.where(mono.pdb_whole_qtmscore.fillna(0) >= mono.pdb_crop_qtmscore.fillna(0),
                                            mono.pdb_whole_description, mono.pdb_crop_description)
    mono["pdb_fold_class_keyword"] = [kw(x) if m else "no_pdb_match" for x, m in zip(mono.pdb_best_description, mono.pdb_match)]
    # AFDB (UniProt protein names of the best hit): ordered keyword classes, beta-sandwich domain names first
    KA = [("uncharacterised_or_DUF", r"uncharacteri|hypothetical|duf\d|unknown function"),
          ("Ig-like_or_beta-sandwich_domain", r"new3|cardb|pkd|ig-like|big[_ -]|fibronectin|carboxypeptidase regulatory|cadherin|invasin|p/homo b|peptidase c-terminal archaeal|carbohydrate[- ]binding|cbm|cohesin|dockerin|lectin|beta[- ]sandwich|ob domain|prealbumin|transthyretin|emp24|bsupi|intracellular proteinase inhibitor"),
          ("pilin_or_archaellin", r"pilin|pilus|fimbri|archaell|flagell|class iii signal"),
          ("S-layer_or_cell-surface_anchor", r"s-layer|surface layer|cell surface|surface glycoprotein|pgf-cterm|lpxtg|cell wall anchor|sortase|archaeosortase"),
          ("lipoprotein_named_only", r"^lipoprotein|lipoprotein$"),
          ("protease_or_peptidase", r"proteas|peptidas|subtilis|thermolysin|proteinase|matrixin"),
          ("glycoside_hydrolase_or_carbohydrate_active", r"glycos|glucan|chitin|cellul|xylan|amylase|galactos|mannos|sialid|pectat|beta helix"),
          ("solute_binding_or_transport", r"binding protein|periplasmic|abc transporter|substrate-binding|solute|transporter|permease|porin"),
          ("other_enzyme", r"ase\b|ase domain|reductase|oxidase|cytochrome|cupredoxin|plastocyanin|blue copper"),
          ("repeat_or_other_named", r".")]
    def kwa(s_):
        s_ = re.sub(r" os=.*", "", str(s_).lower())
        for k, pat in KA:
            if re.search(pat, s_): return k
        return "repeat_or_other_named"
    mono["afdb_class_keyword"] = [kwa(x) if m else "no_afdb_match" for x, m in zip(mono.afdb_whole_description, mono.afdb_match)]
    q = mono[mono.in_981 & mono.qualified]
    S["foldseek_qualified_981"] = {
        "n_qualified": int(len(q)), "pdb_match_whole": int(q.pdb_match_whole.sum()),
        "pdb_match_crop": int(q.pdb_match_crop.sum()), "pdb_match_either": int(q.pdb_match.sum()),
        "afdb_match_E_lt_1e-3": int(q.afdb_match.sum()),
        "pdb_or_afdb": int((q.pdb_match | q.afdb_match).sum()), "neither": int((~q.pdb_match & ~q.afdb_match).sum()),
        "pdb_only": int((q.pdb_match & ~q.afdb_match).sum()), "afdb_only": int((~q.pdb_match & q.afdb_match).sum()),
        "pdb_fold_class_keyword_counts": q.pdb_fold_class_keyword.value_counts().to_dict(),
        "afdb_class_keyword_counts": q.afdb_class_keyword.value_counts().to_dict(),
        "queries_at_pdb_max_seqs_cap_1000_whole": int((pd.read_csv(R / "orix_foldseek_top/pdb_whole.nhits.tsv", sep="\t", header=None)[1] >= 1000).sum()),
        "note": "tables hold the five best hits per query (by qTM for PDB, by E-value for AFDB)"}
    nq = mono[mono.in_981 & mono.modelled & ~mono.qualified]
    S["foldseek_not_qualified_981_exploratory"] = {
        "n": int(len(nq)), "n_with_confident_region_crop": int((nq.n_confident_regions > 0).sum()),
        "pdb_match_crop": int(nq.pdb_match_crop.sum()), "pdb_match_whole": int(nq.pdb_match_whole.sum())}
    top = []
    for nm, d_, key, asc in (("pdb_whole", pw, "qtmscore", False), ("pdb_crop", pc, "qtmscore", False),
                             ("afdb_whole", aw, "evalue", True), ("afdb_crop", ac, "evalue", True)):
        t = d_.sort_values(key, ascending=asc).groupby("query").head(5).copy()
        t.insert(0, "search", nm); top.append(t)
    top = pd.concat(top, ignore_index=True)
    ml = mono[mono.label.notna()].set_index("label")
    top.insert(1, "family_rep", top.label.map(ml.family_rep)); top.insert(2, "qualified", top.label.map(ml.qualified))
    top.drop(columns=["theader"]).to_csv(O / "foldseek_top_hits.tsv", sep="\t", index=False)
# earlier models: only the 47-family audit table on the laptop lists per-gene earlier models (full-length
# sequences); matched by exact gene identifier of the representative, otherwise not guessed.
ev = Path("/home/leechuck/Documents/papers/empty-quarter/mechanism-20261109/exp4_surface/results/followup_v2_review/exp4_47_member_evidence.tsv")
if ev.exists():
    e = pd.read_csv(ev, sep="\t", usecols=["member_gene", "has_colabfold_model", "colabfold_plddt", "has_esmfold_model", "esmfold_plddt"])
    e = e.drop_duplicates("member_gene").rename(columns={"member_gene": "rep_gene", "colabfold_plddt": "earlier_colabfold_plddt_full_length",
                                                          "esmfold_plddt": "earlier_esmfold_plddt_full_length"})
    mono = mono.merge(e, on="rep_gene", how="left")
    mt_ = mono[mono.in_981 & mono.modelled & (mono.has_colabfold_model.astype(str) == "True")]
    me_ = mono[mono.in_981 & mono.modelled & (mono.has_esmfold_model.astype(str) == "True")]
    S["earlier_model_match"] = {
        "source": str(ev), "n_981_representatives_in_audit_table": int(mono[mono.in_981].rep_gene.isin(e.rep_gene).sum()),
        "n_with_earlier_colabfold_model": int(len(mt_)), "n_with_earlier_esmfold_model": int(len(me_)),
        "colabfold_pairs": mt_[["family_rep", "rep_gene", "earlier_colabfold_plddt_full_length", "mean_plddt"]].to_dict("records"),
        "esmfold_pairs": me_[["family_rep", "rep_gene", "earlier_esmfold_plddt_full_length", "mean_plddt"]].to_dict("records"),
        "note": "the 9,644-model ColabFold rescue table is not on the laptop; no cohort-level change is computed"}
mono.to_csv(O / "per_family_monomer_metrics.tsv", sep="\t", index=False)

# ---------------- complexes ----------------
pairs = pd.read_csv(R / "pair_alignment_depth_own.tsv", sep="\t")
cm = cat("complex_*_metrics.shard_*.tsv")
if len(cm):
    cm = cm[cm.status == "ok"].drop_duplicates("id").rename(columns={"id": "pair_id"})
    cx = pairs.merge(cm.drop(columns=["len_A", "len_B"]), on="pair_id", how="left")
    cx["nan_output"] = cx.iptm.isna() & (cx.status_y.eq("ok") if "status_y" in cx else False)
    cx["modelled"] = cx.iptm.notna()
    real, ctrl = cx[(cx.type == "real") & cx.modelled].copy(), cx[(cx.type == "control") & cx.modelled].copy()
    thr = {"iptm_p95": float(ctrl.iptm.quantile(0.95)), "iptm_p99": float(ctrl.iptm.quantile(0.99)), "iptm_max": float(ctrl.iptm.max()),
           "iptm_mean_models_p95": float(ctrl.iptm_mean_models.quantile(0.95)),
           "iptm_min_models_p95": float(ctrl.iptm_min_models.quantile(0.95)),
           "ranking_confidence_p95": float(ctrl.ranking_confidence.quantile(0.95)),
           "interface_pae_p05": float(ctrl.interface_pae.quantile(0.05)),
           "min_interchain_pae_p05": float(ctrl.min_interchain_pae.quantile(0.05)),
           "n_ca_contacts_p95": float(ctrl.n_ca_contacts_8A.quantile(0.95))}
    for d_ in (real, ctrl):
        d_["above_ctrl_p95_iptm"] = d_.iptm > thr["iptm_p95"]
        d_["above_ctrl_p99_iptm"] = d_.iptm > thr["iptm_p99"]
        d_["above_ctrl_max_iptm"] = d_.iptm > thr["iptm_max"]
        d_["above_ctrl_p95_mean_iptm"] = d_.iptm_mean_models > thr["iptm_mean_models_p95"]
        d_["above_ctrl_p95_min_iptm"] = d_.iptm_min_models > thr["iptm_min_models_p95"]
        d_["below_ctrl_p05_interface_pae"] = d_.interface_pae < thr["interface_pae_p05"]
        d_["strict_support"] = (d_.above_ctrl_p95_iptm & d_.above_ctrl_p95_mean_iptm & d_.above_ctrl_p95_min_iptm
                                & d_.below_ctrl_p05_interface_pae.fillna(False))
        d_["ctrl_percentile_iptm"] = [(ctrl.iptm < x).mean() for x in d_.iptm]
    from scipy.stats import mannwhitneyu, wilcoxon
    mw = mannwhitneyu(real.iptm, ctrl.iptm, alternative="greater")
    mt = ctrl[["matched_real_pair", "iptm"]].merge(real[["pair_id", "iptm"]], left_on="matched_real_pair", right_on="pair_id",
                                                    suffixes=("_ctrl", "_real"))
    wx = wilcoxon(mt.iptm_real, mt.iptm_ctrl, alternative="greater") if len(mt) > 5 else (np.nan, np.nan)
    n = len(real)
    S["complex"] = {
        "n_real_pairs_defined": int((pairs.type == "real").sum()), "n_real_modelled": int(n),
        "n_control_defined": int((pairs.type == "control").sum()), "n_control_modelled": int(len(ctrl)),
        "not_modelled_real": sorted(set(pairs[pairs.type == "real"].pair_id) - set(real.pair_id)),
        "not_modelled_real_longer_than_1800_residues": sorted(pairs[(pairs.type == "real") & (pairs.len_A + pairs.len_B > 1800)].pair_id),
        "not_modelled_real_numerical_failure_nan_models": sorted(cx[cx.nan_output].pair_id),
        "real_by_msa_stratum": {k: {"n": int(len(v)), "above_ctrl_p95_iptm": int(v.above_ctrl_p95_iptm.sum()), "strict": int(v.strict_support.sum())} for k, v in real.groupby("prediction_set")},
        "control_by_msa_stratum": {k: {"n": int(len(v)), "iptm_p95": float(v.iptm.quantile(.95)), "iptm_ge_0.5": int((v.iptm >= 0.5).sum()), "strict": int(v.strict_support.sum())} for k, v in ctrl.groupby("prediction_set")},
        "strict_by_partner_class": real[real.strict_support].partner_class_eq_members.value_counts().to_dict(),
        "strict_by_partner_signalp": real[real.strict_support].partner_sp_class.value_counts().to_dict(),
        "real_by_partner_signalp": real.partner_sp_class.value_counts().to_dict(),
        "strict_focal_families": int(real[real.strict_support].focal_family.nunique()),
        "strict_median_paired_depth": float(real[real.strict_support].n_paired_nr95.median()),
        "all_real_median_paired_depth": float(real.n_paired_nr95.median()),
        "control_thresholds": thr,
        "median_iptm_real": float(real.iptm.median()), "median_iptm_control": float(ctrl.iptm.median()),
        "mean_iptm_real": float(real.iptm.mean()), "mean_iptm_control": float(ctrl.iptm.mean()),
        "mannwhitney_real_gt_control_p": float(mw[1]),
        "matched_pairs_n": int(len(mt)), "matched_pairs_real_minus_control_median": float((mt.iptm_real - mt.iptm_ctrl).median()),
        "matched_wilcoxon_p": float(wx[1]),
        "real_above_ctrl_p95_iptm": int(real.above_ctrl_p95_iptm.sum()), "expected_if_null_5pct": round(0.05 * n, 1),
        "real_above_ctrl_p99_iptm": int(real.above_ctrl_p99_iptm.sum()),
        "real_above_ctrl_max_iptm": int(real.above_ctrl_max_iptm.sum()),
        "real_strict_support_all_four_criteria": int(real.strict_support.sum()),
        "control_strict_support_all_four_criteria": int(ctrl.strict_support.sum()),
        "real_iptm_ge_0.5": int((real.iptm >= 0.5).sum()), "control_iptm_ge_0.5": int((ctrl.iptm >= 0.5).sum()),
        "real_iptm_ge_0.6": int((real.iptm >= 0.6).sum()), "control_iptm_ge_0.6": int((ctrl.iptm >= 0.6).sum()),
        "real_iptm_ge_0.8": int((real.iptm >= 0.8).sum()), "control_iptm_ge_0.8": int((ctrl.iptm >= 0.8).sum()),
        "above_p95_by_partner_class": real.groupby("partner_class_eq_members").above_ctrl_p95_iptm.agg(["size", "sum"]).to_dict("index"),
        "focal_families_with_pair_above_ctrl_p95": int(real[real.above_ctrl_p95_iptm].focal_family.nunique()),
    }
    pf = R / "queries_pfam.domtbl"
    if pf.exists():
        rows = [l.split() for l in open(pf) if not l.startswith("#")]
        pfam = pd.DataFrame([(r[0], r[3], r[4], float(r[12])) for r in rows], columns=["qid", "pfam", "acc", "ievalue"])
        pfs = pfam.sort_values("ievalue").groupby("qid").apply(lambda x: ";".join(dict.fromkeys(x.pfam + "(" + x.acc + ")")), include_groups=False)
        real["focal_pfam"] = real.qid_A.map(pfs); real["partner_pfam"] = real.qid_B.map(pfs)
        ctrl["focal_pfam"] = ctrl.qid_A.map(pfs); ctrl["partner_pfam"] = ctrl.qid_B.map(pfs)
    real.sort_values("iptm", ascending=False).to_csv(O / "complex_real_pair_metrics.tsv", sep="\t", index=False)
    ctrl.sort_values("iptm", ascending=False).to_csv(O / "complex_control_pair_metrics.tsv", sep="\t", index=False)
    S["complex"]["paper_modules"] = real[real.source != "conserved_adjacent"][
        ["pair_id", "source", "focal_gene", "partner_gene", "len_A", "len_B", "n_paired_nr95", "iptm", "ptm",
         "iptm_all_models", "interface_pae", "interface_plddt", "n_ca_contacts_8A", "ctrl_percentile_iptm",
         "above_ctrl_p95_iptm", "strict_support"]].to_dict("records")
(O / "summary.json").write_text(json.dumps(S, indent=1, default=str))
print(json.dumps(S, indent=1, default=str))
