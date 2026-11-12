#!/usr/bin/env python3
"""Stage 2: per-protein evidence -> per-family export route, anchor class, coarse role class; summaries; checks.

All calls are PREDICTIONS (SignalP 6.0 class, TMbed topology, HMM hits); role classes are rule-based labels,
not biochemical functions. The rules are in RULES below and are written to results/rules.txt.
"""
import argparse, gzip, hashlib, json, os, re, sys, subprocess
from pathlib import Path
import numpy as np
import pandas as pd

RULES = """
RULES (family level; EQ members only; 'SP+' = members with SignalP class != OTHER)

Export route (route_call)
  Plurality SignalP 6.0 class among SP+ EQ members: SP -> Sec_SPI, LIPO -> Sec_SPII_lipoprotein,
  TAT -> Tat_SPI, TATLIPO -> Tat_SPII_lipoprotein, PILIN -> pilin_like_classIII.
  route_agreement = fraction of SP+ EQ members with the plurality class. Exact ties -> route_call = tie.
  External agreement: fraction of SP+ external members of the same family with the EQ plurality class
  (ext_route_agreement_with_eq) and whether the external plurality class equals the EQ one.

Mature region
  Starts after the SignalP cleavage site; if SignalP reported no site, after the TMbed signal-peptide label;
  if neither, after the median cleavage position of the SignalP class (mature_start_source records which).

Topology (TMbed; labels reused from the 2026-10-30 run for byte-identical sequences, TMbed 1.0.2 CPU run
  for the rest; evaluable = SP+ members with a TMbed prediction)
  tm_helices_mature = TMbed helices whose midpoint lies in the mature region.
  cterm_tm = exactly one mature helix and it ends <= 30 residues from the C-terminus
             (evaluated only on members whose 3' end is complete according to the gene caller).
  tm_helix_in_signal_region & no TMbed signal label = TMbed calls an N-terminal TM helix where SignalP calls
             a cleavable signal peptide (possible uncleaved signal anchor / SignalP false positive).

C-terminal sorting signals (HMMER 3.4 hmmsearch, full-length EQ proteins, model GA cutoff, else TC)
  PGF-CTERM (TIGR04126, Pfam PGF-CTERM), PEF-CTERM (TIGR03024, Pfam PEF-CTERM_ARCH) and the other short
  (<= 50 match states) *-CTERM / sorting-signal models listed in sorting_models_used.tsv; a hit counts only
  if its alignment ends <= 10 residues from the C-terminus.
  Archaeosortase family per genome: TIGR04125 (ArtA), TIGR04144 (ArtB), TIGR03762 (ArtC), TIGR04175 (ArtD),
  TIGR04124 (ArtE), TIGR04293 (ArtF), TIGR04297, TIGR04178 (exosortase/archaeosortase family) and Pfam
  Exosortase_EpsH, at the model GA cutoff, all proteins of each EQ genome.
  Relaxed variant: domain i-Evalue < 1e-3.

Anchor class (anchor_call; first rule that fires; 'majority' = fraction > 0.5 of evaluable members)
  1 multi_TM                 majority of evaluable members have >= 2 mature TM helices
  2 Cterm_sorting_signal     majority of C-complete SP+ members carry a C-terminal *-CTERM sorting-signal hit
  3 Cterm_TM                 majority of C-complete evaluable members have cterm_tm
  4 pilin_Nterm_hydrophobic  route_call = pilin_like_classIII
  5 single_TM_not_Cterm      majority of evaluable members have exactly 1 mature TM helix
  6 lipid_Nterm_lipobox      route_call is a lipoprotein route (SignalP LIPO/TATLIPO)
  7 none_predicted           majority of evaluable members have 0 mature TM helices
  8 unresolved               otherwise (including no evaluable member)

Role class (role_class; first rule that fires)
  1 pilin_like                          route_call = pilin_like_classIII
  2 multipass_membrane_with_SP          anchor_call = multi_TM (possible SignalP false positive)
  3 Cterm_anchored_surface              anchor_call in {Cterm_sorting_signal, Cterm_TM}
  4 single_pass_membrane_with_SP        anchor_call = single_TM_not_Cterm
  5 Nterm_signal_anchor_candidate       majority of evaluable members: TMbed gives no signal-peptide label AND
                                        a TM helix in the SignalP signal region (SignalP/TMbed conflict)
  6 lipoprotein_membrane_anchored       route lipoprotein (Sec or Tat), anchor_call = lipid_Nterm_lipobox
  7 Tat_substrate_soluble               route_call = Tat_SPI, anchor_call = none_predicted
  8 secreted_soluble                    route_call = Sec_SPI, anchor_call = none_predicted
  9 unresolved                          otherwise (route tie, no topology, mixed evidence)

Domain content (HMMER 3.4 hmmsearch on full-length EQ member proteins, Z = number of member proteins;
  only hits whose envelope midpoint lies in the mature region are kept)
  Pfam strict: sequence score >= GA1 and domain score >= GA2 (the Pfam gathering thresholds).
  Pfam relaxed: domain i-Evalue < 1e-3 (includes the strict hits).
  Overlapping hits on one protein are resolved greedily by domain score (envelope overlap > 50% of the
  shorter hit is dropped) before listing 'top' Pfam families; category flags use all hits passing the cutoff.
  dbCAN (HMMdb V14): domain i-Evalue < 1e-15 and HMM coverage > 0.35 (dbCAN defaults).
  MEROPS pepunit (DIAMOND blastp --very-sensitive): best hit with E <= 1e-10 covering >= 50% of the
  peptidase unit; MEROPS contains non-peptidase homologues, so a hit is homology to a peptidase family.
  Domain categories (see domain_category_pfam.tsv for the Pfam families assigned to each):
    peptidase          Pfam name/description matches a peptidase pattern, excluding inhibitors, propeptides,
                       non-catalytic accessory domains (e.g. PPC, PepX_C), signal peptidases, S33 hydrolases
    sbp                Pfam clans PBP, Periplas_BP, or SBP name patterns (ABC/TRAP/TTT substrate-binding)
    adhesin_fold       Pfam clans E-set, Ig, Beta_propeller, Pec_lyase or PKD/Big/cadherin/CARDB name patterns
                       (fold classes common in adhesins/S-layer proteins; not evidence of adhesion)
    pilin_archaellin   Class_IIIsignal, archaeal flagellin/archaellin, pilin families, TIGR02537
  Family flag *_maj: > 0.5 of EQ members carry the category (strict cutoffs); *_any: >= 1 member.
"""

ROUTE = {"SP": "Sec_SPI", "LIPO": "Sec_SPII_lipoprotein", "TAT": "Tat_SPI", "TATLIPO": "Tat_SPII_lipoprotein",
         "PILIN": "pilin_like_classIII"}
SPCL = ["SP", "LIPO", "TAT", "TATLIPO", "PILIN"]
ORDERS = ["o__JACPAO01", "o__JACQPN01", "o__RBG-16-68-12", "o__SW-10-69-26"]
DOMCOLS = ["target", "tacc", "tlen", "model", "macc", "mlen", "full_e", "full_score", "full_bias", "idx", "of",
           "c_e", "i_e", "dom_score", "dom_bias", "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to", "acc"]
PEPT_INC = re.compile(r"peptidase|protease|proteinase|proteolytic|subtilase|subtilisin|trypsin|collagenase|thermolysin|"
                      r"astacin|reprolysin|sedolisin|metallopeptid|aminopeptid|carboxypeptid|endopeptid|"
                      r"oligopeptidase|pepsin|papain|caspase|gingipain|amidohydro.*pept|kumamolisin|lysostaphin", re.I)
PEPT_EXC = re.compile(r"inhibitor|propeptide|pro-peptide|prodomain|pro-domain|PepSY|peptidase.associated|"
                      r"protease.associated|chaperone|ATPase|regulat|activator|substrate|N-terminal domain of|"
                      r"C-terminal domain of|^PA$|signal peptide peptidase.*non|ClpS|anti-|pre-peptidase|non-catalytic|inactive|"
                      r"activation domain|N-terminal beta-propeller|Vault protein|C-terminus Type|S33|transpeptidase|"
                      r"Signal peptidase|N-terminal domain|C-terminal domain$", re.I)
SBP_NAME = re.compile(r"^(SBP_bac|Peripla_BP|ABC_sub_bind|Lipoprotein_9|NMT1|DctP|TctC|Phosphonate-bd|OpuAC|ZnuA|Bmp|"
                      r"PBP_like|Sugar-bind|ANF_receptor|LysR_substrate|TAXI)", re.I)
ADH_NAME = re.compile(r"PKD|^Big_|Cadherin|CARDB|^fn3|^DUF11$|He_PIG|Invasin|Beta_helix|^NHL|^PQQ|^Kelch|^WD40|"
                      r"^TolB|YVTN|SdrD_B|Laminin_G|Calx-beta|^Ig_|^I-set|^TIG|^Y_Y_Y|^VCBS|^Cna_B|S-layer|S_layer", re.I)
PIL_NAME = re.compile(r"Class_IIIsignal|Flagellin_arch|Arch_flagellin|^Pilin|pilin|archaell", re.I)
CTERM_SIG = re.compile(r"PGF-CTERM|PGF_CTERM|PEF-CTERM|PEF_CTERM|PEFG-CTERM|PIP-CTERM|PIP_CTERM|VPxxxP|VPXXXP|"
                       r"PEP-CTERM|PEP_exosort|IPTL-CTERM|VPDSG|VPEID|PID_CTERM|MYXO-CTERM|GlyGly|VPLPA|LPXTG", re.I)
SORTASE = re.compile(r"archaeosortase|exosortase|^sortase|sortase$|Exosortase_EpsH|thauma_sortase|: sortase", re.I)


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def read_domtbl(p):
    d = pd.read_csv(p, sep=r"\s+", header=None, usecols=range(22), names=DOMCOLS, engine="c",
                    dtype={"target": str, "model": str, "macc": str, "tacc": str})
    return d


def resolve_overlaps(d):
    """Greedy per-protein overlap resolution by domain score (envelope overlap > 50% of the shorter)."""
    d = d.sort_values(["target", "dom_score"], ascending=[True, False])
    keep = np.zeros(len(d), bool)
    t = d.target.to_numpy(); s = d.env_from.to_numpy(); e = d.env_to.to_numpy()
    i = 0; n = len(d)
    while i < n:
        j = i; acc = []
        while j < n and t[j] == t[i]:
            ok = True
            for (a, b) in acc:
                ov = min(b, e[j]) - max(a, s[j]) + 1
                if ov > 0.5 * min(b - a + 1, e[j] - s[j] + 1):
                    ok = False; break
            if ok:
                acc.append((s[j], e[j])); keep[j] = True
            j += 1
        i = j
    return d[keep]


def top_str(counts, n, k=5):
    return ";".join(f"{a}:{c / n:.2f}" for a, c in counts.most_common(k)) if n else ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--clans", type=Path, required=True)
    ap.add_argument("--draws", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20261111)
    ap.add_argument("--receipt-extra", default="{}")
    a = ap.parse_args()
    R, O = a.root, a.out
    W = O / "work"; RES = O / "results"; RES.mkdir(exist_ok=True)
    A2 = R / "analysis2_families/results"
    checks = {}
    (RES / "rules.txt").write_text(RULES)

    # ------------------------------------------------------------------ inputs
    ft = pd.read_csv(A2 / "eq_family_table_id30_c80_all_genes.tsv.gz", sep="\t")
    wide981 = pd.read_csv(A2 / "dark_exported_widespread_families_id30_c80.tsv", sep="\t")
    cc = pd.read_csv(A2 / "class_counts.tsv", sep="\t")
    cc = cc[(cc.identity == "id30_c80") & (cc.gene_set == "all_genes")].set_index("class")
    mem = pd.read_csv(W / "members.tsv.gz", sep="\t", dtype={"gene_id": str, "family_rep": str})
    pf = pd.read_csv(W / "eq_protein_features.tsv.gz", sep="\t")
    gm = pd.read_csv(R / "metadata/harmonized_genome_manifest.tsv", sep="\t", dtype=str)
    ani = pd.read_csv(R / "source-metadata/ani_groups/ani_95_groups_1957.tsv", sep="\t", dtype=str)
    gm = gm.merge(ani[["genome_id", "ani95_group_id"]], on="genome_id", how="left")
    gm["block"] = np.where(gm.tax_genus.fillna("g__").isin(["g__", "", "missing", "none"]),
                           "ani::" + gm.ani95_group_id, gm.tax_genus)
    eqg = gm[gm.cohort == "EQ_Thermoplasmatota_237"].set_index("genome_id")
    assert len(eqg) == 237
    exp = ft[ft["class"].isin(["dark_exported", "nondark_exported"])].set_index("family_rep")

    # ------------------------------------------------------------------ Pfam
    pm = pd.read_csv(W / "pfam_models.tsv", sep="\t", header=None, names=["model", "macc", "mlen", "ga1", "ga2", "desc"])
    cl = pd.read_csv(a.clans, sep="\t", header=None, names=["pf", "clan", "clan_name", "name", "cdesc"])
    pm["pf"] = pm.macc.str.split(".").str[0]
    pm = pm.merge(cl[["pf", "clan", "clan_name"]], on="pf", how="left")
    pm["clan_name"] = pm.clan_name.fillna("")
    txt = pm.model + " " + pm.desc.fillna("")
    cat = pd.Series("", index=pm.index)
    is_pept = txt.str.contains(PEPT_INC) & ~txt.str.contains(PEPT_EXC)
    is_sbp = pm.clan_name.isin(["PBP", "Periplas_BP"]) | pm.model.str.contains(SBP_NAME)
    is_pil = txt.str.contains(PIL_NAME)
    is_adh = (pm.clan_name.isin(["E-set", "Ig", "Beta_propeller", "Pec_lyase"]) | pm.model.str.contains(ADH_NAME)) & ~is_pil
    pm["cat_peptidase"] = is_pept; pm["cat_sbp"] = is_sbp & ~is_pept
    pm["cat_adhesin_fold"] = is_adh & ~is_pept & ~pm.cat_sbp; pm["cat_pilin_archaellin"] = is_pil
    CATS = ["peptidase", "sbp", "adhesin_fold", "pilin_archaellin"]
    checks["pfam_models"] = int(len(pm)); checks["pfam_models_with_clan"] = int((pm.clan_name != "").sum())
    checks["pfam_category_model_counts"] = {c: int(pm["cat_" + c].sum()) for c in CATS}

    mstart = pf.set_index("gene_id").mature_start

    def mature_only(x, lo="env_from", hi="env_to", key="target"):
        """keep hits whose midpoint lies in the mature region (searches were run on full-length proteins)"""
        return x[((x[lo] + x[hi]) / 2 - 1) >= x[key].map(mstart)]
    d = read_domtbl(W / "pfam_vs_members.domtbl.gz")
    d = d.merge(pm[["model", "ga1", "ga2"]], on="model", how="left")
    assert d.ga1.notna().all()
    n0 = int(((d.full_score >= d.ga1) & (d.dom_score >= d.ga2)).sum())
    d = mature_only(d)
    checks["pfam_strict_hits_dropped_in_signal_region"] = n0 - int(((d.full_score >= d.ga1) & (d.dom_score >= d.ga2)).sum())
    d["strict"] = (d.full_score >= d.ga1) & (d.dom_score >= d.ga2)
    d["relaxed"] = d.strict | (d.i_e < 1e-3)
    checks["pfam_domtbl_rows"] = int(len(d)); checks["pfam_strict_rows"] = int(d.strict.sum())
    checks["pfam_relaxed_rows"] = int(d.relaxed.sum())
    d = d[d.relaxed].copy()
    pmi = pm.set_index("model")
    for c in CATS:
        d["cat_" + c] = d.model.map(pmi["cat_" + c])
    ds = resolve_overlaps(d[d.strict]); dr = resolve_overlaps(d)
    prot = pf.set_index("gene_id")
    prot["pfam_strict"] = ds.sort_values(["target", "env_from"]).groupby("target").model.agg(";".join).reindex(prot.index).fillna("")
    prot["pfam_relaxed"] = dr.sort_values(["target", "env_from"]).groupby("target").model.agg(";".join).reindex(prot.index).fillna("")
    for c in CATS:
        for lab, dd in (("strict", d[d.strict]), ("relaxed", d)):
            hit = dd[dd["cat_" + c]].groupby("target").model.agg(lambda x: ";".join(sorted(set(x))))
            prot[f"{c}_{lab}"] = hit.reindex(prot.index).fillna("")
    used = d.groupby("model").agg(n_proteins_relaxed=("target", "nunique"),
                                  n_proteins_strict=("target", lambda x: x[d.loc[x.index, "strict"]].nunique()))
    catt = pm.set_index("model").join(used, how="inner")
    catt = catt[catt[["cat_" + c for c in CATS]].any(axis=1)]
    catt.reset_index()[["model", "macc", "clan_name", "desc"] + ["cat_" + c for c in CATS] +
                       ["n_proteins_strict", "n_proteins_relaxed"]].to_csv(RES / "domain_category_pfam.tsv", sep="\t", index=False)

    # ------------------------------------------------------------------ dbCAN / MEROPS
    dc = mature_only(read_domtbl(W / "dbcan_vs_members.domtbl.gz"))
    dc = dc[(dc.i_e < 1e-15) & ((dc.hmm_to - dc.hmm_from + 1) / dc.mlen > 0.35)].copy()
    dc["fam"] = dc.model.str.replace(".hmm", "", regex=False)
    prot["dbcan"] = dc.groupby("target").fam.agg(lambda x: ";".join(sorted(set(x)))).reindex(prot.index).fillna("")
    mcols = ["q", "s", "pident", "length", "qlen", "slen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle"]
    try:
        me = pd.read_csv(W / "merops_pepunit.m8.gz", sep="\t", header=None, names=mcols)
    except pd.errors.EmptyDataError:
        me = pd.DataFrame(columns=mcols)
    checks["merops_raw_rows"] = int(len(me))
    if len(me):
        me = mature_only(me, "qstart", "qend", "q")
    me = me.sort_values("bitscore", ascending=False).drop_duplicates("q")
    me["ucov"] = (me.send - me.sstart + 1) / me.slen
    me["mfam"] = me.stitle.str.extract(r"#([A-Z]\d+[A-Z]?)#")[0]
    me["mcode"] = me.stitle.str.extract(r"\[([A-Z]\d+\.[A-Za-z0-9]+)\]")[0]
    mes = me[(me.evalue <= 1e-10) & (me.ucov >= 0.5)].set_index("q")
    prot["merops_family"] = mes.mfam.reindex(prot.index).fillna("")
    prot["merops_code"] = mes.mcode.reindex(prot.index).fillna("")
    prot["merops_best_evalue"] = me.set_index("q").evalue.reindex(prot.index)
    prot["merops_best_title"] = me.set_index("q").stitle.reindex(prot.index).fillna("").str.slice(0, 120)

    # ------------------------------------------------------------------ sorting signals / sortases (all EQ proteins)
    sm = pd.read_csv(W / "sorting_models.tsv", sep="\t", header=None,
                     names=["model", "macc", "mlen", "ga1", "ga2", "tc1", "tc2", "desc"])
    sm["c1"] = sm.ga1.fillna(sm.tc1); sm["c2"] = sm.ga2.fillna(sm.tc2)
    sm["label"] = np.where(sm.desc.fillna("-").str.contains(":"), sm.desc.str.split(":").str[0], sm.model)
    stxt = sm.model + " " + sm.desc.fillna("")
    sm["is_cterm_signal"] = (stxt.str.contains(r"CTERM|sorting signal|sorting domain|LPXTG|srtB_sig|srtB_target|MprA_tail|Por_Secre_tail", case=False)
                             & (sm.mlen <= 50) & ~stxt.str.contains(r"repeat|H-X9-DG|prepilin|pilin|flagellin|Class III", case=False))
    sm["is_sortase"] = stxt.str.contains(SORTASE) & ~stxt.str.contains(r"associated|system-|sorting|target|signal", case=False)
    sm["is_pgf_pre_pgf"] = stxt.str.contains(r"PGF_pre_PGF|PGF-pre-PGF", case=False)
    sm["is_pilin"] = stxt.str.contains(r"arch_flag_Nterm|Class_IIIsignal|Flagellin_arch|class III signal|Pilin_N", case=False)
    sm.to_csv(RES / "sorting_models_used.tsv", sep="\t", index=False)
    sd = read_domtbl(W / "sorting_models_vs_eq_all.domtbl.gz")
    sd = sd.merge(sm[["model", "label", "c1", "c2", "is_cterm_signal", "is_sortase", "is_pgf_pre_pgf", "is_pilin"]], on="model")
    sd["strict"] = (sd.full_score >= sd.c1) & (sd.dom_score >= sd.c2)
    sd["relaxed"] = sd.strict | (sd.i_e < 1e-3)
    sd["at_cterm"] = (sd.tlen - sd.ali_to) <= 10
    sd["genome_id"] = sd.target.str.split("::").str[0]
    sd = sd[sd.relaxed]
    cs_s = sd[sd.is_cterm_signal & sd.at_cterm & sd.strict].groupby("target").label.agg(lambda x: ";".join(sorted(set(x))))
    cs_r = sd[sd.is_cterm_signal & sd.at_cterm].groupby("target").label.agg(lambda x: ";".join(sorted(set(x))))
    prot["cterm_sorting_signal_strict"] = cs_s.reindex(prot.index).fillna("")
    prot["cterm_sorting_signal_relaxed"] = cs_r.reindex(prot.index).fillna("")
    prot["pgf_pre_pgf_strict"] = prot.index.isin(set(sd[sd.is_pgf_pre_pgf & sd.strict].target))
    prot["pilin_tigr_pfam_strict"] = prot.index.isin(set(sd[sd.is_pilin & sd.strict].target))
    # per genome
    g = pd.DataFrame(index=eqg.index)
    g["tax_order"] = eqg.tax_order; g["completeness"] = eqg.completeness_primary
    for lab, sub in sd[sd.is_sortase & sd.strict].groupby("label"):
        g["n_sortase_" + lab] = sub.groupby("genome_id").target.nunique().reindex(g.index).fillna(0).astype(int)
    sc = [c for c in g.columns if c.startswith("n_sortase_")]
    arch = [c for c in sc if re.search(r"archaeo|exosort_PGF_TRM|exo_archaeo|art[A-G]|thauma|Exosortase_EpsH", c, re.I)]
    g["n_archaeosortase_family_proteins_strict"] = sd[sd.is_sortase & sd.strict & sd.label.isin([c[10:] for c in arch])].groupby("genome_id").target.nunique().reindex(g.index).fillna(0).astype(int)
    g["n_any_sortase_family_proteins_relaxed"] = sd[sd.is_sortase].groupby("genome_id").target.nunique().reindex(g.index).fillna(0).astype(int)
    pg = sd[sd.is_cterm_signal & sd.at_cterm & sd.strict]
    for lab, sub in pg.groupby("label"):
        g["n_proteins_" + lab] = sub.groupby("genome_id").target.nunique().reindex(g.index).fillna(0).astype(int)
    g["n_proteins_any_cterm_sorting_signal_strict"] = pg.groupby("genome_id").target.nunique().reindex(g.index).fillna(0).astype(int)
    g["n_proteins_PGF_pre_PGF_strict"] = sd[sd.is_pgf_pre_pgf & sd.strict].groupby("genome_id").target.nunique().reindex(g.index).fillna(0).astype(int)
    g.to_csv(RES / "archaeosortase_and_sorting_signals_per_genome.tsv", sep="\t")
    summ_sort = {"n_genomes": 237,
                 "genomes_with_archaeosortase_family_strict": int((g.n_archaeosortase_family_proteins_strict > 0).sum()),
                 "genomes_with_any_sortase_family_relaxed": int((g.n_any_sortase_family_proteins_relaxed > 0).sum()),
                 "genomes_with_any_cterm_sorting_signal_protein_strict": int((g.n_proteins_any_cterm_sorting_signal_strict > 0).sum()),
                 "per_model_genomes": {c: int((g[c] > 0).sum()) for c in g.columns if c.startswith("n_sortase_") or c.startswith("n_proteins_")},
                 "per_model_proteins": {c: int(g[c].sum()) for c in g.columns if c.startswith("n_sortase_") or c.startswith("n_proteins_")},
                 "by_order_genomes_with_archaeosortase": g.groupby("tax_order").n_archaeosortase_family_proteins_strict.agg(lambda x: f"{(x > 0).sum()}/{len(x)}").to_dict()}
    checks["sortase_summary"] = summ_sort

    # ------------------------------------------------------------------ per-protein derived flags
    prot["sp_pos"] = prot.signalp_class != "OTHER"
    prot["topo_eval"] = prot.sp_pos & prot.tmbed_available
    prot["tm0"] = prot.topo_eval & (prot.tm_helices_mature == 0)
    prot["tm1"] = prot.topo_eval & (prot.tm_helices_mature == 1)
    prot["tm_ge2"] = prot.topo_eval & (prot.tm_helices_mature >= 2)
    prot["c_eval"] = prot.sp_pos & prot.c_complete.astype(bool)
    prot["c_topo_eval"] = prot.c_eval & prot.tmbed_available
    prot["cterm_tm1"] = prot.c_topo_eval & (prot.tm_helices_mature == 1) & prot.cterm_tm.astype(bool)
    prot["sorting_strict"] = prot.c_eval & (prot.cterm_sorting_signal_strict != "")
    prot["sorting_relaxed"] = prot.c_eval & (prot.cterm_sorting_signal_relaxed != "")
    prot["motif_heur"] = prot.c_eval & prot.cterm_PGF_PEF_motif_heuristic.astype(bool)
    prot["tmbed_signal"] = prot.topo_eval & (prot.tmbed_signal_len > 0)
    prot["nterm_anchor_conflict"] = prot.topo_eval & (prot.tmbed_signal_len == 0) & (prot.tm_helix_in_signal_region > 0)
    prot["st_rich"] = prot.st_max_window30 >= 0.4
    prot["any_pfam_strict"] = prot.pfam_strict != ""; prot["any_pfam_relaxed"] = prot.pfam_relaxed != ""
    prot["has_merops"] = prot.merops_family != ""; prot["has_dbcan"] = prot.dbcan != ""
    for c in CATS:
        prot[f"has_{c}_strict"] = prot[f"{c}_strict"] != ""
        prot[f"has_{c}_relaxed"] = prot[f"{c}_relaxed"] != ""
    prot["has_pilin_archaellin_strict"] |= prot.pilin_tigr_pfam_strict
    prot["has_pilin_archaellin_relaxed"] |= prot.pilin_tigr_pfam_strict
    prot["has_peptidase_evidence_any"] = prot.has_peptidase_strict | prot.has_merops
    prot.to_csv(RES / "eq_protein_evidence.tsv.gz", sep="\t")

    # ------------------------------------------------------------------ family aggregation
    gb = prot.groupby("family_rep", sort=False)
    F = pd.DataFrame(index=exp.index)
    F["class"] = exp["class"]
    for c in ["n_eq_proteins", "n_eq_dark", "n_eq_signal_peptide", "n_eq_genomes", "n_eq_blocks", "median_eq_aa_length"] + \
             ["n_eq_" + o for o in ORDERS] + ["n_eq_unassigned_order", "cluster_member_in_any_external", "search_hit_in_any_external"]:
        F[c] = exp[c]
    F["n_eq_members_recount"] = gb.size().reindex(F.index)
    F["n_eq_sp_pos_recount"] = gb.sp_pos.sum().reindex(F.index)
    cls_counts = prot[prot.sp_pos].groupby(["family_rep", "signalp_class"]).size().unstack(fill_value=0).reindex(F.index).fillna(0).astype(int)
    for c in SPCL:
        F["n_eq_" + c] = cls_counts[c] if c in cls_counts else 0
    cm = F[["n_eq_" + c for c in SPCL]].to_numpy()
    top = cm.max(1); ntop = (cm == top[:, None]).sum(1)
    major = np.array(SPCL)[cm.argmax(1)]
    F["route_major_signalp_class"] = np.where((top > 0) & (ntop == 1), major, "tie")
    F["route_call"] = F.route_major_signalp_class.map(ROUTE).fillna("tie")
    F["route_agreement"] = np.where(cm.sum(1) > 0, top / np.maximum(cm.sum(1), 1), np.nan).round(4)
    # external members
    ext = mem[~mem.is_eq]
    F["n_ext_members"] = ext.groupby("family_rep").size().reindex(F.index).fillna(0).astype(int)
    exs = ext[ext.signalp_class != "OTHER"]
    ec = exs.groupby(["family_rep", "signalp_class"]).size().unstack(fill_value=0).reindex(F.index).fillna(0).astype(int)
    for c in SPCL:
        if c not in ec:
            ec[c] = 0
    ec = ec[SPCL]
    F["n_ext_sp_pos"] = ec.sum(1)
    F["ext_frac_sp_pos"] = (F.n_ext_sp_pos / F.n_ext_members.where(F.n_ext_members > 0)).round(4)
    em = ec.to_numpy(); etop = em.max(1); entop = (em == etop[:, None]).sum(1)
    F["ext_route_major_signalp_class"] = np.where(etop == 0, "none", np.where(entop == 1, np.array(SPCL)[em.argmax(1)], "tie"))
    idx = pd.Series(F.route_major_signalp_class).map({c: i for i, c in enumerate(SPCL)})
    same = np.array([em[i, int(j)] if j == j else np.nan for i, j in enumerate(idx.to_numpy())], float)
    F["ext_route_agreement_with_eq"] = np.where(F.n_ext_sp_pos > 0, same / np.maximum(F.n_ext_sp_pos, 1), np.nan).round(4)
    F["ext_major_equals_eq"] = np.where((F.n_ext_sp_pos > 0) & (F.route_major_signalp_class != "tie"),
                                        (F.ext_route_major_signalp_class == F.route_major_signalp_class).astype(float), np.nan)

    def frac(num, den):
        n = gb[num].sum().reindex(F.index); dn = gb[den].sum().reindex(F.index)
        return (n / dn.where(dn > 0)).round(4), dn.astype(int)
    F["frac_tm0"], F["n_topology_evaluable"] = frac("tm0", "topo_eval")
    F["frac_tm1"], _ = frac("tm1", "topo_eval")
    F["frac_tm_ge2"], _ = frac("tm_ge2", "topo_eval")
    F["frac_tmbed_signal_label"], _ = frac("tmbed_signal", "topo_eval")
    F["frac_nterm_anchor_conflict"], _ = frac("nterm_anchor_conflict", "topo_eval")
    F["frac_cterm_tm"], F["n_cterm_topology_evaluable"] = frac("cterm_tm1", "c_topo_eval")
    F["frac_cterm_sorting_signal_strict"], F["n_cterm_evaluable"] = frac("sorting_strict", "c_eval")
    F["frac_cterm_sorting_signal_relaxed"], _ = frac("sorting_relaxed", "c_eval")
    F["frac_cterm_PGF_PEF_motif_heuristic"], _ = frac("motif_heur", "c_eval")
    F["cterm_sorting_signal_models"] = prot[prot.sorting_strict].groupby("family_rep").cterm_sorting_signal_strict.agg(
        lambda x: ";".join(sorted(set(";".join(x).split(";"))))).reindex(F.index).fillna("")
    F["frac_pgf_pre_pgf"] = gb.pgf_pre_pgf_strict.mean().reindex(F.index).round(4)
    spp = prot[prot.sp_pos].groupby("family_rep")
    F["median_tm_helices_mature"] = prot[prot.topo_eval].groupby("family_rep").tm_helices_mature.median().reindex(F.index)
    F["median_mature_length"] = spp.mature_length.median().reindex(F.index)
    F["median_mature_pI"] = spp.mature_pI.median().reindex(F.index).round(2)
    F["median_sequons_per_100aa"] = spp.sequons_per_100aa.median().reindex(F.index).round(3)
    F["median_st_fraction"] = spp.st_fraction.median().reindex(F.index).round(4)
    F["frac_st_rich_window"] = spp.st_rich.mean().reindex(F.index).round(4)
    F["median_lowcomplexity_fraction"] = spp.lowcomplexity_fraction.median().reindex(F.index).round(4)
    F["frac_mature_start_from_signalp_cs"] = spp.mature_start_source.agg(lambda x: (x == "signalp_cs").mean()).reindex(F.index).round(4)
    for c in ["any_pfam_strict", "any_pfam_relaxed", "has_merops", "has_dbcan", "has_peptidase_evidence_any"] + \
             [f"has_{k}_{l}" for k in CATS for l in ("strict", "relaxed")]:
        F["frac_" + c.replace("has_", "")] = gb[c].mean().reindex(F.index).round(4)
    from collections import Counter

    def top_of(col, k=5):
        out = {}
        for fam, x in prot.loc[prot[col] != "", ["family_rep", col]].groupby("family_rep")[col]:
            cn = Counter()
            for s in x:
                cn.update(set(s.split(";")))
            out[fam] = cn
        return out
    nmem = F.n_eq_members_recount
    tops = {}
    for col in ["pfam_strict", "pfam_relaxed", "merops_family", "dbcan", "peptidase_strict", "peptidase_relaxed"]:
        tops[col] = top_of(col)
        F["top_" + col] = [top_str(tops[col][f], nmem[f]) if f in tops[col] else "" for f in F.index]
    # family Pfam set = strict Pfam families present in > 50% of EQ members
    F["pfam_majority_set"] = [";".join(sorted(p for p, c in tops["pfam_strict"][f].items() if c > 0.5 * nmem[f]))
                              if f in tops["pfam_strict"] else "" for f in F.index]
    for k in CATS:
        F[f"{k}_maj"] = F[f"frac_{k}_strict"] > 0.5
        F[f"{k}_any"] = F[f"frac_{k}_strict"] > 0
        F[f"{k}_relaxed_maj"] = F[f"frac_{k}_relaxed"] > 0.5
        F[f"{k}_relaxed_any"] = F[f"frac_{k}_relaxed"] > 0
    F["merops_maj"] = F.frac_merops > 0.5; F["merops_any"] = F.frac_merops > 0
    F["dbcan_maj"] = F.frac_dbcan > 0.5; F["dbcan_any"] = F.frac_dbcan > 0
    F["peptidase_evidence_maj"] = F.frac_peptidase_evidence_any > 0.5
    F["peptidase_evidence_any"] = F.frac_peptidase_evidence_any > 0
    F["any_pfam_strict_maj"] = F.frac_any_pfam_strict > 0.5
    F["any_pfam_relaxed_maj"] = F.frac_any_pfam_relaxed > 0.5

    # ------------------------------------------------------------------ anchor + role rules
    lipo = F.route_call.isin(["Sec_SPII_lipoprotein", "Tat_SPII_lipoprotein"])
    conds = [F.frac_tm_ge2 > 0.5, F.frac_cterm_sorting_signal_strict > 0.5, F.frac_cterm_tm > 0.5,
             F.route_call == "pilin_like_classIII", F.frac_tm1 > 0.5, lipo, F.frac_tm0 > 0.5]
    labs = ["multi_TM", "Cterm_sorting_signal", "Cterm_TM", "pilin_Nterm_hydrophobic", "single_TM_not_Cterm",
            "lipid_Nterm_lipobox", "none_predicted"]
    F["anchor_call"] = np.select(conds, labs, default="unresolved")
    rc = [F.route_call == "pilin_like_classIII", F.anchor_call == "multi_TM",
          F.anchor_call.isin(["Cterm_sorting_signal", "Cterm_TM"]), F.anchor_call == "single_TM_not_Cterm",
          F.frac_nterm_anchor_conflict > 0.5, lipo & (F.anchor_call == "lipid_Nterm_lipobox"),
          (F.route_call == "Tat_SPI") & (F.anchor_call == "none_predicted"),
          (F.route_call == "Sec_SPI") & (F.anchor_call == "none_predicted")]
    rl = ["pilin_like", "multipass_membrane_with_SP", "Cterm_anchored_surface", "single_pass_membrane_with_SP",
          "Nterm_signal_anchor_candidate", "lipoprotein_membrane_anchored", "Tat_substrate_soluble", "secreted_soluble"]
    F["role_class"] = np.select(rc, rl, default="unresolved")
    F["role_rule"] = np.select(rc, list(range(1, 9)), default=9)
    F["ge2_genomes"] = F.n_eq_genomes >= 2
    F["widespread"] = (F.n_eq_genomes >= 10) & (F.n_eq_blocks >= 3)
    F["no_external_homologue"] = (F.cluster_member_in_any_external == 0) & (F.search_hit_in_any_external == 0)
    w = F.loc[wide981.family_rep]

    # ------------------------------------------------------------------ checks
    dk = F["class"] == "dark_exported"; nd = ~dk
    checks["n_dark_exported"] = [int(dk.sum()), int(cc.n_families["dark_exported"])]
    checks["n_dark_exported_ge2"] = [int((dk & F.ge2_genomes).sum()), int(cc.n_families_ge2_eq_genomes["dark_exported"])]
    checks["n_dark_exported_widespread"] = [int((dk & F.widespread).sum()), int(cc.n_families_ge10genomes_ge3blocks["dark_exported"])]
    checks["n_nondark_exported"] = [int(nd.sum()), int(cc.n_families["nondark_exported"])]
    checks["n_nondark_exported_ge2"] = [int((nd & F.ge2_genomes).sum()), int(cc.n_families_ge2_eq_genomes["nondark_exported"])]
    checks["n_nondark_exported_widespread"] = [int((nd & F.widespread).sum()), int(cc.n_families_ge10genomes_ge3blocks["nondark_exported"])]
    checks["widespread_file_rows"] = int(len(wide981))
    checks["widespread_file_set_equals_recomputed"] = bool(set(wide981.family_rep) == set(F.index[dk & F.widespread]))
    checks["all_981_have_role"] = bool(len(w) == 981 and w.role_class.notna().all() and w.route_call.notna().all())
    checks["member_recount_matches_family_table"] = bool((F.n_eq_members_recount == F.n_eq_proteins).all())
    checks["sp_recount_matches_family_table"] = bool((F.n_eq_sp_pos_recount == F.n_eq_signal_peptide).all())
    checks["all_exported_majority_sp"] = bool((2 * F.n_eq_sp_pos_recount > F.n_eq_members_recount).all())
    checks["eq_protein_rows"] = [int(len(prot)), int(F.n_eq_proteins.sum())]
    checks["tmbed_coverage_sp_pos"] = [int(prot.topo_eval.sum()), int(prot.sp_pos.sum())]
    checks["widespread_no_external"] = {"cluster0_and_search0": int(w.no_external_homologue.sum()),
                                        "search0": int((w.search_hit_in_any_external == 0).sum()),
                                        "cluster0": int((w.cluster_member_in_any_external == 0).sum())}
    checks["mature_start_source_sp_pos"] = prot[prot.sp_pos].mature_start_source.value_counts().to_dict()

    # ------------------------------------------------------------------ summaries
    sets = {"dark_widespread": dk & F.widespread, "dark_ge2_genomes": dk & F.ge2_genomes, "dark_all": dk,
            "nondark_widespread": nd & F.widespread, "nondark_ge2_genomes": nd & F.ge2_genomes, "nondark_all": nd,
            "dark_widespread_no_external_homologue": dk & F.widespread & F.no_external_homologue,
            "dark_widespread_with_external_homologue": dk & F.widespread & ~F.no_external_homologue}
    for base in ["dark_widespread", "dark_ge2_genomes", "nondark_widespread", "nondark_ge2_genomes"]:
        for o in ORDERS:
            sets[f"{base}|present_in_{o}"] = sets[base] & (F["n_eq_" + o] > 0)
    flagvars = [f"{k}_{s}" for k in CATS for s in ("maj", "any", "relaxed_maj", "relaxed_any")] + \
               ["merops_maj", "merops_any", "dbcan_maj", "dbcan_any", "peptidase_evidence_maj", "peptidase_evidence_any",
                "any_pfam_strict_maj", "any_pfam_relaxed_maj"]
    F["Cterm_sorting_signal_any_relaxed"] = F.frac_cterm_sorting_signal_relaxed > 0
    F["pgf_pre_pgf_any"] = F.frac_pgf_pre_pgf > 0
    F["st_rich_majority"] = F.frac_st_rich_window > 0.5
    F["tmbed_signal_label_majority"] = F.frac_tmbed_signal_label > 0.5
    F["ext_route_major_equals_eq"] = F.ext_major_equals_eq.map({1.0: "same", 0.0: "different"}).fillna("no_external_SP_member_or_tie")
    flagvars += ["Cterm_sorting_signal_any_relaxed", "pgf_pre_pgf_any", "st_rich_majority", "tmbed_signal_label_majority"]
    rows = []
    for sn, m in sets.items():
        sub = F[m]; n = len(sub)
        for var in ["route_call", "anchor_call", "role_class", "ext_route_major_equals_eq"]:
            for lev, c in sub[var].value_counts().items():
                rows.append({"set": sn, "variable": var, "level": lev, "n": int(c), "n_total": n, "prop": round(c / n, 4) if n else np.nan})
        for var in flagvars:
            c = int(sub[var].sum())
            rows.append({"set": sn, "variable": var, "level": "True", "n": c, "n_total": n, "prop": round(c / n, 4) if n else np.nan})
        for var in ["route_agreement", "ext_route_agreement_with_eq", "median_mature_length", "median_mature_pI",
                    "median_sequons_per_100aa", "median_st_fraction", "median_lowcomplexity_fraction"]:
            x = sub[var].dropna()
            rows.append({"set": sn, "variable": var, "level": "median_over_families", "n": int(len(x)), "n_total": n,
                         "prop": round(float(x.median()), 4) if len(x) else np.nan})
    S = pd.DataFrame(rows)
    S.to_csv(RES / "summary_class_distributions.tsv", sep="\t", index=False)
    # route x role and role x domain-category cross tabs
    xr = []
    for sn in ["dark_widespread", "dark_ge2_genomes", "nondark_widespread", "nondark_ge2_genomes"]:
        sub = F[sets[sn]]
        for (r1, r2), c in sub.groupby(["route_call", "role_class"]).size().items():
            xr.append({"set": sn, "table": "route_x_role", "row": r1, "col": r2, "n": int(c)})
        for (r1, r2), c in sub.groupby(["route_call", "anchor_call"]).size().items():
            xr.append({"set": sn, "table": "route_x_anchor", "row": r1, "col": r2, "n": int(c)})
        for var in ["peptidase_maj", "peptidase_any", "peptidase_relaxed_any", "peptidase_evidence_any", "sbp_maj", "sbp_any",
                    "adhesin_fold_maj", "adhesin_fold_any", "pilin_archaellin_maj", "pilin_archaellin_any", "dbcan_any", "merops_any"]:
            for by in ["role_class", "route_call"]:
                for r1, c in sub[sub[var]].groupby(by).size().items():
                    xr.append({"set": sn, "table": f"{var}_by_{by}", "row": r1, "col": var, "n": int(c)})
    pd.DataFrame(xr).to_csv(RES / "summary_crosstabs.tsv", sep="\t", index=False)
    # Pfam families per route / role (family counted when the Pfam family is in its majority set)
    desc = pm.set_index("model").desc.to_dict(); clan = pm.set_index("model").clan_name.to_dict()
    pr = []
    for sn in ["nondark_widespread", "nondark_ge2_genomes", "dark_widespread", "dark_ge2_genomes"]:
        sub = F[sets[sn]]
        for by in ["route_call", "role_class"]:
            for lev, s2 in sub.groupby(by):
                cn = Counter()
                for x in s2.pfam_majority_set:
                    if x:
                        cn.update(x.split(";"))
                nn = len(s2); nnone = int((s2.pfam_majority_set == "").sum())
                pr.append({"set": sn, "by": by, "level": lev, "rank": 0, "pfam": "(no Pfam family in >50% of members)", "clan": "",
                           "description": "", "n_families": nnone, "n_families_in_level": nn, "prop": round(nnone / nn, 4)})
                for i, (p, c) in enumerate(cn.most_common(25), 1):
                    pr.append({"set": sn, "by": by, "level": lev, "rank": i, "pfam": p, "clan": clan.get(p, ""),
                               "description": desc.get(p, ""), "n_families": c, "n_families_in_level": nn, "prop": round(c / nn, 4)})
    pd.DataFrame(pr).to_csv(RES / "pfam_by_route_and_role.tsv", sep="\t", index=False)
    # peptidase / SBP / pilin / adhesin detail for the dark sets (auditable list)
    det = F[(dk & F.ge2_genomes) & (F.peptidase_relaxed_any | F.merops_any | F.sbp_relaxed_any | F.pilin_archaellin_relaxed_any |
                                    F.adhesin_fold_relaxed_any | F.dbcan_any | (F.route_call == "pilin_like_classIII"))]
    det[["class", "widespread", "no_external_homologue", "n_eq_proteins", "n_eq_genomes", "n_eq_blocks", "route_call", "route_agreement",
         "anchor_call", "role_class", "median_mature_length", "frac_peptidase_strict", "frac_peptidase_relaxed", "top_peptidase_strict",
         "top_peptidase_relaxed", "frac_merops", "top_merops_family", "frac_sbp_strict", "frac_sbp_relaxed", "frac_adhesin_fold_strict",
         "frac_adhesin_fold_relaxed", "frac_pilin_archaellin_strict", "frac_pilin_archaellin_relaxed", "frac_dbcan", "top_dbcan",
         "top_pfam_strict", "top_pfam_relaxed"]].to_csv(RES / "dark_exported_ge2_domain_category_families.tsv", sep="\t")

    F.to_csv(RES / "family_localisation_exported_families.tsv.gz", sep="\t")
    F.loc[wide981.family_rep].to_csv(RES / "family_localisation_dark_exported_widespread_981.tsv", sep="\t")
    # ------------------------------------------------------------------ block bootstrap (family set and calls fixed)
    # Weight of family f in a draw = sum_b m_b * (share of the family's EQ members that lie in block b); m = multinomial
    # resampling of EQ blocks. With m = 1 every weight is 1, so the point estimate is the plain proportion.
    prot["block"] = prot.genome_id.map(eqg.block)
    assert prot.block.notna().all()
    bcodes, bnames = pd.factorize(prot.block)
    fcodes = pd.Series(np.arange(len(F)), index=F.index).reindex(prot.family_rep).to_numpy()
    M = np.zeros((len(F), len(bnames)))
    np.add.at(M, (fcodes, bcodes), 1.0)
    M /= M.sum(1, keepdims=True)
    rng = np.random.default_rng(a.seed)
    k = len(bnames)
    draws = rng.multinomial(k, np.ones(k) / k, size=a.draws).T.astype(float)     # blocks x draws
    Wt = M @ draws                                                               # families x draws
    brow = []
    for pair, (sa, sb) in {"widespread": ("dark_widespread", "nondark_widespread"),
                           "ge2_genomes": ("dark_ge2_genomes", "nondark_ge2_genomes")}.items():
        ma, mb = sets[sa].to_numpy(), sets[sb].to_numpy()
        for var in ["route_call", "anchor_call", "role_class"]:
            for lev in sorted(set(F[var][ma | mb])):
                x = (F[var] == lev).to_numpy().astype(float)
                pa, pb = x[ma].mean(), x[mb].mean()
                ba = (Wt[ma] * x[ma, None]).sum(0) / Wt[ma].sum(0)
                bb = (Wt[mb] * x[mb, None]).sum(0) / Wt[mb].sum(0)
                dd = ba - bb
                lo, hi = np.quantile(dd, [0.025, 0.975])
                brow.append({"comparison": pair, "variable": var, "level": lev, "n_dark": int(x[ma].sum()), "N_dark": int(ma.sum()),
                             "prop_dark": round(pa, 4), "n_nondark": int(x[mb].sum()), "N_nondark": int(mb.sum()),
                             "prop_nondark": round(pb, 4), "diff_dark_minus_nondark": round(pa - pb, 4),
                             "boot_lo95": round(lo, 4), "boot_hi95": round(hi, 4), "boot_median": round(float(np.median(dd)), 4),
                             "interval_contains_point": bool(lo <= pa - pb <= hi)})
    B = pd.DataFrame(brow)
    B.to_csv(RES / "dark_vs_nondark_block_bootstrap.tsv", sep="\t", index=False)
    checks["bootstrap_n_blocks"] = int(k); checks["bootstrap_draws"] = a.draws
    checks["bootstrap_all_intervals_contain_point"] = bool(B.interval_contains_point.all())
    checks["bootstrap_n_intervals_not_containing_point"] = int((~B.interval_contains_point).sum())

    # ------------------------------------------------------------------ spot check against raw SignalP rows
    rs = np.random.default_rng(a.seed + 1)
    pick = list(rs.choice(wide981.family_rep.to_numpy(), 10, replace=False))
    famfile = pd.read_csv(R / "analysis2_families/clusters/families_id30_c80.tsv", sep="\t", header=None, names=["rep", "gene_id"], dtype=str)
    pm_ = famfile[famfile.rep.isin(pick)]
    want = set(pm_.gene_id)
    raw = []
    with gzip.open(R / "signalp/harmonized_signalp6_predictions.tsv.gz", "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            gidx = line.split("\t", 1)[0]
            if gidx in want:
                raw.append(line.rstrip("\n").split("\t"))
    raw = pd.DataFrame(raw, columns=hdr).merge(pm_, on="gene_id")
    raw["is_eq"] = raw.genome_id.isin(set(eqg.index))
    raw.to_csv(RES / "spotcheck_raw_signalp_rows.tsv.gz", sep="\t", index=False)
    sc = []
    for f in pick:
        r = raw[(raw.rep == f) & raw.is_eq]
        vc = r[r.signalp_class != "OTHER"].signalp_class.value_counts()
        rawmaj = vc.index[0] if len(vc) and (len(vc) == 1 or vc.iloc[0] > vc.iloc[1]) else "tie"
        sc.append({"family_rep": f, "raw_n_eq": len(r), "raw_n_eq_sp_any_true": int((r.sp_any == "True").sum()),
                   "raw_class_counts": ";".join(f"{a_}:{b_}" for a_, b_ in r.signalp_class.value_counts().items()),
                   "raw_majority": rawmaj, "raw_agreement": round(vc.iloc[0] / vc.sum(), 4) if len(vc) else np.nan,
                   "table_n_eq": int(F.n_eq_proteins[f]), "table_n_sp": int(F.n_eq_signal_peptide[f]),
                   "call_major": F.route_major_signalp_class[f], "call_route": F.route_call[f], "call_agreement": float(F.route_agreement[f]),
                   "anchor_call": F.anchor_call[f], "role_class": F.role_class[f],
                   "match": bool(len(r) == F.n_eq_proteins[f] and (r.sp_any == "True").sum() == F.n_eq_signal_peptide[f]
                                 and rawmaj == F.route_major_signalp_class[f]
                                 and abs((vc.iloc[0] / vc.sum()) - F.route_agreement[f]) < 1e-3)})
    sc = pd.DataFrame(sc)
    sc.to_csv(RES / "spotcheck_10_families.tsv", sep="\t", index=False)
    checks["spotcheck_10_all_match"] = bool(sc.match.all())

    # ------------------------------------------------------------------ receipt
    inputs = ["genecall/harmonized_proteins.faa", "genecall/harmonized_gene_manifest.tsv.gz", "metadata/harmonized_genome_manifest.tsv",
              "mmseqs/harmonized_uniref50_per_gene.tsv.gz", "signalp/harmonized_signalp6_predictions.tsv.gz",
              "source-metadata/ani_groups/ani_95_groups_1957.tsv",
              "analysis2_families/clusters/families_id30_c80.tsv", "analysis2_families/results/dark_exported_widespread_families_id30_c80.tsv",
              "analysis2_families/results/eq_family_table_id30_c80_all_genes.tsv.gz", "analysis2_families/results/class_counts.tsv"]
    extra = json.loads(a.receipt_extra)
    failed = [k_ for k_, v in checks.items() if (isinstance(v, bool) and not v) or
              (isinstance(v, list) and len(v) == 2 and k_.startswith("n_") and v[0] != v[1])]
    rec = {"status": "COMPLETE" if not failed else "COMPLETE_WITH_FAILED_CHECKS", "failed_checks": failed,
           "seed": a.seed, "bootstrap_draws": a.draws, "block": "named genus, else ANI95 group (EQ genomes only)",
           "checks": checks, "prep_info": json.load(open(W / "prep_info.json")),
           "features_info": json.load(open(W / "features_info.json")),
           "inputs_sha256": {p: sha256(R / p) for p in inputs},
           "databases": {
               "pfam": {"path": "/project/k10018/sw/db/pfam-ibex/pfam/Pfam-A.hmm", "n_models": int(len(pm)),
                        "sha256_recorded_in_stage_file": "a78a42d6faf265b6bfca59e8f062d06fae6083ce2c6e335d7b381f20b82b7903",
                        "note": "release not labelled in the file; newest model date 2025; Pfam-A.clans.tsv from EBI current_release (file dated 2026-09-18) used for clan names"},
               "pfam_clans": {"path": str(a.clans), "sha256": sha256(a.clans)},
               "tigrfam": "TIGRFAMs 15.0 (InterProScan 5.61-93.0 data dir), subset listed in results/sorting_models_used.tsv",
               "dbcan": {"path": "/project/k10018/sw/db/eq-localisation-20261111/dbCAN-HMMdb-V14.txt",
                         "sha256": sha256("/project/k10018/sw/db/eq-localisation-20261111/dbCAN-HMMdb-V14.txt")},
               "merops": {"path": "/project/k10018/sw/db/eq-localisation-20261111/pepunit.lib", "note": "EBI current_release, file dated 2023-02-22",
                          "sha256": sha256("/project/k10018/sw/db/eq-localisation-20261111/pepunit.lib")},
               "tmbed_reused_predictions": "/project/k10018/work/eq-annot-20261030/results/sp_tm_analysis/tmbed/chunk_*.pred.gz (reused for byte-identical sequences only; run 2026-10-30 on a GPU system; TMbed version not recorded on Shaheen)",
               "prot_t5": "Rostlab/prot_t5_xl_half_uniref50-enc (HF cache /project/k10018/sw/db/eq-localisation-20261111/hf_home)"},
           "tools": {"tmbed": "1.0.2 (github BernhoferM/TMbed tag v1.0.2; env /project/k10018/sw/envs/tmbed; torch 2.15.0+cpu, transformers 4.46.3; CPU, out-format 0)", "hmmer": "3.4 (hmmer.sif)", "diamond": "2.2.8 (diamond.sif)", "signalp": "6.0h fast, organism other (reused table)",
                     "python": sys.version.split()[0], "pandas": pd.__version__, "numpy": np.__version__},
           "parameters": {"hmmsearch": "default reporting thresholds, -Z = number of target sequences, --noali; cutoffs applied post hoc (see rules.txt)",
                          "diamond": "blastp --very-sensitive -e 1e-3 -k 5; post hoc best hit E<=1e-10, unit coverage>=0.5",
                          "cterm_tm_max_tail": 30, "cterm_signal_max_distance_from_end": 10, "st_rich_window": "30 aa, >=40% S/T",
                          "lowcomplexity": "12-aa windows with Shannon entropy <= 2.2 bits", "pI": "Bjellqvist pK set"},
           **extra}
    json.dump(rec, open(O / "receipt.json", "w"), indent=1, default=str)
    print(json.dumps(checks, indent=1, default=str))
    print("FAILED_CHECKS", failed)
    print("COMBINE_DONE")


if __name__ == "__main__":
    main()
