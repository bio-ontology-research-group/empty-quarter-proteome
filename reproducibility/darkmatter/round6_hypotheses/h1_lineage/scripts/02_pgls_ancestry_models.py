#!/usr/bin/env python3
"""02_pgls_ancestry_models.py - Parsimonious Phylogenetic Generalized Least Squares (PGLS)
modeling of within-Thermoplasmatota darkness, conditioning on quality and phylogenetic covariance.

Primary Estimand:
Partial association of external divergence (nearest AAI) and shared gene coverage (frac_rbh_desert)
with dark fraction, controlling for completeness, contamination, log10(length), and log10(N50).

Models & Sensitivities:
1. Primary PGLS model across all 237 desert Thermoplasmatota representatives.
2. Pagel's lambda sensitivity: fixed lambda in {0.0, 0.25, 0.50, 0.75, 1.0} and profiled lambda_ML.
3. Dominant-family omission sensitivity: omitting f__JAJPHT01 (n=140).
4. Restricted population: within c__SW-10-69-26 (n=185).
5. Alternative caller quality proxy sensitivity: adding CheckM2 frac_complete (n=232).
"""

import os
import sys
import gzip
import numpy as np
import pandas as pd
from scipy import linalg
from scipy.optimize import minimize_scalar
from scipy.stats import t as t_dist

def fit_pgls(X, y, C, lam=1.0):
    """Fit PGLS model via Cholesky decomposition of V(lambda) = lambda * C + (1 - lambda) * I.
    
    Returns dict with beta, se, t_val, p_val, sigma2, loglik, aic, bic.
    """
    n, p = X.shape
    lam = np.clip(lam, 0.0, 1.0)
    V = lam * C + (1.0 - lam) * np.eye(n)
    # Add small regularization if needed for inversion
    V += 1e-7 * np.eye(n)
    
    try:
        L = linalg.cholesky(V, lower=True)
    except linalg.LinAlgError:
        # Fallback to eigendecomposition if not strictly positive definite
        evals, evecs = linalg.eigh(V)
        evals = np.maximum(evals, 1e-7)
        V = evecs @ np.diag(evals) @ evecs.T
        L = linalg.cholesky(V, lower=True)
        
    log_det_V = 2.0 * np.sum(np.log(np.diag(L)))
    
    # Transform data: \tilde{y} = L^{-1} y, \tilde{X} = L^{-1} X
    y_tilde = linalg.solve_triangular(L, y, lower=True)
    X_tilde = linalg.solve_triangular(L, X, lower=True)
    
    # OLS on transformed data
    XtX = X_tilde.T @ X_tilde
    Xty = X_tilde.T @ y_tilde
    
    beta = linalg.solve(XtX, Xty)
    residuals = y_tilde - X_tilde @ beta
    rss = np.sum(residuals**2)
    
    sigma2 = rss / n
    # Covariance of beta
    XtX_inv = linalg.inv(XtX)
    var_beta = sigma2 * XtX_inv
    se = np.sqrt(np.maximum(np.diag(var_beta), 1e-12))
    
    t_val = beta / se
    df_resid = n - p
    p_val = 2.0 * (1.0 - t_dist.cdf(np.abs(t_val), df=df_resid))
    
    # Log-likelihood
    loglik = -0.5 * n * np.log(2.0 * np.pi) - 0.5 * n * np.log(sigma2) - 0.5 * log_det_V - 0.5 * n
    aic = 2.0 * (p + 1) - 2.0 * loglik
    bic = np.log(n) * (p + 1) - 2.0 * loglik
    
    # R-squared relative to phylogenetic mean
    ones = np.ones((n, 1))
    ones_tilde = linalg.solve_triangular(L, ones, lower=True)
    mu_null = (ones_tilde.T @ y_tilde) / (ones_tilde.T @ ones_tilde)
    rss_null = np.sum((y_tilde - ones_tilde * mu_null)**2)
    r2 = max(0.0, 1.0 - rss / rss_null)
    
    return {
        "beta": beta,
        "se": se,
        "t_val": t_val,
        "p_val": p_val,
        "sigma2": sigma2,
        "loglik": loglik,
        "aic": aic,
        "bic": bic,
        "r2": r2,
        "df_resid": df_resid,
        "lambda": lam
    }

def profile_lambda(X, y, C):
    """Find maximum likelihood estimate of Pagel's lambda."""
    def neg_loglik(lam):
        res = fit_pgls(X, y, C, lam=lam)
        return -res["loglik"]
    
    opt = minimize_scalar(neg_loglik, bounds=(0.0, 1.0), method="bounded")
    best_lam = float(opt.x)
    return best_lam

def main():
    W_IBEX = "/ibex/scratch/projects/c2014/rob/metagenomics"
    if os.path.exists(W_IBEX):
        res_dir = os.path.join(W_IBEX, "06_darkmatter", "p6_hypotheses_clean", "h1_lineage", "results")
    else:
        here = os.path.dirname(os.path.abspath(__file__))
        res_dir = os.path.abspath(os.path.join(here, "..", "results"))

    audit_tsv = os.path.join(res_dir, "denominator_caller_audit.tsv")
    corr_tsv = os.path.join(res_dir, "phylo_correlation.tsv.gz")

    if not os.path.exists(audit_tsv) or not os.path.exists(corr_tsv):
        sys.exit("Error: audit table or correlation matrix missing.")

    print(f"Loading data from {audit_tsv} and {corr_tsv}...")
    df = pd.read_csv(audit_tsv, sep="\t")
    corr_df = pd.read_csv(corr_tsv, sep="\t", compression="gzip", index_col=0)

    # Align genomes
    common_genomes = sorted(list(set(df["genome"]).intersection(set(corr_df.index))))
    df = df[df["genome"].isin(common_genomes)].sort_values("genome").reset_index(drop=True)
    C_full = corr_df.loc[df["genome"], df["genome"]].values
    print(f"Aligned {len(df)} genomes with phylogenetic correlation matrix ({C_full.shape})")

    # Add log covariates
    df["log_len"] = np.log10(df["length"])
    df["log_n50"] = np.log10(df["N50"])
    df["frac_rbh_desert"] = df["frac_rbh_of_desert_aai"]
    df["fam_clean"] = df["family_r232"].fillna("unassigned").astype(str)

    # Define primary model variables
    # dark_fraction ~ completeness + contamination + log_len + log_n50 + aai + frac_rbh_desert
    covars = ["completeness", "contamination", "log_len", "log_n50", "aai", "frac_rbh_desert"]
    var_names = ["Intercept"] + covars

    # Filter complete cases for primary variables
    valid_mask = df[covars + ["dark_fraction"]].notna().all(axis=1)
    df_prim = df[valid_mask].copy().reset_index(drop=True)
    C_prim = C_full[valid_mask.values, :][:, valid_mask.values]
    n_prim = len(df_prim)
    print(f"Primary analysis population: n = {n_prim} representatives")

    y_prim = df_prim["dark_fraction"].values
    X_prim = np.column_stack([np.ones(n_prim)] + [df_prim[c].values for c in covars])

    # 1. Profile lambda
    opt_lambda = profile_lambda(X_prim, y_prim, C_prim)
    print(f"\nProfiled Maximum Likelihood Pagel's lambda: {opt_lambda:.4f}")

    # Evaluate grid of lambda values
    lam_grid = [0.0, 0.25, 0.50, 0.75, 1.0, opt_lambda]
    grid_records = []
    for lam in sorted(set(lam_grid)):
        fit = fit_pgls(X_prim, y_prim, C_prim, lam=lam)
        grid_records.append({
            "lambda": lam,
            "is_mle": abs(lam - opt_lambda) < 1e-4,
            "loglik": fit["loglik"],
            "aic": fit["aic"],
            "bic": fit["bic"],
            "r2": fit["r2"],
            "sigma2": fit["sigma2"],
            "coef_aai": fit["beta"][var_names.index("aai")],
            "se_aai": fit["se"][var_names.index("aai")],
            "p_aai": fit["p_val"][var_names.index("aai")],
            "coef_frac_rbh": fit["beta"][var_names.index("frac_rbh_desert")],
            "se_frac_rbh": fit["se"][var_names.index("frac_rbh_desert")],
            "p_frac_rbh": fit["p_val"][var_names.index("frac_rbh_desert")],
            "coef_log_len": fit["beta"][var_names.index("log_len")],
            "p_log_len": fit["p_val"][var_names.index("log_len")],
            "coef_log_n50": fit["beta"][var_names.index("log_n50")],
            "p_log_n50": fit["p_val"][var_names.index("log_n50")]
        })
    grid_df = pd.DataFrame(grid_records)
    grid_tsv = os.path.join(res_dir, "pgls_lambda_profile.tsv")
    grid_df.to_csv(grid_tsv, sep="\t", index=False, float_format="%.5f")
    print(f"Saved lambda profile table to {grid_tsv}")

    # 2. Detailed Primary Model Output (at optimal lambda and OLS lambda=0)
    fit_opt = fit_pgls(X_prim, y_prim, C_prim, lam=opt_lambda)
    fit_ols = fit_pgls(X_prim, y_prim, C_prim, lam=0.0)

    model_rows = []
    for i, vname in enumerate(var_names):
        model_rows.append({
            "variable": vname,
            "pgls_coef": fit_opt["beta"][i],
            "pgls_se": fit_opt["se"][i],
            "pgls_t": fit_opt["t_val"][i],
            "pgls_pval": fit_opt["p_val"][i],
            "pgls_ci95_low": fit_opt["beta"][i] - 1.96 * fit_opt["se"][i],
            "pgls_ci95_high": fit_opt["beta"][i] + 1.96 * fit_opt["se"][i],
            "ols_coef": fit_ols["beta"][i],
            "ols_se": fit_ols["se"][i],
            "ols_pval": fit_ols["p_val"][i]
        })
    model_df = pd.DataFrame(model_rows)
    model_tsv = os.path.join(res_dir, "pgls_model_summary.tsv")
    model_df.to_csv(model_tsv, sep="\t", index=False, float_format="%.6f")
    print(f"Saved primary PGLS model table to {model_tsv}")

    # 3. Sensitivity Analyses
    sens_records = []

    # A. Primary Model Baseline
    sens_records.append({
        "cohort": "All Thermoplasmatota (Primary)",
        "n": n_prim,
        "lambda": opt_lambda,
        "r2": fit_opt["r2"],
        "coef_aai": fit_opt["beta"][var_names.index("aai")],
        "se_aai": fit_opt["se"][var_names.index("aai")],
        "p_aai": fit_opt["p_val"][var_names.index("aai")],
        "coef_frac_rbh": fit_opt["beta"][var_names.index("frac_rbh_desert")],
        "se_frac_rbh": fit_opt["se"][var_names.index("frac_rbh_desert")],
        "p_frac_rbh": fit_opt["p_val"][var_names.index("frac_rbh_desert")],
        "coef_log_len": fit_opt["beta"][var_names.index("log_len")],
        "p_log_len": fit_opt["p_val"][var_names.index("log_len")],
        "coef_log_n50": fit_opt["beta"][var_names.index("log_n50")],
        "p_log_n50": fit_opt["p_val"][var_names.index("log_n50")]
    })

    # B. Population restriction: within c__SW-10-69-26 only
    sw_mask = (df_prim["cls"] == "SW-10-69-26").values
    df_sw = df_prim[sw_mask].copy().reset_index(drop=True)
    C_sw = C_prim[sw_mask, :][:, sw_mask]
    n_sw = len(df_sw)
    y_sw = df_sw["dark_fraction"].values
    X_sw = np.column_stack([np.ones(n_sw)] + [df_sw[c].values for c in covars])
    opt_lam_sw = profile_lambda(X_sw, y_sw, C_sw)
    fit_sw = fit_pgls(X_sw, y_sw, C_sw, lam=opt_lam_sw)
    sens_records.append({
        "cohort": "c__SW-10-69-26 Only",
        "n": n_sw,
        "lambda": opt_lam_sw,
        "r2": fit_sw["r2"],
        "coef_aai": fit_sw["beta"][var_names.index("aai")],
        "se_aai": fit_sw["se"][var_names.index("aai")],
        "p_aai": fit_sw["p_val"][var_names.index("aai")],
        "coef_frac_rbh": fit_sw["beta"][var_names.index("frac_rbh_desert")],
        "se_frac_rbh": fit_sw["se"][var_names.index("frac_rbh_desert")],
        "p_frac_rbh": fit_sw["p_val"][var_names.index("frac_rbh_desert")],
        "coef_log_len": fit_sw["beta"][var_names.index("log_len")],
        "p_log_len": fit_sw["p_val"][var_names.index("log_len")],
        "coef_log_n50": fit_sw["beta"][var_names.index("log_n50")],
        "p_log_n50": fit_sw["p_val"][var_names.index("log_n50")]
    })

    # C. Dominant-family omission sensitivity: omitting f__JAJPHT01 (n=140)
    dom_mask = (df_prim["fam_clean"] != "JAJPHT01").values
    df_nodom = df_prim[dom_mask].copy().reset_index(drop=True)
    C_nodom = C_prim[dom_mask, :][:, dom_mask]
    n_nodom = len(df_nodom)
    y_nodom = df_nodom["dark_fraction"].values
    X_nodom = np.column_stack([np.ones(n_nodom)] + [df_nodom[c].values for c in covars])
    opt_lam_nodom = profile_lambda(X_nodom, y_nodom, C_nodom)
    fit_nodom = fit_pgls(X_nodom, y_nodom, C_nodom, lam=opt_lam_nodom)
    sens_records.append({
        "cohort": "Dominant Family Omitted (No JAJPHT01)",
        "n": n_nodom,
        "lambda": opt_lam_nodom,
        "r2": fit_nodom["r2"],
        "coef_aai": fit_nodom["beta"][var_names.index("aai")],
        "se_aai": fit_nodom["se"][var_names.index("aai")],
        "p_aai": fit_nodom["p_val"][var_names.index("aai")],
        "coef_frac_rbh": fit_nodom["beta"][var_names.index("frac_rbh_desert")],
        "se_frac_rbh": fit_nodom["se"][var_names.index("frac_rbh_desert")],
        "p_frac_rbh": fit_nodom["p_val"][var_names.index("frac_rbh_desert")],
        "coef_log_len": fit_nodom["beta"][var_names.index("log_len")],
        "p_log_len": fit_nodom["p_val"][var_names.index("log_len")],
        "coef_log_n50": fit_nodom["beta"][var_names.index("log_n50")],
        "p_log_n50": fit_nodom["p_val"][var_names.index("log_n50")]
    })

    # D. Caller Quality Proxy Sensitivity (adding frac_complete)
    if "frac_complete" in df_prim and df_prim["frac_complete"].notna().sum() > 50:
        call_mask = df_prim["frac_complete"].notna().values
        df_call = df_prim[call_mask].copy().reset_index(drop=True)
        C_call = C_prim[call_mask, :][:, call_mask]
        n_call = len(df_call)
        covars_call = covars + ["frac_complete"]
        vnames_call = ["Intercept"] + covars_call
        y_call = df_call["dark_fraction"].values
        X_call = np.column_stack([np.ones(n_call)] + [df_call[c].values for c in covars_call])
        opt_lam_call = profile_lambda(X_call, y_call, C_call)
        fit_call = fit_pgls(X_call, y_call, C_call, lam=opt_lam_call)
        sens_records.append({
            "cohort": "With Caller Complete Proxy (frac_complete)",
            "n": n_call,
            "lambda": opt_lam_call,
            "r2": fit_call["r2"],
            "coef_aai": fit_call["beta"][vnames_call.index("aai")],
            "se_aai": fit_call["se"][vnames_call.index("aai")],
            "p_aai": fit_call["p_val"][vnames_call.index("aai")],
            "coef_frac_rbh": fit_call["beta"][vnames_call.index("frac_rbh_desert")],
            "se_frac_rbh": fit_call["se"][vnames_call.index("frac_rbh_desert")],
            "p_frac_rbh": fit_call["p_val"][vnames_call.index("frac_rbh_desert")],
            "coef_log_len": fit_call["beta"][vnames_call.index("log_len")],
            "p_log_len": fit_call["p_val"][vnames_call.index("log_len")],
            "coef_log_n50": fit_call["beta"][vnames_call.index("log_n50")],
            "p_log_n50": fit_call["p_val"][vnames_call.index("log_n50")]
        })

    sens_df = pd.DataFrame(sens_records)
    sens_tsv = os.path.join(res_dir, "pgls_sensitivity.tsv")
    sens_df.to_csv(sens_tsv, sep="\t", index=False, float_format="%.5f")
    print(f"Saved sensitivity table to {sens_tsv}")

    # Summary text
    rep = []
    rep.append("=" * 80)
    rep.append("PGLS ANCESTRY-AWARE ASSOCIATION REPORT: 237 THERMOPLASMATOTA")
    rep.append("=" * 80)
    rep.append("")
    rep.append("1. PRIMARY MODEL SPECIFICATION AND FIT")
    rep.append("  Response: dark_fraction (N_dark / N_proteins)")
    rep.append(f"  Population: All desert Thermoplasmatota representatives with valid covariates (n = {n_prim})")
    rep.append(f"  Phylogenetic Correlation: ar53 FastTree branch-length covariance")
    rep.append(f"  Profiled Pagel's lambda MLE: {opt_lambda:.4f}")
    rep.append(f"  Model R^2: {fit_opt['r2']:.4f}, Residual variance sigma^2: {fit_opt['sigma2']:.6f}")
    rep.append("")
    rep.append("2. PRIMARY ESTIMAND: EXTERNAL DIVERGENCE & SHARED RBH COVERAGE")
    aai_idx = var_names.index("aai")
    rbh_idx = var_names.index("frac_rbh_desert")
    rep.append(f"  Nearest-relative AAI slope (per 1% AAI): beta = {fit_opt['beta'][aai_idx]:+.5f} (SE = {fit_opt['se'][aai_idx]:.5f}, t = {fit_opt['t_val'][aai_idx]:+.2f}, p = {fit_opt['p_val'][aai_idx]:.4f})")
    rep.append(f"    95% Model-Conditional CI: [{fit_opt['beta'][aai_idx] - 1.96*fit_opt['se'][aai_idx]:+.5f}, {fit_opt['beta'][aai_idx] + 1.96*fit_opt['se'][aai_idx]:+.5f}]")
    rep.append(f"  Shared RBH coverage slope (frac_rbh_desert): beta = {fit_opt['beta'][rbh_idx]:+.5f} (SE = {fit_opt['se'][rbh_idx]:.5f}, t = {fit_opt['t_val'][rbh_idx]:+.2f}, p = {fit_opt['p_val'][rbh_idx]:.4f})")
    rep.append(f"    95% Model-Conditional CI: [{fit_opt['beta'][rbh_idx] - 1.96*fit_opt['se'][rbh_idx]:+.5f}, {fit_opt['beta'][rbh_idx] + 1.96*fit_opt['se'][rbh_idx]:+.5f}]")
    rep.append("")
    rep.append("3. QUALITY & FRAGMENTATION COVARIATE EFFECTS")
    len_idx = var_names.index("log_len")
    n50_idx = var_names.index("log_n50")
    comp_idx = var_names.index("completeness")
    rep.append(f"  Assembly length (log10 length): beta = {fit_opt['beta'][len_idx]:+.5f} (SE = {fit_opt['se'][len_idx]:.5f}, t = {fit_opt['t_val'][len_idx]:+.2f}, p = {fit_opt['p_val'][len_idx]:.4e})")
    rep.append(f"  Contig fragmentation (log10 N50): beta = {fit_opt['beta'][n50_idx]:+.5f} (SE = {fit_opt['se'][n50_idx]:.5f}, t = {fit_opt['t_val'][n50_idx]:+.2f}, p = {fit_opt['p_val'][n50_idx]:.4e})")
    rep.append(f"  CheckM2 completeness: beta = {fit_opt['beta'][comp_idx]:+.5f} (SE = {fit_opt['se'][comp_idx]:.5f}, t = {fit_opt['t_val'][comp_idx]:+.2f}, p = {fit_opt['p_val'][comp_idx]:.4e})")
    rep.append("")
    rep.append("4. SENSITIVITY SUMMARY ACROSS COHORTS AND ASSUMPTIONS")
    for _, r in sens_df.iterrows():
        rep.append(f"  Cohort: {r['cohort']} (n={int(r['n'])}, lambda={r['lambda']:.2f}):")
        rep.append(f"    AAI beta = {r['coef_aai']:+.5f} (p = {r['p_aai']:.4f}), RBH frac beta = {r['coef_frac_rbh']:+.5f} (p = {r['p_frac_rbh']:.4f}), log_len beta = {r['coef_log_len']:+.5f} (p = {r['p_log_len']:.4e})")
    rep.append("=" * 80)

    summary_text = "\n".join(rep)
    print(summary_text)
    out_txt = os.path.join(res_dir, "pgls_ancestry_summary.txt")
    with open(out_txt, "w") as f:
        f.write(summary_text + "\n")
    print(f"Saved summary to {out_txt}")

if __name__ == "__main__":
    main()
