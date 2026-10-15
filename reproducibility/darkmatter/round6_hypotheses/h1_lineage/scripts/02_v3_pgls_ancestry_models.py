#!/usr/bin/env python3
"""02_v3_pgls_ancestry_models.py - Repaired v3 Phylogenetic Generalized Least Squares (PGLS)
modeling of within-Thermoplasmatota darkness across the 237 desert representatives.

Repairs implemented (C1, C2, C3):
1. Fixes the 1D-to-2D array broadcasting bug in whitened GLS R^2 (y_tilde and ones_tilde
   strictly 1D shape (n,), verified against statsmodels.GLS.rsquared to <1e-12).
2. Separates concentrated ML residual variance sigma2_ml = RSS / n (used strictly inside
   profile log-likelihood) from unbiased fixed-covariance residual variance
   sigma2_unb = RSS / (n - p) (used for SEs, two-sided Student-t statistics with df = n - p,
   two-sided p-values, and exact t_{0.975, n-p} 95% confidence intervals).
3. Computes exact 95% profile-likelihood confidence intervals for Pagel's lambda via
   likelihood-ratio root-finding (2 * [loglik(lambda_ML) - loglik(lambda)] <= chi2_{1,0.95})
   and reports model-conditional coefficient intervals alongside values at the 95% lambda CI bounds.
4. Verifies all fixed-lambda fits (lambda = 0.0, 0.5, lambda_ML, and rounded 0.9158) against
   statsmodels.regression.linear_model.GLS.
"""

import os
import sys
import time
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import linalg
from scipy.optimize import minimize_scalar, brentq
from scipy.stats import t as t_dist, chi2


def fit_pgls(X, y, C_mat, lam=1.0, eps=1e-7):
    """Fit PGLS model for fixed Pagel's lambda via Cholesky whitening of V(lambda)."""
    y_1d = np.asarray(y, dtype=np.float64).reshape(-1)
    X_2d = np.asarray(X, dtype=np.float64)
    C_2d = np.asarray(C_mat, dtype=np.float64)

    n, p = X_2d.shape
    assert y_1d.shape == (n,), f"Expected 1D y of shape ({n},), got {y_1d.shape}"
    assert C_2d.shape == (n, n), f"Expected C of shape ({n},{n}), got {C_2d.shape}"
    assert np.max(np.abs(C_2d - C_2d.T)) < 1e-10, "Covariance matrix C_mat is not symmetric"
    assert np.linalg.matrix_rank(X_2d) == p, f"Design matrix X is rank-deficient: rank < {p}"

    lam = float(np.clip(lam, 0.0, 1.0))
    V = lam * C_2d + (1.0 - lam) * np.eye(n, dtype=np.float64) + eps * np.eye(n, dtype=np.float64)

    evals = np.linalg.eigvalsh(V)
    min_eig = float(evals.min())
    assert min_eig > 1e-6, f"V(lambda={lam}) not strictly positive definite: min_eig={min_eig}"

    L = linalg.cholesky(V, lower=True)
    log_det_V = float(2.0 * np.sum(np.log(np.diag(L))))

    y_tilde = linalg.solve_triangular(L, y_1d, lower=True)
    X_tilde = linalg.solve_triangular(L, X_2d, lower=True)
    assert y_tilde.ndim == 1 and y_tilde.shape == (n,), f"Bad y_tilde shape: {y_tilde.shape}"
    assert X_tilde.ndim == 2 and X_tilde.shape == (n, p), f"Bad X_tilde shape: {X_tilde.shape}"

    XtX = X_tilde.T @ X_tilde
    Xty = X_tilde.T @ y_tilde
    beta = linalg.solve(XtX, Xty)
    assert beta.ndim == 1 and beta.shape == (p,), f"Bad beta shape: {beta.shape}"

    residuals = y_tilde - (X_tilde @ beta)
    assert residuals.ndim == 1 and residuals.shape == (n,), f"Bad residuals shape: {residuals.shape}"
    rss = float(np.sum(residuals ** 2))

    df_resid = int(n - p)
    assert df_resid > 0, f"Non-positive residual degrees of freedom: {df_resid}"

    # 1. Concentrated ML variance (strictly for profile log-likelihood)
    sigma2_ml = float(rss / n)
    loglik = float(-0.5 * n * np.log(2.0 * np.pi) - 0.5 * n * np.log(sigma2_ml) - 0.5 * log_det_V - 0.5 * n)
    aic = float(2.0 * (p + 1) - 2.0 * loglik)
    bic = float(np.log(n) * (p + 1) - 2.0 * loglik)

    # 2. Unbiased fixed-covariance residual variance (for SEs, Student-t tests, and 95% CIs)
    sigma2_unb = float(rss / df_resid)
    XtX_inv = linalg.inv(XtX)
    var_beta_unb = sigma2_unb * XtX_inv
    se_unb = np.sqrt(np.maximum(np.diag(var_beta_unb), 1e-18))

    # Also compute ML plug-in SE for audit transparency
    se_ml = np.sqrt(np.maximum(np.diag(sigma2_ml * XtX_inv), 1e-18))

    t_unb = beta / se_unb
    p_unb = 2.0 * t_dist.sf(np.abs(t_unb), df=df_resid)
    t_crit = float(t_dist.ppf(0.975, df=df_resid))
    ci95_low = beta - t_crit * se_unb
    ci95_high = beta + t_crit * se_unb

    assert np.all(ci95_low < beta) and np.all(beta < ci95_high), "CI ordering violated"
    assert np.all((p_unb > 0.0) & (p_unb <= 1.0)), "Invalid p-values detected"

    # 3. Correct 1D whitened GLS R-squared relative to fixed-lambda GLS intercept-only null
    ones_1d = np.ones(n, dtype=np.float64)
    ones_tilde = linalg.solve_triangular(L, ones_1d, lower=True)
    assert ones_tilde.ndim == 1 and ones_tilde.shape == (n,), f"Bad ones_tilde shape: {ones_tilde.shape}"
    mu_null = float(np.dot(ones_tilde, y_tilde) / np.dot(ones_tilde, ones_tilde))
    null_resid = y_tilde - (ones_tilde * mu_null)
    assert null_resid.ndim == 1 and null_resid.shape == (n,), f"Bad null_resid shape: {null_resid.shape}"
    tss_gls = float(np.sum(null_resid ** 2))
    r2_gls = float(1.0 - rss / tss_gls)
    assert 0.0 <= r2_gls <= 1.0, f"Invalid R^2={r2_gls}"

    return {
        "n": n,
        "p": p,
        "df_resid": df_resid,
        "lambda": lam,
        "min_eig_V": min_eig,
        "beta": beta,
        "se_unb": se_unb,
        "se_ml": se_ml,
        "t_unb": t_unb,
        "p_unb": p_unb,
        "t_crit_975": t_crit,
        "ci95_low": ci95_low,
        "ci95_high": ci95_high,
        "rss": rss,
        "tss_gls": tss_gls,
        "r2_gls": r2_gls,
        "sigma2_unb": sigma2_unb,
        "sigma2_ml": sigma2_ml,
        "loglik": loglik,
        "aic": aic,
        "bic": bic,
        "V": V,
    }


def profile_lambda_with_ci(X, y, C_mat, eps=1e-7):
    """Find profile MLE of Pagel's lambda and its 95% profile-likelihood CI."""
    def neg_ll(lam):
        return -fit_pgls(X, y, C_mat, lam=lam, eps=eps)["loglik"]

    opt = minimize_scalar(neg_ll, bounds=(0.0, 1.0), method="bounded", options={"xatol": 1e-8})
    lam_mle = float(opt.x)
    ll_mle = float(-opt.fun)
    cutoff = ll_mle - 0.5 * float(chi2.ppf(0.95, df=1))  # ll_mle - 1.9207294

    # Lower 95% profile CI bound
    ll_0 = fit_pgls(X, y, C_mat, lam=0.0, eps=eps)["loglik"]
    if ll_0 >= cutoff:
        lam_low = 0.0
    else:
        lam_low = float(brentq(lambda l: fit_pgls(X, y, C_mat, lam=l, eps=eps)["loglik"] - cutoff, 0.0, lam_mle, xtol=1e-7))

    # Upper 95% profile CI bound
    ll_1 = fit_pgls(X, y, C_mat, lam=1.0, eps=eps)["loglik"]
    if ll_1 >= cutoff:
        lam_high = 1.0
    else:
        lam_high = float(brentq(lambda l: fit_pgls(X, y, C_mat, lam=l, eps=eps)["loglik"] - cutoff, lam_mle, 1.0, xtol=1e-7))

    assert 0.0 <= lam_low <= lam_mle <= lam_high <= 1.0, "Invalid profile-likelihood CI bounds for lambda"
    lrt_vs_zero = float(2.0 * (ll_mle - ll_0))
    p_lrt_vs_zero = float(chi2.sf(lrt_vs_zero, df=1))
    return lam_mle, lam_low, lam_high, ll_mle, lrt_vs_zero, p_lrt_vs_zero


def verify_against_statsmodels(X, y, C_mat, var_names, lam_list, eps=1e-7):
    """Cross-check fit_pgls against statsmodels.regression.linear_model.GLS."""
    records = []
    for lam in lam_list:
        mine = fit_pgls(X, y, C_mat, lam=lam, eps=eps)
        sm_fit = sm.GLS(y, X, sigma=mine["V"]).fit()

        max_beta_diff = float(np.max(np.abs(mine["beta"] - sm_fit.params)))
        max_se_diff = float(np.max(np.abs(mine["se_unb"] - sm_fit.bse)))
        max_t_diff = float(np.max(np.abs(mine["t_unb"] - sm_fit.tvalues)))
        max_p_diff = float(np.max(np.abs(mine["p_unb"] - sm_fit.pvalues)))
        r2_diff = float(abs(mine["r2_gls"] - sm_fit.rsquared))
        scale_diff = float(abs(mine["sigma2_unb"] - sm_fit.scale))
        ll_diff = float(abs(mine["loglik"] - sm_fit.llf))

        assert max_beta_diff < 1e-10, f"beta mismatch vs statsmodels at lam={lam}: {max_beta_diff}"
        assert max_se_diff < 1e-10, f"SE mismatch vs statsmodels at lam={lam}: {max_se_diff}"
        assert max_t_diff < 1e-9, f"t mismatch vs statsmodels at lam={lam}: {max_t_diff}"
        assert max_p_diff < 1e-10, f"p mismatch vs statsmodels at lam={lam}: {max_p_diff}"
        assert r2_diff < 1e-12, f"R^2 mismatch vs statsmodels at lam={lam}: {r2_diff}"
        assert scale_diff < 1e-12, f"scale mismatch vs statsmodels at lam={lam}: {scale_diff}"
        assert ll_diff < 1e-10, f"loglik mismatch vs statsmodels at lam={lam}: {ll_diff}"

        for idx, vname in enumerate(var_names):
            records.append({
                "lambda": lam,
                "variable": vname,
                "custom_beta": mine["beta"][idx],
                "statsmodels_beta": float(sm_fit.params[idx]),
                "abs_diff_beta": float(abs(mine["beta"][idx] - sm_fit.params[idx])),
                "custom_se_unb": mine["se_unb"][idx],
                "statsmodels_bse": float(sm_fit.bse[idx]),
                "abs_diff_se": float(abs(mine["se_unb"][idx] - sm_fit.bse[idx])),
                "custom_t_unb": mine["t_unb"][idx],
                "statsmodels_t": float(sm_fit.tvalues[idx]),
                "custom_p_unb": mine["p_unb"][idx],
                "statsmodels_p": float(sm_fit.pvalues[idx]),
                "custom_r2_gls": mine["r2_gls"],
                "statsmodels_rsquared": float(sm_fit.rsquared),
                "abs_diff_r2": r2_diff,
                "custom_sigma2_unb": mine["sigma2_unb"],
                "statsmodels_scale": float(sm_fit.scale),
                "custom_loglik": mine["loglik"],
                "statsmodels_llf": float(sm_fit.llf),
            })
    return pd.DataFrame(records)


def main():
    t0 = time.time()
    W_IBEX = "/ibex/scratch/projects/c2014/rob/metagenomics"
    if os.path.exists(W_IBEX):
        res_dir = os.path.join(W_IBEX, "06_darkmatter", "p6_hypotheses_clean", "h1_lineage", "results")
    else:
        here = os.path.dirname(os.path.abspath(__file__))
        res_dir = os.path.abspath(os.path.join(here, "..", "results"))

    audit_tsv = os.path.join(res_dir, "denominator_caller_audit.tsv")
    corr_tsv = os.path.join(res_dir, "phylo_correlation.tsv.gz")
    vcv_tsv = os.path.join(res_dir, "phylo_vcv.tsv.gz")

    for p in (audit_tsv, corr_tsv, vcv_tsv):
        if not os.path.exists(p):
            sys.exit(f"Fatal error: missing input file {p}")

    df = pd.read_csv(audit_tsv, sep="\t").sort_values("genome").reset_index(drop=True)
    corr_df = pd.read_csv(corr_tsv, sep="\t", compression="gzip", index_col=0)
    vcv_df = pd.read_csv(vcv_tsv, sep="\t", compression="gzip", index_col=0)

    assert len(df) == 237, f"Expected 237 rows, found {len(df)}"
    assert df["genome"].nunique() == 237, "Duplicate genome IDs in denominator_caller_audit.tsv"
    assert set(df["genome"]) == set(corr_df.index), "Genome ID mismatch with phylo_correlation.tsv.gz"
    assert set(df["genome"]) == set(vcv_df.index), "Genome ID mismatch with phylo_vcv.tsv.gz"

    C_full = corr_df.loc[df["genome"], df["genome"]].values
    V_raw_full = vcv_df.loc[df["genome"], df["genome"]].values
    V_scaled_full = V_raw_full / float(np.mean(np.diag(V_raw_full)))

    df["log_len"] = np.log10(df["length"])
    df["log_n50"] = np.log10(df["N50"])
    df["frac_rbh_desert"] = df["frac_rbh_of_desert_aai"]
    df["fam_clean"] = df["family_r232"].fillna("unassigned").astype(str)

    covars = ["completeness", "contamination", "log_len", "log_n50", "aai", "frac_rbh_desert"]
    var_names = ["Intercept"] + covars

    assert df[covars + ["dark_fraction"]].notna().all().all(), "Unexpected missing values in primary covariates"
    n_prim = len(df)
    y_prim = df["dark_fraction"].values
    X_prim = np.column_stack([np.ones(n_prim)] + [df[c].values for c in covars])

    # 1. Profile Pagel's lambda and 95% profile-likelihood CI on primary model
    lam_mle, lam_ci_low, lam_ci_high, ll_mle, lrt_zero, p_lrt_zero = profile_lambda_with_ci(X_prim, y_prim, C_full)
    print(f"[02_v3] Primary PGLS (n={n_prim}): lambda_ML = {lam_mle:.6f}, 95% profile-likelihood CI = [{lam_ci_low:.6f}, {lam_ci_high:.6f}], LRT vs lambda=0 = {lrt_zero:.4f} (p = {p_lrt_zero:.4e})")

    # 2. Verify against statsmodels.GLS at lambda = 0.0, 0.5, lambda_ML, and 0.9158
    sm_ver_df = verify_against_statsmodels(X_prim, y_prim, C_full, var_names, [0.0, 0.5, lam_mle, 0.9158])
    sm_ver_tsv = os.path.join(res_dir, "v3_pgls_statsmodels_verification.tsv")
    sm_ver_df.to_csv(sm_ver_tsv, sep="\t", index=False, float_format="%.10e")
    print(f"[02_v3] Verified custom Cholesky GLS against statsmodels.GLS (<1e-10 max diff); saved {sm_ver_tsv}")

    # 3. Evaluate lambda grid including 95% profile CI bounds and MLE
    lam_grid_points = [
        (0.0, "OLS (lambda=0)"),
        (0.25, "Fixed lambda=0.25"),
        (0.50, "Fixed lambda=0.50"),
        (0.75, "Fixed lambda=0.75"),
        (lam_ci_low, "Lambda 95% Profile CI Lower"),
        (lam_mle, "Profile MLE (lambda_ML)"),
        (lam_ci_high, "Lambda 95% Profile CI Upper"),
        (1.0, "Strict Brownian (lambda=1.0)"),
    ]
    lam_grid_points.sort(key=lambda x: x[0])

    aai_idx = var_names.index("aai")
    rbh_idx = var_names.index("frac_rbh_desert")
    len_idx = var_names.index("log_len")
    n50_idx = var_names.index("log_n50")
    comp_idx = var_names.index("completeness")

    grid_rows = []
    for lam_val, status_label in lam_grid_points:
        f = fit_pgls(X_prim, y_prim, C_full, lam=lam_val)
        grid_rows.append({
            "lambda": f"{lam_val:.6f}",
            "grid_role": status_label,
            "is_mle": abs(lam_val - lam_mle) < 1e-7,
            "in_lambda_95ci": (lam_ci_low - 1e-7) <= lam_val <= (lam_ci_high + 1e-7),
            "loglik": f"{f['loglik']:.6f}",
            "aic": f"{f['aic']:.6f}",
            "bic": f"{f['bic']:.6f}",
            "r2_gls": f"{f['r2_gls']:.6f}",
            "rss": f"{f['rss']:.6f}",
            "tss_gls": f"{f['tss_gls']:.6f}",
            "sigma2_unb": f"{f['sigma2_unb']:.8f}",
            "sigma2_ml": f"{f['sigma2_ml']:.8f}",
            "coef_aai": f"{f['beta'][aai_idx]:.6f}",
            "se_unb_aai": f"{f['se_unb'][aai_idx]:.6f}",
            "t_unb_aai": f"{f['t_unb'][aai_idx]:.4f}",
            "p_unb_aai": f"{f['p_unb'][aai_idx]:.6e}",
            "ci95_low_aai": f"{f['ci95_low'][aai_idx]:.6f}",
            "ci95_high_aai": f"{f['ci95_high'][aai_idx]:.6f}",
            "coef_frac_rbh": f"{f['beta'][rbh_idx]:.6f}",
            "se_unb_frac_rbh": f"{f['se_unb'][rbh_idx]:.6f}",
            "t_unb_frac_rbh": f"{f['t_unb'][rbh_idx]:.4f}",
            "p_unb_frac_rbh": f"{f['p_unb'][rbh_idx]:.6e}",
            "ci95_low_frac_rbh": f"{f['ci95_low'][rbh_idx]:.6f}",
            "ci95_high_frac_rbh": f"{f['ci95_high'][rbh_idx]:.6f}",
            "coef_log_len": f"{f['beta'][len_idx]:.6f}",
            "se_unb_log_len": f"{f['se_unb'][len_idx]:.6f}",
            "p_unb_log_len": f"{f['p_unb'][len_idx]:.6e}",
            "coef_log_n50": f"{f['beta'][n50_idx]:.6f}",
            "se_unb_log_n50": f"{f['se_unb'][n50_idx]:.6f}",
            "p_unb_log_n50": f"{f['p_unb'][n50_idx]:.6e}",
        })
    grid_df = pd.DataFrame(grid_rows)
    grid_tsv = os.path.join(res_dir, "v3_pgls_lambda_profile.tsv")
    grid_df.to_csv(grid_tsv, sep="\t", index=False)
    print(f"[02_v3] Saved lambda profile table to {grid_tsv}")

    # 4. Primary Model Summary Table (at lambda_ML and OLS lambda=0)
    fit_opt = fit_pgls(X_prim, y_prim, C_full, lam=lam_mle)
    fit_ols = fit_pgls(X_prim, y_prim, C_full, lam=0.0)

    model_rows = []
    for i, vname in enumerate(var_names):
        model_rows.append({
            "variable": vname,
            "pgls_coef": f"{fit_opt['beta'][i]:.6f}",
            "pgls_se_unb": f"{fit_opt['se_unb'][i]:.6f}",
            "pgls_se_ml": f"{fit_opt['se_ml'][i]:.6f}",
            "pgls_t_unb": f"{fit_opt['t_unb'][i]:.4f}",
            "pgls_pval_unb": f"{fit_opt['p_unb'][i]:.6e}",
            "pgls_ci95_low_unb": f"{fit_opt['ci95_low'][i]:.6f}",
            "pgls_ci95_high_unb": f"{fit_opt['ci95_high'][i]:.6f}",
            "ols_coef": f"{fit_ols['beta'][i]:.6f}",
            "ols_se_unb": f"{fit_ols['se_unb'][i]:.6f}",
            "ols_t_unb": f"{fit_ols['t_unb'][i]:.4f}",
            "ols_pval_unb": f"{fit_ols['p_unb'][i]:.6e}",
            "ols_ci95_low_unb": f"{fit_ols['ci95_low'][i]:.6f}",
            "ols_ci95_high_unb": f"{fit_ols['ci95_high'][i]:.6f}",
        })
    model_df = pd.DataFrame(model_rows)
    model_tsv = os.path.join(res_dir, "v3_pgls_model_summary.tsv")
    model_df.to_csv(model_tsv, sep="\t", index=False)
    print(f"[02_v3] Saved primary PGLS model summary to {model_tsv}")

    # 5. Pre-specified Sensitivity Analyses (S1 through S7)
    specs = [
        ("S1_Primary_All_237", np.ones(n_prim, dtype=bool), C_full, "frac_rbh_of_desert_aai", []),
        ("S2_Unnormalized_Brownian_VCV_237", np.ones(n_prim, dtype=bool), V_scaled_full, "frac_rbh_of_desert_aai", []),
        ("S3_Catalogue_Denominator_RBH_237", np.ones(n_prim, dtype=bool), C_full, "frac_rbh_of_desert_cat", []),
        ("S4_Restricted_c__SW-10-69-26_185", (df["cls"] == "SW-10-69-26").values, C_full, "frac_rbh_of_desert_aai", []),
        ("S5_Dominant_Family_JAJPHT01_Omitted_97", (df["fam_clean"] != "JAJPHT01").values, C_full, "frac_rbh_of_desert_aai", []),
        ("S6_Same_Class_Comparators_SW_plus_E2_236", (df["cls"] != "Thermoplasmata").values, C_full, "frac_rbh_of_desert_aai", []),
        ("S7_With_Caller_Complete_Proxy_232", df["frac_complete"].notna().values, C_full, "frac_rbh_of_desert_aai", ["frac_complete"]),
    ]

    sens_rows = []
    for label, mask, base_mat, rbh_col, extra_cols in specs:
        sub = df[mask].copy().reset_index(drop=True)
        C_sub = base_mat[mask, :][:, mask]
        n_s = len(sub)
        cols_s = ["completeness", "contamination", "log_len", "log_n50", "aai", rbh_col] + extra_cols
        vnames_s = ["Intercept"] + cols_s
        X_s = np.column_stack([np.ones(n_s)] + [sub[c].values for c in cols_s])
        y_s = sub["dark_fraction"].values

        lam_s, lam_lo_s, lam_hi_s, ll_s, lrt_s, p_lrt_s = profile_lambda_with_ci(X_s, y_s, C_sub)
        f_s = fit_pgls(X_s, y_s, C_sub, lam=lam_s)

        ia = vnames_s.index("aai")
        ir = vnames_s.index(rbh_col)
        il = vnames_s.index("log_len")
        in50 = vnames_s.index("log_n50")

        sens_rows.append({
            "model_id": label,
            "n": n_s,
            "p": f_s["p"],
            "df_resid": f_s["df_resid"],
            "lambda_ml": f"{lam_s:.6f}",
            "lambda_ci95_low": f"{lam_lo_s:.6f}",
            "lambda_ci95_high": f"{lam_hi_s:.6f}",
            "r2_gls": f"{f_s['r2_gls']:.6f}",
            "sigma2_unb": f"{f_s['sigma2_unb']:.8f}",
            "sigma2_ml": f"{f_s['sigma2_ml']:.8f}",
            "loglik": f"{f_s['loglik']:.6f}",
            "aic": f"{f_s['aic']:.6f}",
            "coef_aai": f"{f_s['beta'][ia]:.6f}",
            "se_unb_aai": f"{f_s['se_unb'][ia]:.6f}",
            "t_unb_aai": f"{f_s['t_unb'][ia]:.4f}",
            "p_unb_aai": f"{f_s['p_unb'][ia]:.6e}",
            "ci95_low_aai": f"{f_s['ci95_low'][ia]:.6f}",
            "ci95_high_aai": f"{f_s['ci95_high'][ia]:.6f}",
            "coef_frac_rbh": f"{f_s['beta'][ir]:.6f}",
            "se_unb_frac_rbh": f"{f_s['se_unb'][ir]:.6f}",
            "t_unb_frac_rbh": f"{f_s['t_unb'][ir]:.4f}",
            "p_unb_frac_rbh": f"{f_s['p_unb'][ir]:.6e}",
            "ci95_low_frac_rbh": f"{f_s['ci95_low'][ir]:.6f}",
            "ci95_high_frac_rbh": f"{f_s['ci95_high'][ir]:.6f}",
            "coef_log_len": f"{f_s['beta'][il]:.6f}",
            "se_unb_log_len": f"{f_s['se_unb'][il]:.6f}",
            "p_unb_log_len": f"{f_s['p_unb'][il]:.6e}",
            "coef_log_n50": f"{f_s['beta'][in50]:.6f}",
            "se_unb_log_n50": f"{f_s['se_unb'][in50]:.6f}",
            "p_unb_log_n50": f"{f_s['p_unb'][in50]:.6e}",
        })

    sens_df = pd.DataFrame(sens_rows)
    sens_tsv = os.path.join(res_dir, "v3_pgls_sensitivity.tsv")
    sens_df.to_csv(sens_tsv, sep="\t", index=False)
    print(f"[02_v3] Saved sensitivity table to {sens_tsv}")

    # 6. Human-readable summary report
    rep = [
        "=" * 88,
        "V3 PGLS ANCESTRY-AWARE ASSOCIATION REPORT (REPAIRED R^2, UNBIASED SE, LAMBDA CI)",
        "=" * 88,
        "",
        "1. PRIMARY MODEL SPECIFICATION AND FIT (S1_Primary_All_237)",
        f"  Population: n = {n_prim} desert Thermoplasmatota species representatives (p = {fit_opt['p']}, df_resid = {fit_opt['df_resid']})",
        f"  Response: dark_fraction (N_dark / N_proteins_cat)",
        f"  Profiled Pagel's lambda MLE: {lam_mle:.6f} (95% Profile-Likelihood CI: [{lam_ci_low:.6f}, {lam_ci_high:.6f}])",
        f"  LRT vs lambda=0 (OLS): 2*Delta_logL = {lrt_zero:.4f} (p = {p_lrt_zero:.4e})",
        f"  Log-likelihood: {fit_opt['loglik']:.6f}, AIC: {fit_opt['aic']:.6f}, BIC: {fit_opt['bic']:.6f}",
        f"  Fixed-covariance whitened GLS R^2 (at lambda_ML): {fit_opt['r2_gls']:.6f} (RSS = {fit_opt['rss']:.6f}, TSS_GLS = {fit_opt['tss_gls']:.6f})",
        f"  Unbiased residual variance sigma^2_unb = RSS/(n-p): {fit_opt['sigma2_unb']:.8f} (ML sigma^2_ML = RSS/n: {fit_opt['sigma2_ml']:.8f})",
        f"  Student-t critical value t_{{0.975, 230}}: {fit_opt['t_crit_975']:.6f}",
        "",
        "2. PRIMARY TWO-SIDED ESTIMANDS (UNBIASED FIXED-LAMBDA SEs & STUDENT-t INTERVALS)",
        f"  aai slope: beta = {fit_opt['beta'][aai_idx]:+.6f}, SE_unb = {fit_opt['se_unb'][aai_idx]:.6f} (SE_ML = {fit_opt['se_ml'][aai_idx]:.6f}), t = {fit_opt['t_unb'][aai_idx]:+.4f}, two-sided p = {fit_opt['p_unb'][aai_idx]:.6e}",
        f"    95% Conditional Student-t CI (at lambda_ML): [{fit_opt['ci95_low'][aai_idx]:+.6f}, {fit_opt['ci95_high'][aai_idx]:+.6f}]",
        f"  frac_rbh_desert slope: beta = {fit_opt['beta'][rbh_idx]:+.6f}, SE_unb = {fit_opt['se_unb'][rbh_idx]:.6f} (SE_ML = {fit_opt['se_ml'][rbh_idx]:.6f}), t = {fit_opt['t_unb'][rbh_idx]:+.4f}, two-sided p = {fit_opt['p_unb'][rbh_idx]:.6e}",
        f"    95% Conditional Student-t CI (at lambda_ML): [{fit_opt['ci95_low'][rbh_idx]:+.6f}, {fit_opt['ci95_high'][rbh_idx]:+.6f}]",
        f"  log10(length) slope: beta = {fit_opt['beta'][len_idx]:+.6f}, SE_unb = {fit_opt['se_unb'][len_idx]:.6f}, t = {fit_opt['t_unb'][len_idx]:+.4f}, two-sided p = {fit_opt['p_unb'][len_idx]:.6e}",
        f"  log10(N50) slope: beta = {fit_opt['beta'][n50_idx]:+.6f}, SE_unb = {fit_opt['se_unb'][n50_idx]:.6f}, t = {fit_opt['t_unb'][n50_idx]:+.4f}, two-sided p = {fit_opt['p_unb'][n50_idx]:.6e}",
        f"  completeness slope: beta = {fit_opt['beta'][comp_idx]:+.6f}, SE_unb = {fit_opt['se_unb'][comp_idx]:.6f}, t = {fit_opt['t_unb'][comp_idx]:+.4f}, two-sided p = {fit_opt['p_unb'][comp_idx]:.6e}",
        "",
        "3. SENSITIVITY SUMMARY ACROSS COHORTS AND COVARIANCE/DENOMINATOR SPECIFICATIONS",
    ]
    for r in sens_rows:
        rep.append(
            f"  {r['model_id']} (n={r['n']}, lambda_ML={r['lambda_ml']} [{r['lambda_ci95_low']}, {r['lambda_ci95_high']}], R^2_GLS={r['r2_gls']}):"
        )
        rep.append(
            f"    aai: beta={r['coef_aai']} (SE_unb={r['se_unb_aai']}, p={r['p_unb_aai']}), frac_rbh: beta={r['coef_frac_rbh']} (SE_unb={r['se_unb_frac_rbh']}, p={r['p_unb_frac_rbh']})"
        )
    rep.append("=" * 88)

    summary_txt = os.path.join(res_dir, "v3_pgls_ancestry_summary.txt")
    with open(summary_txt, "w") as f:
        f.write("\n".join(rep) + "\n")
    print("\n".join(rep))
    print(f"[02_v3] Completed in {time.time()-t0:.2f}s")


if __name__ == "__main__":
    main()
