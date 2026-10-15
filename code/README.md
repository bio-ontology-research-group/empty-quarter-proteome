# code/ — analysis scripts and figure generators

## Layout
- **`figures/`** — LIVE figure generators, current basis (5,438 reps / 106,042 families / R232).
  One file per figure. Each script asserts every plotted number against the registered
  (`notes/CATALOGUE_NUMBERS.md`) values before drawing, reads workbench `*_clean` results (or repo
  `results/`), and writes to repo-root `figures/` (design options to `figures/options/`).
  `paper/figures/` is the compile-time mirror: copy new or updated PDFs there;
  `manuscript/figures/` syncs at ship time via `reproducibility/sync_manuscript_snapshot.sh`.
- everything else at this level — PI-era analysis scripts, kept as provenance. Never edited here.

## Live generators (code/figures/)
| script | output | section |
|---|---|---|
| make_bgc_amr_figures.py | fig_bgc_product, fig_amr_class | 2.4 |
| make_darkfam_env_figures.py | fig_darkmatter_pcoa, fig_darkmatter_permanova | 2.5 |
| make_specificity_figures.py | fig_specificity_controls | 3.1 |
| make_structural_calibration_figure.py | fig_structural_calibration (+2 options) | 1.3 |
| make_darkfam_reality_figure.py | fig_darkfam_reality (+1 option) | 1.4 |
| (moved) make_lifeonair_figure.py | now `lifeonair/code/` (P3) | n/a |
| make_refbias_figure.py | fig_refbias_ceiling | 2.1 |
| make_stress_figure.py | fig_stress_modes | 2.3 |
| make_darkfilter_option_figure.py | options/fig_darkfam_filter_opt | 2.5 option |

## Dead-basis figure scripts — do NOT rerun (outputs superseded)
- `make_global_specificity_figure.py`, `make_darkmatter_figures.py` — 23,382-orthogroup basis
  (already flagged in `notes/PLAN.md` and `notes/PAPER_STRUCTURE.md`).
- `make_darkmatter_dnds_figure.py` — its output prints 22,767 orthogroups / median 0.391, the
  pre-2026-08-14 basis; superseded by `figures/make_darkfam_reality_figure.py`.
- `make_picrust2_bias_figure.py` — 15,469-MAG / 436-sample snapshot; superseded by
  `figures/make_refbias_figure.py`.

Palette across all paper figures: dark/focal set `#4477aa`, control/comparator `#cc6677`
(CVD-validated; see the docstring of `figures/make_specificity_figures.py`), grey `#999999` for
de-emphasis. Colour follows the entity, never the panel.
