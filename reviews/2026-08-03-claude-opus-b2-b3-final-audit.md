# Claude Code Opus final audit: dark-matter B2--B3

Date: 2026-08-03

Reviewer: Claude Code, model Opus, high effort, read-only plan mode. Codex
independently checked each finding, made the accepted corrections, regenerated
the artifacts, and reran the executable checks.

## Audit sequence

1. Opus rejected the first revision because LaTeX resolved a stale copy of
   `fig_darkmatter_dnds.pdf` under `paper/figures`; the rendered manuscript
   therefore still displayed the withdrawn non-dark comparator. The current
   dark-only figure was synchronized and enlarged to full text width.
2. Opus rejected the next rendered artifact because the B3 occupancy panel
   claimed that "most" dark families were site-restricted, although only 8.5%
   occur at three or fewer sites and 45.6% occur at 30 or more sites. The panel
   now reports the observed occupancy distribution without that claim.
3. The post-fix rendered-artifact audit returned **VERDICT: PASS**, with no
   blocker or major issue.

## PASS evidence

- No withdrawn dN/dS comparator survives in `paper/main.tex`, the active
  claims table, figure code, or the rendered manuscript. The retained raw
  `control` rows remain checksum-bound historical data and are not computed or
  plotted.
- B2 consistently reports 22,767 evaluable recurrent dark orthogroups and
  medians 0.391/0.418/0.449/0.489/0.499 at member floors
  3/5/10/20/30 (n=22,767/8,769/2,479/671/302).
- The dN/dS panel explicitly labels omega=1 as a reference rather than a
  calibrated null, discloses the nine values above the displayed range, and is
  legible at full manuscript width.
- The B3 occupancy panel reports 2.3% single-site and 45.6% in at least 30 of
  60 sites. These percentages are computed from the plotted array rather than
  hardcoded.
- The environmental-structure and occupancy/global-endemism composites each
  have one unambiguous figure label.
- `reproducibility/verify.py --quiet` passes 54 computed claims with zero
  failures; `b2.dark.n` is now computed through the same floor-3 path used by
  the figure.
- The three withdrawn legacy scripts fail closed unless
  `ALLOW_WITHDRAWN_DNDS_PIPELINE=1` is set, and their historical Git/blob and
  data hashes are recorded in `reproducibility/dnds_legacy_provenance.tsv`.

## Final artifact controls

- `paper/main.pdf`: 16 pages; SHA-256
  `836e1c7c50738369cd35dd36e04fcba8c5d5724f2b00915669b4153d7e6987d9`
- `figures/fig_darkmatter_dnds.pdf` and the manuscript copy: SHA-256
  `c8094fa8edd480ffaf00d2cb49aad1584255c8a27ecb92c1c958e76d0db37ecb`
- `figures/fig_darkmatter_endemism_sites.pdf` and the manuscript copy:
  SHA-256
  `e2b2b3595cfc18d5fcc5305dbbb1e86050d8f488a3db8f30ce25f9db9e052d8b`

The superseded report counter and historical unclustered p-value were also
corrected after the PASS audit. Opus's remaining comments concerned cosmetic
label size in other B3 panels and did not affect the audited claims.
