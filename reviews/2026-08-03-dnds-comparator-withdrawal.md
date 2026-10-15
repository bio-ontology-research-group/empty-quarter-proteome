# Dark-family dN/dS comparator: consolidated withdrawal decision

Date: 2026-08-03

Status: comparator withdrawn from the manuscript; historical analysis retained
only as a bounded descriptive screen.

## Independent checks and decision

The 2026-07-29 Codex review and a fresh, read-only Claude Code Opus audit both
reproduced the historical table. The dark arm contains 22,767 evaluable
families from 430 representative source samples (median family dN/dS 0.391).
The historical control contains 21,579 families from only 28 source samples
(median 0.286), with 25 source samples shared between arms. Restricting the
dark arm to those 28 samples leaves 570 families with median 0.302, reducing
the raw median gap from 0.105 to 0.016.

Opus identified the deeper, decisive problem. Darkness was assigned as no
UniRef50 hit or less than 30% identity. Homology detectability is itself
related to long-run evolutionary rate, so the arm definition is associated
with the outcome being compared. Matching can repair the 28-versus-430 sample
imbalance but cannot turn this contrast into a valid test of the
random/erroneous-ORF null. The control was also incorrectly called
"annotated"; non-darkness does not establish a functional annotation.

The manuscript therefore no longer reports or interprets the 0.391-versus-
0.286 contrast. It retains only the conditional statement that recurrent,
alignable dark families yielded a median pairwise NG86 estimate of 0.391. The
median is stable under increasing total-family-size floors:

| minimum total members | evaluable families | median family dN/dS |
|---:|---:|---:|
| 3 | 22,767 | 0.391 |
| 5 | 8,769 | 0.418 |
| 10 | 2,479 | 0.449 |
| 20 | 671 | 0.489 |
| 30 | 302 | 0.499 |

This argues against but does not by itself refute an erroneous-ORF
explanation. A decisive replacement would calibrate the complete estimator on
omega=1 codon simulations and on recurrent ORFs called from intergenic and
dinucleotide-preserving shuffled sequence, processed through the identical
pipeline. Those nulls are an optional strengthening analysis, not a condition
for the bounded wording now used.

## Additional audit findings

The legacy pipeline is not approved for new inference because it also:

- appended to member FASTAs without cleaning the output tree, so rerunning
  could silently duplicate sequences;
- aggregated whatever array parts happened to exist without asserting the
  expected part count;
- sampled controls and sequence pairs lexicographically rather than with a
  recorded random design;
- reported medians of pairwise ratios, with silent filtering at zero or
  saturated dS, rather than a calibrated family-level estimator;
- did not record retained alignment length, partial-ORF status, or dropout
  reasons; and
- used pooled tests without sample, MAG, or taxonomic clustering.

The three historical scripts now fail by default and require the explicit
environment acknowledgement `ALLOW_WITHDRAWN_DNDS_PIPELINE=1`. Their exact
pre-guard versions remain recoverable from Git commit
`ecef3a3ea620fb80f51f16d10c22f857e4a4a59a`; hashes are frozen in
`reproducibility/dnds_legacy_provenance.tsv`.

## Manuscript disposition

- Removed "size-matched annotated control 0.29" from B2 and its caption.
- Replaced "refutes" with bounded, conditional wording.
- Removed the unclustered endemic-versus-wide p-value from B2; the medians are
  explicitly descriptive pending clustered inference.
- Replaced the old figure with a dark-only distribution and family-size
  sensitivity panel.
- Added every printed family-size count and median to the executable claims
  verifier.
