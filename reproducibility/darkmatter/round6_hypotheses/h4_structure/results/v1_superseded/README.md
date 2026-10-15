# SUPERSEDED: V1 Structure Pilot (2026-10-13)

This initial pilot (Slurm Job 53700367) is superseded per Root review (ROOT_AUDIT_NOTES.md).

### Reasons for Supersedure:
1. **Whole-Chain Overlap vs Region TM**: V1 mapped whole-chain hit spans onto predicted segments rather than cropping confident regions and re-searching them directly with Foldseek. Whole-chain hit overlap does not provide region-specific TM scores.
2. **Biological / Mechanistic Overreach**: V1 inappropriately labeled confidence-thresholded segments "bona fide novel domains", claimed disordered-tail and linker-reorientation mechanisms without direct testing, and formulated an ungrounded "Confirmed Novel Structural Domains" go criterion.
3. **Control Calibration Limits**: V1 labeled controls "homology-excluded" when only source 4-letter PDB identifiers were removed. Sequence-homology exclusion (e.g. sequence identity thresholding) was not implemented in the pilot.
4. **Equality Inferences**: V1 did not properly distinguish lack of rejection from equivalence.

V2 replaces this with actual outcome-independent cropped-region searches, explicit sequence-identity exclusion sensitivities, and a reference-conditional recovery estimand without mechanistic overreach.
