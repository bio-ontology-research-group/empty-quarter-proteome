# RQ-A1 — How much of the EQ proteome is dark, and is it organised into families/folds?

**Status:** ✅ answered

## Answer
**A historical 51.9 M-protein snapshot contains 2.84 M UniRef50-dark representatives and 334,810 recurrent sequence-dark families.** Its 68,808 pre-correction structural candidates are not fold calls. Under the current full-coordinate confidence/length gate, 608/1,875 are qTM dual no-hits and 97/1,875 are alignment-TM dual no-hits, against error floors of 101/3,013 and 68/3,013 in known-fold controls.

**Evidence / reproduce:** `reproducibility/verify.py` (claims `b1.*`, PROVENANCE); source counts in the workbench darkmatter pipeline.

**References:** pavlopoulos2023darkmatter; vanni2022agnostos; durairaj2023folds; yeo2025afesm; rodriguezdelrio2024fesnov (keys in `manuscript/references.bib`).
