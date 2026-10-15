# GitHub Issue #578: Foldseek for pairwise structural alignment

**Repository:** steineggerlab/foldseek  
**URL:** https://github.com/steineggerlab/foldseek/issues/578  
**Author:** mohammadian-reza  
**Date:** March 22, 2026  
**Status:** Open / Unresolved  

## Issue Description

I am working on a template-based docking study where I benchmark pairwise structural alignment methods for protein interfaces among 78 different clusters in 4 different scenarios. I have been testing Foldseek for aligning interface structures against full dimer structures, but I encountered some results that I find difficult to interpret.

With easy-search and `--alignment-type 1`, I sometimes obtain suspicious outputs, including TM-scores greater than 1, alignment lengths that seem inconsistent with the actual structures, and superpositions that do not appear correct.

### Command

```bash
foldseek easy-search Q.pdb T.pdb out.tsv tmp --alignment-type 1 -a --format-output "query,target,qlen,tlen,alnlen,qtmscore,ttmscore,alntmscore,lddt"
```

### Examples Reported

1. **3KX9_P_U_P.pdb vs 2JD6_3_Q_Q.pdb:**
   - Foldseek-TMalign: RMSD 0.197, L_align 13, TM-score 1.042
   - TMalign: RMSD 2.61, L_align 19, TM-score 0.46941

2. **1LWF_A_B_A.pdb vs 1BQN_A_B_A.pdb:**
   - Foldseek-TMalign: RMSD 1.621, L_align 544, TM-score 0.931
   - TMalign: RMSD 1.8, L_align 134, TM-score 0.79824

These differences make me unsure how to interpret the Foldseek TM-align-based output, especially when the reported alignment length spans nearly the full query or the TM-score exceeds 1.

Is Foldseek methodologically suitable for benchmarking pairwise interface alignment in this kind of study, especially for interface-vs-global dimer comparisons? If so, which Foldseek mode or workflow is recommended for this use case? If not, is there any recommendation for how Foldseek could be evaluated fairly in such a benchmark?

## Audit Analysis & Resolution

The reporter observes TM-score > 1.0 in `foldseek easy-search` with `--alignment-type 1`. Specifically in case 1, `L_align = 13` and `TM-score = 1.042`.

Source code inspection of Foldseek (commit `941cd33ff0771cd2e3f144e3293e22a2b87e9fda`) in `src/strucclustutils/structureconvertalis.cpp` lines 895–900 reveals the exact underlying cause:

```cpp
case LocalParameters::OUTFMT_ALNTMSCORE:
    tmres = tmaligner->computeTMscore(targetCaData, &targetCaData[res.dbLen], &targetCaData[res.dbLen+res.dbLen], res.dbLen,
                                      res.qStartPos, res.dbStartPos, tmpBt,
                                      std::min(res.qEndPos - res.qStartPos, res.dbEndPos - res.dbStartPos));
    result.append(SSTR(tmres.tmscore));
    break;
```

Foldseek passes `std::min(res.qEndPos - res.qStartPos, res.dbEndPos - res.dbStartPos)` as the normalization length parameter (`normalizationLen`) for `alntmscore`.

In Foldseek, `res.qStartPos` and `res.qEndPos` are 0-indexed inclusive residue indices (`0 <= qStartPos <= qEndPos < qLen`), formatted to 1-indexed output as `res.qStartPos + 1` and `res.qEndPos + 1`. The span length of aligned residues is `(qEndPos + 1) - (qStartPos + 1) + 1 = qEndPos - qStartPos + 1`.

Because line 898 computes `res.qEndPos - res.qStartPos` WITHOUT `+ 1`, the normalization length is strictly off-by-one (`Lnorm = alnlen - 1`). When the alignment is of length $L = 13$, Foldseek normalizes by $Lnorm = 13 - 1 = 12$. The TM-score sum over 13 residues divided by 12 produces a score scaled up by $\approx 13/12 = 1.0833$. For a high-geometry match, this causes the reported `alntmscore` to exceed 1.0 (e.g., $0.962 \times 13 / 12 = 1.042$).

By contrast, `qtmscore` normalizes by `res.qLen` and `ttmscore` normalizes by `res.dbLen`, both of which are always $\ge$ the number of aligned residues, preventing `qtmscore` and `ttmscore` from ever exceeding 1.0.
