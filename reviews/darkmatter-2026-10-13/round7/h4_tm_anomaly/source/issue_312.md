# GitHub Issue #312: `alntmscore` output is wrong

**Repository:** steineggerlab/foldseek  
**URL:** https://github.com/steineggerlab/foldseek/issues/312  
**Author:** ekiefl  
**Participants:** Martin Steinegger (martin-steinegger)  
**Date:** July 2024  
**Status:** Closed / Resolved in upstream format parser  

## Issue Description

### Expected Behavior
I expect `alntmscore` to be different from `ttmscore`.

### Current Behavior
All `alntmscore` values are equal to `ttmscore`.

### Steps to Reproduce
```bash
rm -rf tmp
foldseek easy-search ./structures/ ./structures/ output.txt tmp -a --exhaustive-search --format-mode 4 --format-output query,target,qlen,tlen,alnlen,qtmscore,ttmscore,alntmscore,cigar,qseq,tseq
```

The user tested six combinations of the `--format-output` specification and found that when `qtmscore,ttmscore,alntmscore` was requested, `alntmscore` was assigned the value of `ttmscore`. When `ttmscore,qtmscore,alntmscore` was requested, `alntmscore` was assigned the value of `qtmscore`.

### Author (Martin Steinegger) Response

> "How do you compute the alntmscore? I checked the code we do use `std::min(static_cast<unsigned int>(res.backtrace.size()), std::min(res.dbLen, res.qLen)))` as normalization factor. Normalizing by `res.backtrace.size()` might be better, maybe we should change...
> 
> Okay I think I know what is going on. If you print out `qtmscore,ttmscore,alntmscore` then the `alntmscore = ttmscore`.
> If you print `ttmscore,qtmscore,alntmscore` then `alntmscore = qtmscore`;
> However, if you print `alntmscore,qtmscore,ttmscore` then it should work. This is a bug, which I will fix..."

## Relevance to H4 Audit

1. **Format Order in H4 Pipeline:**
   The H4 pipeline uses the format:
   ```
   "query,target,alntmscore,qtmscore,ttmscore,prob,lddt,evalue,bits,fident,pident,qcov,tcov,qlen,tlen,alnlen,qstart,qend,tstart,tend,cigar"
   ```
   Because `alntmscore` appears BEFORE `qtmscore` and `ttmscore`, it does **not** suffer from the field-overwrite bug described in Issue #312.
2. **Normalization Discussion:**
   Martin Steinegger's comment confirms that `alntmscore` was intended to normalize by alignment length (`res.backtrace.size()` or alignment span). However, in `structureconvertalis.cpp`, the implementation passed `std::min(res.qEndPos - res.qStartPos, res.dbEndPos - res.dbStartPos)` which was missing the `+ 1` required for 0-indexed inclusive coordinates. This connects directly to the off-by-one denominator bug observed in Issue #578 and H4.
