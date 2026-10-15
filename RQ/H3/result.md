# RQ-H3 — Is there taxonomic novelty (new deep lineages)?

**Status:** ✅ answered (bounded at high ranks; lower ranks open)

## Answer
**No current phylum- or class-gap candidate, but substantial lower-rank novelty.** GTDB-Tk 2.7.2/R232 assigns all 5,438 ANI95 HQ/MQ representatives through class across 51 named phyla. First assignment gaps occur for 11 representatives at order, 142 at family, 1,343 at genus and 3,755 at species. Ten order gaps entered de novo phylogenies and none advanced to proposal; one 46.49%-MSA candidate remains unadjudicated and in recovery.

**Evidence / reproduce:** `reproducibility/verify.py` (`scale.*`); checksum-pinned R232 census and phylogeny artifacts under `reproducibility/`.

**References:** parks2018gtdb; chaumeil2022gtdbtk (keys in `manuscript/references.bib`).
