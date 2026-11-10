# Manuscript snapshot

This buildable snapshot belongs to the proteome reproducibility package (`~/Public/software/empty-quarter-proteome/`). Editing checkout: `~/Documents/papers/empty-quarter-proteome/`. Source (`main.tex`, `supplementary_notes.tex`), authors (`authors.tex`), bibliography (`references.bib`), and compiled `main.pdf` match that checkout byte-for-byte (`reproducibility/manuscript_snapshot.sha256`).

- Synchronise from the editing checkout: `bash reproducibility/sync_manuscript_snapshot.sh ~/Documents/papers/empty-quarter-proteome`
- Regenerate the four reframed main figures + panel TSVs and run tests: `bash reproducibility/reframed_20261108/run_reframed.sh`
- Build `main.pdf`: `bash reproducibility/build_manuscript.sh`

The old `paper/` submodule is not part of this package. Reconcile edits with the dedicated checkout before refreshing this snapshot.
