# Boltz-2 ligand validation — confident-function exemplars

Boltz-2 (single-sequence mode, `msa: empty`; ORIX H200) on the two confident-function families with
an unambiguous cognate ligand class: the dark **lectin** (58PRr3, structural homolog *Photorhabdus*
lectin 6fhy) against a carbohydrate panel, and the dark **ABC transporter** (S11PRr1, NosDFY) against
ATP. Complex + binding-affinity prediction.

| target | ligand | confidence | ligand_iptm | complex_pLDDT | binding_prob | affinity (log, lower=tighter) |
|---|---|---|---|---|---|---|
| lectin 58PRr3 | **GlcNAc (NAG)** | 0.65 | **0.55** | 0.68 | **0.48–0.59** | 1.25 |
| lectin 58PRr3 | mannose (MAN) | 0.67 | 0.47 | 0.72 | 0.26 | 1.54 |
| lectin 58PRr3 | galactose (GAL) | 0.64 | 0.41 | 0.70 | 0.28 | 2.22 |
| ABC S11PRr1 | ATP | 0.55 | 0.34 | 0.61 | 0.29 | 0.20 |

## Reading (honest)
- **The dark lectin most confidently binds N-acetylglucosamine (GlcNAc):** best ligand-interface
  confidence (ligand_iptm 0.55) and highest binder probability among the sugars — concordant with its
  *Photorhabdus*-lectin fold. This is a concrete, testable function hypothesis for a protein that is
  **both sequence-dark (UniRef50-novel) and structure-dark-adjacent** — reached purely by ESMFold →
  FoldSeek → Boltz-2.
- **Confidence is modest, not decisive** (iptm ~0.4–0.55, binder prob <0.6): a hypothesis-strengthening
  signal, not proof. The ABC-transporter × ATP prediction is weak/inconclusive.
- **Interface-localized selection could not be confirmed:** of 34 nominal positive sites in 58PRr3,
  0 fall within the predicted 3-residue NAG-contact site (≈chance). But the lectin family has only 4
  members, so its per-site FUBAR is underpowered and the 34 sites are unreliable. This exposes the
  core limitation of the whole structure×selection test: **the function-assignable families are shallow
  (3–4 members), exactly the regime where per-site selection has no power** — while the deep families
  (>=8 members, real per-site signal) mostly lack a confident structural homolog. The two requirements
  (assignable function vs powered per-site test) anti-correlate in this dataset.

## Bottom line
ESMFold→FoldSeek→Boltz-2 turns a sequence-dark, positively-selected protein into a falsifiable
function + ligand hypothesis (dark lectin → GlcNAc binder). It does NOT, in this dataset, let us place
the diversifying residues on the binding interface, because functionally-decodable dark families are
too small for residue-level selection inference.
