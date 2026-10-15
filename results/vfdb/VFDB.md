# Virulence-factor census (VFDB setB) — dual-use framing

mmseqs of 5438 dRep-rep proteomes (20.5M proteins) vs VFDB setB (29,804 VF proteins). 541,508 proteins have a VF hit (id>=40, cov>=70, e<=1e-10). **This count is inflated by dual-use genes** — most 'virulence factors' in VFDB are secretion/uptake/surface functions ubiquitous in soil bacteria. Identity tiers below show how the signal collapses with stringency.

MAGs with >=1 VF hit: id>=40 5435/5438, id>=70 3769, id>=90 202.

## By VFDB functional category (proteins)

| category | id>=40 | id>=50 | id>=70 | id>=90 | MAGs(>=40) |
|---|---|---|---|---|---|
| Immune modulation | 140845 | 50582 | 3872 | 205 | 5406 |
| Nutritional/Metabolic factor | 140169 | 45892 | 3539 | 257 | 5304 |
| Adherence | 56257 | 24988 | 4751 | 256 | 5375 |
| Regulation | 50183 | 15405 | 2646 | 256 | 4776 |
| Effector delivery system | 46310 | 16352 | 1936 | 202 | 5244 |
| Motility | 43508 | 16617 | 2025 | 251 | 4819 |
| Stress survival | 27679 | 17476 | 1347 | 57 | 5180 |
| Biofilm | 12238 | 3445 | 439 | 87 | 2995 |
| Exoenzyme | 8048 | 4581 | 145 | 3 | 4597 |
| Exotoxin | 6260 | 1447 | 40 | 16 | 3570 |
| Others | 3773 | 2406 | 631 | 45 | 2192 |
| Post-translational modification | 2756 | 476 | 0 | 0 | 1866 |
| Antimicrobial activity/Competitive advantage | 2724 | 993 | 25 | 2 | 1358 |
| Invasion | 758 | 59 | 9 | 6 | 672 |

## Genuine toxins/exotoxins (id>=50): 3020 proteins

The pathogenesis-specific slice (exotoxins, hemolysins, etc.) — the part that actually speaks to 'toxins', not competition/nutrition:

| MAG | %id | VF class | phylum |
|---|---|---|---|
| F35PRr1__self_SemiBin_714 | 100.0 | Hemolysin III |  |
| 30PRr3__single_SemiBin_4_sub | 100.0 | Hemolysin III |  |
| T1PRr2__self_SemiBin_31 | 100.0 | Hemolysin III |  |
| 64PRr2__self_SemiBin_167 | 100.0 | Diphtheria toxin repressor DtxR |  |
| F4Sr2__single_SemiBin_30_sub | 99.8 | Phospholipase D |  |
| 30PRr3__single_SemiBin_4_sub | 99.5 | Nhe |  |
| 30PRr3__single_SemiBin_4_sub | 99.0 | Hemolysin II |  |
| V53PRr2__self_SemiBin_704 | 99.0 | Hemolysin III |  |
| 30PRr3__single_SemiBin_4_sub | 99.0 | CytK |  |
| 30PRr3__single_SemiBin_4_sub | 98.9 | Nhe |  |
| 30PRr3__single_SemiBin_4_sub | 98.6 | ALO |  |
| V53PRr2__self_SemiBin_82 | 98.1 | Hemolysin III |  |
| V35PRr2__self_SemiBin_99 | 97.6 | Hemolysin III |  |
| 30PRr1__single_SemiBin_73 | 97.5 | Brevilysin (BVL) |  |
| F28Dr1__self_SemiBin_285 | 96.7 | Hemolysin III |  |
