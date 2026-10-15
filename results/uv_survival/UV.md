# UV-survival / DNA-damage-repair census (15457/15469 MAGs carry >=1 marker)

> **SUPERSEDED (2026-08-13).** These tables are the HISTORICAL 15,469-single-sample-MAG basis and the
> per-MAG compartment test is pseudoreplicated. Current numbers (5,438 catalogue, rep_ko_v2,
> site-paired design): workbench `M/06_metabolism/uv_survival_clean/` + `notes/CORRECTIONS.md`
> 2026-08-13 entry. Tables kept unmodified because `reproducibility/verify.py` (`a2.*`) still binds
> to them; rebinding is pending PI coordination.

EggNOG KO census (job 47736827) over a 28-KO UV-damage-response panel. The Empty Quarter genome is **UV/desiccation-hardened across the board** — repair is near-universal, not a specialist trait.

## Prevalence by functional group (MAGs with >=1 marker)

| group | MAGs | % catalog |
|---|---|---|
| base_excision_repair | 15398 | 99.5% |
| nucleotide_excision_repair | 15382 | 99.4% |
| SOS_translesion | 14570 | 94.2% |
| dna_protection_oxidative | 14215 | 91.9% |
| dsb_nhej | 12122 | 78.4% |
| photoprotective_pigment | 9050 | 58.5% |
| mismatch_repair | 7548 | 48.8% |
| photoreactivation | 2710 | 17.5% |

## Compartment test — does repair track sun exposure?

Surface = UV-exposed top layer; Deep = shielded subsurface. Mean distinct UV-repair KOs and UV-specific marker prevalence per compartment:

| compartment | MAGs | mean UV KOs | % photolyase | % NHEJ | % carotenoid |
|---|---|---|---|---|---|
| Surface | 2241 | 16.01 | 22.4% | 73.0% | 67.6% |
| Deep | 3889 | 15.20 | 13.1% | 81.1% | 55.1% |
| PlantRoot | 9163 | 16.29 | 17.8% | 78.6% | 57.6% |

Surface vs Deep, UV-repair load (one-sided MWU, Surface>Deep): U p=1.12e-14. Kruskal-Wallis across compartments: p=8.42e-49.

## Phylum breakdown — UV-diagnostic markers

| marker | phylum | MAGs |
|---|---|---|
| photolyase | Actinomycetota | 1085 |
| photolyase | Pseudomonadota | 968 |
| photolyase | Bacteroidota | 156 |
| photolyase | Halobacteriota | 117 |
| NHEJ | Actinomycetota | 6305 |
| NHEJ | Pseudomonadota | 2868 |
| NHEJ | Bacteroidota | 445 |
| NHEJ | Gemmatimonadota | 419 |
| carotenoid | Actinomycetota | 3562 |
| carotenoid | Pseudomonadota | 2914 |
| carotenoid | Bacteroidota | 549 |
| carotenoid | Myxococcota | 338 |

## Reading

- **Excision/recombinational repair is the floor, not the ceiling**: NER (uvrABCD) and BER (mutM/nth/mutY) are in ~99% of MAGs — every desert genome is pre-loaded for DNA-damage repair.
- **The signal is mechanistically coherent — defenses partition by what they need.** The LIGHT-dependent UV defenses concentrate where the sun hits: photolyase Surface 22% vs Deep 13%, carotenoid 68% vs 55%, and total UV-repair load is higher at the surface (MWU p<1e-13). **NHEJ runs the other way** (Deep 81% > Surface 73%): double-strand-break repair is driven by desiccation/dormancy, not light, so it does not follow the sun gradient. Two distinct stressors, two distinct spatial signatures.
- **NHEJ (ku+ligD) is unusually common** (~78%) for environmental bacteria — the hallmark of genomes built to rejoin breaks accumulated while desiccated and metabolically dormant.
- Caveat: KO-based (EggNOG) — spore photoproduct lyase (splB), 6-4 photolyase and MAA/scytonemin screening are not cleanly KO-resolvable and would need an HMM follow-up.
