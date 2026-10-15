# Independent source audit

Flash syntheses and original manifests are preserved as model output, not authoritative verification. Use verified_sources.tsv and verified/ first. DOI verification establishes source identity, not support for a scientific claim.

## Excluded mislabeled downloads

- `Vanni_2022_AGNOSTOS_fulltext.xml`: embedded title/DOI is an unrelated paper; retained in quarantine, excluded from scientific evidence.
- `RodriguezdelRio_2024_FESNov_fulltext.xml`: embedded title/DOI is an unrelated paper; retained in quarantine, excluded from scientific evidence.
- `Coelho_2022_GMGC_Nature_preprint.pdf`: embedded title/DOI is an unrelated paper; retained in quarantine, excluded from scientific evidence.
- `Coelho_2022_GMGC_fulltext.xml`: embedded title/DOI is an unrelated paper; retained in quarantine, excluded from scientific evidence.
- `Karatzas_2026_proteinfamilies_GigaScience_preprint.pdf`: embedded title/DOI is an unrelated paper; retained in quarantine, excluded from scientific evidence.
- `Schmidt_2024_SPIRE_NAR_preprint.pdf`: embedded title/DOI is an unrelated paper; retained in quarantine, excluded from scientific evidence.

## Metadata and interpretation corrections

- Orr and Bateman, DOI 10.1093/bioadv/vbag160, was published in 2026, not 2024 (Europe PMC XML publication date 10 June 2026). The original Flash filename/year is wrong.
- Flash structural synthesis misstates the 4.29-million cluster census as derived from the 5,438 refined representatives. It actually uses membership in the 37,949 unrefined bins. It also treats local structural matches as proof of domain architecture; that inference requires reviewer scrutiny.
- The structural PDF integrity scan flagged the TED phrase “should likely be disregarded.” Page 10 shows this refers to low-confidence predicted folds, not instructions to a reviewer. Other flags examined in this batch are tiny figure labels or publisher/preprint margin text. No targeted instruction identified in this batch.
- Some source PDFs are preprints rather than the final DOI-associated article. Preserve version identity and use DOI-verified journal XML when available.

Full-text identity audit is ongoing while the remaining searches complete.

Additional unrelated full-text downloads quarantined after embedded-DOI checks:
- `Shulgina_2021_Codetta_alternative_genetic_codes_eLife.xml`
- `Borges_2022_stop_codon_recoding_phages_Nature.xml`
- `Ivanova_2014_stop_codon_reassignments_metagenomics_Science.xml`
- `Tang_2023_pangenomic_assessment_MAGs_BriefBioinform.xml`

The Penn GUIDANCE XML has DOI 10.1093/nar/gkq443 (the web-server article), whereas the Flash manifest requests 10.1093/molbev/msq066 (the method paper). Both are relevant but must not be cited interchangeably.
The existing bibliography DOI checker reports one title warning only: Edgell, Arabian Deserts, omits its subtitle in DOI metadata; this is not evidence of a wrong book.

Corrected two unresolvable Flash DOI leads using publisher pages: Borges is Nature Microbiology 10.1038/s41564-022-01128-6; Rinke is ISME Journal 10.1038/s41396-018-0282-y, with correction 10.1038/s41396-019-0556-z. Original incorrect identifiers remain in the raw search manifests for audit.

The supposed Tang 2023 pangenomic-assessment reference is unverified: DOI 10.1093/bib/bbac484 resolves to an unrelated DNA-storage error-correction article. Exclude this DOI and its attributed claims; no correct replacement has been established.

The Austin/SCRuB PDF is also mislabeled (its first-page title concerns YAP/collagen in uterine leiomyosarcoma). Quarantined; the actual SCRuB paper still needs a verified full-text copy.

Further unrelated downloads excluded after first-page title or XML DOI inspection:
- `04_lineage_accessory/Hellewell_2024_CELEBRIMBOR_core_accessory_MAGs_Bioinformatics.pdf`
- `04_lineage_accessory/Majernik_2025_PanSweep_misassembly_pangenome_mSphere.pdf`
- `05_ecology_design/02_Liu_2022_SqueegeeContamination.pdf`
- `05_ecology_design/10_Monzon_2022_AntiFamAlphaFold.pdf`
- `05_ecology_design/07_SantosJunior_2024_AMPSpherePeptides.xml`
- `05_ecology_design/10_Monzon_2022_AntiFamAlphaFold.xml`

Correct full-text XML remains available for CELEBRIMBOR, PanSweep, Squeegee and SCRuB. The Monzon DOI in the raw manifest is itself incorrect (vbac083 is Metadensity); exclude that DOI as evidence about spurious proteins.

## Round-2 search qualifications (13 October 2026)

The root checked embedded DOI identities for the 27 newly available round-2 XML
papers; these matched the requested manifest DOIs. This verifies identity only.
The Flash synthesis remains an unreviewed set of leads: its causal claims about
all metagenomic Prodigal predictions and exact shadow-ORF exclusion thresholds
must be checked against the cited papers before use. A high-confidence model is
not by itself a validation that an ORF is translated.

Flash overstates journal policy in its final 'Key Author Mandates'. The saved
Nature Portfolio reporting-standards page requires access to central custom code
and a Code Availability statement; it describes DOI repository deposition and
open-source licensing with different recommendation/requirement wording.
Predicted ESMFold models must not be assumed eligible for mandatory wwPDB or
AlphaFold Database deposition. Nature Communications' chemical-characterisation
page encourages predicted-structure sharing in appropriate repositories:
https://www.nature.com/ncomms/submit/chemical-characterisation . The author-facing
submission checklist must use the applicable primary policy, not Flash's summary.
No repository deposition or publication was performed.

## Round 3
Root verified all 32 retained round-3 XMLs against the primary front-matter DOI (`round3_identity_checks.json`). Read the Mugal, AntiFam, Orr, tRNAscan-SE2, skani and FastTree abstracts and relevant method passages before citation. Publisher HTML is retained separately. Flash recommendations are search guidance, not independently established facts or journal requirements.

### Important corrections to Flash round-3 synthesis
The downloaded source identities passed; the synthesis contains errors and is NOT an analysis result. Admission requires three samples and TWO SITES, not two expeditions. The review's proposed Foldseek prefilter mechanism remains untested. Culture percentages quoted from round2 are superseded by explicit base alignments (516/1,105 dark-unshared; post hoc20-pair subset496/940). AntiFam screening has now completed with v8.0, not v7.0; a negative AntiFam result cannot guarantee a real protein. Gene-caller effects are supported, but the entire sharing difference has not been causally assigned to them. “Radiation” is not a tested evolutionary claim. Bibliography inclusion follows source inspection, not this synthesis.

### Final primary-source follow-up
Nine additional free PMC HTML articles passed exact `citation_doi` checks, and the retained Yeo v1 PDF passed the hidden-text scan and first-page DOI/title check. The download ledger is `papers/additional_free/retrievals.json`. These complement68canonical DOI-verified XML fulltexts. A current search of the exact Yeo title found the bioRxiv v2 page (posted16March2026; https://www.biorxiv.org/content/10.1101/2025.04.23.650224v2) and another preprint, but no verified journal version. The bibliography therefore retains the preprint designation; absence of a search hit is not proof that no journal version exists.
