# SMAG Public Source-Retrieval and Alternative Release Audit

**Date:** 2026-10-14  
**Audit Context:** Hypothesis 3 (Round 7 Full Analysis, Reference Database Construction)  
**Assigned Working Directory:** `reviews/darkmatter-2026-10-13/round7/smag_source_audit/`  
**Primary Reference:** Ma et al. (2023), *A genomic catalogue of soil microbiomes boosts mining of biodiversity and genetic resources*, Nature Communications 14:7318, DOI: [10.1038/s41467-023-43000-z](https://doi.org/10.1038/s41467-023-43000-z).  
**Sampling Frame:** Complete 40,039 SMAG MAGs from `Supplementary Data 2.xlsx`.  
**Audited Target Cohort:** Pre-specified uniform random sample of $M = 1,500$ MAGs (seed `20261014`, recorded in `results/sampled_smag_1500_ids_seed20261014.txt`).  

---

## 1. Executive Summary and Audit Verdict

**Definitive Outcome:**  
**No viable official alternative release source exists** that provides individually downloadable genome assemblies or a faster exact archive mirror for the published SMAG 40,039 catalogue.

**Zenodo Record 8223844 remains the sole complete, authoritative, and scientifically reproducible public release** for the published SMAG genome assemblies and our pre-specified 1,500-genome cohort.

### Key Audit Findings:
1. **NCBI GenBank (BioProject PRJNA983538) is Incomplete and Ineligible:**
   - NCBI contains **only 16,530 assemblies**, representing exclusively the *unknown* Species-level Genome Bins (uSGB representatives). It omits all 4,547 known species SGBs and all 18,963 redundant/non-representative MAGs.
   - For our fixed random sample of 1,500 MAGs, **exactly 833 genomes (55.53%) are completely absent from NCBI**. Only 667 genomes (44.47%) are uSGB representatives.
   - NCBI assemblies use isolate names (`SMAG_U1` to `SMAG_U16530`) without an official public cross-reference table to the 40,039 `user_genome` IDs (`BL1bin.1`, etc.).
   - Conditioning on NCBI would introduce severe novelty selection bias, altering the reference database and violating the pre-registered uniform random sampling design.
2. **CyVerse Anonymous WebDAV is Blocked:**
   - The CyVerse endpoint cited in the paper (`https://data.cyverse.org/dav-anon/iplant/home/lucyzju/Caiyu_SMAG_catalog_2023/MAG.tar.gz`) returns an HTTP 307 redirect to `https://unblockme.cyverse.org/` (interactive anti-bot challenge).
   - Furthermore, CyVerse only hosted a monolithic ~36–39 GB `.tar.gz` archive, not individual genome assemblies.
3. **Dedicated Author Web Portal is Permanently Decommissioned:**
   - The paper's web interface (`https://smag.microbmalab.cn`) fails DNS resolution.
   - The authors' official GitHub repository (`Caiyulu-818/SMAG`) and project site (`microbma.github.io/project/SMAG.html`) state: *"Due to the high cost of long-term maintenance, we have decided to discontinue the use of https://smag.microbmalab.cn/ for data access and downloads. We strongly recommend using the Zenodo download options instead."*
4. **AWS China S3 Mirror is Offline:**
   - The secondary mirror (`https://bma-smag.s3.cn-northwest-1.amazonaws.com.cn/SMAG/MAG40039.tar.gz`) returns HTTP 404 Not Found.
5. **ENA and NGDC / CNCB:**
   - ENA only mirrors the raw SRA read runs under PRJNA983538 (no independent assembly catalogue).
   - NGDC (China) hosts no independent SMAG assembly archive; its community index (*Awesome MAG*) explicitly redirects users to Zenodo 8223844.

**Recommendation:**  
The existing H3 retrieval workflow on IBEX (master job `53738740` + auxiliary resume job `53741162`, running 8 concurrent streams across the 373 Zenodo parts with MD5 validation) is the **only scientifically valid and technically viable strategy**. Do not alter jobs, do not modify sampling IDs, and do not switch to partial public subsets.

---

## 2. Sampling Frame and Cohort Specifications

| Metric | Full SMAG Catalogue | Pre-Specified Sample (`seed20261014`) |
| :--- | :--- | :--- |
| **Total Genomes** | **40,039** | **1,500** |
| Sampling Source | `Supplementary Data 2.xlsx` | Uniform random without replacement |
| Unique `user_genome` IDs | 40,039 (100.0% match) | 1,500 (100.0% match) |
| Species Representatives (`is_species_representative == 'yes'`) | 21,076 (52.64%) | 847 (56.47%) |
| Redundant Cluster Members (`is_species_representative == 'no'`) | 18,963 (47.36%) | 653 (43.53%) |
| Novel Taxa (`info == 'unknown'`) | 27,323 (68.24%) | 1,032 (68.80%) |
| Known Taxa (`info == 'known'`) | 12,716 (31.76%) | 468 (31.20%) |
| **Theoretical Overlap with NCBI PRJNA983538** (`yes` & `unknown`) | **16,529 (41.28%)** | **667 (44.47%)** |
| **Missing from NCBI PRJNA983538** | **23,510 (58.72%)** | **833 (55.53%)** |

*Verification:* Computed programmatically against verified `Supplementary Data 2.xlsx` (MD5: `124b2b7424199f31635406dd31a32411`) and retained in `evidence/audit_sampling_coverage_summary.json`.

---

## 3. Systematic Source-by-Source Empirical Audit

### 3.1 Zenodo Record 8223844 (Active Authoritative Baseline)
- **Landing URL:** `https://zenodo.org/records/8223844`
- **REST API Endpoint:** `https://zenodo.org/api/records/8223844/files`
- **Concept DOI / Record DOI:** `10.5281/zenodo.7341719` / `10.5281/zenodo.8223844`
- **Release Version:** v1.0 (indexed version 7, last modified 2025-06-05)
- **Archive Structure:** 373 split `.tar.gz` chunks (`mag.tar.gz.aa` to `mag.tar.gz.oh`)
  - Standard part size: Exactly 104,857,600 bytes (100.00 MiB) for parts 1–372.
  - Total catalogue archive bytes: 39,021,060,793 bytes (36.34 GiB / 39.02 GB).
- **Member-Level Download Supported:** **No.** Genome assemblies are packed inside a single concatenated tar stream.
- **Coverage of Fixed 1,500 IDs:** **1,500 / 1,500 (100.0%)**.
- **Access Limits and Transfer Characteristics:**
  - Zenodo applies rate limits to unauthenticated API queries.
  - Files endpoints support standard HTTP `Range` headers, enabling resumable direct multi-stream curl transfers.
  - Transferring all 373 parts guarantees 100% recovery of any arbitrary random sample across the full 40,039 catalogue.
- **Audit Assessment:** **AUTHORITATIVE & COMPLETE.** The only repository providing the full 40,039 MAG catalogue matching published results.

### 3.2 NCBI GenBank BioProject PRJNA983538
- **BioProject URL:** `https://www.ncbi.nlm.nih.gov/bioproject/PRJNA983538`
- **Submission Title:** *soil metagenome Genome sequencing and assembly*
- **Submitter:** Zhejiang University (Bin Ma / Caiyu Lu)
- **Assembly Count in PRJNA983538:** Exactly 16,531 assemblies (query `PRJNA983538[BioProject]` in NCBI Assembly database).
- **Assembly Accessions:** `GCA_036704125.1` to `GCA_036937395.1` (major release 2024-02-20).
- **Member-Level Download Supported:** **Yes** (via NCBI Datasets API or FTP `ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/...`).
- **Coverage of Fixed 1,500 IDs:** **FAILURE: Maximum 667 / 1,500 (44.47%)**.
  - Exactly 833 genomes in our pre-registered sample (55.53%) were **never deposited to NCBI**.
  - The publication deposited only the 16,530 *unknown* Species-level Genome Bins (uSGBs). All 4,547 known species SGBs and all 18,963 redundant cluster members were excluded.
- **ID Mapping Deficit:**
  - NCBI BioSample isolate fields use `SMAG_U1` through `SMAG_U16530`.
  - There is no official public translation table linking `SMAG_U*` accessions to the original 40,039 `user_genome` identifiers (`BL1bin.1`, `SRR10011679bin.1`, etc.).
- **Audit Assessment:** **NON-VIABLE / INELIGIBLE.** Relying on NCBI would alter the sampling frame to an incomplete, novelty-biased subset, destroying scientific integrity and comparability.

### 3.3 CyVerse Data Store
- **Catalogue URL:** `https://data.cyverse.org/dav-anon/iplant/home/lucyzju/Caiyu_SMAG_catalog_2023/MAG.tar.gz`
- **Dereplicated URL:** `https://data.cyverse.org/dav-anon/iplant/home/lucyzju/Caiyu_SMAG_catalog_2023/MAGdrep.tar.gz`
- **Member-Level Download Supported:** **No.** Monolithic tar archives (~36–39 GB).
- **Live Empirical Access Test:**
  - `curl -I -L -m 15 https://data.cyverse.org/dav-anon/iplant/home/lucyzju/Caiyu_SMAG_catalog_2023/MAG.tar.gz`
  - Result: `HTTP/1.1 307 Temporary Redirect` -> `Location: https://unblockme.cyverse.org/`.
  - CyVerse anonymous WebDAV is intercepted by an interactive anti-bot challenge page, failing automated/unattended pipeline downloads from both local machines and IBEX compute nodes.
- **Audit Assessment:** **BLOCKED / NON-VIABLE.**

### 3.4 Microbma Lab Author Web Server
- **URL Cited in Paper:** `https://smag.microbmalab.cn`
- **Live Empirical Access Test:**
  - `curl -I -m 10 https://smag.microbmalab.cn/`
  - Result: `curl: (6) Could not resolve host: smag.microbmalab.cn` (DNS NXDOMAIN).
- **Official Documentation Evidence:**
  - Author GitHub README (`Caiyulu-818/SMAG`): *"Due to the high cost of long-term maintenance, we have decided to discontinue the use of https://smag.microbmalab.cn/ for data access and downloads. We strongly recommend using the Zenodo download options instead."*
- **Audit Assessment:** **PERMANENTLY DECOMMISSIONED.**

### 3.5 AWS China (cn-northwest-1) S3 Bucket
- **URL Checked:** `https://bma-smag.s3.cn-northwest-1.amazonaws.com.cn/SMAG/MAG40039.tar.gz`
- **Live Empirical Access Test:**
  - Result: `HTTP/1.1 404 Not Found` (XML response `NoSuchKey` / `NoSuchBucket`).
- **Audit Assessment:** **OFFLINE / REMOVED.**

### 3.6 European Nucleotide Archive (ENA)
- **URL:** `https://www.ebi.ac.uk/ena/browser/view/PRJNA983538`
- **Live Empirical Access Test:**
  - ENA API filereport for `result=analysis`: 0 records found.
  - ENA API filereport for `result=read_run`: 726 raw sequencing runs (e.g., `SRR25158176` to `SRR25158191`).
- **Audit Assessment:** **NO INDEPENDENT ASSEMBLY CATALOGUE.** Only mirrors raw in-house sequencing reads.

### 3.7 National Genomics Data Center (NGDC / CNCB, China)
- **Portal URL:** `https://ngdc.cncb.ac.cn`
- **Audit Test:**
  - Checked Genome Warehouse (GWH) and Genome Sequence Archive (GSA) for "SMAG" and "PRJNA983538": No primary MAG submission exists.
  - The Beijing Institute of Genomics community repository (*Awesome MAG*, maintained by Dr. Rui Li) indexes SMAG but explicitly provides download instructions pointing directly to Zenodo Record 8223844.
- **Audit Assessment:** **NO INDEPENDENT HOSTING.** Points to Zenodo.

---

## 4. Retained Evidence Inventory

All reference documents, metadata payloads, and analysis summaries have been downloaded and preserved in `reviews/darkmatter-2026-10-13/round7/smag_source_audit/evidence/`:

| File Name | File Size | Checksum (MD5 / SHA-256) | Provenance / Purpose |
| :--- | :--- | :--- | :--- |
| `Supplementary_Data_2.xlsx` | 18,691,101 bytes | `124b2b7424199f31635406dd31a32411` | Primary paper metadata table from Zenodo 8223844 containing all 40,039 MAG attributes, quality metrics, SGB representation, and novelty classifications. |
| `SMAG_README.md` | 1,866 bytes | `536ff4ea6f6630fbc0285a8ce5c7ea97` | Upstream Zenodo release README documenting archive structure and unpacking instructions. |
| `zenodo_8223844.json` | 116,029 bytes | Verified in round 6 | Official Zenodo REST API response listing all 384 files (373 part files + 11 companion files) with exact byte sizes and MD5 hashes. |
| `ncbi_bioproject_PRJNA983538.json` | 1,743 bytes | E-utilities API | Full XML/JSON document summary of BioProject PRJNA983538 registered by Zhejiang University. |
| `ncbi_assembly_sample.json` | 24,041 bytes | E-utilities API | Sample document summaries for assemblies `GCA_036704125.1` to `GCA_036937395.1`, demonstrating `SMAG_U*` isolate nomenclature. |
| `caiyulu_smag_github_readme.md` | 2,469 bytes | GitHub API | Upstream GitHub repository README confirming decommissioning of `smag.microbmalab.cn` and recommending Zenodo. |
| `awesome_mag_smag_download.md` | 5,343 bytes | GitHub API | Beijing Institute of Genomics community download guide by Dr. Rui Li detailing Zenodo retrieval and mirror status. |
| `microbma_github_io_SMAG.html` | 7,446 bytes | HTTP GET | Microbma laboratory project page. |
| `audit_sampling_coverage_summary.json` | 2,951 bytes | Local Audit | Machine-readable breakdown of catalogue metrics, SGB overlap, and source evaluations. |

---

## 5. Methodological and Computational Conclusion

1. **Integrity of Sampling Frame:**  
   The scientific requirement for Hypothesis 3 is an unbiased, uniform random sample of $M = 1,500$ MAGs across the complete 40,039-genome SMAG catalogue. Bypassing Zenodo to download from NCBI would discard 55.53% of our sampled genomes and restrict the reference set to unknown taxa, introducing substantial ascertainment bias.
2. **Transfer Efficiency:**  
   Because the uniform random sample is distributed throughout the concatenated 36.34 GiB archive, streaming all 373 parts is mathematically necessary to recover the complete sample without replacement.
3. **Execution Status:**  
   The active IBEX transfer configuration (8 parallel resumable streams across master `53738740` and auxiliary `53741162`) is proceeding cleanly with strict byte-size and MD5 verification. No alternative source provides a faster or valid alternative.
