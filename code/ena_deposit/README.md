# ENA deposit scripts (RAK-WGS = PRJEB125153)

Decisions, status and the exact command sequence live in `notes/ENA_DEPOSIT.md`. Work dir (not in git):
`/ibex/scratch/projects/c2014/alelopezv/ena_wgs/`. Credential: `~/.webin_pass` only. Every script is idempotent and
re-runnable; sbatch scripts take an explicit `MODE` (no default) where a step is irreversible.

| phase | script | does |
|---|---|---|
| shared | `ena_common.py` | credential, test/prod service URLs, curl, drop-box fetch, XML submission + receipt parsing |
| reads | `00_build_runs_table.py` | `wgs_runs.tsv` from the raw batches (instrument, kit, insert size, run alias); `--symlinks` builds the FASTQ farm |
| reads | `01_upload_ftp.sbatch` | md5 + FTP upload per run, resumable |
| study | `02_register_study.py` | registers a study by XML (all fields explicit) |
| samples | `03_register_biosamples.py` | `register` missing soil BioSamples (ERC000022) / `cleanup` MODIFY pass; writes `wgs_biosample_map.tsv` |
| reads | `04_generate_manifests.py` | one Webin-CLI reads manifest per run |
| reads | `05_webin_cli.sbatch` | Webin-CLI reads array, `MODE=validate|submit` |
| study | `06_link_umbrella.py` | attaches child studies to an existing umbrella (MODIFY) |
| reads | `07_insert_size.sbatch` | insert size per sample from the binning BAMs |
| reads | `08_submit_runs_xml.py` | experiment + run XML per run for pre-uploaded files |
| reads | `09_fold_accessions.py` | ERX/ERR from ENA's live report into `wgs_runs.tsv` |
| assemblies | `10_build_assembly_table.py` | `--prepare` sample lists · build `wgs_assemblies.tsv` + genome manifests · `--fold` ERZ accessions |
| assemblies | `11_gzip_assemblies.sbatch` | gzip each contig file into the work dir + contig count / total bp |
| assemblies | `12_coverage_bam.sbatch` | `samtools coverage` for assemblies without a binning depth table |
| assemblies, bins, MAGs | `13_webin_cli_genome.sbatch` | Webin-CLI genome array, `MODE=validate|submit`; `MANIFEST_DIR` selects the object set |
| bins | `14_build_bin_table.py` | `wgs_bins.tsv`: the 12,434 HQ+MQ bins with quality, GTDB lineage, rRNA/tRNA, coverage, paths (run under Slurm) |
| bins | `15_resolve_bin_taxa.py` | NCBI taxon per bin from the GTDB lineage: species → "Genus sp." → rank names (GTDB name, then the NCBI names its genomes carry in the r232 metadata) → uncultured; every name checked against ENA's taxonomy service (cached); `ncbi_level` records the step |

| bins, MAGs | `16_register_bin_samples.py` | derived BioSamples by XML: `--tier binned` (ERC000050, derived from the soil sample) / `--tier mag` (ERC000047, derived from the binned sample); `--dryrun` writes XML only; `--service test|prod` |
| bins | `17_gzip_bins.sbatch` | gzip the 12,434 bin FASTAs into the work dir (25 tasks) |
| bins, MAGs | `18_build_bin_manifests.py` | one Webin-CLI genome manifest per bin for a tier (chromosome list for a single-contig MAG); `--fold` puts the ERZ (ENA live list) and, for MAGs, the GCA + sequence accessions (saved accession emails) into `wgs_bins.tsv` |

Tables: `wgs_runs.tsv` (554 runs), `wgs_biosample_map.tsv` (551 samples), `wgs_assemblies.tsv` (549 assemblies),
`wgs_bins.tsv` (12,434 bins; the representatives with contamination < 10 form the MAG tier).
