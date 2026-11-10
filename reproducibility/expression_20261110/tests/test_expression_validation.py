#!/usr/bin/env python3
"""
Negative-path and regression unit test suite for Analysis 4 (Karaoz et al. 2026
external metaproteomic expression validation).

Tests all mandatory controls from:
- ROOT-PSM-SCORES.txt (upstream MS-GF:PepQValue search_engine_score[3] and MS-GF:QValue
  search_engine_score[6] metadata verification; no reconstructed peptide FDR; discordance
  between opt_global PSM q and engine QValue rejected; NaN/q>1/missing scores rejected;
  opt_global_valid != VALID excluded; no silent fallback to 11 in Pyrodigal translation_table)
- ROOT-TRANSLATION-CORRECTION.txt (table 4 vs table 11 internal TGA stop regression;
  canonical sequence normalization preserving internal stops and invalid characters)
- ROOT-STAGE1-REVIEW.txt (alias collisions, unmapped IDs, short-domain rejection,
  coordinate/evalue bounds)
- ROOT-DOWNLOAD-REVIEW.txt (mzTab >100KB HTML rejection, trailing newline check,
  PSM/PRT row count match, PSH required columns, SHA-256 tamper detection)
- ROOT-ANALYSIS-CONTRACT.txt Rule 7 (wrong IDs, same counts with replaced sequence hash,
  decoy hits, shared peptides across MAGs, I/L isobaric ambiguity against whole background,
  duplicate biological replicates / n=15 depth pseudoreplication / 2016 DNA mismatch,
  missing PSM/Pep q-values, and non-estimable absent search-database entries).
"""
import csv
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest

_THIS_DIR = Path(__file__).resolve().parent
_PACKAGE_DIR = _THIS_DIR.parent
sys.path.insert(0, str(_PACKAGE_DIR / "code"))
sys.path.insert(0, str(_THIS_DIR))

from download_massive_mztab import validate_mztab_file
from expression_validation_core import (
    aggregate_unique_spectra_for_entity,
    benjamini_hochberg_fdr,
    build_unique_alias_map,
    canonicalize_and_validate_protein_seq,
    classify_peptide_specificity_against_background,
    classify_primary_warming_estimability,
    classify_protein_and_family_detection,
    compute_paired_warming_effect_n3_plots,
    evaluate_homology_alignment,
    filter_and_audit_mztab_psms,
    lookup_strict,
    parse_and_validate_flank_evidence_rows,
    parse_mztab_score_metadata,
    stream_multipattern_peptide_search_records,
    translate_pyrodigal_gene_consistent,
    validate_karaoz_experimental_design,
    verify_sequence_hash_join,
)

REAL_MTD_LINES = [
    "MTD\tmzTab-version\t1.0.0",
    "MTD\tpsm_search_engine_score[1]\t[MS, MS:1002049, MS-GF:RawScore, ]",
    "MTD\tpsm_search_engine_score[2]\t[MS, MS:1002050, MS-GF:DeNovoScore, ]",
    "MTD\tpsm_search_engine_score[3]\t[MS, MS:1002055, MS-GF:PepQValue, ]",
    "MTD\tpsm_search_engine_score[4]\t[MS, MS:1002053, MS-GF:EValue, ]",
    "MTD\tpsm_search_engine_score[5]\t[MS, MS:1002052, MS-GF:SpecEValue, ]",
    "MTD\tpsm_search_engine_score[6]\t[MS, MS:1002054, MS-GF:QValue, ]",
]


class FakeTable4Gene:
    """
    Simulates a gene called under genetic code 4 (where TGA encodes Trp 'W',
    whereas under genetic code 11 TGA is a stop codon '*').
    """
    def __init__(self, dna_seq: str, called_table: int = 4):
        self.dna_seq = dna_seq.upper()
        self.translation_table = called_table

    def translate(self, translation_table=None) -> str:
        tbl = self.translation_table if translation_table is None else int(translation_table)
        codons = [self.dna_seq[i:i + 3] for i in range(0, len(self.dna_seq) - 2, 3)]
        aa = []
        for c in codons:
            if c in ("TAA", "TAG"):
                aa.append("*")
            elif c == "TGA":
                aa.append("W" if tbl == 4 else "*")
            elif c == "ATG":
                aa.append("M")
            else:
                aa.append("A")
        return "".join(aa)


class TestTranslationAndCanonicalizationCorrection(unittest.TestCase):
    """Tests for ROOT-TRANSLATION-CORRECTION.txt & ROOT-PSM-SCORES.txt Pyrodigal provenance."""

    def test_table4_internal_tga_uses_calling_code_and_rejects_forced_table11(self):
        gene_t4 = FakeTable4Gene("ATGGCTTGAGCTTAA", called_table=4)

        with self.assertRaises(ValueError) as ctx:
            translate_pyrodigal_gene_consistent(gene_t4, requested_translation_table=11)
        self.assertIn("translation_table=4", str(ctx.exception))

        res = translate_pyrodigal_gene_consistent(gene_t4)
        self.assertEqual(res["called_translation_table"], 4)
        self.assertEqual(res["canonical_seq"], "MAWA")
        self.assertTrue(res["is_valid"])
        self.assertEqual(res["n_internal_stops"], 0)
        self.assertTrue(res["had_terminal_stop"])

    def test_missing_translation_table_attribute_raises_no_silent_11_fallback(self):
        class GeneWithoutTableAttr:
            def translate(self):
                return "MAAAA*"

        with self.assertRaises(AttributeError) as ctx:
            translate_pyrodigal_gene_consistent(GeneWithoutTableAttr())
        self.assertIn("translation_table", str(ctx.exception))

    def test_real_pyrodigal_api_when_available(self):
        try:
            import pyrodigal
        except ImportError:
            return
        # Train a single-genome GeneFinder with translation_table=4 on a synthetic sequence containing internal TGA
        coding_t4 = "ATGGCT" + ("GCTTGA" * 45) + "GCTTAA"
        contig = ("A" * 120) + coding_t4 + ("A" * 120)
        gf4 = pyrodigal.GeneFinder(meta=False)
        gf4.train(contig.encode("ascii") * 50, translation_table=4)
        genes = gf4.find_genes(contig.encode("ascii"))
        self.assertGreaterEqual(len(genes), 1)
        g = genes[0]
        self.assertEqual(g.translation_table, 4)
        res = translate_pyrodigal_gene_consistent(g)
        self.assertEqual(res["called_translation_table"], 4)
        self.assertTrue(res["is_valid"])
        self.assertEqual(res["n_internal_stops"], 0)
        with self.assertRaises(ValueError):
            translate_pyrodigal_gene_consistent(g, requested_translation_table=11)

    def test_internal_stop_is_never_deleted_during_canonicalization(self):
        corrupted_aa = "MA*A*"
        res = canonicalize_and_validate_protein_seq(corrupted_aa)
        self.assertEqual(res["canonical_seq"], "MA*A")
        self.assertNotEqual(res["canonical_seq"], "MAA")
        self.assertFalse(res["is_valid"])
        self.assertEqual(res["validity_status"], "invalid_internal_stop_codon")
        self.assertEqual(res["n_internal_stops"], 1)

        res_deleted = canonicalize_and_validate_protein_seq("MAA*")
        self.assertNotEqual(res["canonical_sha256"], res_deleted["canonical_sha256"])
        self.assertFalse(
            verify_sequence_hash_join("orf1", corrupted_aa, "orf2", "MAA*", is_raw_sequence=True)
        )

    def test_gaps_and_punctuation_flagged_invalid_not_silently_repaired(self):
        gapped = "MKT-LVG*A*"
        res = canonicalize_and_validate_protein_seq(gapped)
        self.assertFalse(res["is_valid"])
        self.assertIn("-", res["canonical_seq"])
        self.assertIn("*", res["canonical_seq"])


class TestStage1ReviewControls(unittest.TestCase):
    """Tests for ROOT-STAGE1-REVIEW.txt."""

    def test_alias_collision_raises_instead_of_overwriting(self):
        records = [
            {"MAG": "metabat.bin_101.fa", "score": 1},
            {"MAG": "maxbin.bin_101.fa", "score": 2},
        ]

        def bad_normalizer(s: str) -> str:
            for p in ("metabat.", "maxbin."):
                if s.startswith(p):
                    s = s[len(p):]
            return s.replace(".fa", "")

        with self.assertRaises(ValueError) as ctx:
            build_unique_alias_map(records, "MAG", bad_normalizer, strict=True)
        self.assertIn("collision", str(ctx.exception).lower())

    def test_unmapped_query_or_target_id_raises_keyerror(self):
        query_meta = {"gene_A": {"family_id": "fam_1"}}
        with self.assertRaises(KeyError):
            lookup_strict(query_meta, "unexpected_gene_B", "DIAMOND query lookup")

    def test_short_transmembrane_or_domain_hit_rejected_even_with_low_evalue(self):
        res = evaluate_homology_alignment(
            pident=65.0,
            qstart=10,
            qend=44,
            qlen=350,
            sstart=20,
            send=54,
            slen=340,
            evalue=1e-12,
            bitscore=95.0,
        )
        self.assertFalse(res["passes_strict_homology_gate"])
        self.assertFalse(res["passes_moderate_homology_gate"])
        self.assertEqual(res["homology_gate_status"], "rejected_short_or_generic_domain_only")

    def test_invalid_alignment_bounds_and_nan_evalue_rejected(self):
        with self.assertRaises(ValueError):
            evaluate_homology_alignment(50.0, 1, 400, 300, 1, 200, 250, 1e-10, 100.0)
        with self.assertRaises(ValueError):
            evaluate_homology_alignment(50.0, 1, 200, 300, 1, 200, 250, float("nan"), 100.0)


class TestDownloadReviewControls(unittest.TestCase):
    """Tests for ROOT-DOWNLOAD-REVIEW.txt."""

    def test_large_html_200_error_page_over_100kb_rejected(self):
        with tempfile.TemporaryDirectory(dir=os.path.dirname(os.path.abspath(__file__))) as tmpd:
            html_file = os.path.join(tmpd, "fake.mzTab")
            with open(html_file, "w", encoding="utf-8") as f:
                f.write("<!DOCTYPE html><html><body>Rate limit exceeded" + ("A" * 150000) + "</body></html>\n")
            ok, reason, _ = validate_mztab_file(html_file, 10, 5, "sample_1")
            self.assertFalse(ok)
            self.assertEqual(reason, "html_error_payload_not_mztab")

    def test_truncated_mztab_and_row_count_mismatch_rejected(self):
        with tempfile.TemporaryDirectory(dir=os.path.dirname(os.path.abspath(__file__))) as tmpd:
            mz_path = os.path.join(tmpd, "test.mzTab")
            valid_content = (
                "MTD\tmzTab-version\t1.0.0\n"
                "MTD\tms_run[1]-location\tfile://sample_1.mzML\n"
                "PRH\taccession\tdescription\n"
                "PRT\tprot1\tdesc1\n"
                "PSH\tsequence\tPSM_ID\taccession\topt_global_cv_MS:1002217_decoy_peptide\topt_global_cv_MS:1002354_PSM-level_q-value\n"
                "PSM\tACDEFGHIK\t1\tprot1\t0\t0.001\n"
                "PSM\tLMNPQRSTV\t2\tprot1\t0\t0.002\n"
            )
            with open(mz_path, "w", encoding="utf-8") as f:
                f.write(valid_content.rstrip("\n"))
            ok, reason, _ = validate_mztab_file(mz_path, 2, 1, "sample_1")
            self.assertFalse(ok)
            self.assertEqual(reason, "missing_trailing_newline_truncated")

            with open(mz_path, "w", encoding="utf-8") as f:
                f.write(valid_content)
            ok, reason, stats = validate_mztab_file(mz_path, 3, 1, "sample_1")
            self.assertFalse(ok)
            self.assertIn("psm_row_count_mismatch", reason)

            ok, reason, stats = validate_mztab_file(mz_path, 2, 1, "sample_1")
            self.assertTrue(ok)
            recorded_sha = stats["sha256"]

            with open(mz_path, "w", encoding="utf-8") as f:
                f.write(valid_content.replace("ACDEFGHIK", "ACDEFGHIR"))
            ok2, reason2, _ = validate_mztab_file(mz_path, 2, 1, "sample_1", previous_sha256=recorded_sha)
            self.assertFalse(ok2)
            self.assertIn("sha256_mismatch_with_previously_verified", reason2)


class TestAnalysisContractRule7AndRootPsmScores(unittest.TestCase):
    """Tests all mandatory negative paths in ROOT-ANALYSIS-CONTRACT.txt Rule 7 and ROOT-PSM-SCORES.txt."""

    def test_1_and_2_same_counts_with_replaced_id_or_corrupted_sequence_hash_rejected(self):
        seq_orig = "MKTLLLTLVVVTIVCLDLGYT*"
        seq_corrupted = "MKTLLLTLVVVTIVCLDLGYA*"
        self.assertFalse(
            verify_sequence_hash_join("Ga0577_1", seq_orig, "Ga0577_1", seq_corrupted, is_raw_sequence=True)
        )

    def test_3_decoy_hits_excluded_and_inconsistent_xxx_decoy_rejected(self):
        rows = [
            {
                "PSM_ID": "1",
                "sequence": "PEPTIDEAK",
                "accession": "XXX_Ga0577_decoy",
                "search_engine_score[3]": "0.0005",
                "search_engine_score[5]": "1e-12",
                "search_engine_score[6]": "0.0001",
                "opt_global_cv_MS:1002217_decoy_peptide": "1",
                "opt_global_cv_MS:1002354_PSM-level_q-value": "0.0001",
                "opt_global_valid": "VALID",
            },
            {
                "PSM_ID": "2",
                "sequence": "TARGETPEPK",
                "accession": "Ga0577_target",
                "search_engine_score[3]": "0.002",
                "search_engine_score[5]": "1e-14",
                "search_engine_score[6]": "0.001",
                "opt_global_cv_MS:1002217_decoy_peptide": "0",
                "opt_global_cv_MS:1002354_PSM-level_q-value": "0.001",
                "opt_global_valid": "VALID",
            },
        ]
        res = filter_and_audit_mztab_psms(rows, mtd_lines=REAL_MTD_LINES, max_psm_q=0.01, max_peptide_q=0.01)
        self.assertEqual(len(res["qualifying_target_psms"]), 1)
        self.assertEqual(res["qualifying_target_psms"][0]["accession"], "Ga0577_target")

        bad_rows = [
            {
                "PSM_ID": "1",
                "sequence": "PEPTIDEAK",
                "accession": "XXX_Ga0577_decoy",
                "search_engine_score[3]": "0.0005",
                "search_engine_score[5]": "1e-12",
                "search_engine_score[6]": "0.0001",
                "opt_global_cv_MS:1002217_decoy_peptide": "0",
                "opt_global_cv_MS:1002354_PSM-level_q-value": "0.0001",
            }
        ]
        with self.assertRaises(ValueError):
            filter_and_audit_mztab_psms(bad_rows, mtd_lines=REAL_MTD_LINES)

    def test_upstream_msgf_pepqvalue_enforced_not_reconstructed_from_psmq(self):
        # Row 1 has passing PSM QValue (0.004) and 0 decoys in table, so a naive reconstructed
        # TDC peptide FDR would give 0.0, but upstream MS-GF:PepQValue (search_engine_score[3]) is 0.025 (>0.01)!
        rows = [
            {
                "PSM_ID": "1",
                "sequence": "LOWPSMQBUTHIGHPEPQK",
                "accession": "Ga0577_target1",
                "search_engine_score[3]": "0.025",  # PepQValue > 0.01 -> MUST BE EXCLUDED
                "search_engine_score[5]": "1e-11",
                "search_engine_score[6]": "0.004",  # PSM QValue <= 0.01
                "opt_global_cv_MS:1002217_decoy_peptide": "0",
                "opt_global_cv_MS:1002354_PSM-level_q-value": "0.004",
                "opt_global_valid": "VALID",
            },
            {
                "PSM_ID": "2",
                "sequence": "BOTHPASSINGPEPK",
                "accession": "Ga0577_target2",
                "search_engine_score[3]": "0.008",  # PepQValue <= 0.01 -> RETAINED
                "search_engine_score[5]": "1e-13",
                "search_engine_score[6]": "0.002",
                "opt_global_cv_MS:1002217_decoy_peptide": "0",
                "opt_global_cv_MS:1002354_PSM-level_q-value": "0.002",
                "opt_global_valid": "VALID",
            },
        ]
        res = filter_and_audit_mztab_psms(rows, mtd_lines=REAL_MTD_LINES, max_psm_q=0.01, max_peptide_q=0.01)
        self.assertEqual(len(res["qualifying_target_psms"]), 1)
        self.assertEqual(res["qualifying_target_psms"][0]["clean_peptide"], "BOTHPASSINGPEPK")

    def test_discordant_psm_q_and_out_of_range_or_invalid_rows_rejected(self):
        # 1. Discordance between opt_global PSM q-value and search_engine_score[6] (MS-GF:QValue)
        discordant_rows = [
            {
                "PSM_ID": "1",
                "sequence": "TARGETPEPK",
                "accession": "Ga0577_target",
                "search_engine_score[3]": "0.002",
                "search_engine_score[5]": "1e-12",
                "search_engine_score[6]": "0.050",  # Discordant with 0.001!
                "opt_global_cv_MS:1002217_decoy_peptide": "0",
                "opt_global_cv_MS:1002354_PSM-level_q-value": "0.001",
            }
        ]
        with self.assertRaises(ValueError) as ctx:
            filter_and_audit_mztab_psms(discordant_rows, mtd_lines=REAL_MTD_LINES)
        self.assertIn("discordant PSM q-values", str(ctx.exception))

        # 2. PepQValue > 1.0 or NaN rejected
        for bad_q in ("1.25", "nan", "-0.01", "null"):
            bad_q_rows = [
                {
                    "PSM_ID": "1",
                    "sequence": "TARGETPEPK",
                    "accession": "Ga0577_target",
                    "search_engine_score[3]": bad_q,
                    "search_engine_score[5]": "1e-12",
                    "search_engine_score[6]": "0.001",
                    "opt_global_cv_MS:1002217_decoy_peptide": "0",
                    "opt_global_cv_MS:1002354_PSM-level_q-value": "0.001",
                }
            ]
            with self.assertRaises(ValueError):
                filter_and_audit_mztab_psms(bad_q_rows, mtd_lines=REAL_MTD_LINES)

        # 3. Missing MS:1002055 (PepQValue) in MTD metadata rejected
        incomplete_mtd = [line for line in REAL_MTD_LINES if "MS:1002055" not in line]
        with self.assertRaises(ValueError) as ctx_mtd:
            parse_mztab_score_metadata(incomplete_mtd)
        self.assertIn("pep_q_col", str(ctx_mtd.exception))

        # 4. opt_global_valid != VALID excluded and audited
        invalid_flag_rows = [
            {
                "PSM_ID": "99",
                "sequence": "TARGETPEPK",
                "accession": "Ga0577_target",
                "search_engine_score[3]": "0.001",
                "search_engine_score[5]": "1e-12",
                "search_engine_score[6]": "0.001",
                "opt_global_cv_MS:1002217_decoy_peptide": "0",
                "opt_global_cv_MS:1002354_PSM-level_q-value": "0.001",
                "opt_global_valid": "INVALID",
                "opt_global_invalid_reason": "BelowThreshold",
            }
        ]
        res_inv = filter_and_audit_mztab_psms(invalid_flag_rows, mtd_lines=REAL_MTD_LINES)
        self.assertEqual(res_inv["n_excluded_invalid_rows"], 1)
        self.assertEqual(res_inv["n_qualifying_target_psms"], 0)

    def test_4_shared_peptide_across_mags_cannot_assign_single_mag(self):
        bg_proteins = {
            "orf_mag1": "MKTAAAPEPTIDESHAREDKVVV",
            "orf_mag2": "MGGGGGPEPTIDESHAREDKLLL",
        }
        prot_to_mags = {
            "orf_mag1": {"MAG_Thermo_1"},
            "orf_mag2": {"MAG_Bacteria_99"},
        }
        prot_to_fams = {
            "orf_mag1": {"fam_dark_1"},
            "orf_mag2": {"fam_other_2"},
        }
        spec = classify_peptide_specificity_against_background(
            "PEPTIDESHAREDK",
            {"orf_mag1"},
            bg_proteins,
            prot_to_mags,
            prot_to_fams,
        )
        self.assertFalse(spec["is_unambiguous_for_single_mag"])
        self.assertFalse(spec["is_unambiguous_for_family"])
        self.assertEqual(spec["specificity_status"], "ambiguous_shared_across_multiple_mags_or_families")

    def test_5_isobaric_il_collision_with_background_flagged_ambiguous(self):
        bg_proteins = {
            "target_orf": "MMMMAVILPEPTIDEKGGG",
            "bg_competitor_orf": "CCCCAVLLPEPTIDEKTTT",
        }
        prot_to_mags = {
            "target_orf": {"MAG_E2_JACPAO01"},
            "bg_competitor_orf": {"MAG_Acidobacteriota_1"},
        }
        spec = classify_peptide_specificity_against_background(
            "AVILPEPTIDEK",
            {"target_orf"},
            bg_proteins,
            prot_to_mags,
        )
        self.assertTrue(spec["has_isobaric_il_collision"])
        self.assertIn("bg_competitor_orf", spec["isobaric_il_collision_proteins"])
        self.assertFalse(spec["is_unambiguous_for_single_mag"])
        self.assertEqual(spec["specificity_status"], "ambiguous_isobaric_IL_collision_with_background")

    def test_6_pseudoreplication_n15_and_duplicate_replicates_and_2016_dna_rejected(self):
        depths = ["0-10", "10-30", "30-45", "45-60", "60-80"]
        samples = []
        for p in ("1", "2", "3"):
            for d in depths:
                samples.append({
                    "sample_id": f"P{p}_C_{d}",
                    "plot_pair_id": p,
                    "treatment": "control",
                    "depth_cm": d,
                    "sampling_year": 2018,
                    "paired_dna_year": 2018,
                    "value": 1.0,
                })
                samples.append({
                    "sample_id": f"P{p}_H_{d}",
                    "plot_pair_id": p,
                    "treatment": "heated",
                    "depth_cm": d,
                    "sampling_year": 2018,
                    "paired_dna_year": 2018,
                    "value": 2.0,
                })

        with self.assertRaises(ValueError) as ctx:
            compute_paired_warming_effect_n3_plots(samples, claimed_independent_n=15)
        self.assertIn("Pseudoreplication error", str(ctx.exception))

        res = compute_paired_warming_effect_n3_plots(samples)
        self.assertEqual(res["n_independent_plot_pairs"], 3)
        self.assertEqual(res["n_nested_depth_pairs"], 15)
        self.assertAlmostEqual(res["exact_whole_plot_sign_flip_two_sided_p"], 0.25)
        self.assertAlmostEqual(res["min_achievable_two_sided_p_for_n3"], 0.25)

        dup_samples = samples + [{
            "sample_id": "P1_C_0-10_dup",
            "plot_pair_id": "1",
            "treatment": "control",
            "depth_cm": "0-10",
            "sampling_year": 2018,
            "paired_dna_year": 2018,
            "value": 1.0,
        }]
        with self.assertRaises(ValueError):
            validate_karaoz_experimental_design(dup_samples)

        bad_year_samples = [dict(x) for x in samples]
        bad_year_samples[0]["paired_dna_year"] = 2016
        with self.assertRaises(ValueError) as ctx2:
            validate_karaoz_experimental_design(bad_year_samples)
        self.assertIn("2014/2016 metagenomes must not be pooled", str(ctx2.exception))

    def test_7_missing_psm_qvalue_rejected(self):
        rows_missing_q = [
            {
                "PSM_ID": "1",
                "sequence": "TARGETPEPK",
                "accession": "Ga0577_target",
                "search_engine_score[3]": "0.001",
                "search_engine_score[5]": "1e-12",
                "search_engine_score[6]": "0.001",
                "opt_global_cv_MS:1002217_decoy_peptide": "0",
                "opt_global_cv_MS:1002354_PSM-level_q-value": "null",
            }
        ]
        with self.assertRaises(ValueError) as ctx:
            filter_and_audit_mztab_psms(rows_missing_q, mtd_lines=REAL_MTD_LINES)
        self.assertIn("missing per-hit", str(ctx.exception))

    def test_8_absent_search_db_classified_nonestimable_never_negative(self):
        s1 = classify_protein_and_family_detection(False, False, 0, 0, 0)
        self.assertEqual(s1, "nonestimable_no_homolog_in_study_genomes")

        # When search_database_exclusion_verified is False (or unverified full search background),
        # absence from observed EMSL does NOT prove absent from search DB (ROOT-FINAL-CLAIM-REVIEW.txt line 2):
        s2_unverified = classify_protein_and_family_detection(True, False, 0, 0, 0, full_background_verified=False)
        self.assertEqual(s2_unverified, "nonestimable_search_database_inclusion_unverified")

        # Only when explicit complete search FASTA proves exclusion is it classified as absent from search DB:
        s2_verified = classify_protein_and_family_detection(True, False, 0, 0, 0, search_database_exclusion_verified=True)
        self.assertEqual(s2_verified, "nonestimable_absent_from_search_database")

        s3 = classify_protein_and_family_detection(True, True, 1, 0, 5, full_background_verified=True)
        self.assertEqual(s3, "provisional_single_unambiguous_peptide_q01")
        s4 = classify_protein_and_family_detection(True, True, 2, 0, 5, full_background_verified=True)
        self.assertEqual(s4, "confirmed_detected_ge2_unambiguous_peptides_q01")

        qmap = benjamini_hochberg_fdr({"lead_1": None, "fam_a": 0.25, "fam_b": 0.50})
        self.assertIsNone(qmap["lead_1"])
        self.assertAlmostEqual(qmap["fam_a"], 0.50)


class TestBackgroundAndRuntimeRepairControls(unittest.TestCase):
    """
    Tests required by ROOT-BACKGROUND-REVIEW.txt and ROOT-RUNTIME-REPAIR.txt:
    1. Same peptide present only in the broader background (not in EMSL 122k subset) is flagged as shared/ambiguous.
    2. Two accessions/aliases for a single spectrum count ONCE via unique spectrum union, never summed to 2.
    3. Missing full search background (full_background_verified=False) is never classified as
       confirmed_detected_ge2_unambiguous_peptides_q01.
    4. Ambiguous family/paralog peptides (including two paralogs in different families within the same MAG)
       cannot drive candidate-family warming abundance.
    5. Fixture with actual 27-column exp4_47_flank_evidence.tsv header validates both Primary Leads
       and all 5 Partner Families after join, and raises on missing header columns or missing partner records.
    """

    def test_1_peptide_present_only_in_broader_background_flagged_shared(self):
        # Peptide is unique in EMSL 122k subset ('emsl_thermo_orf1'), but a broader background
        # protein ('unbinned_bg_orf_999999') contains an exact or I/L-isobaric match.
        records = [
            ("emsl_thermo_orf1", "MKTAAAPEPTIDEUNIQEMSLKVVV"),
            ("unbinned_bg_orf_999999", "MGGGGGPEPTLDEUNLQEMSLKLLL"),  # I/L isobaric match in broader background!
        ]
        stream_res = stream_multipattern_peptide_search_records(records, ["PEPTIDEUNIQEMSLK"])
        self.assertEqual(stream_res["PEPTIDEUNIQEMSLK"]["exact_matches"], {"emsl_thermo_orf1"})
        self.assertEqual(
            stream_res["PEPTIDEUNIQEMSLK"]["il_matches"],
            {"emsl_thermo_orf1", "unbinned_bg_orf_999999"},
        )

        prot_to_mags = {
            "emsl_thermo_orf1": {"MAG_Thermo_1"},
            "unbinned_bg_orf_999999": {"unbinned_or_non_target_search_db_entry"},
        }
        spec = classify_peptide_specificity_against_background(
            "PEPTIDEUNIQEMSLK",
            {"emsl_thermo_orf1"},
            background_proteins={},
            protein_to_mags=prot_to_mags,
            precomputed_exact_matches=stream_res["PEPTIDEUNIQEMSLK"]["exact_matches"],
            precomputed_il_matches=stream_res["PEPTIDEUNIQEMSLK"]["il_matches"],
        )
        self.assertFalse(spec["is_unambiguous_for_single_mag"])
        self.assertFalse(spec["is_unambiguous_for_family"])
        self.assertEqual(spec["specificity_status"], "ambiguous_isobaric_IL_collision_with_background")

    def test_2_two_accessions_for_one_spectrum_count_once_not_summed(self):
        # Same spectrum (spectra_ref='ms_run[1]:index=4412') matched to two duplicate/alias accessions
        # for the same target entity in sample P1_H_0-10. Must count as 1 spectrum, never 2!
        psm_rows = [
            {
                "PSM_ID": "101",
                "spectra_ref": "ms_run[1]:index=4412",
                "accession": "Ga0577_10001_1_200",
                "clean_peptide": "UNIQTHERMOPEPK",
            },
            {
                "PSM_ID": "101",
                "spectra_ref": "ms_run[1]:index=4412",
                "accession": "Ga0577_20002_1_200",  # Duplicate/alias protein for same spectrum
                "clean_peptide": "UNIQTHERMOPEPK",
            },
            {
                "PSM_ID": "205",
                "spectra_ref": "ms_run[1]:index=8819",
                "accession": "Ga0577_10001_1_200",
                "clean_peptide": "UNIQTHERMOPEPK",
            },
        ]
        unique_specs = aggregate_unique_spectra_for_entity(
            psm_rows,
            target_accessions={"Ga0577_10001_1_200", "Ga0577_20002_1_200"},
            allowed_clean_peptides={"UNIQTHERMOPEPK"},
        )
        self.assertEqual(len(unique_specs), 2)
        self.assertEqual(unique_specs, {"ms_run[1]:index=4412", "ms_run[1]:index=8819"})

    def test_3_missing_full_background_never_called_confirmed_unambiguous(self):
        status_unverified = classify_protein_and_family_detection(
            has_karaoz_homolog=True,
            in_search_database=True,
            n_unambiguous_il_peptides_q01=3,
            n_family_shared_il_peptides_q01=0,
            reported_spectral_count=12,
            full_background_verified=False,
        )
        self.assertNotEqual(status_unverified, "confirmed_detected_ge2_unambiguous_peptides_q01")
        self.assertEqual(
            status_unverified,
            "provisional_observed_subset_ge2_peptides_q01_full_background_unverified_primary_uniqueness_nonestimable",
        )

    def test_4_ambiguous_family_or_within_mag_paralog_peptides_excluded_from_family_warming(self):
        # Two paralogs in the SAME MAG ('MAG_Thermo_1'), but belonging to different families
        # ('fam_candidate_A' vs 'fam_paralog_B' or unassigned).
        bg_proteins = {
            "orf_paralog_A": "MKTAAASHAREDPARALOGPEPKVVV",
            "orf_paralog_B": "MGGGGGSHAREDPARALOGPEPKLLL",
            "orf_uniq_A": "MKTAAAONLYINFAMAPEPKVVV",
        }
        prot_to_mags = {
            "orf_paralog_A": {"MAG_Thermo_1"},
            "orf_paralog_B": {"MAG_Thermo_1"},
            "orf_uniq_A": {"MAG_Thermo_1"},
        }
        prot_to_fams = {
            "orf_paralog_A": {"fam_candidate_A"},
            "orf_paralog_B": {"fam_paralog_B"},
            "orf_uniq_A": {"fam_candidate_A"},
        }
        spec_shared = classify_peptide_specificity_against_background(
            "SHAREDPARALOGPEPK",
            {"orf_paralog_A"},
            bg_proteins,
            prot_to_mags,
            prot_to_fams,
        )
        # Unambiguous for the single MAG, but NOT unambiguous for the family!
        self.assertTrue(spec_shared["is_unambiguous_for_single_mag"])
        self.assertFalse(spec_shared["is_unambiguous_for_family"])
        self.assertEqual(
            spec_shared["specificity_status"],
            "ambiguous_shared_across_multiple_paralog_families_within_mag",
        )

        spec_uniq = classify_peptide_specificity_against_background(
            "ONLYINFAMAPEPK",
            {"orf_uniq_A"},
            bg_proteins,
            prot_to_mags,
            prot_to_fams,
        )
        self.assertTrue(spec_uniq["is_unambiguous_for_family"])

        # Only peptides with is_unambiguous_for_family=True are allowed to drive family warming abundance:
        allowed_fam_peptides = {
            pep for pep, s in [("SHAREDPARALOGPEPK", spec_shared), ("ONLYINFAMAPEPK", spec_uniq)]
            if s["is_unambiguous_for_family"]
        }
        self.assertEqual(allowed_fam_peptides, {"ONLYINFAMAPEPK"})

        psm_rows = [
            {
                "PSM_ID": "1",
                "spectra_ref": "scan=10",
                "accession": "orf_paralog_A",
                "clean_peptide": "SHAREDPARALOGPEPK",  # Must be excluded from family count!
            },
            {
                "PSM_ID": "2",
                "spectra_ref": "scan=20",
                "accession": "orf_uniq_A",
                "clean_peptide": "ONLYINFAMAPEPK",
            },
        ]
        fam_specs = aggregate_unique_spectra_for_entity(
            psm_rows,
            target_accessions={"orf_paralog_A", "orf_uniq_A"},
            allowed_clean_peptides=allowed_fam_peptides,
        )
        self.assertEqual(fam_specs, {"scan=20"})

    def test_5_real_exp4_47_flank_evidence_schema_fixture_and_partner_validation(self):
        import csv
        import io

        primary_leads = {
            "22PRr1::22PRr1_2531955_3": "primary_lead_1_PilJ_archaellin_module",
            "V60PRr2::V60PRr2_1314979_4": "primary_lead_2_Slayer_coat_module",
        }
        partner_families = {
            "F37Dr2::F37Dr2_752038_3": "partner_Lead1_PilJ_archaellin_PF04964",
            "S20PRr1::S20PRr1_1478121_8": "partner_Lead1_archaellin_PF04964",
            "F22PRr2::F22PRr2_1079523_133": "partner_Lead2_Slayer_csg_PF10731",
            "38Dr3::38Dr3_563471_1": "partner_Lead2_Slayer_csg_PF18204",
            "25PRr1::25PRr1_1174463_140": "partner_Lead2_Slayer_csg_PF04711",
        }
        # Exact 27-column header from Shaheen exp4_47_flank_evidence.tsv
        header = (
            "family_id\tfocal_member_gene\tfocal_member_is_verified_interior\tfocal_genome\t"
            "focal_contig\trelative_gene_offset\tflanking_gene\tflanking_family_id\t"
            "flanking_is_dark\tflanking_signalp_class\tflanking_functional_description\t"
            "flanking_specific_guilds\tannotation_source\tannotation_profile_or_hit_id\t"
            "annotation_confidence_or_evalue\tgenomic_strand\tfocal_strand\tflanking_strand\t"
            "same_strand\tfocal_start\tfocal_end\tflanking_start\tflanking_end\t"
            "intergenic_distance_bp\tcoordinate_provenance\tannotated_guild_support\t"
            "machinery_experimental_verification_status\n"
        )
        rows_txt = (
            "22PRr1::22PRr1_2531955_3\t1PRr2__self_SemiBin_81::1PRr2_468829_3\tTrue\t1PRr2__self_SemiBin_81\t1PRr2_468829\t-1\t1PRr2__self_SemiBin_81::1PRr2_468829_2\tF37Dr2::F37Dr2_752038_3\tTrue\tPILIN\tPfam:PF04964\tsurface_assembly_pilin_archaellin\tPfam_profile_match\tPF04964\t1e-10\tfocal:-1;flank:-1\t-1\t-1\tTrue\t2394\t3638\t1918\t2388\t5\tverified\tTrue\tunverified\n"
            "22PRr1::22PRr1_2531955_3\t51Dr3__self_SemiBin_1002::51Dr3_586870_7\tTrue\t51Dr3__self_SemiBin_1002\t51Dr3_586870\t-1\t51Dr3__self_SemiBin_1002::51Dr3_586870_6\tS20PRr1::S20PRr1_1478121_8\tFalse\tOTHER\tPfam:PF04964\tsurface_assembly_pilin_archaellin\tPfam_profile_match\tPF04964\t1e-12\tfocal:-1;flank:-1\t-1\t-1\tTrue\t11125\t12297\t10637\t11128\t-4\tverified\tTrue\tunverified\n"
            "V60PRr2::V60PRr2_1314979_4\tF22PRr2__single_SemiBin_97::F22PRr2_1079523_132\tTrue\tF22PRr2__single_SemiBin_97\tF22PRr2_1079523\t1\tF22PRr2__single_SemiBin_97::F22PRr2_1079523_133\tF22PRr2::F22PRr2_1079523_133\tFalse\tSP\tPfam:PF10731\tcell_coat_or_envelope\tPfam_profile_match\tPF10731\t1e-15\tfocal:1;flank:1\t1\t1\tTrue\t100\t500\t510\t900\t10\tverified\tTrue\tunverified\n"
            "V60PRr2::V60PRr2_1314979_4\t36PRr1__self_SemiBin_11::36PRr1_1538782_6\tTrue\t36PRr1__self_SemiBin_11\t36PRr1_1538782\t-1\t36PRr1__self_SemiBin_11::36PRr1_1538782_5\t38Dr3::38Dr3_563471_1\tFalse\tSP\tPfam:PF18204\tcell_coat_or_envelope\tPfam_profile_match\tPF18204\t1e-11\tfocal:-1;flank:-1\t-1\t-1\tTrue\t100\t500\t510\t900\t10\tverified\tTrue\tunverified\n"
            "V60PRr2::V60PRr2_1314979_4\t25PRr1__self_SemiBin_1693::25PRr1_1174463_141\tTrue\t25PRr1__self_SemiBin_1693\t25PRr1_1174463\t-1\t25PRr1__self_SemiBin_1693::25PRr1_1174463_140\t25PRr1::25PRr1_1174463_140\tFalse\tSP\tPfam:PF04711\tcell_coat_or_envelope\tPfam_profile_match\tPF04711\t1e-14\tfocal:-1;flank:-1\t-1\t-1\tTrue\t100\t500\t510\t900\t10\tverified\tTrue\tunverified\n"
        )
        reader = list(csv.DictReader(io.StringIO(header + rows_txt), delimiter="\t"))
        res = parse_and_validate_flank_evidence_rows(reader, primary_leads, partner_families)
        self.assertEqual(len(res["found_leads"]), 2)
        self.assertEqual(len(res["found_partners"]), 5)

        # If a partner family is missing after join, must raise ValueError (never silently skip)
        with self.assertRaises(ValueError):
            parse_and_validate_flank_evidence_rows(reader[:3], primary_leads, partner_families)

        # If wrong header is supplied, must raise KeyError
        bad_row = [{"family_id": "22PRr1::22PRr1_2531955_3", "flank_family_id": "wrong"}]
        with self.assertRaises(KeyError):
            parse_and_validate_flank_evidence_rows(bad_row, primary_leads, partner_families)


class TestFinalClaimReviewControls(unittest.TestCase):
    """
    Tests required by ROOT-FINAL-CLAIM-REVIEW.txt:
    1. Distinguishes unverified search database inclusion ('nonestimable_search_database_inclusion_unverified')
       from true verified exclusion ('nonestimable_absent_from_search_database') and absence of homologs
       in study genomes ('nonestimable_no_homolog_in_study_genomes').
    2. Enforces that primary family/clade warming estimability is strictly marked
       'nonestimable_full_search_background_unavailable' when full search FASTAs are missing,
       and that retained statistics are observed-subset provisional sensitivity diagnostics only.
    3. Verifies zero handling in paired warming differences without arbitrary 1e-6 pseudocounts
       producing artifactual log2FC ~24-30 values: pairs involving zero are reported with explicit
       detection states ('both_zero', 'heated_only_zero_control', 'control_only_zero_heated'),
       descriptive log2 ratio is None/NA for zero pairs, and primary differences use true (heated - control).
    """

    def test_1_unverified_search_database_inclusion_vs_explicit_exclusion(self):
        # Case A: Family has no homolog in study MAGs -> nonestimable_no_homolog_in_study_genomes
        s_no_homolog = classify_protein_and_family_detection(
            has_karaoz_homolog=False,
            in_search_database=False,
            n_unambiguous_il_peptides_q01=0,
            n_family_shared_il_peptides_q01=0,
            reported_spectral_count=0,
            full_background_verified=False,
        )
        self.assertEqual(s_no_homolog, "nonestimable_no_homolog_in_study_genomes")

        # Case B: Family has homolog in study MAGs, but no match in observed EMSL subset and full search FASTA missing:
        # Absence from observed EMSL does NOT prove absence from search DB!
        s_unverified = classify_protein_and_family_detection(
            has_karaoz_homolog=True,
            in_search_database=False,
            n_unambiguous_il_peptides_q01=0,
            n_family_shared_il_peptides_q01=0,
            reported_spectral_count=0,
            full_background_verified=False,
            search_database_exclusion_verified=False,
        )
        self.assertEqual(s_unverified, "nonestimable_search_database_inclusion_unverified")
        self.assertNotEqual(s_unverified, "nonestimable_absent_from_search_database")

        # Case C: Family has homolog in study MAGs, and explicit complete search FASTA proves absence:
        s_verified_absent = classify_protein_and_family_detection(
            has_karaoz_homolog=True,
            in_search_database=False,
            n_unambiguous_il_peptides_q01=0,
            n_family_shared_il_peptides_q01=0,
            reported_spectral_count=0,
            full_background_verified=True,
            search_database_exclusion_verified=True,
        )
        self.assertEqual(s_verified_absent, "nonestimable_absent_from_search_database")

    def test_2_primary_warming_estimability_status(self):
        # Even when a family or clade has >=2 qualifying peptides in the observed subset,
        # if full_background_verified=False, primary warming is nonestimable_full_search_background_unavailable:
        w_status = classify_primary_warming_estimability(
            has_karaoz_homolog=True,
            in_search_database=True,
            has_qualifying_unambiguous_peptides_in_observed_bg=True,
            reported_spectral_count=10,
            full_background_verified=False,
        )
        self.assertEqual(w_status, "nonestimable_full_search_background_unavailable")

        # When full background is verified, primary warming is estimable:
        w_status_estimable = classify_primary_warming_estimability(
            has_karaoz_homolog=True,
            in_search_database=True,
            has_qualifying_unambiguous_peptides_in_observed_bg=True,
            reported_spectral_count=10,
            full_background_verified=True,
            search_database_exclusion_verified=True,
        )
        self.assertEqual(w_status_estimable, "estimable_confirmed_full_search_background")

        # When no homolog exists:
        w_status_no_hom = classify_primary_warming_estimability(
            has_karaoz_homolog=False,
            in_search_database=False,
            has_qualifying_unambiguous_peptides_in_observed_bg=False,
            reported_spectral_count=0,
            full_background_verified=False,
        )
        self.assertEqual(w_status_no_hom, "nonestimable_no_homolog_in_study_genomes")

    def test_3_paired_warming_zero_handling_no_arbitrary_pseudocount_logfc(self):
        depths = ["0-10", "10-30", "30-45", "45-60", "60-80"]
        samples = []
        for p in ("1", "2", "3"):
            for d in depths:
                # Plot 1: both positive at 0-10, heated only at 10-30, control only at 30-45, both zero at 45-60 and 60-80
                if p == "1" and d == "0-10":
                    val_c, val_h = 4.0, 8.0
                elif p == "1" and d == "10-30":
                    val_c, val_h = 0.0, 5.0
                elif p == "1" and d == "30-45":
                    val_c, val_h = 6.0, 0.0
                elif p == "2" and d == "0-10":
                    val_c, val_h = 2.0, 6.0
                else:
                    val_c, val_h = 0.0, 0.0

                samples.append({
                    "sample_id": f"P{p}_C_{d}",
                    "plot_pair_id": p,
                    "treatment": "control",
                    "depth_cm": d,
                    "sampling_year": 2018,
                    "paired_dna_year": 2018,
                    "value": val_c,
                })
                samples.append({
                    "sample_id": f"P{p}_H_{d}",
                    "plot_pair_id": p,
                    "treatment": "heated",
                    "depth_cm": d,
                    "sampling_year": 2018,
                    "paired_dna_year": 2018,
                    "value": val_h,
                })

        res = compute_paired_warming_effect_n3_plots(samples)
        dp = res["depth_profile"]

        # 0-10 on Plot 1: both positive -> descriptive log2 ratio is log2(8/4) = 1.0
        self.assertEqual(dp["0-10"]["plot1_detection_state"], "both_positive")
        self.assertAlmostEqual(dp["0-10"]["plot1_descriptive_log2_ratio_both_positive_only"], 1.0)
        self.assertEqual(dp["0-10"]["plot1_delta"], 4.0)

        # 10-30 on Plot 1: heated only (0 vs 5) -> log2 ratio must be None, NOT ~24
        self.assertEqual(dp["10-30"]["plot1_detection_state"], "heated_only_zero_control")
        self.assertIsNone(dp["10-30"]["plot1_descriptive_log2_ratio_both_positive_only"])
        self.assertEqual(dp["10-30"]["plot1_delta"], 5.0)

        # 30-45 on Plot 1: control only (6 vs 0) -> log2 ratio must be None, NOT ~-24
        self.assertEqual(dp["30-45"]["plot1_detection_state"], "control_only_zero_heated")
        self.assertIsNone(dp["30-45"]["plot1_descriptive_log2_ratio_both_positive_only"])
        self.assertEqual(dp["30-45"]["plot1_delta"], -6.0)

        # 45-60 on Plot 1: both zero -> log2 ratio must be None, delta must be 0.0
        self.assertEqual(dp["45-60"]["plot1_detection_state"], "both_zero")
        self.assertIsNone(dp["45-60"]["plot1_descriptive_log2_ratio_both_positive_only"])
        self.assertEqual(dp["45-60"]["plot1_delta"], 0.0)

        # Summary must NOT contain arbitrary mean log2FC over changing subsets
        self.assertNotIn("overall_mean_log2fc_heated_vs_control", res)
        # Summary must report true mean paired difference (heated - control)
        self.assertIn("overall_mean_delta_heated_minus_control", res)
        self.assertEqual(res["n_depth_pairs_both_positive"], 2)
        self.assertEqual(res["n_depth_pairs_heated_only_zero_control"], 1)
        self.assertEqual(res["n_depth_pairs_control_only_zero_heated"], 1)
        self.assertEqual(res["n_depth_pairs_both_zero"], 11)

    def test_4_actual_caller_paths_author_subset_and_conservative_both_inclusion_unverified_when_unmatched(self):
        """
        Regression test required by ROOT-SUBSET-LABEL-FINAL.txt:
        In both protein and family caller paths, when full search background FASTA is unavailable
        (full_background_verified=False), an unmatched sequence (in_search_database=False) must be
        classified as 'nonestimable_search_database_inclusion_unverified' in BOTH conservative AND
        dataset_reported_observed_subset tiers. Observed subset detection grades where matches exist
        are preserved.
        """
        full_background_verified = False
        has_homolog = True
        in_search_db_confirmed = False

        # Caller path 1: Conservative detection tier
        cons_status = classify_protein_and_family_detection(
            has_karaoz_homolog=has_homolog,
            in_search_database=in_search_db_confirmed,
            n_unambiguous_il_peptides_q01=0,
            n_family_shared_il_peptides_q01=0,
            reported_spectral_count=0,
            full_background_verified=full_background_verified,
            search_database_exclusion_verified=full_background_verified,
        )
        self.assertEqual(cons_status, "nonestimable_search_database_inclusion_unverified")

        # Caller path 2: Author / dataset-reported observed subset tier
        author_status = classify_protein_and_family_detection(
            has_karaoz_homolog=has_homolog,
            in_search_database=in_search_db_confirmed,
            n_unambiguous_il_peptides_q01=0,
            n_family_shared_il_peptides_q01=0,
            reported_spectral_count=0,
            full_background_verified=True,  # Observed-subset reporting mode
            search_database_exclusion_verified=full_background_verified,  # Full search FASTA is unavailable
        )
        self.assertEqual(author_status, "nonestimable_search_database_inclusion_unverified")
        self.assertNotEqual(author_status, "nonestimable_absent_from_search_database")

        # Caller path 3: When observed matches DO exist, observed subset detection grade is preserved
        author_detected_status = classify_protein_and_family_detection(
            has_karaoz_homolog=True,
            in_search_database=True,
            n_unambiguous_il_peptides_q01=2,
            n_family_shared_il_peptides_q01=0,
            reported_spectral_count=5,
            full_background_verified=True,
            search_database_exclusion_verified=full_background_verified,
        )
        self.assertEqual(author_detected_status, "confirmed_detected_ge2_unambiguous_peptides_q01")


class TestFrozenStage1AndStage2Results(unittest.TestCase):
    """Verifies the frozen Stage 1 v2 and Stage 2 results packaged in reproducibility/expression_20261110."""

    def test_stage1_v2_and_stage2_frozen_summaries_and_sha256_manifest(self):
        s1_path = _PACKAGE_DIR / "results" / "stage1_v2_homology_eligibility" / "stage1_v2_homology_eligibility_summary.json"
        s2_path = _PACKAGE_DIR / "results" / "stage2_expression_validation" / "stage2_expression_validation_summary.json"
        self.assertTrue(s1_path.is_file())
        self.assertTrue(s2_path.is_file())

        s1 = json.loads(s1_path.read_text(encoding="utf-8"))
        s2 = json.loads(s2_path.read_text(encoding="utf-8"))

        self.assertEqual(s1["counts"]["total_derep_mags"], 2321)
        self.assertEqual(s1["counts"]["derep_archaeal_mags"], 113)
        self.assertEqual(s1["counts"]["derep_thermoplasmatota_mags"], 38)
        self.assertEqual(s1["counts"]["predicted_orfs_in_113_archaeal_mags_v2"], 274006)
        self.assertEqual(s1["counts"]["emsl50964_total_valid_detected_proteins"], 122370)
        self.assertEqual(s1["counts"]["emsl50964_proteins_linked_to_archaeal_mags"], 333)
        self.assertEqual(s1["counts"]["emsl50964_proteins_linked_to_thermoplasmatota_mags"], 78)
        self.assertEqual(s1["counts"]["total_verified_query_families"], 1514)

        self.assertEqual(s2["runtime_verified_audit_counts"]["n_mztab_files_audited"], 30)
        self.assertEqual(s2["runtime_verified_audit_counts"]["total_discordant_psm_q_rows"], 0)
        self.assertFalse(s2["runtime_verified_audit_counts"]["full_15split_composite_search_fasta_deposited"])
        self.assertEqual(
            s2["primary_clade_warming_estimability_status"],
            "nonestimable_full_search_background_unavailable",
        )
        self.assertEqual(
            s2["eq_family_conservative_detection_summary_by_cohort_role"]["primary_lead_1_pilj_archaellum_context"],
            {"nonestimable_search_database_inclusion_unverified": 1},
        )
        self.assertEqual(
            s2["eq_family_conservative_detection_summary_by_cohort_role"]["primary_lead_2_slayer_coat_context"],
            {"nonestimable_no_homolog_in_study_genomes": 1},
        )

        manifest_path = _PACKAGE_DIR / "manifests" / "sha256_manifest.tsv"
        self.assertTrue(manifest_path.is_file())
        with manifest_path.open("r", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        self.assertGreaterEqual(len(rows), 20)
        for row in rows:
            p = _PACKAGE_DIR / row["relative_path"]
            self.assertTrue(p.is_file(), f"Missing file: {row['relative_path']}")
            self.assertEqual(p.stat().st_size, int(row["size_bytes"]))
            self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(), row["sha256"])


if __name__ == "__main__":
    unittest.main()
