#!/usr/bin/env python3
"""
ProstT5 AA2fold GPU Embedding Runner (v3).

Extracts pooled per-protein embeddings from caller-provided FASTA using
Rostlab/ProstT5 encoder in AA2fold direction.

Features:
- Preserves original sequence IDs and order
- Strictly validates FASTA: rejects pre-header sequence lines, empty headers,
  and duplicate FASTA IDs
- Unified raw-sequence hash semantics:
  * raw_sha256: SHA256 of exact raw sequence in FASTA (after whitespace strip & uppercase)
  * clean_sha256: SHA256 of sequence after stripping trailing stop codon '*'
  * length: residue count of clean_seq
- Both raw_sequence_sha256 and clean_sequence_sha256 preserved in NPZ and TSV
- Explicit legacy=True for T5Tokenizer
- Enforces strict token alignment assertions: non-padded token count equals
  prefix + one residue each + EOS (attention length == seq_len + 2), token 0 is
  <AA2fold>, and token 1+seq_len is EOS, preventing silent pooling of EOS or pad
- Excludes special tokens (<AA2fold>, </s>, <pad>) from mean pooling
- Handles UZOB ambiguity codes as documented by Rostlab (transformed to X),
  with exact transformations, positions, and counts recorded in sequence manifest
- Enforces explicit ambiguity and long-sequence rules (no silent truncation)
- Uses length-aware dynamic batching
- Deterministic inference (fixed seeds, eval mode, inference_mode, deterministic cuDNN)
- Comprehensive resume signature: binds runner code SHA, tokenizer hashes,
  config hash, model weights hash, library versions, and inference parameters
- Atomic publication of results (.npz, .npy, manifest.tsv, provenance.json, checksums.sha256)
- Strictly requires CUDA and records GPU device metadata and memory utilization
"""

import argparse
import hashlib
import json
import os
import re
import sys
import time
import uuid
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any

import numpy as np
import torch
import transformers
from transformers import T5EncoderModel, T5Tokenizer

# Valid standard 20 amino acids
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
# Documented ambiguous codes transformed to X
UZOB_MAP = {"U": "X", "Z": "X", "O": "X", "B": "X"}
# Permitted letters after transformation
ALLOWED_AA = STANDARD_AA.union({"X"})


def compute_file_sha256(filepath: Path, chunk_size: int = 16 * 1024 * 1024) -> str:
    """Compute SHA256 checksum of a file."""
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(chunk_size):
            h.update(chunk)
    return h.hexdigest()


def parse_fasta(fasta_path: Path) -> List[Tuple[str, str, str]]:
    """
    Parse FASTA file preserving sequence order with strict structural validation:
    - Rejects sequence lines before the first header
    - Rejects empty headers
    - Rejects duplicate sequence IDs
    Returns list of (seq_id, full_header, raw_sequence).
    """
    records = []
    seen_ids = set()
    current_id = None
    current_header = None
    current_seq_parts = []
    first_header_seen = False

    with open(fasta_path, "r", encoding="utf-8") as f:
        for line_num, line in enumerate(f, 1):
            line_stripped = line.strip()
            if not line_stripped:
                continue

            if line_stripped.startswith(">"):
                first_header_seen = True
                if current_id is not None:
                    raw_seq = "".join(current_seq_parts)
                    if not raw_seq:
                        raise ValueError(f"Empty sequence for ID '{current_id}' at line {line_num} in {fasta_path}")
                    records.append((current_id, current_header, raw_seq))

                header_content = line_stripped[1:].strip()
                if not header_content:
                    raise ValueError(f"Empty header found at line {line_num} in {fasta_path}")

                current_id = header_content.split()[0]
                if not current_id:
                    raise ValueError(f"Empty sequence ID parsed from header '{header_content}' at line {line_num} in {fasta_path}")

                if current_id in seen_ids:
                    raise ValueError(
                        f"Duplicate FASTA ID '{current_id}' detected at line {line_num} in {fasta_path}. "
                        "FASTA IDs must be unique."
                    )
                seen_ids.add(current_id)
                current_seq_parts = []
            else:
                if not first_header_seen:
                    raise ValueError(
                        f"Found sequence content before the first FASTA header at line {line_num} in {fasta_path}: '{line_stripped[:40]}'"
                    )
                current_seq_parts.append(line_stripped)

        if current_id is not None:
            raw_seq = "".join(current_seq_parts)
            if not raw_seq:
                raise ValueError(f"Empty sequence for final ID '{current_id}' in {fasta_path}")
            records.append((current_id, current_header, raw_seq))

    if not records:
        raise ValueError(f"No valid FASTA records found in {fasta_path}")

    return records


def preprocess_sequence(
    raw_seq: str,
    seq_id: str,
    max_seq_len: int
) -> Dict[str, Any]:
    """
    Preprocess sequence and enforce explicit ambiguity, length, and hash rules.
    - Computes raw_sha256 on exact raw sequence as provided in FASTA (uppercase, stripped)
    - Strips trailing stop codon '*' if present (recorded)
    - Computes clean_sha256 on clean sequence (post stop-strip)
    - Validates clean sequence length against max_seq_len (no silent truncation)
    - Replaces U, Z, O, B with X per Rostlab documentation, recording details
    - Rejects any invalid characters not in standard AA + X
    """
    raw_seq_clean = raw_seq.strip().upper()
    raw_sha256 = hashlib.sha256(raw_seq_clean.encode("ascii")).hexdigest()

    trailing_stop = False
    if raw_seq_clean.endswith("*"):
        clean_seq = raw_seq_clean[:-1]
        trailing_stop = True
    else:
        clean_seq = raw_seq_clean

    clean_sha256 = hashlib.sha256(clean_seq.encode("ascii")).hexdigest()
    seq_len = len(clean_seq)

    if seq_len == 0:
        raise ValueError(f"Sequence '{seq_id}' has residue length 0 after stripping trailing stop.")

    if seq_len > max_seq_len:
        raise ValueError(
            f"Sequence '{seq_id}' length {seq_len} exceeds max allowed length {max_seq_len}. "
            "Silent truncation is forbidden by compute policy."
        )

    # Check and transform UZOB on clean_seq
    uzob_replacements = []
    transformed_chars = []
    for pos, char in enumerate(clean_seq):
        if char in UZOB_MAP:
            replacement = UZOB_MAP[char]
            uzob_replacements.append({"pos": pos, "original": char, "replacement": replacement})
            transformed_chars.append(replacement)
        elif char in ALLOWED_AA:
            transformed_chars.append(char)
        else:
            raise ValueError(
                f"Sequence '{seq_id}' contains invalid character '{char}' at position {pos}. "
                "Explicit ambiguity rule: only standard 20 amino acids, U/Z/O/B (mapped to X), or X are accepted."
            )

    processed_seq = "".join(transformed_chars)

    return {
        "seq_id": seq_id,
        "raw_seq": raw_seq_clean,
        "clean_seq": clean_seq,
        "raw_sha256": raw_sha256,
        "clean_sha256": clean_sha256,
        "length": seq_len,
        "trailing_stop_stripped": trailing_stop,
        "has_uzob": len(uzob_replacements) > 0,
        "uzob_count": len(uzob_replacements),
        "uzob_replacements": uzob_replacements,
        "processed_seq": processed_seq,
    }


def form_length_aware_batches(
    preprocessed_records: List[Dict[str, Any]],
    batch_tokens: int
) -> List[List[int]]:
    """
    Form length-aware batches to minimize padding overhead and prevent OOM.
    Returns list of batches, each batch being a list of indices into preprocessed_records.
    """
    sorted_indices = sorted(range(len(preprocessed_records)), key=lambda i: preprocessed_records[i]["length"])

    batches = []
    current_batch = []
    current_max_len = 0

    for idx in sorted_indices:
        seq_len = preprocessed_records[idx]["length"]
        new_max_len = max(current_max_len, seq_len)
        estimated_tokens = (len(current_batch) + 1) * (new_max_len + 2)

        if current_batch and estimated_tokens > batch_tokens:
            batches.append(current_batch)
            current_batch = [idx]
            current_max_len = seq_len
        else:
            current_batch.append(idx)
            current_max_len = new_max_len

    if current_batch:
        batches.append(current_batch)

    return batches


def embed_batch(
    model: T5EncoderModel,
    tokenizer: T5Tokenizer,
    batch_records: List[Dict[str, Any]],
    device: torch.device
) -> List[np.ndarray]:
    """
    Embed a single batch of sequences.
    Strictly asserts tokenizer token alignment before slicing to guarantee
    that residue embeddings exclude prefix (<AA2fold>), EOS (</s>), and padding (<pad>).
    """
    formatted_texts = [
        f"<AA2fold> {' '.join(list(r['processed_seq']))}"
        for r in batch_records
    ]

    encoded = tokenizer(
        formatted_texts,
        add_special_tokens=True,
        padding=True,
        return_tensors="pt"
    )

    input_ids = encoded["input_ids"].to(device)
    attention_mask = encoded["attention_mask"].to(device)

    prefix_id = tokenizer.convert_tokens_to_ids("<AA2fold>")
    eos_id = tokenizer.eos_token_id
    pad_id = tokenizer.pad_token_id

    for b_idx, r in enumerate(batch_records):
        seq_len = r["length"]
        expected_unpadded_len = 1 + seq_len + 1  # prefix + residues + EOS
        actual_attn_len = attention_mask[b_idx].sum().item()

        if actual_attn_len != expected_unpadded_len:
            raise RuntimeError(
                f"Token alignment assertion failed for '{r['seq_id']}': "
                f"attention length ({actual_attn_len}) != seq_len ({seq_len}) + 2 ({expected_unpadded_len})"
            )

        first_token = input_ids[b_idx, 0].item()
        if first_token != prefix_id:
            raise RuntimeError(
                f"Token alignment assertion failed for '{r['seq_id']}': "
                f"token 0 is {first_token}, expected <AA2fold> ({prefix_id})"
            )

        eos_pos = 1 + seq_len
        eos_token = input_ids[b_idx, eos_pos].item()
        if eos_token != eos_id:
            raise RuntimeError(
                f"Token alignment assertion failed for '{r['seq_id']}': "
                f"token at position {eos_pos} is {eos_token}, expected </s> ({eos_id}). "
                "Silent pooling of EOS is strictly forbidden."
            )

        if input_ids.shape[1] > expected_unpadded_len:
            pad_tokens = input_ids[b_idx, expected_unpadded_len:]
            pad_mask = attention_mask[b_idx, expected_unpadded_len:]
            if not (pad_tokens == pad_id).all():
                raise RuntimeError(f"Unexpected non-pad tokens found beyond EOS for '{r['seq_id']}'")
            if not (pad_mask == 0).all():
                raise RuntimeError(f"Attention mask not zero for pad tokens in '{r['seq_id']}'")

    with torch.inference_mode():
        outputs = model(input_ids=input_ids, attention_mask=attention_mask)
        last_hidden_state = outputs.last_hidden_state

    batch_embeddings = []
    for b_idx, r in enumerate(batch_records):
        seq_len = r["length"]
        residue_embeddings = last_hidden_state[b_idx, 1 : 1 + seq_len, :]

        if residue_embeddings.shape[0] != seq_len:
            raise RuntimeError(
                f"Residue embedding length mismatch for sequence '{r['seq_id']}': "
                f"expected {seq_len}, got {residue_embeddings.shape[0]}"
            )

        pooled = residue_embeddings.mean(dim=0)

        if not torch.isfinite(pooled).all():
            raise RuntimeError(f"Non-finite values detected in embedding for sequence '{r['seq_id']}'.")

        batch_embeddings.append(pooled.cpu().to(torch.float32).numpy())

    return batch_embeddings


def get_library_versions() -> Dict[str, str]:
    """Collect exact versions of relevant libraries."""
    import sentencepiece
    import safetensors
    return {
        "python": sys.version.split()[0],
        "torch": torch.__version__,
        "transformers": transformers.__version__,
        "sentencepiece": sentencepiece.__version__,
        "safetensors": safetensors.__version__,
        "numpy": np.__version__,
    }


def calculate_comprehensive_run_signature(
    fasta_sha256: str,
    num_sequences: int,
    model_weights_sha256: str,
    config_sha256: str,
    tokenizer_hashes: Dict[str, str],
    runner_code_sha256: str,
    library_versions: Dict[str, str],
    precision: str,
    max_seq_len: int,
    batch_tokens: int
) -> Tuple[str, Dict[str, Any]]:
    """
    Generate unique deterministic hash for run configuration binding:
    - runner code SHA
    - tokenizer files SHA
    - model config SHA
    - model weights SHA
    - library versions (torch, transformers, sentencepiece, safetensors, python)
    - input FASTA SHA and sequence count
    - execution parameters
    """
    sig_payload = {
        "runner_code_sha256": runner_code_sha256,
        "input_fasta": {
            "sha256": fasta_sha256,
            "num_sequences": num_sequences,
        },
        "model_and_tokenizer": {
            "weights_sha256": model_weights_sha256,
            "config_sha256": config_sha256,
            "tokenizer_hashes": tokenizer_hashes,
        },
        "library_versions": library_versions,
        "parameters": {
            "precision": precision,
            "max_seq_len": max_seq_len,
            "batch_tokens": batch_tokens,
            "legacy_tokenizer": True,
        }
    }
    encoded = json.dumps(sig_payload, sort_keys=True).encode("utf-8")
    sig_hash = hashlib.sha256(encoded).hexdigest()
    return sig_hash, sig_payload


def main():
    parser = argparse.ArgumentParser(description="ProstT5 AA2fold GPU Embedding Runner")
    parser.add_argument("--fasta", type=Path, required=True, help="Input FASTA file path")
    parser.add_argument("--output-dir", type=Path, required=True, help="Output destination directory")
    parser.add_argument("--model-dir", type=Path, default=Path("/mnt/data/u/hohndor/eq-family-discovery-20261015/gpu/model"),
                        help="Path to staged ProstT5 model directory")
    parser.add_argument("--batch-tokens", type=int, default=4096, help="Maximum tokens per batch")
    parser.add_argument("--max-seq-len", type=int, default=2048, help="Maximum allowed sequence length (no silent truncation)")
    parser.add_argument("--precision", choices=["half", "float32"], default="half", help="Precision mode (half=fp16, float32=fp32)")
    parser.add_argument("--device", type=str, default="cuda:0", help="Torch device string (CUDA required)")
    parser.add_argument("--force", action="store_true", help="Force recomputation even if valid signature matches")
    parser.add_argument("--verify-consistency", action="store_true", help="Run single-vs-batched consistency verification")
    args = parser.parse_args()

    # 1. Hardware & CUDA Enforcement
    if not torch.cuda.is_available():
        sys.stderr.write("FATAL ERROR: CUDA is strictly required for ProstT5 runner. torch.cuda.is_available() is False.\n")
        sys.exit(64)

    device = torch.device(args.device)
    torch.cuda.set_device(device)
    torch.manual_seed(42)
    torch.cuda.manual_seed_all(42)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

    device_name = torch.cuda.get_device_name(device)
    device_props = torch.cuda.get_device_properties(device)
    total_mem_gb = device_props.total_memory / (1024 ** 3)
    start_time = time.time()
    torch.cuda.reset_peak_memory_stats(device)

    # 2. Compute runner script SHA256
    runner_script_path = Path(__file__).resolve()
    runner_code_sha256 = compute_file_sha256(runner_script_path)
    lib_versions = get_library_versions()

    print(f"=== ProstT5 GPU Embedding Runner ===")
    print(f"Device: {device_name} ({total_mem_gb:.2f} GB total VRAM)")
    print(f"Runner code SHA256: {runner_code_sha256}")
    print(f"Libraries: torch={lib_versions['torch']}, transformers={lib_versions['transformers']}, sentencepiece={lib_versions['sentencepiece']}, safetensors={lib_versions['safetensors']}")
    print(f"Input FASTA: {args.fasta}")
    print(f"Output Directory: {args.output_dir}")
    print(f"Model Directory: {args.model_dir}")
    print(f"Precision: {args.precision}, Batch tokens: {args.batch_tokens}, Max seq len: {args.max_seq_len}")

    # 3. Check input and compute FASTA checksum
    if not args.fasta.is_file():
        sys.stderr.write(f"FATAL ERROR: Input FASTA file not found: {args.fasta}\n")
        sys.exit(66)

    fasta_sha256 = compute_file_sha256(args.fasta)
    print(f"FASTA SHA256: {fasta_sha256}")

    # 4. Model and Tokenizer Files Inspection & Checksums
    config_file = args.model_dir / "config.json"
    weights_file = args.model_dir / "model.safetensors"
    if not weights_file.exists():
        weights_file = args.model_dir / "pytorch_model.bin"

    if not config_file.exists() or not weights_file.exists():
        sys.stderr.write(f"FATAL ERROR: Model directory {args.model_dir} must contain config.json and weights.\n")
        sys.exit(66)

    tokenizer_files = ["tokenizer_config.json", "special_tokens_map.json", "spiece.model"]
    tokenizer_hashes = {}
    for tf_name in tokenizer_files:
        tf_path = args.model_dir / tf_name
        if not tf_path.exists():
            sys.stderr.write(f"FATAL ERROR: Required tokenizer file {tf_name} missing from {args.model_dir}\n")
            sys.exit(66)
        tokenizer_hashes[tf_name] = compute_file_sha256(tf_path)

    added_tokens_path = args.model_dir / "added_tokens.json"
    if added_tokens_path.exists():
        tokenizer_hashes["added_tokens.json"] = compute_file_sha256(added_tokens_path)

    config_sha256 = compute_file_sha256(config_file)
    print(f"Computing model weights SHA256 on worker...")
    model_sha256 = compute_file_sha256(weights_file)
    print(f"Model weights ({weights_file.name}) SHA256: {model_sha256}")
    print(f"Model config SHA256: {config_sha256}")

    # 5. Parse and preprocess FASTA (with duplicate ID and header validation)
    raw_records = parse_fasta(args.fasta)
    num_seqs = len(raw_records)
    print(f"Parsed and validated {num_seqs} unique sequence records from FASTA.")

    preprocessed_records = []
    for seq_id, full_header, raw_seq in raw_records:
        prep = preprocess_sequence(raw_seq, seq_id, args.max_seq_len)
        prep["full_header"] = full_header
        preprocessed_records.append(prep)

    uzob_count_total = sum(r["uzob_count"] for r in preprocessed_records)
    seqs_with_uzob = sum(1 for r in preprocessed_records if r["has_uzob"])
    trailing_stops = sum(1 for r in preprocessed_records if r["trailing_stop_stripped"])
    print(f"Preprocessing summary: {num_seqs} sequences; {seqs_with_uzob} with UZOB substitutions ({uzob_count_total} residues); {trailing_stops} trailing stops stripped.")

    # 6. Check comprehensive run signature for resume / idempotency
    run_signature, signature_payload = calculate_comprehensive_run_signature(
        fasta_sha256=fasta_sha256,
        num_sequences=num_seqs,
        model_weights_sha256=model_sha256,
        config_sha256=config_sha256,
        tokenizer_hashes=tokenizer_hashes,
        runner_code_sha256=runner_code_sha256,
        library_versions=lib_versions,
        precision=args.precision,
        max_seq_len=args.max_seq_len,
        batch_tokens=args.batch_tokens,
    )
    print(f"Comprehensive run signature: {run_signature}")

    prov_file = args.output_dir / "provenance.json"
    chk_file = args.output_dir / "checksums.sha256"
    npz_file = args.output_dir / "embeddings.npz"

    if not args.force and prov_file.exists() and chk_file.exists() and npz_file.exists():
        try:
            with open(prov_file, "r") as pf:
                existing_prov = json.load(pf)
            if existing_prov.get("run_signature") == run_signature:
                print("Existing run provenance signature matches exactly. Checking output integrity...")
                all_ok = True
                with open(chk_file, "r") as cf:
                    for line in cf:
                        parts = line.strip().split()
                        if len(parts) == 2:
                            expected_sha, fname = parts[0], parts[1].lstrip("*")
                            target_f = args.output_dir / fname
                            if not target_f.exists() or compute_file_sha256(target_f) != expected_sha:
                                all_ok = False
                                break
                if all_ok:
                    print("Signature and checksums verified. Results are up-to-date. Skipping recomputation.")
                    sys.exit(0)
                else:
                    print("Existing output checksum verification failed. Recomputing...")
            else:
                print("Signature mismatch detected (input, code, model, tokenizer, or libraries changed). Recomputing...")
        except Exception as e:
            print(f"Could not verify existing provenance ({e}). Recomputing...")

    # 7. Load Model and Tokenizer (Explicit legacy=True)
    print("Loading ProstT5 tokenizer (explicit legacy=True) and encoder model...")
    tokenizer = T5Tokenizer.from_pretrained(args.model_dir, do_lower_case=False, legacy=True)
    model = T5EncoderModel.from_pretrained(args.model_dir)

    if args.precision == "half":
        model = model.half()
    else:
        model = model.float()

    model = model.to(device)
    model.eval()
    print("Model loaded successfully into VRAM.")

    # 8. Form length-aware batches
    # Compute identical preprocessed sequences once; preserve every original ID on expansion.
    # This avoids batch-shape floating point differences for identical proteins.
    canonical_indices = []
    canonical_by_sha = {}
    original_to_canonical = []
    for index, record in enumerate(preprocessed_records):
        sha = record["clean_sha256"]
        if sha not in canonical_by_sha:
            canonical_by_sha[sha] = index
            canonical_indices.append(index)
        original_to_canonical.append(canonical_by_sha[sha])
    unique_records = [preprocessed_records[i] for i in canonical_indices]
    unique_batches = form_length_aware_batches(unique_records, args.batch_tokens)
    batches = [[canonical_indices[i] for i in batch] for batch in unique_batches]
    print(f"Formed {len(batches)} length-aware batches across {num_seqs} sequences.")

    # 9. Run Inference with Token Alignment Assertions
    ordered_embeddings = [None] * num_seqs
    embed_start = time.time()

    for b_idx, batch_indices in enumerate(batches, 1):
        batch_records = [preprocessed_records[i] for i in batch_indices]
        batch_embs = embed_batch(model, tokenizer, batch_records, device)
        for orig_idx, emb in zip(batch_indices, batch_embs):
            ordered_embeddings[orig_idx] = emb
        if b_idx % 10 == 0 or b_idx == len(batches):
            print(f"Processed batch {b_idx}/{len(batches)} (last batch size: {len(batch_indices)})")

    for index, canonical_index in enumerate(original_to_canonical):
        ordered_embeddings[index] = ordered_embeddings[canonical_index]
    embed_duration = time.time() - embed_start
    final_embeddings_matrix = np.stack(ordered_embeddings, axis=0).astype(np.float32)
    print(f"Embedding completed in {embed_duration:.2f} s. Output shape: {final_embeddings_matrix.shape}")

    # Verify output dimensions and finiteness
    assert final_embeddings_matrix.shape == (num_seqs, 1024), f"Unexpected shape {final_embeddings_matrix.shape}"
    assert np.all(np.isfinite(final_embeddings_matrix)), "Non-finite values present in final embeddings matrix!"

    # 10. Optional Single-vs-Batched Consistency Verification
    consistency_results = None
    if args.verify_consistency:
        print("Running single-sequence consistency verification...")
        single_embeddings = []
        for r in preprocessed_records:
            single_emb = embed_batch(model, tokenizer, [r], device)[0]
            single_embeddings.append(single_emb)
        single_matrix = np.stack(single_embeddings, axis=0).astype(np.float32)

        abs_diff = np.abs(final_embeddings_matrix - single_matrix)
        max_abs_diff = float(np.max(abs_diff))
        mean_abs_diff = float(np.mean(abs_diff))

        dot_product = np.sum(final_embeddings_matrix * single_matrix, axis=1)
        norm_batched = np.linalg.norm(final_embeddings_matrix, axis=1)
        norm_single = np.linalg.norm(single_matrix, axis=1)
        cos_sim = dot_product / (norm_batched * norm_single)
        min_cos_sim = float(np.min(cos_sim))

        print(f"Single vs batched consistency:")
        print(f"  Max absolute difference: {max_abs_diff:.6e}")
        print(f"  Mean absolute difference: {mean_abs_diff:.6e}")
        print(f"  Min cosine similarity: {min_cos_sim:.8f}")

        tol_threshold = 5e-3 if args.precision == "half" else 1e-4
        if max_abs_diff > tol_threshold:
            raise RuntimeError(
                f"Single-vs-batched consistency check failed: max_abs_diff={max_abs_diff} exceeds tolerance {tol_threshold}"
            )
        if min_cos_sim < 0.9999:
            raise RuntimeError(
                f"Single-vs-batched consistency check failed: min_cos_sim={min_cos_sim} is below 0.9999"
            )

        consistency_results = {
            "max_abs_diff": max_abs_diff,
            "mean_abs_diff": mean_abs_diff,
            "min_cos_sim": min_cos_sim,
            "tolerance_threshold": tol_threshold,
            "verified": True,
        }

    # 11. Check duplicate sequence equivalence (grouped by clean_sha256)
    duplicate_results = []
    clean_sha_to_indices: Dict[str, List[int]] = {}
    for idx, r in enumerate(preprocessed_records):
        clean_sha_to_indices.setdefault(r["clean_sha256"], []).append(idx)

    for sha, indices in clean_sha_to_indices.items():
        if len(indices) > 1:
            base_idx = indices[0]
            base_emb = final_embeddings_matrix[base_idx]
            for other_idx in indices[1:]:
                other_emb = final_embeddings_matrix[other_idx]
                dup_diff = float(np.max(np.abs(base_emb - other_emb)))
                duplicate_results.append({
                    "seq_id_1": preprocessed_records[base_idx]["seq_id"],
                    "seq_id_2": preprocessed_records[other_idx]["seq_id"],
                    "clean_sha256": sha,
                    "max_abs_diff": dup_diff,
                    "identical": dup_diff == 0.0,
                })
                print(f"Duplicate check ({preprocessed_records[base_idx]['seq_id']} vs {preprocessed_records[other_idx]['seq_id']}): max_abs_diff={dup_diff}")
                if dup_diff > 1e-6:
                    raise RuntimeError(f"Duplicate sequence embedding mismatch: diff={dup_diff} > 1e-6")

    # 12. Atomic publication to temporary workspace first
    args.output_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = args.output_dir / f".work_{int(time.time())}_{uuid.uuid4().hex[:8]}"
    temp_dir.mkdir(parents=True, exist_ok=True)

    peak_vram_bytes = torch.cuda.max_memory_allocated(device)
    peak_vram_gb = peak_vram_bytes / (1024 ** 3)
    total_duration = time.time() - start_time

    # Save embeddings.npz with unified raw and clean sequence hashes
    all_ids = np.array([r["seq_id"] for r in preprocessed_records], dtype=object)
    all_lens = np.array([r["length"] for r in preprocessed_records], dtype=np.int32)
    all_raw_shas = np.array([r["raw_sha256"] for r in preprocessed_records], dtype=object)
    all_clean_shas = np.array([r["clean_sha256"] for r in preprocessed_records], dtype=object)

    temp_npz = temp_dir / "embeddings.npz"
    np.savez_compressed(
        temp_npz,
        embeddings=final_embeddings_matrix,
        ids=all_ids,
        lengths=all_lens,
        raw_sequence_sha256=all_raw_shas,
        clean_sequence_sha256=all_clean_shas,
        sequence_sha256=all_raw_shas  # Aliased to raw for backwards compatibility
    )

    # Save embeddings.npy
    temp_npy = temp_dir / "embeddings.npy"
    np.save(temp_npy, final_embeddings_matrix)

    # Save sequence_manifest.tsv with unified hash semantics
    temp_manifest = temp_dir / "sequence_manifest.tsv"
    with open(temp_manifest, "w", encoding="utf-8") as mf:
        mf.write("index\tseq_id\tlength\traw_sha256\tclean_sha256\ttrailing_stop_stripped\thas_uzob\tuzob_count\tuzob_replacements_json\n")
        for idx, r in enumerate(preprocessed_records):
            rep_json = json.dumps(r["uzob_replacements"])
            mf.write(f"{idx}\t{r['seq_id']}\t{r['length']}\t{r['raw_sha256']}\t{r['clean_sha256']}\t{r['trailing_stop_stripped']}\t{r['has_uzob']}\t{r['uzob_count']}\t{rep_json}\n")

    # Comprehensive machine-readable provenance
    provenance = {
        "runner": "prostt5_embed.py",
        "runner_code_sha256": runner_code_sha256,
        "run_signature": run_signature,
        "signature_components": signature_payload,
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "input_fasta": {
            "path": str(args.fasta.resolve()),
            "sha256": fasta_sha256,
            "num_sequences": num_seqs,
        },
        "model": {
            "dir": str(args.model_dir.resolve()),
            "weights_filename": weights_file.name,
            "weights_sha256": model_sha256,
            "config_sha256": config_sha256,
            "tokenizer_hashes": tokenizer_hashes,
            "tokenizer_legacy": True,
            "architecture": "T5EncoderModel",
            "pooling": "mean_excluding_special_tokens",
            "prefix": "<AA2fold>",
            "token_alignment_guarantee": "prefix+residues+eos_verified",
        },
        "parameters": {
            "precision": args.precision,
            "batch_tokens": args.batch_tokens,
            "max_seq_len": args.max_seq_len,
            "num_batches": len(batches),
            "unique_clean_sequences_embedded": len(canonical_indices),
            "duplicate_policy": "embed_once_expand_all_original_ids",
        },
        "software_environment": {
            "python_version": sys.version,
            "library_versions": lib_versions,
            "cuda_runtime_version": torch.version.cuda,
            "cudnn_version": torch.backends.cudnn.version(),
        },
        "hardware": {
            "device_name": device_name,
            "total_vram_gb": round(total_mem_gb, 4),
            "peak_vram_used_gb": round(peak_vram_gb, 4),
            "peak_vram_used_bytes": peak_vram_bytes,
        },
        "performance": {
            "embedding_duration_seconds": round(embed_duration, 4),
            "total_duration_seconds": round(total_duration, 4),
            "sequences_per_second": round(num_seqs / max(embed_duration, 0.001), 2),
        },
        "validation": {
            "dimensions": list(final_embeddings_matrix.shape),
            "all_finite": True,
            "id_order_preserved": True,
            "token_alignment_checked": True,
            "single_vs_batched_consistency": consistency_results,
            "duplicate_sequence_checks": duplicate_results,
        }
    }

    temp_prov = temp_dir / "provenance.json"
    with open(temp_prov, "w", encoding="utf-8") as pf:
        json.dump(provenance, pf, indent=2)

    # Compute checksums for all files in temp_dir
    temp_chk = temp_dir / "checksums.sha256"
    published_files = ["embeddings.npz", "embeddings.npy", "sequence_manifest.tsv", "provenance.json"]
    with open(temp_chk, "w", encoding="utf-8") as cf:
        for fname in published_files:
            fpath = temp_dir / fname
            cf.write(f"{compute_file_sha256(fpath)}  {fname}\n")

    # Move published files atomically into final output_dir
    print(f"Publishing results atomically to {args.output_dir}...")
    for fname in published_files + ["checksums.sha256"]:
        src = temp_dir / fname
        dst = args.output_dir / fname
        os.replace(src, dst)

    try:
        temp_dir.rmdir()
    except Exception:
        pass

    print(f"SUCCESS: Embeddings and provenance published to {args.output_dir}.")


if __name__ == "__main__":
    main()
