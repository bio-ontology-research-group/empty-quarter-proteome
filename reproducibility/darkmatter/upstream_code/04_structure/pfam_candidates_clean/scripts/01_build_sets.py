#!/usr/bin/env python
"""Build the three protein sets for the Pfam scan of the structural candidates.

Sets
  K_members  : every member protein of the 1,031 candidate fold families K.
               per_rep_structures.tsv is keyed on the MODELLED representative,
               which is fam_reps.tsv `chosen_rep`; the cluster key in
               clusters.tsv is fam_reps.tsv `family`, and the two differ for the
               165 REUSE families. Membership is therefore resolved
               chosen_rep -> family -> members, never chosen_rep -> members.
  QnotK_reps : the 1,807 representatives that are confidence-qualified but not
               dual-novel (confidence_qualified == 1, conf_qualified_dual_novel == 0)
  controls   : the confidence-qualified reference controls (PDB chain crops),
               confidence_qualified == 1 in per_control_structures.tsv

Writes id lists, one FASTA per set and a combined FASTA for a single hmmsearch run.
Read-only on every input. Idempotent: rerunning overwrites its own outputs.
"""
import gzip
import os
import sys

W = "/ibex/scratch/projects/c2014/rob/metagenomics"
OUT = os.path.join(W, "04_structure/pfam_candidates_clean")

PER_REP = os.path.join(W, "04_structure/00_structure_clean/05_validate/per_rep_structures.tsv")
PER_CTRL = os.path.join(W, "04_structure/00_structure_clean/05_validate/per_control_structures.tsv")
FAM_REPS = os.path.join(W, "04_structure/00_structure_clean/05_validate/fam_reps.tsv")
CLUSTERS = os.path.join(W, "04_structure/00_structure_clean/04_cluster/clusters.tsv")
DARK_FAA = os.path.join(W, "06_darkmatter/darkfam_basis_clean/results/dark_reps.faa")
CTRL_FAA = os.path.join(
    W, "rebuild_20260729/structural_calibration/preparation/controls/controls_blinded.faa"
)

RES = os.path.join(OUT, "results")
SEQ = os.path.join(OUT, "sequences")
os.makedirs(SEQ, exist_ok=True)


def read_flagged(path, want):
    """Return the set of ids whose flag columns match `want` (dict col -> value)."""
    keep = set()
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(hdr)}
        for col in want:
            if col not in idx:
                sys.exit("missing column %s in %s" % (col, path))
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if all(f[idx[c]] == v for c, v in want.items()):
                keep.add(f[0])
    return keep


def read_fasta(path, wanted, opener=open):
    """Stream a FASTA, yield (id, seq) for ids in `wanted`. Id = first whitespace token."""
    cur, buf = None, []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None and cur in wanted:
                    yield cur, "".join(buf)
                cur = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if cur is not None and cur in wanted:
        yield cur, "".join(buf)


def write_fasta(path, records):
    with open(path, "w") as fh:
        for sid, seq in records:
            fh.write(">%s\n" % sid)
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")


def main():
    k_reps = read_flagged(PER_REP, {"conf_qualified_dual_novel": "1"})
    q_reps = read_flagged(PER_REP, {"confidence_qualified": "1"})
    qnotk = q_reps - k_reps
    print("K reps            %d" % len(k_reps))
    print("Q reps            %d" % len(q_reps))
    print("Q minus K reps    %d" % len(qnotk))
    assert len(k_reps) == 1031, len(k_reps)
    assert len(q_reps) == 2838, len(q_reps)
    assert len(qnotk) == 1807, len(qnotk)

    # chosen_rep -> family key, and the declared family size
    fam_of = {}
    declared = {}
    with open(FAM_REPS) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rep = f[ix["chosen_rep"]]
            if rep in k_reps:
                fam_of[rep] = f[ix["family"]]
                declared[rep] = int(f[ix["n_members"]])
    print("K reps mapped     %d" % len(fam_of))
    assert len(fam_of) == len(k_reps)
    print("declared members  %d" % sum(declared.values()))
    fam_to_rep = {v: k for k, v in fam_of.items()}
    assert len(fam_to_rep) == len(fam_of)

    # membership: family key -> members, restricted to the K families
    k_members = {}
    n_lines = 0
    with open(CLUSTERS) as fh:
        for line in fh:
            n_lines += 1
            fam, mem = line.rstrip("\n").split("\t")[:2]
            if fam in fam_to_rep:
                k_members.setdefault(fam_to_rep[fam], []).append(mem)
    print("clusters.tsv rows %d" % n_lines)
    print("K families seen   %d" % len(k_members))
    assert len(k_members) == len(k_reps)
    for rep, mm in k_members.items():
        assert len(mm) == declared[rep], (rep, len(mm), declared[rep])
    all_k_members = set()
    for mm in k_members.values():
        all_k_members.update(mm)
    print("K member proteins %d" % len(all_k_members))

    ctrl_ids = read_flagged(PER_CTRL, {"confidence_qualified": "1"})
    print("controls (conf)   %d" % len(ctrl_ids))
    assert len(ctrl_ids) == 3013, len(ctrl_ids)

    overlap = all_k_members & qnotk
    assert not overlap, "K members and Q-not-K reps overlap: %d" % len(overlap)

    # membership table
    with open(os.path.join(RES, "k_family_members.tsv"), "w") as fh:
        fh.write("modelled_rep\tfamily_key\tmember\n")
        for rep in sorted(k_members):
            for mem in sorted(k_members[rep]):
                fh.write("%s\t%s\t%s\n" % (rep, fam_of[rep], mem))

    for name, ids in (
        ("k_reps", k_reps),
        ("k_members", all_k_members),
        ("qnotk_reps", qnotk),
        ("controls", ctrl_ids),
    ):
        with open(os.path.join(RES, "ids_%s.txt" % name), "w") as fh:
            for i in sorted(ids):
                fh.write(i + "\n")

    # sequences: one streaming pass over the 4.29 M dark representatives
    need = all_k_members | qnotk
    got = {}
    opener = gzip.open if DARK_FAA.endswith(".gz") else open
    for sid, seq in read_fasta(DARK_FAA, need, opener):
        got[sid] = seq
    print("dark_reps.faa hits %d of %d requested" % (len(got), len(need)))

    miss_k = sorted(all_k_members - set(got))
    miss_q = sorted(qnotk - set(got))
    with open(os.path.join(RES, "missing_sequences.txt"), "w") as fh:
        for i in miss_k:
            fh.write("k_members\t%s\n" % i)
        for i in miss_q:
            fh.write("qnotk_reps\t%s\n" % i)
    print("missing K members %d, missing Q-not-K reps %d" % (len(miss_k), len(miss_q)))

    ctrl_seqs = dict(read_fasta(CTRL_FAA, ctrl_ids))
    print("control seqs      %d of %d" % (len(ctrl_seqs), len(ctrl_ids)))

    write_fasta(
        os.path.join(SEQ, "k_members.faa"),
        sorted((i, got[i]) for i in all_k_members if i in got),
    )
    write_fasta(
        os.path.join(SEQ, "qnotk_reps.faa"), sorted((i, got[i]) for i in qnotk if i in got)
    )
    write_fasta(os.path.join(SEQ, "controls.faa"), sorted(ctrl_seqs.items()))

    # combined target file: the three sets are disjoint by construction
    combined = os.path.join(SEQ, "all_sets.faa")
    with open(combined, "w") as out:
        for part in ("k_members.faa", "qnotk_reps.faa", "controls.faa"):
            with open(os.path.join(SEQ, part)) as fh:
                out.write(fh.read())

    with open(os.path.join(RES, "set_membership.tsv"), "w") as fh:
        fh.write("protein\tset\tis_family_rep\tlength\n")
        for i in sorted(all_k_members):
            if i in got:
                fh.write(
                    "%s\tK_members\t%d\t%d\n" % (i, 1 if i in k_reps else 0, len(got[i]))
                )
        for i in sorted(qnotk):
            if i in got:
                fh.write("%s\tQnotK_reps\t1\t%d\n" % (i, len(got[i])))
        for i in sorted(ctrl_seqs):
            fh.write("%s\tcontrols\t1\t%d\n" % (i, len(ctrl_seqs[i])))

    n_comb = sum(1 for line in open(combined) if line.startswith(">"))
    print("combined FASTA    %d sequences -> %s" % (n_comb, combined))


if __name__ == "__main__":
    main()
