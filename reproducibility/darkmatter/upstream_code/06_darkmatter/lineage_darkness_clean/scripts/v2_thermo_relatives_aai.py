#!/usr/bin/env python3
"""v2_thermo_relatives_aai - reciprocal-best-hit AAI for the desert Thermoplasmatota + comparators.

No CompareM or EzAAI exists in any environment here, so AAI is computed directly, in the CompareM
style: proteins predicted with Pyrodigal 3.7.1 (single-genome mode, trained per genome) on every
genome alike (MGnify's own .faa files are not used), one DIAMOND 2.1.16 blastp all-against-all
(--very-sensitive, e-value 1e-5, identity >= 30 %, query and subject coverage >= 70 %, all targets
kept), best hit per (query protein, target genome) by bitscore, reciprocal best hits between each
genome pair, AAI = mean identity over the reciprocal pairs (averaged over both directions).

Subcommands (called from v2_thermo_relatives.sbatch; run on compute nodes, never on login):
  proteins <threads>           tmp/thermo_relatives/faa/all.faa, genome_index.tsv, chunk_*.faa
  reduce <chunk>               stdin = DIAMOND outfmt 6 (qseqid sseqid pident bitscore);
                               writes tmp/thermo_relatives/hits/chunk_<chunk>.npz
  aai                          results/v2_thermo_relatives_aai.tsv.gz (all genome pairs with RBH)

Protein ids are "<genome index>_<n>", so a protein's genome is its prefix and its integer key is
index * 100000 + n (no genome here has 100,000 proteins; asserted).
"""
import os, sys, csv, gzip, glob
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
P = os.path.abspath(os.path.join(HERE, ".."))
RES = os.path.join(P, "results")
TMP = os.path.join(P, "tmp", "thermo_relatives")
FAA = os.path.join(TMP, "faa")
HITS = os.path.join(TMP, "hits")
N_CHUNKS = 20
K = 100000


def read_fasta(path):
    op = gzip.open if path.endswith(".gz") else open
    name, seq = None, []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if name:
        yield name, "".join(seq)


def _predict(args):
    import pyrodigal
    gi, path = args
    seqs = [s.upper() for _, s in read_fasta(path)]
    gf = pyrodigal.GeneFinder(meta=False)
    gf.train(*[s.encode() for s in seqs])
    out, n = [], 0
    for s in seqs:
        for g in gf.find_genes(s.encode()):
            n += 1
            out.append(f">{gi}_{n}\n{g.translate(include_stop=False)}\n")
    assert n < K, (path, n)
    return gi, n, "".join(out)


def proteins(threads):
    from multiprocessing import Pool
    os.makedirs(FAA, exist_ok=True)
    rows = list(csv.DictReader(open(os.path.join(TMP, "genome_paths.tsv")), delimiter="\t"))
    jobs = [(i, r["path"]) for i, r in enumerate(rows)]
    with Pool(threads) as pool:
        res = sorted(pool.map(_predict, jobs, chunksize=4))
    with open(os.path.join(FAA, "all.faa"), "w") as o:
        for _, _, txt in res:
            o.write(txt)
    per = int(np.ceil(len(res) / N_CHUNKS))
    for c in range(N_CHUNKS):
        with open(os.path.join(FAA, f"chunk_{c:02d}.faa"), "w") as o:
            for _, _, txt in res[c * per:(c + 1) * per]:
                o.write(txt)
    with open(os.path.join(TMP, "genome_index.tsv"), "w") as o:
        o.write("gidx\tgenome_id\tset\tclass\tn_proteins\n")
        for (gi, n, _), r in zip(res, rows):
            o.write(f"{gi}\t{r['genome_id']}\t{r['set']}\t{r['class']}\t{n}\n")
    print(f"genomes {len(res)}, proteins {sum(n for _, n, _ in res)}, chunks {N_CHUNKS}")


def reduce(chunk):
    """Stream DIAMOND hits (grouped by query) -> best hit per (query, target genome)."""
    os.makedirs(HITS, exist_ok=True)
    Q, S, I = [], [], []
    cur, best, seen = None, {}, set()

    def flush():
        for sg, (bs, s, pid) in best.items():
            Q.append(cur); S.append(s); I.append(pid)

    n = 0
    for line in sys.stdin:
        q, s, pid, bs = line.split("\t")
        qg, qn = q.split("_"); sg, sn = s.split("_")
        if qg == sg:
            continue
        qk = int(qg) * K + int(qn)
        if qk != cur:
            if cur is not None:
                flush()
            assert qk not in seen, f"query {q} not contiguous in DIAMOND output"
            seen.add(qk)
            cur, best = qk, {}
        bs = float(bs)
        if sg not in best or bs > best[sg][0]:
            best[sg] = (bs, int(sg) * K + int(sn), float(pid))
        n += 1
    if cur is not None:
        flush()
    np.savez(os.path.join(HITS, f"chunk_{chunk}.npz"), q=np.array(Q, dtype=np.int64),
             s=np.array(S, dtype=np.int64), pid=np.array(I, dtype=np.float32))
    print(f"chunk {chunk}: hit lines {n}, best-per-genome rows {len(Q)}")


def aai():
    files = sorted(glob.glob(os.path.join(HITS, "chunk_*.npz")))
    assert len(files) == N_CHUNKS, f"{len(files)} hit chunks, expected {N_CHUNKS}"
    q = np.concatenate([np.load(f)["q"] for f in files])
    s = np.concatenate([np.load(f)["s"] for f in files])
    pid = np.concatenate([np.load(f)["pid"] for f in files])
    print(f"best-hit rows {len(q)}", flush=True)
    M = np.int64(K * 1000)
    fwd = q * M + s
    rev = np.sort(s * M + q)
    pos = np.searchsorted(rev, fwd)
    pos[pos == len(rev)] = 0
    rbh = rev[pos] == fwd
    qg, sg = (q // K)[rbh], (s // K)[rbh]
    G =int(max(q.max(), s.max()) // K) + 1
    pair = qg * G + sg
    cnt = np.bincount(pair, minlength=G * G).reshape(G, G)
    tot = np.bincount(pair, weights=pid[rbh], minlength=G * G).reshape(G, G)
    gi = list(csv.DictReader(open(os.path.join(TMP, "genome_index.tsv")), delimiter="\t"))
    assert len(gi) == G
    names = [r["genome_id"] for r in gi]
    nprot = [int(r["n_proteins"]) for r in gi]
    n_out = 0
    with gzip.open(os.path.join(RES, "v2_thermo_relatives_aai.tsv.gz"), "wt") as o:
        o.write("genome_a\tgenome_b\taai\tn_rbh\tn_prot_a\tn_prot_b\tfrac_rbh_of_smaller\n")
        for a in range(G):
            for b in range(a + 1, G):
                n = cnt[a, b] + cnt[b, a]
                if n == 0:
                    continue
                v = (tot[a, b] + tot[b, a]) / n
                nr = int(round(n / 2))
                o.write(f"{names[a]}\t{names[b]}\t{v:.3f}\t{nr}\t{nprot[a]}\t{nprot[b]}\t"
                        f"{nr / min(nprot[a], nprot[b]):.4f}\n")
                n_out += 1
    print(f"RBH rows {int(rbh.sum())}, genome pairs written {n_out}")


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "proteins":
        proteins(int(sys.argv[2]))
    elif cmd == "reduce":
        reduce(sys.argv[2])
    elif cmd == "aai":
        aai()
