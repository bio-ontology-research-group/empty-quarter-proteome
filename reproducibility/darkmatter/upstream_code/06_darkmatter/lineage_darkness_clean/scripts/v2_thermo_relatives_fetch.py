#!/usr/bin/env python3
"""v2_thermo_relatives_fetch - comparator genomes for the 237 desert Thermoplasmatota.

Question (with v2_thermo_relatives.py): how many of the 237 desert Thermoplasmatota species have
any known relative, and how close. This script only assembles the comparator set and the genome
path list; the compute is v2_thermo_relatives.sbatch and v2_thermo_relatives_gtdbtk.sbatch.

Comparators
  (a) every GTDB r232 species representative in c__SW-10-69-26 (21) and c__E2 (197), accessions
      from the r232 ar53 reps taxonomy the tree script uses (GB_/RS_ prefix stripped). FASTAs from
      the NCBI FTP (genomes/all/GCA|GCF/nnn/nnn/nnn/<acc>_<asm>/), md5 checked against the
      assembly's md5checksums.txt.
  (b) every MGnify genome whose taxon lineage is in either class (genomes API, taxon_lineage
      filter), FASTA from the genome's downloads endpoint. MGnify publishes no md5; the local md5 is
      recorded. The MGnify lineage comes from the catalogue's own GTDB release, which is why the
      GTDB-Tk r232 step re-classifies these.
  (c) SMAG (Ma et al. 2023, Zenodo 8223844) is NOT fetched: the record holds one genome archive
      split into ~370 parts of 100 MB (mag.tar.gz.aa ... .oi), no per-genome files.

Idempotent: a file already on disk whose md5 matches is not downloaded again.
Downloads go to tmp/thermo_relatives/ and are deleted once the results exist.

Outputs
  results/v2_thermo_relatives_comparators.tsv   manifest: accession, source, catalogue, biome,
                                                lineage, ncbi accession, bytes, md5, md5 check
  tmp/thermo_relatives/genome_paths.tsv         genome_id, set, class, path (237 + comparators)
"""
import os, re, csv, gzip, json, hashlib, time, urllib.request
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
P = os.path.abspath(os.path.join(HERE, ".."))
RES = os.path.join(P, "results")
TMP = os.path.join(P, "tmp", "thermo_relatives")
W = "/ibex/scratch/projects/c2014/rob/metagenomics"
TAX = f"{W}/rebuild_20260729/references/gtdbtk_r232/taxonomy/ar53_taxonomy_r232_reps.tsv"
DREP = f"{W}/05_refine/_drep/out/dereplicated_genomes"
ECOL = os.path.join(RES, "v2_thermo_ecology_by_rep.tsv")
CLASSES = ["SW-10-69-26", "E2"]
NCBI = "https://ftp.ncbi.nlm.nih.gov/genomes/all"
MGN = "https://www.ebi.ac.uk/metagenomics/api/v1"
EXPECT_GTDB = {"SW-10-69-26": 21, "E2": 197}


def get(url, tries=4):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read()
        except Exception as e:
            if i == tries - 1:
                raise
            time.sleep(5 * (i + 1))


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def rank(lin, letter):
    for part in lin.split(";"):
        if part.startswith(letter + "__"):
            return part[3:]
    return ""


# ------------------------------------------------------------------ (a) GTDB r232 species reps
def fetch_ncbi(row):
    acc = row["accession"]
    pre, num = acc[:3], acc[4:].split(".")[0]
    base = f"{NCBI}/{pre}/{num[0:3]}/{num[3:6]}/{num[6:9]}/"
    out = os.path.join(TMP, "gtdb", acc + ".fna.gz")
    try:
        listing = get(base).decode()
        dirs = sorted(set(re.findall(r'href="(%s_[^"/]+)/?"' % re.escape(acc), listing)))
        if not dirs:
            row.update(status="not_on_ncbi_ftp", path="", bytes="", md5="", md5_ok="")
            return row
        d = dirs[0]
        sums = get(f"{base}{d}/md5checksums.txt").decode()
        want = {l.split()[1].lstrip("./"): l.split()[0] for l in sums.splitlines() if l.strip()}
        fname = f"{d}_genomic.fna.gz"
        if not (os.path.exists(out) and md5(out) == want.get(fname)):
            data = get(f"{base}{d}/{fname}")
            with open(out + ".part", "wb") as o:
                o.write(data)
            os.replace(out + ".part", out)
        m = md5(out)
        row.update(asm_dir=d, status="ok", path=out, bytes=os.path.getsize(out), md5=m,
                   md5_ok="yes" if m == want.get(fname) else "NO")
    except Exception as e:
        row.update(status=f"error:{type(e).__name__}", path="", bytes="", md5="", md5_ok="")
    return row


def gtdb_rows():
    rows = []
    with open(TAX) as fh:
        for line in fh:
            gid, lin = line.rstrip("\n").split("\t")[:2]
            c = rank(lin, "c")
            if c in CLASSES:
                rows.append(dict(genome_id=gid[3:], accession=gid[3:], source="GTDB_r232_rep",
                                 catalogue="GTDB r232", biome="", lineage=lin, cls=c,
                                 ncbi_accession=gid[3:], geographic_origin=""))
    for c, n in EXPECT_GTDB.items():
        got = sum(r["cls"] == c for r in rows)
        assert got == n, f"{c}: {got} r232 reps, expected {n}"
    return rows


# ------------------------------------------------------------------ (b) MGnify genomes
def mgnify_rows():
    rows = []
    for c in CLASSES:
        url = f"{MGN}/genomes?taxon_lineage=d__Archaea;p__Thermoplasmatota;c__{c}&page_size=100"
        while url:
            d = json.loads(get(url))
            for x in d["data"]:
                a, rel = x["attributes"], x["relationships"]
                rows.append(dict(
                    genome_id=a["accession"], accession=a["accession"], source="MGnify",
                    catalogue=(rel.get("catalogue", {}).get("data") or {}).get("id", ""),
                    biome=(rel.get("biome", {}).get("data") or {}).get("id", ""),
                    lineage=a["taxon-lineage"], cls=rank(a["taxon-lineage"], "c"),
                    ncbi_accession=a.get("ncbi-genome-accession") or "",
                    geographic_origin=a.get("geographic-origin") or "",
                    completeness=a.get("completeness"), contamination=a.get("contamination"),
                    mgnify_type=a.get("type")))
            url = d["links"].get("next")
    return rows


def fetch_mgnify(row):
    acc = row["accession"]
    out = os.path.join(TMP, "mgnify", acc + ".fna")
    try:
        if not (os.path.exists(out) and os.path.getsize(out) > 0):
            data = get(f"{MGN}/genomes/{acc}/downloads/{acc}.fna")
            assert data.startswith(b">"), "not a FASTA"
            with open(out + ".part", "wb") as o:
                o.write(data)
            os.replace(out + ".part", out)
        row.update(status="ok", path=out, bytes=os.path.getsize(out), md5=md5(out),
                   md5_ok="no_published_md5")
    except Exception as e:
        row.update(status=f"error:{type(e).__name__}", path="", bytes="", md5="", md5_ok="")
    return row


def main():
    for sub in ("gtdb", "mgnify"):
        os.makedirs(os.path.join(TMP, sub), exist_ok=True)

    g = gtdb_rows()
    print(f"GTDB r232 reps: {len(g)} ({EXPECT_GTDB})", flush=True)
    with ThreadPoolExecutor(8) as ex:
        g = list(ex.map(fetch_ncbi, g))
    m = mgnify_rows()
    print(f"MGnify genomes: {len(m)} "
          f"({ {c: sum(r['cls'] == c for r in m) for c in CLASSES} })", flush=True)
    with ThreadPoolExecutor(6) as ex:
        m = list(ex.map(fetch_mgnify, m))

    # MGnify genomes that are the same NCBI assembly as a GTDB r232 rep (version-agnostic)
    reps = {r["accession"].split(".")[0][4:] for r in g}
    for r in m:
        core = r["ncbi_accession"].split(".")[0][4:] if r["ncbi_accession"] else ""
        r["same_assembly_as_gtdb_rep"] = "yes" if core and core in reps else "no"

    cols = ["genome_id", "source", "catalogue", "biome", "cls", "lineage", "ncbi_accession",
            "geographic_origin", "completeness", "contamination", "mgnify_type",
            "same_assembly_as_gtdb_rep", "asm_dir", "status", "bytes", "md5", "md5_ok"]
    with open(os.path.join(RES, "v2_thermo_relatives_comparators.tsv"), "w") as o:
        w = csv.DictWriter(o, cols, delimiter="\t", extrasaction="ignore", restval="")
        w.writeheader()
        for r in g + m:
            w.writerow({**r, "cls": r["cls"]})

    # genome path list: desert reps first
    desert = []
    with open(ECOL) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["phylum"].replace("p__", "") == "Thermoplasmatota":
                f = os.path.join(DREP, r["genome"] + ".fa")
                assert os.path.exists(f), f
                desert.append((r["genome"], "desert", r["class"], f))
    assert len(desert) == 237, len(desert)
    with open(os.path.join(TMP, "genome_paths.tsv"), "w") as o:
        o.write("genome_id\tset\tclass\tpath\n")
        for t in desert:
            o.write("\t".join(t) + "\n")
        for r in g + m:
            if r["status"] == "ok":
                o.write(f"{r['genome_id']}\t{r['source']}\t{r['cls']}\t{r['path']}\n")

    for name, rows in (("GTDB", g), ("MGnify", m)):
        bad = [(r["genome_id"], r["status"]) for r in rows if r["status"] != "ok"]
        mism = [r["genome_id"] for r in rows if r.get("md5_ok") == "NO"]
        print(f"{name}: ok {len(rows) - len(bad)} / {len(rows)}; failed {bad}; md5 mismatch {mism}")
    print("MGnify same assembly as a GTDB rep:",
          sum(r["same_assembly_as_gtdb_rep"] == "yes" for r in m))


if __name__ == "__main__":
    main()
