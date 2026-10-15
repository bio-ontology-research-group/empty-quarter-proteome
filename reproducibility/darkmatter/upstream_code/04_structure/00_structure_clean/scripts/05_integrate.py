import sys, json, csv, math, os
V, PIREPS, MANIFEST = sys.argv[1:4]

def norm_plddt(x):
    return x / 100.0 if x > 1.5 else x

# rep metadata: chosen_rep -> (family, n_members, source)
meta = {}
with open(f"{V}/fam_reps.tsv") as fh:
    r = csv.DictReader(fh, delimiter="\t")
    for row in r:
        meta[row["chosen_rep"]] = (row["family"], int(row["n_members"]), row["source"])

# pLDDT + length per query
info = {}   # id -> (length, plddt)
for fn in sorted(os.listdir(f"{V}/models")):
    if fn.startswith("plddt_shard_") and fn.endswith(".tsv"):
        with open(f"{V}/models/{fn}") as fh:
            next(fh)
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if p[2] == "OOM": continue
                info[p[0]] = (int(p[1]), norm_plddt(float(p[2])))
def from_model(path):
    """Length + mean pLDDT straight from the model: ESMFold writes per-residue pLDDT into the
    B-factor column, so CA count = length and mean CA B-factor = mean pLDDT. Used for the reused
    models, which keeps every rep measured by one instrument instead of trusting a second table
    (224 of our 3,451 reused reps have empty length/mean_plddt in the PI's per-candidate table)."""
    n, tot = 0, 0.0
    with open(path) as fh:
        for line in fh:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                n += 1
                tot += float(line[60:66])
    return (n, norm_plddt(tot / n)) if n else None

for q in meta:
    if q not in info:
        mp = f"{V}/models/reuse/{q}.pdb"
        if os.path.exists(mp):
            m = from_model(mp)
            if m:
                info[q] = m

# last-resort fallback for anything still unmeasured: the PI's table, populated rows only
with open(PIREPS) as fh:
    r = csv.DictReader(fh, delimiter="\t")
    for row in r:
        if (row["rep"] in meta and row["rep"] not in info
                and row["length"] and row["mean_plddt"]):
            info[row["rep"]] = (int(row["length"]), norm_plddt(float(row["mean_plddt"])))
controls = {}
with open(MANIFEST) as fh:
    r = csv.DictReader(fh, delimiter="\t")
    for row in r:
        if row["fold_status"] == "ok":
            controls[row["blind_id"]] = (int(row["fold_length"]), norm_plddt(float(row["fold_mean_plddt"])))

def best_scores(path):
    best = {}   # query -> [qtm, ttm, alntm]
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            q = p[0][:-4] if p[0].endswith(".pdb") else p[0]
            alntm, qtm, ttm = float(p[2]), float(p[3]), float(p[4])
            b = best.setdefault(q, [0.0, 0.0, 0.0])
            b[0] = max(b[0], qtm); b[1] = max(b[1], ttm); b[2] = max(b[2], alntm)
    return best

pdb = best_scores(f"{V}/screen/pdb_hits.raw.tsv")
afdb = best_scores(f"{V}/screen/afdb_hits.raw.tsv")
Z = [0.0, 0.0, 0.0]

def wilson(k, n, z=1.959964):
    if n == 0: return (0.0, 0.0)
    p = k / n; d = 1 + z*z/n
    c = (p + z*z/(2*n)) / d; h = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / d
    return (c - h, c + h)

def classify(ids, source_info, out_tsv):
    rows, stats = [], dict(searched=0, conf_q=0, dual_qtm=0, conf_dual_qtm=0, dual_alntm=0, conf_dual_alntm=0)
    for q in sorted(ids):
        L, pl = source_info[q]
        pq, aq = pdb.get(q, Z), afdb.get(q, Z)
        conf = (L >= 100 and pl >= 0.70)
        dual_q = (pq[0] < 0.5 and aq[0] < 0.5)
        dual_a = (pq[2] < 0.5 and aq[2] < 0.5)
        stats["searched"] += 1
        stats["conf_q"] += conf; stats["dual_qtm"] += dual_q
        stats["conf_dual_qtm"] += (conf and dual_q)
        stats["dual_alntm"] += dual_a; stats["conf_dual_alntm"] += (conf and dual_a)
        rows.append([q, L, f"{pl:.4f}", f"{pq[0]:.4f}", f"{pq[1]:.4f}", f"{pq[2]:.4f}",
                     f"{aq[0]:.4f}", f"{aq[1]:.4f}", f"{aq[2]:.4f}", int(conf), int(dual_q), int(conf and dual_q)])
    with open(out_tsv, "w") as fh:
        fh.write("id\tlength\tmean_plddt\tbest_pdb_qtm\tbest_pdb_ttm\tbest_pdb_alntm\tbest_afdb_qtm\tbest_afdb_ttm\tbest_afdb_alntm\tconfidence_qualified\tdual_novel_qtm\tconf_qualified_dual_novel\n")
        for r_ in rows: fh.write("\t".join(map(str, r_)) + "\n")
    return stats

rep_ids = [q for q in meta if q in info]
missing = [q for q in meta if q not in info]
rs = classify(rep_ids, info, f"{V}/per_rep_structures.tsv")
cs = classify(controls, controls, f"{V}/per_control_structures.tsv")
flo, fhi = wilson(cs["conf_dual_qtm"], cs["conf_q"])
summary = {
    "admitted_families": len(meta),
    "reps_searched": rs["searched"], "reps_missing_model": len(missing),
    "reps_confidence_qualified_Q": rs["conf_q"],
    "reps_dual_novel_qtm_ungated": rs["dual_qtm"],
    "reps_conf_qualified_dual_novel_K": rs["conf_dual_qtm"],
    "K_over_Q_pct": round(100.0 * rs["conf_dual_qtm"] / rs["conf_q"], 2) if rs["conf_q"] else None,
    "sensitivity_alntm_conf_qualified_dual_novel": rs["conf_dual_alntm"],
    "controls_searched": cs["searched"], "controls_confidence_qualified": cs["conf_q"],
    "control_false_dual_novel_conf_qualified": cs["conf_dual_qtm"],
    "control_floor_pct": round(100.0 * cs["conf_dual_qtm"] / cs["conf_q"], 2) if cs["conf_q"] else None,
    "control_floor_wilson95_pct": [round(100*flo, 2), round(100*fhi, 2)],
    "control_floor_alntm": cs["conf_dual_alntm"],
    "gates": {"length_min": 100, "mean_plddt_min": 0.70, "tm_threshold": 0.5, "primary_norm": "qtm"},
}
json.dump(summary, open(f"{V}/summary_05c.json", "w"), indent=1)
with open(f"{V}/SUMMARY_05c.txt", "w") as fh:
    for k, v in summary.items(): fh.write(f"{k}\t{v}\n")
print(json.dumps(summary, indent=1))
if missing: print(f"WARNING: {len(missing)} reps had no model/plddt (first: {missing[:3]})")
