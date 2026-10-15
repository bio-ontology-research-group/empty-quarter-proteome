#!/usr/bin/env python3
"""Step 6 -- within-species ecotype test: is the SAME species differently equipped by compartment?

327 species clusters have member genomes recovered from BOTH Surface and Deep samples (77 with
HQ members in both). For each such species: fraction of Surface members carrying the trait vs
fraction of Deep members, Wilcoxon signed-rank across species. A positive result means gene
gain/loss WITHIN species tracks compartment -- sub-species structuring, a stronger claim than
community sorting. Completeness confound: MQ members under-call genes, so (a) HQ-only variant is
reported, (b) surface-vs-deep member completeness balance is checked.

Traits: the three focal repair traits (photolyase K01669/K01859, NHEJ K10979/K01971,
carotenoid K02291/K10027). Member KOs read from per-bin eggNOG annotations (resolver includes
the missing-token fix from repko_clean).

Out: results/ecotype_per_species.tsv , results/ecotype_test.tsv
"""
import csv, os, collections, statistics
from scipy.stats import wilcoxon

M="/ibex/scratch/projects/c2014/rob/metagenomics"
D=f"{M}/05_refine/_drep/out/data_tables"
RES=f"{M}/06_metabolism/uv_survival_clean/results"
TREE="/ibex/scratch/projects/c2014/EmptyQuarter_Data/processed/sequencing-results/sites"
TRAITS={"photolyase":{"K01669","K01859"},"nhej":{"K10979","K01971"},"carotenoid":{"K02291","K10027"}}

meta={r["sample"]:r for r in csv.DictReader(open(f"{M}/metadata/sample_metadata.tsv"),delimiter="\t")}
gi={r["genome"]:(float(r["completeness"]),float(r["contamination"]))
    for r in csv.DictReader(open(f"{D}/genomeInfo.csv"))}

def ann_path(g):
    rep=g[:-3] if g.endswith(".fa") else g
    s=rep.split("__")[0]; m=meta.get(s)
    if "__" not in rep or not m: return None
    rest=rep.split("__",1)[1]
    if rest.startswith("single_"): annos,local=["single_sample"],rest[len("single_"):]
    elif rest.startswith("self_"): annos,local=["self_supervised"],rest[len("self_"):]
    else: annos,local=["self_supervised","single_sample"],rest
    base=f"{TREE}/{m['trip']}/site{m['site']}/{s}/wgs/annotation"
    for a in annos:
        for lc in [local]+([local[:-4]] if local.endswith("_sub") else []):
            p=f"{base}/{a}/eggnog/{lc}/{lc}.emapper.annotations"
            if os.path.exists(p): return p
    return None

def panel_kos(path):
    want=set().union(*TRAITS.values()); kos=set(); hdr=None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#query"): hdr=line.rstrip("\n").lstrip("#").split("\t"); continue
            if line.startswith("#") or not hdr: continue
            row=dict(zip(hdr,line.rstrip("\n").split("\t")))
            ko=row.get("KEGG_ko","-")
            if ko and ko not in ("-",""):
                for k in ko.split(","):
                    k=k[3:] if k.startswith("ko:") else k
                    if k in want: kos.add(k)
    return kos

clusters=collections.defaultdict(list)
for r in csv.DictReader(open(f"{D}/Cdb.csv")):
    g=r["genome"]; s=g.split("__")[0]; c=meta.get(s,{}).get("compartment","?")
    clusters[r["secondary_cluster"]].append((g,c))
span={cl:mem for cl,mem in clusters.items()
      if {"Surface","Deep"}<={c for _,c in mem}}
todo=[(cl,g,c) for cl,mem in span.items() for g,c in mem if c in ("Surface","Deep")]
print(f"[ecotype] species spanning Surface+Deep: {len(span)}; member bins to annotate: {len(todo)}",flush=True)

flags={}   # genome -> (trait flags dict, completeness) ; None if unresolvable
for i,(cl,g,c) in enumerate(todo):
    p=ann_path(g)
    if p:
        kos=panel_kos(p)
        flags[g]={t:int(bool(kos&ks)) for t,ks in TRAITS.items()}
    else:
        flags[g]=None
    if i%200==0: print(f"[ecotype] {i}/{len(todo)}",flush=True)
unres=sum(1 for v in flags.values() if v is None)
print(f"[ecotype] unresolved annotations: {unres}/{len(todo)}",flush=True)

rows=[]
for cl,mem in span.items():
    for hq_only in (False,True):
        grp={"Surface":[],"Deep":[]}
        comp={"Surface":[],"Deep":[]}
        for g,c in mem:
            if c not in grp or flags.get(g) is None: continue
            cp,ct=gi.get(g[:-3] if g.endswith(".fa") else g, gi.get(g,(0,100)))
            if hq_only and not (cp>=90 and ct<=5): continue
            grp[c].append(flags[g]); comp[c].append(cp)
        if not grp["Surface"] or not grp["Deep"]: continue
        row=dict(cluster=cl, hq_only=int(hq_only),
                 n_surf=len(grp["Surface"]), n_deep=len(grp["Deep"]),
                 comp_surf=statistics.median(comp["Surface"]), comp_deep=statistics.median(comp["Deep"]))
        for t in TRAITS:
            row[f"{t}_surf"]=sum(f[t] for f in grp["Surface"])/len(grp["Surface"])
            row[f"{t}_deep"]=sum(f[t] for f in grp["Deep"])/len(grp["Deep"])
        rows.append(row)

cols=list(rows[0].keys()) if rows else []
with open(f"{RES}/ecotype_per_species.tsv","w") as o:
    o.write("\t".join(cols)+"\n")
    for r in rows: o.write("\t".join(f"{r[c]:.3f}" if isinstance(r[c],float) else str(r[c]) for c in cols)+"\n")

with open(f"{RES}/ecotype_test.tsv","w") as o:
    o.write("trait\thq_only\tn_species\tn_diff_nonzero\tmean_surf\tmean_deep\tdelta\tp_wilcoxon\tmedian_comp_surf\tmedian_comp_deep\n")
    for hq in (0,1):
        sub=[r for r in rows if r["hq_only"]==hq]
        if not sub: continue
        cs=statistics.median([r["comp_surf"] for r in sub]); cd=statistics.median([r["comp_deep"] for r in sub])
        for t in TRAITS:
            a=[r[f"{t}_surf"] for r in sub]; b=[r[f"{t}_deep"] for r in sub]
            d=[x-y for x,y in zip(a,b)]; nz=sum(1 for x in d if x!=0)
            try: pv=wilcoxon(a,b).pvalue
            except ValueError: pv=float("nan")
            o.write(f"{t}\t{hq}\t{len(sub)}\t{nz}\t{statistics.mean(a):.4f}\t{statistics.mean(b):.4f}\t{statistics.mean(d):+.4f}\t{pv:.4g}\t{cs:.1f}\t{cd:.1f}\n")
            print(f"[ecotype] hq={hq} {t:11s} n={len(sub):3d} nonzero-diff={nz:3d} surf={statistics.mean(a):.3f} deep={statistics.mean(b):.3f} p={pv:.3g}",flush=True)
print("[ecotype] done",flush=True)
