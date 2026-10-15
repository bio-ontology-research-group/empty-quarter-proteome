#!/usr/bin/env python3
"""Phase B prep: build the per-sample environment matrix for the dark-matter ecology analysis.
Joins (all local TSVs, no parquet):
  - sample_metadata.tsv  -> sample, site, compartment, trip, trip_num
  - geodata/trip<N>_geodata.tsv (keyed by Site) -> lat, lon, AnnualMeanTemp, AnnualTotalPrecip
  - geochemistry/xrf_lab_table_all_trips.tsv (SampleID == sample) -> recorded elemental values; units unresolved
Output: env_matrix.tsv (one row per sample with complete env) + coverage report.
Run locally: python3 reproducibility/darkmatter/rebuild_environment.py
"""
import os, csv, sys

from pathlib import Path
HERE = Path(__file__).resolve().parent
WB = str(HERE / 'environment_inputs')
META = str(HERE / 'environment_inputs/metadata/sample_metadata.tsv')
GEO = str(HERE / 'environment_inputs/geodata')
XRF = str(HERE / 'environment_inputs/geochemistry/xrf_lab_table_all_trips.tsv')
OUT = str(HERE / 'environment_rebuild')
os.makedirs(OUT, exist_ok=True)

# elemental columns to keep (elements, not the oxide duplicates, to limit collinearity)
ELEMENTS = ["Al","Ba","Br","Ca","Ce","Cl","Co","Cr","Cs","Cu","Eu","Fe","Ga","Gd","K","La",
            "Mg","Mn","Mo","Na","Nd","Ni","P","Pr","S","Sc","Si","Sm","Sr","Ti","V","Zn","Zr"]

meta = {r["sample"]: r for r in csv.DictReader(open(META), delimiter="\t")}

# geodata per trip, keyed by site
geo = {}  # (trip_num, site) -> dict
for f in sorted(os.listdir(GEO)):
    if not f.endswith("_geodata.tsv"):
        continue
    tn = f.replace("trip", "").replace("_geodata.tsv", "")
    for r in csv.DictReader(open(os.path.join(GEO, f)), delimiter="\t"):
        geo[(tn, r["Site"])] = r

xrf = {r["SampleID"]: r for r in csv.DictReader(open(XRF), delimiter="\t")}

def fnum(x):
    try: return float(x)
    except (TypeError, ValueError): return None

rows, miss_geo, miss_xrf = [], [], []
for s, m in meta.items():
    tn, site = m.get("trip_num", ""), m.get("site", "")
    g = geo.get((tn, site)) or geo.get((tn, str(site)))
    x = xrf.get(s)
    if not g:
        miss_geo.append(s); continue          # climate/space required; XRF optional
    if not x:
        miss_xrf.append(s)
    lat, lon = fnum(g.get("Latitude")), fnum(g.get("Longitude"))
    temp, precip = fnum(g.get("AnnualMeanTemp")), fnum(g.get("AnnualTotalPrecip"))
    row = {"sample": s, "site": site, "compartment": m.get("compartment", ""),
           "trip": m.get("trip", ""), "trip_num": tn, "has_xrf": 1 if x else 0,
           "lat": lat, "lon": lon, "temp": temp, "precip": precip}
    # aridity proxy: lower precip = more arid; keep raw precip + a simple index precip/(temp+10)
    row["aridity_index"] = (precip / (temp + 10.0)) if (precip is not None and temp is not None) else None
    for e in ELEMENTS:
        row[e] = fnum(x.get(e)) if x else None
    # salinity proxy = Na + Cl + S (elemental %)
    sal = [row.get(k) for k in ("Na", "Cl", "S")]
    row["salinity_proxy"] = (sum(v for v in sal if v is not None)
                             if (x and any(v is not None for v in sal)) else None)
    rows.append(row)

cols = (["sample","site","compartment","trip","trip_num","has_xrf","lat","lon","temp","precip",
         "aridity_index","salinity_proxy"] + ELEMENTS)
with open(os.path.join(OUT, "env_matrix.tsv"), "w") as o:
    o.write("\t".join(cols) + "\n")
    for r in rows:
        o.write("\t".join("" if r.get(c) is None else str(r.get(c)) for c in cols) + "\n")

print(f"env_matrix.tsv: {len(rows)}/{len(meta)} samples with complete geo+XRF")
print(f"  missing geodata: {len(miss_geo)}  {miss_geo[:8]}")
print(f"  missing XRF:     {len(miss_xrf)}  {miss_xrf[:8]}")
# quick gradient sanity: V and precip ranges
for k in ("V","precip","aridity_index","salinity_proxy"):
    vals = sorted(r[k] for r in rows if r.get(k) is not None)
    if vals:
        print(f"  {k}: n={len(vals)} min={vals[0]:.4g} median={vals[len(vals)//2]:.4g} max={vals[-1]:.4g}")
