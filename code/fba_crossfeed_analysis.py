#!/usr/bin/env python3
"""Guild-as-feeder test from FBA cross-feeding potential (replaces SMETANA detailed, which needs a
commercial MILP solver unavailable on this cluster). Joins per-MAG donor/receiver potential to
is_guild + CheckM2 completeness; tests whether the atmospheric guild are net metabolite donors, with
the completeness-sensitivity control. Also folds in the SMETANA global MRO (resource overlap) already
computed. Result feeds outline Part C2.

In:  results/smetana/fba_crossfeed/{mag_donor_receiver.tsv, community_summary.tsv}  (pulled from IBEX)
     results/smetana/smetana_out/global/*_global.tsv  (MRO)
     results/smetana/pilot_mags.tsv ; workbench .../checkm2_all_samples.tsv
Out: results/smetana/{CROSSFEED.md, guild_feeder_test.tsv}
"""
import os, csv, glob
from statistics import median, mean

HERE = os.path.dirname(os.path.abspath(__file__))
SM = os.path.join(HERE, "..", "results", "smetana")
FX = os.path.join(SM, "fba_crossfeed")
CK = os.path.join(HERE, "..", "workbench", "results", "soil_gtdbtk", "checkm2_all_samples.tsv")

def sp(x):
    try: return float(x)
    except (TypeError, ValueError): return float("nan")

def mwu(a, b):
    import math
    a=[x for x in a if x==x]; b=[x for x in b if x==x]
    if len(a)<3 or len(b)<3: return float("nan"), float("nan")
    allv=sorted([(v,0) for v in a]+[(v,1) for v in b]); ranks={}; i=0; vals=[v for v,_ in allv]
    while i<len(vals):
        j=i
        while j+1<len(vals) and vals[j+1]==vals[i]: j+=1
        for k in range(i,j+1): ranks[k]=(i+j)/2.0+1
        i=j+1
    Ra=sum(ranks[k] for k,(_,g) in enumerate(allv) if g==0); na,nb=len(a),len(b)
    Ua=Ra-na*(na+1)/2.0; mu=na*nb/2.0; sd=(na*nb*(na+nb+1)/12.0)**0.5
    if sd==0: return Ua, float("nan")
    z=(Ua-mu)/sd; p_one=0.5*math.erfc(z/2**0.5)  # a>b
    return Ua, p_one

def spearman(xs, ys):
    pts=[(x,y) for x,y in zip(xs,ys) if x==x and y==y]; n=len(pts)
    if n<4: return float("nan"), n
    def rank(v):
        order=sorted(range(len(v)), key=lambda i:v[i]); r=[0.0]*len(v); i=0
        while i<len(v):
            j=i
            while j+1<len(v) and v[order[j+1]]==v[order[i]]: j+=1
            for k in range(i,j+1): r[order[k]]=(i+j)/2.0+1
            i=j+1
        return r
    rx,ry=rank([p[0] for p in pts]),rank([p[1] for p in pts]); mx,my=mean(rx),mean(ry)
    num=sum((a-mx)*(b-my) for a,b in zip(rx,ry))
    den=(sum((a-mx)**2 for a in rx)*sum((b-my)**2 for b in ry))**0.5
    return (num/den if den else float("nan")), n

# meta
guild={}; sampof={}
for r in csv.DictReader(open(os.path.join(SM,"pilot_mags.tsv")), delimiter="\t"):
    guild[r["mag_id"]]=int(r["is_guild"]); sampof[r["mag_id"]]=r["sample"]
compl={}
for r in csv.DictReader(open(CK), delimiter="\t"):
    compl[f"{r['sample']}__{r['Name']}"]=sp(r["Completeness"])
# MRO
mro={}
for f in glob.glob(os.path.join(SM,"smetana_out","global","*_global.tsv")):
    for r in csv.DictReader(open(f), delimiter="\t"):
        mro[os.path.basename(f).replace("_global.tsv","")]=sp(r.get("mro"))

# donor/receiver
rows=[]
for r in csv.DictReader(open(os.path.join(FX,"mag_donor_receiver.tsv")), delimiter="\t"):
    mag=r["mag"]; do=int(r["donor_organic"]); re=int(r["receiver_organic"])
    tot=do+re
    rows.append(dict(mag=mag, sample=r["sample"], is_guild=guild.get(mag,0),
                     completeness=compl.get(mag,float("nan")), donor=do, receiver=re,
                     donor_frac=(do/tot if tot>0 else float("nan"))))

g_d=[r["donor"] for r in rows if r["is_guild"]]; ng_d=[r["donor"] for r in rows if not r["is_guild"]]
g_df=[r["donor_frac"] for r in rows if r["is_guild"]]; ng_df=[r["donor_frac"] for r in rows if not r["is_guild"]]
U1,p1=mwu(g_d,ng_d); U2,p2=mwu(g_df,ng_df)
rho_rc,n_rc=spearman([r["completeness"] for r in rows],[r["receiver"] for r in rows])
rho_dc,_=spearman([r["completeness"] for r in rows],[r["donor"] for r in rows])

with open(os.path.join(SM,"guild_feeder_test.tsv"),"w") as o:
    o.write("mag\tsample\tis_guild\tcompleteness\tdonor_organic\treceiver_organic\tdonor_frac\n")
    for r in sorted(rows,key=lambda r:-r["donor"]):
        o.write(f"{r['mag']}\t{r['sample']}\t{r['is_guild']}\t{r['completeness']:.1f}\t{r['donor']}\t"
                f"{r['receiver']}\t{r['donor_frac'] if r['donor_frac']!=r['donor_frac'] else round(r['donor_frac'],3)}\n")

mro_vals=[v for v in mro.values() if v==v]
ml=["# Community metabolic interdependence in the Empty Quarter (pilot, 16 communities)\n",
    "GEMs: gapseq on HQ single-sample MAGs. **No commercial MILP solver is available on this cluster "
    "(CPLEX/Gurobi access-restricted), so SMETANA *detailed* mode (MILP) is intractable with GLPK at "
    "these community sizes.** We therefore report (i) SMETANA *global* MRO (computed) and (ii) an "
    "LP-based FBA cross-feeding-potential analysis for the guild-as-feeder test.\n",
    "## (1) Resource overlap (SMETANA global MRO)\n",
    f"MRO median {median(mro_vals):.2f} (range {min(mro_vals):.2f}-{max(mro_vals):.2f}) across "
    f"{len(mro_vals)} communities — moderate-to-high metabolic resource overlap (a competition signal "
    "typical of oligotrophic resource sharing). MIP returned n/a (non-interacting baseline infeasible "
    "for communities containing non-growing gap-filled models).\n",
    "## (2) Are the atmospheric guild the metabolic donors? (FBA cross-feeding potential)\n",
    f"Donor potential = # (receiver,metabolite) organic cross-feeds a MAG can supply on a complete "
    f"environment. Guild MAGs (n={len(g_d)}) vs non-guild (n={len(ng_d)}):",
    f"- donor potential: guild median {median(g_d):.1f} vs non-guild {median(ng_d):.1f}; one-sided MWU "
    f"(guild>non-guild) p={p1:.3f}." if g_d and ng_d else "- insufficient data.",
    f"- donor FRACTION: guild {median([x for x in g_df if x==x] or [float('nan')]):.3f} vs non-guild "
    f"{median([x for x in ng_df if x==x] or [float('nan')]):.3f}; one-sided MWU p={p2:.3f}.",
    "\n## CONTROL — is it a MAG-incompleteness artefact?\n",
    f"receiver potential vs completeness Spearman rho={rho_rc:.2f} (n={n_rc}); donor vs completeness "
    f"rho={rho_dc:.2f}. All MAGs are HQ (>=90%), and |rho|<~0.3 would support a biological (not "
    "completeness-driven) signal.",
    "\n## Caveats",
    "- Complete-environment FBA => an UPPER BOUND on secretion/cross-feeding potential; a defined "
    "desert-minimal medium would tighten it.",
    "- Potential (capability) cross-feeding, not realised flux or MILP-scored SMETANA cross-feeding.",
    "- GEMs see only decodable metabolism (dark proteome invisible); pilot = 16/547 samples.",
    "- Commercial-solver access would enable full SMETANA detailed + scale-up (flagged for revision)."]
open(os.path.join(SM,"CROSSFEED.md"),"w").write("\n".join(x for x in ml if x)+"\n")

print(f"MAGs={len(rows)}; guild donor median={median(g_d):.1f} vs non-guild={median(ng_d):.1f} (MWU p={p1:.3f})" if g_d and ng_d else "no data")
print(f"MRO median={median(mro_vals):.2f}; completeness control receiver~compl rho={rho_rc:.2f}")
