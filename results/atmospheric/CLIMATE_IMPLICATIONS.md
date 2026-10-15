# Climate-change implications of the Empty Quarter atmospheric trace-gas guild

*Framing for the discussion section. We measure genetic **capacity** (enzyme class, host taxon, metabolic
context), not in-situ **flux** — every claim below is an implication/hypothesis to be tempered with that
caveat and ideally followed by gas-flux incubations. But the genomic census constrains parameters that
global models currently guess.*

## 1. Soil microbes are the dominant sink for atmospheric H₂ — and we just enumerated them in the biggest sand desert
Molecular H₂ (~0.5 ppmv) is removed from the atmosphere mainly by **high-affinity [NiFe]-hydrogenases in
soil bacteria** (groups 1h and 1l) — the largest term in the global H₂ budget (soil uptake ≈ 60–80 Tg yr⁻¹,
~70–80% of total removal). Our data show **3,386 MAGs (22%)** of the EQ community carry exactly these
enzymes, **dominated by Group 1l (5,292)** over the classic mesic-soil Group-1h "Huc" (986). Global H₂-sink
parameterizations are built largely on temperate/boreal/agricultural soils; **drylands (~35–45% of land
area, and expanding) are under-represented.** A dense high-affinity guild in the world's largest sand
desert implies the **dryland contribution to the global H₂ sink is larger and more taxonomically distinct
than current budgets assume.**

## 2. H₂ is an indirect greenhouse gas — so the sink is climate-relevant
H₂ scavenges tropospheric **OH radicals**, extending the lifetime of **methane**, and raises tropospheric
ozone and stratospheric water vapour. Recent assessments give H₂ a **GWP₁₀₀ ≈ 11 (±5)**. The size and
resilience of the **biological soil sink is the key uncertainty** in how rising H₂ translates to warming.

## 3. The hydrogen economy makes this urgent
A global hydrogen economy will **leak H₂** (production, storage, transport; estimates ~1–10%), pushing
atmospheric H₂ upward and partially offsetting the climate benefit of displacing fossil fuels. The **soil
high-affinity hydrogenase reservoir is the natural buffer.** Characterising that reservoir in deserts —
guild size, enzyme class (1h vs 1l kinetics/moisture optima), and host taxa — directly informs projections
of how much extra anthropogenic H₂ the biosphere can absorb.

## 4. CO follows the same logic
Soils are a net sink for atmospheric **CO** via form-I CO dehydrogenase; CO also perturbs OH/methane
chemistry. **3,828 MAGs (25%)** carry bona fide form-I CODH, and **1,204 MAGs co-encode both** CO and H₂
oxidation — a guild built for scavenging multiple atmospheric reductants for survival energy.

## 5. Atmospheric chemosynthesis → a photosynthesis-independent carbon sink in drylands
In Antarctic desert soils, atmospheric H₂/CO oxidation **powers CO₂ fixation** ("atmospheric chemosynthesis";
Ji et al. 2017, *Nature*) — primary production decoupled from light and water. If EQ guild members couple
trace-gas oxidation to carbon fixation (testable from the gapseq models + RuBisCO/CBB markers we already
have — `../metabolism_rbcL_by_phylum.tsv`), the Empty Quarter may host **trace-gas-driven primary production**,
a candidate contributor to the under-explained **dryland carbon sink**.

## 6. Feedbacks are two-way — and that is the point
- **Buffering (negative feedback):** warming + aridification **expand drylands**, potentially expanding the
  range of this trace-gas guild → a larger global H₂/CO sink.
- **Suppression:** high-affinity oxidation needs thin water films; **extreme desiccation suppresses activity**,
  which could shrink the sink. Temperature sensitivity of desert-clade (esp. Group 1l) hydrogenases is
  largely **unknown**.
- Net direction is genuinely uncertain — and a genome-resolved guild census in the most extreme dryland is
  exactly the input needed to resolve it.

## What this dataset uniquely provides to climate science
A **genome-resolved inventory** of the high-affinity H₂/CO oxidisers — enzyme class (1h vs 1l), host taxonomy,
abundance (via 16S), and metabolic context (via gapseq FBA) — for a biome essentially **absent from current
global trace-gas models**. It replaces a parameterized guess with measured guild structure for the largest
sand desert on Earth.

## Honest limits / next experiments
Genomic capacity ≠ measured flux. The decisive follow-ups: (i) **gas-chromatography soil incubations** to
measure actual atmospheric H₂/CO uptake rates and their moisture/temperature response; (ii) metatranscriptomics
to confirm in-situ expression; (iii) kinetic characterisation of representative Group-1l enzymes (K_m, O₂
tolerance). Our genomic census tells modellers *who* and *how many*; the flux work tells them *how fast*.

## Key references
- Greening et al. 2016, *ISME J* — hydrogenase classification (HydDB).
- Søndergaard et al. 2016, *Sci Rep* 6:34212 — HydDB.
- Ji et al. 2017, *Nature* 552:400 — atmospheric chemosynthesis in Antarctic desert soils.
- Cordero et al. 2019, *ISME J* — atmospheric CO oxidation supports microbial survival.
- Bay et al. 2021, *Nature Microbiology* / Greening & Grinter 2022, *Nat Rev Microbiol* — trace-gas survival.
- Paulot et al. / Ocko & Hamburg 2022 — H₂ indirect GWP and hydrogen-economy leakage.
