# Multi-catalogue EQ-specificity of the dark proteome

Fraction of the 23,382 UniRef50-dark family representatives with **no detected homolog** in each named catalogue under the same criterion (DIAMOND blastp `--sensitive`, e-value <1e-3). Higher means more database-relative EQ specificity under this search rule. Catalogue depth is shown because shallow catalogues inflate no-hit fractions.

| Catalogue | Type | # target proteins | EQ dark-family no-hits |
|-----------|------|-------------------:|------------------------:|
| MGnify90 (all biomes) | global | 718M | **36.6%** (8,552) |
| GEM (52k MAGs) | global | 123M | **47.4%** (11,075) |
| ALL deserts (union) | desert | 147M | **31.9%** (7,469) |
| Mojave | desert | 54M | **51.9%** (12,146) |
| Antarctic cold-desert | desert | 49M | **55.9%** (13,064) |
| arid soils | desert | 31M | **41.7%** (9,743) |
| trace-gas-autotrophy desert | desert | 7M | **62.8%** (14,693) |
| Atacama | desert | 6M | **89.1%** (20,838) |

The no-hit fraction remains substantial in the deepest named comparisons — MGnify90 (718M protein clusters, all biomes), GEM (123M proteins from 52k MAGs), and the union of five desert catalogues (147M proteins). This is a database-relative result, not evidence that homologs are geographically absent, and it cannot distinguish true restriction from incomplete reference sampling. The historical no-hit/hit dN/dS contrast is descriptive and cannot validate biogeographic restriction.

