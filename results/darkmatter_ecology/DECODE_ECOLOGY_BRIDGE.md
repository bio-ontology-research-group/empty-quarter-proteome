# Phase E — decoding<->ecology bridge
families with a decoded function: 14227
families matching a desert-adaptive category: 1842

## (1) desert-adaptive functional composition of the decoded dark slice
  aromatic_scav           454 families
  c1_oligotrophy          443 families
  metal_handling          394 families
  osmolyte                263 families
  desiccation_dormancy    221 families
  oxidative_stress         79 families

## (2) is each category over-represented among environmentally-STRUCTURED dark families?
   (structured=9425, cosmopolitan=13342)
  osmolyte             struct=  113 cosmo=  145  OR=1.10  p=4.46e-01  
  oxidative_stress     struct=   30 cosmo=   49  OR=0.87  p=5.69e-01  
  desiccation_dormancy struct=   93 cosmo=  118  OR=1.12  p=4.40e-01  
  c1_oligotrophy       struct=  174 cosmo=  255  OR=0.97  p=7.29e-01  
  metal_handling       struct=  157 cosmo=  228  OR=0.97  p=8.35e-01  
  aromatic_scav        struct=  170 cosmo=  271  OR=0.89  p=2.23e-01  

## (3) does category total abundance track the aridity gradient? (Spearman vs aridity index; low index = more arid)
  osmolyte             n_fam=  263  Spearman r=+0.108 p=2.11e-02  decreases toward arid
  oxidative_stress     n_fam=   79  Spearman r=+0.082 p=8.16e-02  
  desiccation_dormancy n_fam=  221  Spearman r=+0.092 p=4.87e-02  decreases toward arid
  c1_oligotrophy       n_fam=  443  Spearman r=+0.106 p=2.40e-02  decreases toward arid
  metal_handling       n_fam=  394  Spearman r=+0.089 p=5.80e-02  
  aromatic_scav        n_fam=  454  Spearman r=+0.092 p=4.90e-02  decreases toward arid
