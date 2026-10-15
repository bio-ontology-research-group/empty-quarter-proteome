import pandas as pd, numpy as np
d=pd.read_csv('culture_nucleotide/per_gene.tsv',sep='\t',dtype={'partial':str})
e=d[d.eligible==1]
def rate(x): return f"{x.dna_supported.sum()}/{len(x)}={x.dna_supported.mean():.3f}"
for (dk,sh),g in e.groupby(['dark','protein_shared']):
    print('dark',dk,'shared',sh,rate(g),
      'cov>=.9&id>=.85',f"{((g.best_coverage>=.9)&(g.best_identity>=.85)).mean():.3f}",
      'cov>=.5&id>=.9',f"{((g.best_coverage>=.5)&(g.best_identity>=.9)).mean():.3f}",
      'nohit',f"{(g.best_coverage==0).mean():.3f}",
      'complete',f"{(g.partial=='00').mean():.3f}",'median nt',g.nt.median())
sh=e[(e.protein_shared==1)&(e.dna_supported==0)]
print('shared but unsupported: best cov median',sh.best_coverage.median(),'best id median',sh.best_identity.median())
print(pd.cut(sh.best_identity,[0,.8,.85,.9,.95,1.01]).value_counts().sort_index())
print(pd.cut(sh.best_coverage,[-.01,0,.5,.9,1.01]).value_counts().sort_index())
# per pair
pp=e.groupby('rep').apply(lambda g: pd.Series({'ann_sh':g[(g.dark==0)&(g.protein_shared==1)].dna_supported.mean(),'dark_un':g[(g.dark==1)&(g.protein_shared==0)].dna_supported.mean(),'ann_un':g[(g.dark==0)&(g.protein_shared==0)].dna_supported.mean(),'n_dark_un':((g.dark==1)&(g.protein_shared==0)).sum()}))
print(pp.round(3).sort_values('ann_sh').to_string())
# complete-call dark
c=e[(e.dark==1)&(e.partial=='00')]
for sh_,g in c.groupby('protein_shared'): print('dark complete shared',sh_,rate(g))
c=e[(e.dark==1)&(e.partial!='00')]
for sh_,g in c.groupby('protein_shared'): print('dark partial shared',sh_,rate(g))
good=pp[pp.ann_sh>=0.9].index
g=e[e.rep.isin(good)]
print('pairs with ann_sh>=0.9:',len(good))
for (dk,sh),x in g.groupby(['dark','protein_shared']): print(' dark',dk,'shared',sh,rate(x))
