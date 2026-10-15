import pandas as pd, numpy as np
D='/home/leechuck/Documents/papers/empty-quarter-mags/reproducibility/darkmatter/'
d=pd.read_csv(D+'inputs/06_darkmatter/darkfam_basis_clean/results/selection/dnds.tsv.gz',sep='\t')
anti=set(open(D+'antifam/flagged_families.txt').read().split())
f=pd.read_csv(D+'shadow_overlap/family_flags.tsv',sep='\t');ov=set(f.loc[f.opposite_nondark_overlap_members>0,'family_id'])
dk=d[d['set']=='dark']
for lab,s in [('antifam',anti),('overlap',ov),('union',anti|ov)]:
    g=dk[dk.family_id.isin(s)].median_dNdS.astype(float)
    print(lab,len(s),'evaluable',len(g),'median',round(g.median(),3),'frac<1',round((g<1).mean(),3), 'IQR',g.quantile([.25,.75]).round(3).tolist())
print('all dark',len(dk),dk.median_dNdS.median())
# family sizes of flagged: members flagged fraction
print(f[f.opposite_nondark_overlap_members>0][['members','opposite_nondark_overlap_members']].describe().T)
c=pd.read_csv(D+'culture_nucleotide_cigar/per_gene.tsv',sep='\t',dtype={'partial':str})
l=pd.read_csv(D+'culture_locus/per_gene.tsv',sep='\t')
m=l.merge(c[['protein','nt','partial']],on='protein')
m['complete']=m.partial=='00'
print(m.groupby(['dark','classification']).agg(n=('nt','size'),median_nt=('nt','median'),complete=('complete','mean')).round(3))
e=c[c.eligible==1]
e['complete']=e.partial=='00'
print(e.groupby(['dark','protein_shared','dna_supported']).agg(n=('nt','size'),median_nt=('nt','median'),complete=('complete','mean')).round(3))
# DNA-level sharing combined
for dd in [0,1]:
    g=e[e.dark==dd]; s=((g.protein_shared==1)|(g.dna_supported==1)).sum()
    print('dark',dd,'eligible',len(g),'protein_shared',g.protein_shared.sum(),round(g.protein_shared.mean(),3),'protein or DNA',s,round(s/len(g),3))
# no-CDS share among DNA-present calls
for dd in [0,1]:
    g=e[e.dark==dd]; dnapresent=((g.protein_shared==1)|(g.dna_supported==1)).sum()
    nocds=((m.dark==dd)&(m.classification=='no_major_CDS_overlap')).sum()
    print('dark',dd,'no culture CDS',nocds,'of DNA-or-protein present',dnapresent,round(nocds/dnapresent,4))
