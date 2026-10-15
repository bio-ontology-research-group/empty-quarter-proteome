import pandas as pd, gzip
R='/home/leechuck/Documents/papers/empty-quarter-mags/reproducibility/darkmatter/inputs/04_structure/00_structure_clean/05_validate/'
g=pd.read_csv('gated_max.tsv',sep='\t',header=None,names=['ref','id','n','mq','mql','ma','mal','mg','mg2'])
w=g.pivot_table(index='id',columns='ref',values=['ma','mg','mg2','mq','mal','mql'])
w.columns=[a+'_'+b for a,b in w.columns]; w=w.fillna(0)
for arm,fn in [('dark','per_rep_structures.tsv.gz'),('ctrl','per_control_structures.tsv.gz')]:
    t=pd.read_csv(R+fn,sep='\t'); t=t[t.confidence_qualified==1].set_index('id').join(w,how='left').fillna(0)
    ma=t[['ma_pdb','ma_afdb']].max(1); mg=t[['mg_pdb','mg_afdb']].max(1); mg2=t[['mg2_pdb','mg2_afdb']].max(1)
    a_tab=t[['best_pdb_alntm','best_afdb_alntm']].max(1)
    print(arm,len(t),'alnTM<0.5 ungated',(ma<.5).sum(),'(table',(a_tab<.5).sum(),') gated aln>=50',(mg<.5).sum(),'gated aln>=100&qcov>=.5',(mg2<.5).sum())
    K=t[t[['best_pdb_qtm','best_afdb_qtm']].max(1)<.5]
    mgK=K[['mg_pdb','mg_afdb']].max(1); mg2K=K[['mg2_pdb','mg2_afdb']].max(1)
    print('  among qTM<0.5:',len(K),'aln>=50 alnTM>=.5:',(mgK>=.5).sum(),' aln>=100&qcov>=.5 alnTM>=.5:',(mg2K>=.5).sum())
    # alignment length of best alnTM hit
    mal=t[['mal_pdb','mal_afdb']].max(1)
    print('  best-alnTM-hit aln length median (qTM<.5 set):', K[['mal_pdb','mal_afdb']].max(1).median(), 'query len median', K.length.median())
