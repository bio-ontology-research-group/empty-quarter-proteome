from common import *
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=N/'results';O=N/'figures';O.mkdir(exist_ok=True)
p=pd.read_csv(R/'ecology/paired_tests.tsv',sep='\t');a=p[(p.variant=='primary')&(p.contrast=='trip4-trip3')&(p.compartment=='PlantRoot')].set_index('family');names=['HSP20','HSP70','DnaJ','GrpE','Cpn60','GroES','HSP90','ClpP'];a=a.loc[names];b=pd.read_csv(R/'ecology/HS4_environment_models.tsv',sep='\t');b=b[(b.visit_fixed_effect)&(b.term=='temp10')].set_index('family').loc[names]
fig,ax=plt.subplots(1,2,figsize=(10,5),sharey=True,layout='constrained');y=np.arange(8)
for axis,d,col,title,label in [(ax[0],a,'mean','August − February, root-adjacent soil\n32 matched sites','Difference in abundance-weighted copies'),(ax[1],b,'coefficient','Collection-temperature association\nSite, visit, compartment, humidity and hour adjusted','Difference per 10 °C')]:
 axis.axvline(0,color='#888888',lw=1);axis.errorbar(d[col],y,xerr=np.vstack([d[col]-d.lo,d.hi-d[col]]),fmt='o',color='#235b79',ecolor='#6096ac',capsize=3);axis.set_title(title,fontsize=11);axis.set_xlabel(label,fontsize=10);axis.set_yticks(y,names);axis.spines[['top','right']].set_visible(False);axis.grid(axis='x',alpha=.15)
ax[0].invert_yaxis();fig.suptitle('Encoded protein-maintenance capacity',fontsize=14);fig.savefig(O/'HS1_HS4_effects.pdf');fig.savefig(O/'HS1_HS4_effects.png',dpi=180);a.to_csv(O/'HS1_plot_data.tsv',sep='\t');b.to_csv(O/'HS4_plot_data.tsv',sep='\t');print('DONE figures')
