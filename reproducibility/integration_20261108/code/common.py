from pathlib import Path
import ast,re,json
import numpy as np,pandas as pd
R=Path('/scratch/project/k10018/eq-proteome-integration-20261106');E=Path('/scratch/project/k10018/eq-ecology-mag-20261105');A=Path('/scratch/project/k10018/eq-archaeal-function-20261106')
M=E/'inputs/ibex/rebuild_20260729/picrust2_current_20260803/mag_ko'
def read(p):return pd.read_csv(p,sep='\t')
def save(d,n):d.to_csv(R/'results'/n,sep='\t',index=False)
def evaluate(expr,values,n):
 def ev(a):
  if isinstance(a,ast.Constant) and isinstance(a.value,str):return values.get(a.value,np.full(n,np.nan))
  if isinstance(a,ast.BinOp) and isinstance(a.op,(ast.BitAnd,ast.BitOr)):
   x,y=ev(a.left),ev(a.right);missing=np.isnan(x)|np.isnan(y)
   if isinstance(a.op,ast.BitAnd):return np.where((x==0)|(y==0),0,np.where(missing,np.nan,1))
   return np.where((x==1)|(y==1),1,np.where(missing,np.nan,0))
  raise ValueError(ast.dump(a))
 return ev(ast.parse(expr,mode='eval').body)
assert evaluate("'a' & 'b'",{'a':np.array([0.,1.,1.]),'b':np.array([np.nan,np.nan,1.])},3).tolist()[0]==0
assert np.isnan(evaluate("'a' & 'b'",{'a':np.array([1.])},1)[0])
assert evaluate("'a' | 'b'",{'a':np.array([1.])},1)[0]==1
