from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path('/home/leechuck/Documents/papers/empty-quarter-mags')
DM = ROOT / 'reproducibility/darkmatter'
INP = DM / 'inputs'


def find_inp(name):
  matches = list(INP.rglob(name))
  if not matches:
    raise FileNotFoundError(name)
  return matches[0]


print('=== 1. DNDS & SELECTION TABLES ===')
dnds = pd.read_csv(find_inp('dnds.tsv.gz'), sep='\t')
dnds_sub = pd.read_csv(find_inp('dnds_dark_subset.tsv.gz'), sep='\t')
dnds_null = pd.read_csv(find_inp('dnds_null_frame1.tsv.gz'), sep='\t')
print('dnds rows:', len(dnds), 'cols:', list(dnds.columns))
for grp, df_g in dnds.groupby('set'):
  v = df_g['median_dNdS'].dropna()
  ds = df_g['mean_dS'].dropna()
  print(
      f'{grp}: n={len(v)}, median_dNdS median={v.median():.4f}'
      f' (IQR={np.percentile(v,25):.4f}-{np.percentile(v,75):.4f}),'
      f' <1={(v<1).mean()*100:.2f}%, <0.5={(v<0.5).mean()*100:.2f}%, mean_dS'
      f' median={ds.median():.4f}'
      f' (IQR={np.percentile(ds,25):.4f}-{np.percentile(ds,75):.4f})'
  )

vs = dnds_sub['median_dNdS'].dropna()
print(
    f'dnds_dark_subset: n={len(vs)}, median={vs.median():.4f},'
    f' IQR={np.percentile(vs,25):.4f}-{np.percentile(vs,75):.4f},'
    f' <1={(vs<1).mean()*100:.2f}%'
)
vn = dnds_null['median_dNdS'].dropna()
print(
    f'null frame1: n={len(vn)}, median={vn.median():.4f},'
    f' IQR={np.percentile(vn,25):.4f}-{np.percentile(vn,75):.4f},'
    f' <1={(vn<1).mean()*100:.2f}%'
)

print('\n=== 2. STRUCTURAL TABLES ===')
rep_str = pd.read_csv(find_inp('per_rep_structures.tsv.gz'), sep='\t')
ctl_str = pd.read_csv(find_inp('per_control_structures.tsv.gz'), sep='\t')
print('per_rep_structures:', len(rep_str), 'cols:', list(rep_str.columns))
print('per_control_structures:', len(ctl_str), 'cols:', list(ctl_str.columns))
print(
    'rep pLDDT median:',
    rep_str['mean_plddt'].median(),
    'IQR:',
    np.percentile(rep_str['mean_plddt'], 25),
    np.percentile(rep_str['mean_plddt'], 75),
    '>=70:',
    (rep_str['mean_plddt'] >= 70).sum(),
    (rep_str['mean_plddt'] >= 70).mean() * 100,
)
print(
    'ctl pLDDT median:',
    ctl_str['mean_plddt'].median(),
    'IQR:',
    np.percentile(ctl_str['mean_plddt'], 25),
    np.percentile(ctl_str['mean_plddt'], 75),
    '>=70:',
    (ctl_str['mean_plddt'] >= 70).sum(),
    (ctl_str['mean_plddt'] >= 70).mean() * 100,
)
print(
    'rep confidence_qualified sum:',
    rep_str['confidence_qualified'].sum(),
    'conf_qualified_dual_novel sum:',
    rep_str['conf_qualified_dual_novel'].sum(),
)
print(
    'ctl confidence_qualified sum:',
    ctl_str['confidence_qualified'].sum(),
    'conf_qualified_dual_novel sum:',
    ctl_str['conf_qualified_dual_novel'].sum(),
)

print('\n=== 3. MIMAG TIERS & SAMPLE METADATA ===')
mimag = pd.read_csv(find_inp('mimag_tiers.tsv.gz'), sep='\t')
print('mimag rows:', len(mimag), 'cols:', list(mimag.columns))
for col in mimag.columns:
  if mimag[col].dtype == object or mimag[col].nunique() < 20:
    print(col, mimag[col].value_counts().to_dict())

smeta = pd.read_csv(find_inp('sample_metadata.tsv.gz'), sep='\t')
print('sample_metadata rows:', len(smeta), 'cols:', list(smeta.columns))
for col in ['compartment', 'position', 'zone']:
  if col in smeta.columns:
    print(col, smeta[col].value_counts(dropna=False).to_dict())

print('\n=== 4. ALL DERIVED TSVS UNDER reproducibility/darkmatter ===')
for p in sorted(DM.rglob('*.tsv')):
  rel = p.relative_to(DM)
  if 'source_snapshot' in str(rel) or 'upstream_code' in str(rel):
    continue
  if p.stat().st_size < 15000:
    print(f'\n--- {rel} ---')
    print(p.read_text())
  else:
    print(f'\n--- {rel} (size={p.stat().st_size} bytes, head) ---')
    with open(p) as fh:
      for _ in range(10):
        line = fh.readline()
        if not line:
          break
        print(line.rstrip('\n'))
