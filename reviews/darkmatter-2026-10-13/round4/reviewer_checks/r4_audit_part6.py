from pathlib import Path
import pandas as pd

ROOT = Path('/home/leechuck/Documents/papers/empty-quarter-mags')
DM = ROOT / 'reproducibility/darkmatter'
SNAP = DM / 'source_snapshot'

env = pd.read_csv(DM / 'environment_rebuild/env_matrix.tsv', sep='\t')
xrf = pd.read_csv(
    DM / 'environment_inputs/geochemistry/xrf_lab_table_all_trips.tsv', sep='\t'
)
print('=== XRF & ENV_MATRIX FOR SITE 44 / HIGH FE ===')
print(
    env[env['Fe'] > 2][
        ['sample', 'site', 'compartment', 'trip', 'Fe', 'Si', 'Ca', 'Al']
    ]
)
print(
    xrf[xrf['Fe'] > 2][[
        'SampleID',
        'trip',
        'site',
        'compartment',
        'Fe',
        'Fe2O3',
        'Si',
        'SiO2',
        'Ca',
        'Al',
    ]]
)
print(
    xrf[xrf['site'].astype(str) == '44'][[
        'SampleID',
        'trip',
        'site',
        'compartment',
        'Fe',
        'Fe2O3',
        'Si',
        'SiO2',
        'Ca',
        'Al',
    ]]
)

print('\n=== DIFFERENTIAL FUNCTION CLEAN RESULTS ===')
for p in sorted(
    (SNAP / '06_metabolism/differential_function_clean/results').glob('*')
):
  if p.is_file() and p.stat().st_size < 25000:
    print(f'\n--- {p.name} ---')
    print(p.read_text()[:4000])
