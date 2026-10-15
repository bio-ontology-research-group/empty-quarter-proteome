#!/usr/bin/env python3
"""Read-only checks of round-4 corrections, independent of fixed review summaries."""
from pathlib import Path
import json, tempfile, importlib.util
import numpy as np
import pandas as pd
H=Path(__file__).resolve().parent;O=H/'round4_corrections';n=0

def check(actual,expected,name):
 global n
 assert actual==expected,(name,actual,expected)
 n+=1
r=pd.read_csv(H/'inputs/06_reference_bias/refbias_clean/results/per_sample_comparison.tsv.gz',sep='\t')
check(int(r.n_pic_ko.eq(10543).sum()),419,'PICRUSt reference saturation')
check(float(r.pic_only.median()),3077.5,'predicted only orthologues')
check(round(float(r.jaccard.median()),3),.634,'Jaccard')
m=pd.read_csv(O/'inputs/v2_thermo_metabolism_summary.tsv.gz',sep='\t')
for feature,counts in [('M00012',(139,63)),('ilv_carbon_skeleton',(3,165)),('catalase',(49,153)),('sugar_ABC_msm',(4,132))]:
 for group,want in zip(['thermo_for_halo','halo_matched'],counts):
  a=m[(m.feature==feature)&(m.group==group)].iloc[0]
  check(int(a.n),181,feature+' denominator '+group);check(int(a.n_present),want,feature+' numerator '+group)
check(int(m[(m.feature=='archaellum_fla')&(m.group=='thermo_237')].iloc[0].n_present),236,'archaellum markers')
a=pd.read_csv(O/'results/indicators/plantroot_vs_other_all.tsv',sep='\t')
check(len(a),101897,'root tested');check(int(((a.q<.05)&(a.spearman_r>0)).sum()),20136,'root enriched');check(int(((a.q<.05)&(a.spearman_r<0)).sum()),41046,'root depleted')
cohort=json.loads((O/'results/compartment_cohort.json').read_text());check(cohort['n_samples'],547,'root samples');check(cohort['excluded'],['62PRr2','64Dr3','64PRr2','64Sr2'],'excluded sample identities')
a=pd.read_csv(O/'inputs/gradient_tests.tsv.gz',sep='\t');b=pd.read_csv(O/'results/gradient_tests.tsv',sep='\t');keys=['block','axis','n','trait','metric'];j=a.merge(b,on=keys,validate='one_to_one',suffixes=('_old','_new'))
check(len(j),224,'gradient tests');check(int((j.q_partial_new<.05).sum()),57,'corrected gradient discoveries');check(int(((j.q_partial_old<.05)!=(j.q_partial_new<.05)).sum()),1,'changed FDR decisions');check(float(abs(j.rho_partial_old-j.rho_partial_new).max()),0.0,'unchanged effect estimates')
r=json.loads((O/'results/rrna_parser_audit.json').read_text());check(r['files'],5438,'rRNA files');check(r['changed_calls'],[],'rRNA actual call changes')
# Exercise the corrected parser on a genuine failure mode of the old longest-only rule.
p=O/'scripts/mimag_summarise_corrected.py';spec=importlib.util.spec_from_file_location('rrna',p);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
with tempfile.TemporaryDirectory() as t:
 p=Path(t)/'example.gff';p.write_text('ctg\tbarrnap\trRNA\t1\t900\t.\t+\t.\tName=16S_rRNA\nctg\tbarrnap\trRNA\t1001\t2200\t.\t+\t.\tName=16S_rRNA;partial=true\n')
 check(mod.has_full_rrna(mod.parse_barrnap_gff(p),'16S_rRNA'),True,'complete rRNA survives longer partial hit')
q=pd.read_csv(O/'permanova_blocks_corrected.tsv',sep='\t').set_index('block')
check(int(q.loc['geochem|space[core]','n']),251,'geochem sample size');check(int(q.loc['climate|space[core]','n']),456,'climate sample size')
s=(H.parents[1]/'manuscript/main.tex').read_text()
for phrase in ['20{,}136','41{,}046','181 completeness-matched','419 samples','271{,}000{,}000-sequence anchor','Among 237 Thermoplasmatota representatives, 236 exceeded']:
 check(phrase in s,True,'manuscript: '+phrase)
print(f'PASS: {n} round-4 checks, including key denominators, corrected FDR decisions and rRNA parser regression.')
