from pathlib import Path
p=Path('reviews/darkmatter-2026-10-13/round7/root_checks')
s=(p/'h3_download_resumable.py').read_text().replace("TMP=BASE/'auxiliary_tmp'","TMP=BASE/'middle_tmp'")
s=s.replace("rows=[]\nwith", "# Work outward from the middle while existing workers cover both ends.\nparts=sorted(parts,key=lambda x:abs(sorted(y['key'] for y in parts).index(x['key'])-186))\nrows=[]\nwith")
s=s.replace('max_workers=4','max_workers=8').replace('smag_resumable_verified_parts_manifest.tsv','smag_middle_verified_parts_manifest.tsv')
(p/'h3_download_middle.py').write_text(s)
s=(p/'h3_download_resumable.sbatch').read_text().replace('h3_smag_resume','h3_smag_middle').replace('--cpus-per-task=4','--cpus-per-task=8').replace('--mem=4G','--mem=8G').replace('download_resume_%j','download_middle_%j').replace('h3_download_resumable.py','h3_download_middle.py')
(p/'h3_download_middle.sbatch').write_text(s)
