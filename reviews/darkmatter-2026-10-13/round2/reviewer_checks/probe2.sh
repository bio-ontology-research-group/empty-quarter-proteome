cd /ibex/scratch/projects/c2014/rob/metagenomics/04_structure/00_structure_clean/05_validate/
while read i g; do f=$(readlink -f query_pdb/$i.pdb); echo -e "$i\t$g\t$(basename $(dirname $(dirname $f)))\t$(grep -c '^ATOM' $f)\t$(grep -m1 -v -E '^(ATOM|PARENT)' $f | head -c 60)"; done
