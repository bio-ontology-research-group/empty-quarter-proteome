cd /ibex/scratch/projects/c2014/rob/metagenomics/04_structure/00_structure_clean/05_validate/
while read i g; do
 f=$(readlink -f query_pdb/$i.pdb)
 awk -v i="$i" -v g="$g" '
 function ab(v){return v<0?-v:v}
 /^ATOM/{x=substr($0,31,8)+0;y=substr($0,39,8)+0;z=substr($0,47,8)+0;m=ab(x);if(ab(y)>m)m=ab(y);if(ab(z)>m)m=ab(z);if(m>M)M=m;
  if(substr($0,13,4)==" CA "){n++; if(n>1){d=sqrt((x-px)^2+(y-py)^2+(z-pz)^2); if(d>4.2||d<3.4)bad++}; px=x;py=y;pz=z}}
 END{print i"\t"g"\t"M"\t"n"\t"bad+0}' "$f"
done
