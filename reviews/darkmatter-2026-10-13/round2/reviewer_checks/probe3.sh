cd /ibex/scratch/projects/c2014/rob/metagenomics/04_structure/00_structure_clean/05_validate/
while read i g; do f=$(readlink -f "query_pdb/$i.pdb")
 awk -v i="$i" -v g="$g" '/^ATOM/ && substr($0,13,4)==" CA "{n++;x[n]=substr($0,31,8)+0;y[n]=substr($0,39,8)+0;z[n]=substr($0,47,8)+0}
 END{for(k=1;k<=n;k++){cx+=x[k];cy+=y[k];cz+=z[k]};cx/=n;cy/=n;cz/=n;for(k=1;k<=n;k++)r+=(x[k]-cx)^2+(y[k]-cy)^2+(z[k]-cz)^2;rg=sqrt(r/n);
 h=0;for(k=1;k+3<=n;k++){d=sqrt((x[k]-x[k+3])^2+(y[k]-y[k+3])^2+(z[k]-z[k+3])^2);if(d<6)h++};
 printf "%s\t%s\t%d\t%.1f\t%.3f\t%.2f\n",i,g,n,rg,rg/(2.2*n^0.38),h/(n-3)}' "$f"; done
