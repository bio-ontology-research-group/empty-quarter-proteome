# Input 1: genome_groups.tsv (genome_id, group). Input 2: search hits sorted by query
# (query,target,pident,alnlen,qcov,tcov,evalue,bits). Output: one row per query and target group
# with the number of hits and the best hit (highest bit score).
BEGIN { FS = OFS = "\t" }
NR == FNR { grp[$1] = $2; next }
{
  if ($1 != q) { flush(); q = $1 }
  split($2, t, "::"); g = grp[t[1]]
  if (g == "") { print "unmapped target " $2 > "/dev/stderr"; exit 3 }
  n[g]++; tot++
  if (!(g in bb) || $8 + 0 > bb[g]) { bb[g] = $8 + 0; be[g] = $7; bp[g] = $3; bq[g] = $5; bt[g] = $6; btg[g] = $2 }
}
END { flush() }
function flush(   g) {
  if (q != "") for (g in n) print q, g, n[g], bb[g], be[g], bp[g], bq[g], bt[g], btg[g], tot
  delete n; delete bb; delete be; delete bp; delete bq; delete bt; delete btg; tot = 0
}
