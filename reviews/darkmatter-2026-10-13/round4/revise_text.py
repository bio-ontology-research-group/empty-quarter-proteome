from pathlib import Path
p=Path('manuscript/main.tex');s=p.read_text()
def sub(old,new):
 global s
 assert old in s,old[:100]
 s=s.replace(old,new)
sub('Persistence: with reads from all 551 samples mapped to the member genes, the median family is','With reads from all 551 samples mapped to the member genes, the median family is')
sub('Selection: the median pairwise dN/dS across 97{,}200 evaluable dark families is','The median pairwise dN/dS across 97{,}200 evaluable dark families is')
sub('their axis: their share is 0.55\\% of reads in subsurface samples against 0.17\\% at the surface, higher\nbelow than above at 40 of 45 sites with both positions (Wilcoxon $p = 4 \\times 10^{-7}$, unchanged\nafter dividing by the mapped fraction), and indifferent to plant roots (26 of 53 sites, $p = 0.66$).',
'''their clearest soil-position association: across 45 sites with both depths, the median site-level
share is 0.72\\% of reads below the surface and 0.27\\% at the surface, with the subsurface share
higher at 40 sites (Wilcoxon $p = 4 \\times 10^{-7}$). After division by the mapped fraction,
the subsurface share is higher at 38 of 45 sites ($p = 8 \\times 10^{-7}$).
The root-adjacent contrast has mixed directions across sites (26 of 53, $p = 0.66$).''')
a=s.index('Annotations in the remaining approximately 70\\%');b=s.index('\n\n',a)
s=s[:a]+r'''Annotations in the remaining approximately 70\% provide a provisional metabolic profile.
Across all 237 Thermoplasmatota representatives, the annotation rules identify glycolysis in
82\%, a glyoxylate-cycle module in 79\%, 2-oxoglutarate:ferredoxin oxidoreductase markers in 96\%,
both cytochrome c oxidase subunits in 59\% and archaeal ATPase markers in 91\%.
Branched-chain amino-acid transport markers occur in 79\% and sugar transport markers in 2\%;
serine-protease, endoglucanase and chitinase annotations occur in 96\%, 79\% and 44\%, respectively.
Archaellum-associated markers occur in 236 of 237 representatives (99.6\%) and CRISPR-associated
markers in 79\%. These percentages describe the full lineage set. In the separate comparison
of 181 completeness-matched Thermoplasmatota--Halobacteriota pairs, glyoxylate-module calls occur
in 139 versus 63 genomes (76.8\% versus 34.8\%), branched-chain amino-acid carbon-skeleton
markers in three versus 165 (1.7\% versus 91.2\%), and catalase markers in 49 versus 153
(27.1\% versus 84.5\%). Sugar transport markers occur in four versus 132 (2.2\% versus 72.9\%).
The full 237-genome set includes no detected marker annotations for methanogenesis,
ammonia oxidation or rhodopsin; rhodopsin markers occur in 51 of all 269 Halobacteriota.
These annotations motivate tests of peptide and polysaccharide use in subsurface archaea.
Physiological pathway completeness, auxotrophy and rhodopsin absence require validation with
lineage-sensitive models and experimental measurements. Approximately a third of the
Thermoplasmatota proteins remain without a UniRef50 match.'''+s[b:]
a=s.index('samples with both a genome profile and a 16S profile');b=s.index('\n\n',a)
s=s[:a]+r'''samples with both a genome profile and a 16S profile from the companion survey, PICRUSt2
emits all 10{,}543 orthologues in its reference set for 419 samples (93.5\%) at the presence
threshold of abundance greater than zero. Under this nearly saturated prediction rule,
median recall of genome-encoded orthologues is 86.0\% and median Jaccard similarity is 0.634.
A median of 3{,}077.5 predicted orthologues per sample are absent from that sample's recovered
genomes; incomplete genome recovery limits their interpretation. Spearman $\rho$ on shared
orthologues is 0.839 (median weighted nearest-taxon index 0.159; Fig.~\ref{fig:reference}a).
Across the catalogue, 2{,}517 of the 11{,}374 encoded orthologues (22\%) lie outside the predictor's
emitted set. In 420 of 448 samples, every genome-encoded orthologue inside that set is predicted,
so presence-recall in those samples measures overlap with the reference universe.'''+s[b:]
sub('\\citealp{','\\citep{')
needle='\n\\subsection*{Structural analysis}'
methods=r'''
The metabolic comparison used a separate 1:1 greedy nearest-neighbour match without replacement
on CheckM2 completeness, with a caliper of one percentage point: 181 pairs with Halobacteriota,
216 with other archaea (including Halobacteriota), and 237 with bacteria. We report full-lineage
prevalences separately from paired-subset prevalences. KEGG module calls required at least
two thirds of the top-level steps, each satisfied by any listed alternative orthologue.
Custom marker sets required specified subunit counts: cytochrome c oxidase required both
K02274 and K02275; archaellum markers required two of nine orthologues; CRISPR markers required
one of seven; catalase required one of K03781, K03782 or K07217. The carbon-skeleton and sugar
transport sets required two of five and two of six orthologues, respectively. The complete
orthologue lists and module rules are retained with the metabolic-analysis source code.
'''
sub(needle,methods+needle)
sub('Recall of genome-encoded orthologues and rank correlation on shared orthologues across 448 paired profiles.', 'Presence-recall of genome-encoded orthologues and rank correlation on shared orthologues across 448 paired profiles. PICRUSt2 emits its full orthologue set in 419 samples at the greater-than-zero threshold; recall therefore largely reflects reference overlap.')
sub("at each catalogue's size.", "at the indexed catalogue sizes, except the retained SPIRE expectation uses a 271.0-million-sequence anchor for its 269.5-million-sequence index.")
# Preserve the analytical estimate while making its non-exact SPIRE anchor explicit in Methods.
needle='ordering is descriptive and the content-versus-depth comparison rests on a rarefaction.'
sub(needle,needle+' The SPIRE expectation uses the original 271{,}000{,}000-sequence anchor; the indexed catalogue contains 269{,}496{,}580 sequences.')
s=s.replace('(80.0\\% against 85.7\\%)\n.', '(80.0\\% against 85.7\\%).')
p.write_text(s)
