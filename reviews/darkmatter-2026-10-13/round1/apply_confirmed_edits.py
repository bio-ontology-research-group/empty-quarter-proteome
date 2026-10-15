from pathlib import Path
p = Path('paper/main.tex')
s = p.read_text()
def replace(old, new):
    global s
    assert old in s, old[:100]
    s = s.replace(old, new)
def between(start, end, new):
    global s
    a=s.index(start); b=s.index(end,a)
    s=s[:a]+new+'\n\n'+s[b:]
between('Most microbial proteins match no characterised family,', '% ============================================================================\n\\section*{Introduction}', r'''Reference coverage limits interpretation of microbial gene content. We analysed 551 soil
metagenomes from the Rub' al-Khali across five expeditions and recovered 5{,}438 species-level
genome representatives. Among 43.4 million protein clusters linked to unrefined genome bins,
9.89\% lacked a detected UniRef50 match. A recurrent subset comprised 106{,}042 families containing
536{,}437 cluster representatives, 12.5\% of the dark set. Among 97{,}200 evaluable families,
96.7\% had median pairwise dN/dS below one; read recruitment detected 79.4\% of the admitted
families in all five expeditions. These results support sequence constraint within a recurrent
subset of the reference-absent repertoire. Per-genome darkness varied strongly by lineage:
237 Thermoplasmatota representatives had a median dark fraction of 31.0\%, and this phylum
carried 19.8\% of the single-phylum-attributable dark families. Soil- and genome-derived
catalogues recovered more dark families than expected from random subsets of MGnify90.
The catalogue links reference-absent gene content to genomes, sampling history and environmental
measurements, providing targets for tests of coding status, strain variation and function.''')
replace('\\subsection*{The organisms are ordinary at every rank above order}', '\\subsection*{Genome recovery and taxonomic coverage}')
replace('Every representative was placed to class rank, within 51 named phyla. The first rank without a\nnamed taxon was order', 'Every representative was placed to class rank, within 51 GTDB phylum labels. The first empty\ntaxonomic rank was order')
replace('187 carry a named species', '187 carry a species assignment')
replace('\\subsection*{A tenth of the genome-resolved protein clusters have no reference homolog, and they behave as coding sequence}', '\\subsection*{A recurrent subset of sequence-dark clusters shows frame-dependent sequence constraint}')
between('proteins or of species. Length, not assembly completeness, separates dark from annotated clusters:', 'Persistence: with reads', r'''proteins or of species. Within the genome-resolved dark set, 65.8\% of representatives are
shorter than 100 residues and 56.4\% are partial gene calls. Across the assembly-wide catalogue,
including unbinned sequence, the corresponding dark-set values are 75.7\% and 86.7\%; among all
449 million representatives, including dark clusters, 29.7\% are shorter than 100 residues and
87.5\% are partial. These distributions use different catalogue scopes (Fig.~\ref{fig:dark}b).

We selected recurrent families to examine sequence constraint and distribution. Clustering the
4{,}292{,}446 genome-resolved dark representatives at 50\% identity produced 3{,}658{,}797 clusters:
3{,}353{,}521 singletons, 196{,}450 two-member clusters and 108{,}826 clusters with at least three
members. Requiring at least three samples and two sites retained 106{,}042 families with 536{,}437
members, or 12.5\% of the dark representatives. The other 2{,}784 clusters with at least three
members failed the spread criterion. An additional 782{,}616 singleton or two-member clusters
recurred across the required samples and sites but remained outside the admitted family set.
The median admitted family was assembled in 12 samples at 9 sites. Subsequent family analyses
are conditional on this recurrence-based selection.''')
between('The dark families are measurably less\nconstrained than the matched control.', 'Cultivation is the fourth test.', r'''The pooled distributions differ, and the groups also differ in sequence divergence: median
family mean dS is 0.607 in the dark arm and 0.863 in the control. Comparisons therefore depend
on divergence strata. For families with mean dS at most 0.10, median dN/dS is 0.582 in the dark
arm and 0.321 in the control; for mean dS between 0.5 and 1.0, the medians are 0.435 and 0.452.
The control provides a reference distribution for the estimation pipeline, with divergence and
family-admission differences that remain relevant to interpretation.

For a seeded input subset of 10{,}000 families, shifting member genes by one nucleotide and
masking stop codons produced 9{,}148 evaluable families, with median dN/dS 1.086 and 37.6\%
below one. The corresponding masked annotated-frame analysis yielded 9{,}182 families, with
median 0.427 and 97.0\% below one. Both arms were evaluable for 9{,}029 families, and the shifted
frame gave the higher estimate in 92.4\%. This contrast shows that the estimated constraint
depends on reading frame. Interpreting coding status also requires tests of overlapping gene
calls and of the effects of sequence ascertainment.''')
between('Cultivation is the fourth test.', '\\subsection*{Lineage, not genome quality, sets how much of a genome is unknown}', r'''Genomes recovered from cultures provide a complementary comparison. Long-read assemblies of
722 cultured genomes from saline pans, oases and campground soils of the basin were assessed
with the same UniRef50 search rule. They shared a species-level ANI match with 25 catalogue
representatives, yielding 23 representative--culture groups after assignment to the best-ANI
representative. Excluding the \emph{Escherichia coli} representative detected in the no-template
control left 22 groups. Cultured genomes had a lower dark fraction in 18 of these groups,
with a median difference of $-0.36$ percentage points (Wilcoxon $p = 0.0066$).
The cultured assemblies had a median of two contigs, and 49\% were single-contig assemblies.
Across all cultured genomes, dark proteins had median length 59 residues, compared with 70 for
genome-resolved dark representatives.

In the 22 retained representative--culture pairs, reciprocal protein matches covered 65{,}478
of 81{,}723 UniRef50-matched MAG proteins (80.1\%) and 238 of 1{,}507 dark MAG proteins (15.8\%).
The annotated fraction was more frequently shared in every pair; shared matches had median
identity 99.7\%. Among the 166 complete dark calls of at least 100 residues, 26 had a reciprocal
partner. These comparisons measure sharing between two predicted proteomes. Gene-caller
sensitivity, assembly fragmentation, bin membership and strain variation can each contribute to
unmatched calls. Nucleotide-level comparisons and matched gene prediction would distinguish
these contributions. The cultured comparison covers well-referenced bacteria with low dark
fractions; archaeal lineages require separate validation.''')
replace('\\subsection*{Lineage, not genome quality, sets how much of a genome is unknown}', '\\subsection*{Per-genome darkness varies strongly among lineages}')
replace('\\subsection*{Two poorly sampled archaeal classes carry a fifth of the attributable dark families}', '\\subsection*{Thermoplasmatota carry a large share of the attributable dark families}')
replace('Matching enlarges the contrast, because the unmatched\ncomparators of similar novelty were themselves unusually dark. Rarity, completeness and contiguity\ndo not account for these genomes, and the other archaea of the same soils do not show the pattern.', 'The contrast persists after this quality and abundance matching. In a separate comparison\nrestricted to genomes whose first alphanumeric placeholder occurs at class rank, the 236\nThermoplasmatota representatives have median darkness 0.310, compared with 0.125 in 324 other\nrepresentatives, a 2.5-fold contrast. These descriptive comparisons share phylogenetic dependence\nwithin lineages.')
replace('The two classes are a minor part of the community and a major part of the unknown.', 'Thermoplasmatota contribute a small share of community abundance and a larger share of the\nattributable unknown repertoire.')
replace('\\subsection*{Structural novelty is architectural, not local, against a blinded floor}', '\\subsection*{Whole-protein and local structural matches give different reference-coverage estimates}')
replace('(the 76\\% of dark representatives under 100 residues are\noutside this test)', '(the 65.8\\% of genome-resolved dark representatives under 100 residues are\noutside this test)')
replace('novelty is therefore architectural, below the reach of sequence search.', 'the frequency of reference absence depends on whether the score is normalised by query length\nor aligned length.')
between('Structure comparison recovers relatedness in proteins where the\nstandard sequence method finds nothing,', '\\subsection*{The gap is in the reference set: biome-matched content beats depth}', r'''Structure comparison identifies locally similar regions in sequence-dark proteins. Unmatched
sequence may comprise additional domains, extensions, flexible regions or prediction error;
domain-level analysis would distinguish these possibilities. The local search used a permissive
E-value threshold and no minimum aligned-length gate. The 1{,}031 families are candidates for
whole-protein reference absence under the stated metric. The 3.35\% no-hit rate in known-fold
controls measures retrieval of predicted versions of PDB-derived proteins against references
containing their source structures. Sequence-homology exclusion and matched triage would be
needed to calibrate remote-homology detection in the candidate set.''')
replace('\\subsection*{The gap is in the reference set: biome-matched content beats depth}', '\\subsection*{Reference composition affects recovery of dark families}')
replace('a rarefaction model that reproduces the\nobserved MGnify90 absence to 0.15 percentage points', 'a rarefaction model whose full-catalogue endpoint agrees with the\nobserved MGnify90 absence to 0.15 percentage points')
replace('Per sequence added, references matched to the biome\nbuy more coverage of this proteome than generic depth.', 'Soil- and genome-derived catalogues recover more dark families than the random-MGnify90\nexpectation. Their redundancy, assembly basis and biome composition differ, so the comparison\ncombines these properties. Searches of actual random subsets and matched non-desert soil\ncatalogues would separate their contributions.')
between('The unknown tenth of this soil proteome is the accessory genome', 'The global family catalogues established', r'''This catalogue links reference-absent protein clusters to genomes and repeated sampling of a
hyperarid environment. Evidence is strongest for the recurrent subset: 106{,}042 families hold
12.5\% of genome-resolved dark representatives, and 79.4\% recruit reads in all five expeditions.
Most evaluable families have median pairwise dN/dS below one, and the annotated frame usually
has a lower estimate than the shifted frame. These observations support frame-dependent
sequence constraint within the selected families. Recurrence and recruitment describe the
persistence of homologous DNA; expression measurements and screens for overlapping spurious
calls would further resolve coding status and biological activity.

Darkness also varies strongly among lineages. Quality- and abundance-matched comparisons
retain a large Thermoplasmatota contrast, while a comparison at the same placeholder rank gives
a smaller 2.5-fold difference. Reference proximity and shared ancestry therefore remain relevant
to the lineage pattern. \citet{maeke2025thermoplasmatota} examined the relationship between rarity
and unknown-gene content in other Thermoplasmatota; direct searches against close relatives
would test how much of the desert catalogue's darkness reflects representation in UniRef50.

The cultured-genome comparison identifies lower sharing of dark than annotated protein calls,
even after exclusion of the confirmed \emph{E. coli} contaminant. Its interpretation depends on
gene prediction and genome reconstruction as well as strain variation. Re-calling both partners
with one pipeline and aligning unmatched coding sequences to the partner's genome would provide
a direct test of accessory content. The current control evidence also has limited campaign
coverage: one no-template library was sequenced with expedition four, and the catalogue-level
lineage and environmental summaries still include the flagged representatives. Sensitivity
analyses across these summaries are required to assess their contribution.''')
replace('36.3\\% of confident models without a covering match, 6.6 to 8.1 times a blinded floor, sits inside\nthat consensus once read as whole-protein architecture: novelty below the reach of domain databases,\nrecovered by structure comparison.', '36.3\\% of confident models fall below the whole-query similarity threshold, while local\nsimilarity is detected much more frequently. The known-fold comparison provides a positive-control\nretrieval baseline; domain boundaries and remote-homology controls would refine the biological\ninterpretation of the unmatched regions.')
replace('biome-matched references explain 24 points more of the dark proteome\nthan a random reference slice of the same size, against one point for known proteins. Per sequence\nadded, soil references grow faster in coverage by sampling under-sampled biomes than by deepening\nsampled ones.', 'the desert-derived comparison explains 24 percentage points more of the dark families\nthan expected from a random MGnify90 subset of the same size, against about one point for the\ncontrol. Comparable gains in other genome-derived catalogues motivate a test of reference\ncomposition, assembly basis and redundancy using matched sequence sets.')
replace('and a\nresistome led by $\\beta$-lactam and arsenic determinants in a biome with almost no human input is\npre-anthropogenic', 'and the detected $\\beta$-lactam and arsenic resistance determinants can be considered alongside\nreports of environmental resistance in ancient and remote samples')
replace('Medium-quality genomes\nunder-call gene presence, which biases position contrasts toward the null;', 'Incomplete genomes\nreduce sensitivity to gene presence, with consequences that depend on how genome quality varies\namong the compared communities;')
replace('sequenced on an Illumina NovaSeq 6000 in 150-nt paired-end mode at KAUST.', 'sequenced at KAUST in 150-nt paired-end mode on Illumina NovaSeq 6000 and NovaSeq X instruments.\nThe submission comprises 554 runs for 551 samples: 221 NovaSeq 6000 runs and 333 NovaSeq X runs\nacross four sequencing batches. The sample-to-run manifest records instrument, library kit and\nbatch; three samples have two submitted runs.\\todo{verify the exact duplicate-run merge commands\nand include batch sensitivity in environmental analyses}')
replace('Both push the\ncontrol toward more reference absence and lower dN/dS, so both comparisons in this paper are\nconservative.', 'These differences, together with the divergent dS distributions, affect comparability; their\ndirection and size depend on the endpoint.')
replace('It was run on a seeded random subset of 10{,}000 evaluable dark families (seed 20260920), together\nwith the same subset in the annotated frame as a check that the masking step does not move the\nestimate (0.427 against 0.425 for the full arm).', 'It was run on a seeded subset of 10{,}000 families with at least three indexed coding\nsequences (seed 20260920), together with the same input subset in the annotated frame. Stop\ncodons were masked in both calibration arms. The null, annotated-frame and paired output sets\ncontain 9{,}148, 9{,}182 and 9{,}029 families, respectively.')
replace('Long-read genomes of 722 strains cultured', 'Long-read genomes recovered from 722 cultured assemblies')
replace('compared by Wilcoxon signed-rank; pooled rates over\nstrains are reported only as descriptive.', 'before contaminant exclusion. We excluded representative \\texttt{V42PRr2\\_\\_single\\_SemiBin\\_33},\nidentified as \\emph{E. coli} in the no-template control, leaving 22 groups for the paired rate\ncomparison by Wilcoxon signed-rank. Protein sharing used the best-ANI cultured genome of each\nrepresentative and reciprocal best DIAMOND hits (\\texttt{--very-sensitive -e 1e-5 -k 5}),\nrequiring at least 50\\% identity and 50\\% coverage of both proteins. Sharing percentages pool\nprotein counts across the 22 pairs. The complete-call sensitivity subset contains 166 dark MAG\nproteins of at least 100 residues.')
replace('Taxonomic novelty was coded as the first\nunnamed rank (a placeholder name contains a digit)', 'Placeholder-name depth was coded as the first rank with an empty or alphanumeric label\n(a placeholder name contains a digit), separately from the taxonomy census of empty ranks')
replace('The 14{,}426 refined bins and the 5{,}438 species', 'Successful private submission receipts map all 554 runs to ENA accessions; the public\nfilereport returned zero runs on 13 October 2026. The sample-to-run accession map is retained\nwith the reproducibility records. The 14{,}426 refined bins and the 5{,}438 species')
replace('dark against all representatives: length separates them (76\\% against 30\\% under\n100 residues), completeness does not (86.7\\% against 87.5\\% partial).', 'genome-resolved dark representatives (65.8\\% under 100 residues; 56.4\\% partial),\nwith assembly-wide comparisons shown separately.')
replace('0.425 and 0.388), with the frame-shifted null on a 10{,}000-family subset (median 1.087) and the same\nsubset in the annotated frame (0.427).', '0.425 and 0.388), with 9{,}148 evaluable shifted-frame families (median 1.086) and\n9{,}182 annotated-frame calibration families (0.427) from a 10{,}000-family input subset.')
replace('Dark fraction of 23 soil species', 'Dark fraction of 22 retained species groups')
replace('\\textbf{Lineage sets the unknown fraction of a genome, and two archaeal classes carry a\nfifth of it.}', '\\textbf{Lineage-associated variation in per-genome darkness and the Thermoplasmatota contribution.}')
replace('\\textbf{Structural novelty is architectural, not local, against a blinded floor.}', '\\textbf{Whole-protein and local structural matches under two score normalisations.}')
replace('1{,}031 candidate architectures', '1{,}031 whole-query no-hit candidates')
p.write_text(s)
Path('manuscript/main.tex').write_text(s)
