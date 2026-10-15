from pathlib import Path
p=Path('manuscript/main.tex');s=p.read_text()
changes={
'excluding the single Thermoplasmata representative whose primary comparator belonged to': 'excluding the single Thermoplasmata representative, a class absent from the external comparator\npool, whose primary comparator belonged to',
'Thus a below-one estimate does not distinguish a coding gene from an artifact-flagged family. A flag on any member need not invalidate every member, and overlapping frames can inherit sequence constraints':'Below-one estimates therefore occur in both unflagged and artifact-flagged families. Member-level flags identify particular sequences within a family, and overlapping frames can inherit sequence constraints',
'These checks identify\nknown spurious-family matches and overlaps among existing calls; they do not detect every\nuncalled alternative ORF or establish translation':'These checks identify\nknown spurious-family matches and overlaps among existing calls; detecting uncalled alternative\nORFs and establishing translation require additional sequence analysis and experimental evidence',
'This union includes protein-supported calls that failed the nucleotide threshold; it is not a DNA-only detection rate.':'This union combines nucleotide-threshold support with protein support, including protein-supported calls that failed the nucleotide threshold.',
'This local diagnostic separates possible gene-call\nand protein-match effects; it does not establish full codon correspondence or validate a gene.':'This local diagnostic separates possible gene-call\nand protein-match effects using overlap and phase at the alignment start. Full codon correspondence\nand coding status require additional validation.',
'The culture result therefore describes a broader bacterial\ncall set and cannot be transferred wholesale to the recurrent-family analysis.':'The culture result therefore describes a broader bacterial\ncall set; 47 admitted families provide its direct connection to the recurrent-family analysis.',
'so genome size does not explain darkness.':'so corrected genome size accounts for a small share of the variation.',
'Novelty does not absorb the lineage term: with first empty or digit-bearing rank in the':'The lineage term persists after adjustment for novelty: with first empty or digit-bearing rank in the',
'and flanking annotation opportunity above versus at or below the pooled median. All 353':'and flanking annotation opportunity above versus at or below the pooled median. Annotation\nopportunity was the fraction of the four flanks per evaluable species carrying a mapped KEGG\northologue, pooled across that family\'s selected loci. All 353',
'The larger sulfur count after aggregation can reflect between-site structure and is not evidence that spatial confounding has been removed.':'The larger sulfur count after aggregation can reflect between-site structure; its interpretation remains conditional on spatial and measurement assumptions.',
'PDB\nrelease 250101, MGnify90 and MGnify90 2026\\_07, GEM, SPIRE, and the five desert studies listed in':'PDB\nrelease 250101 for the initial screen and the frozen PDB Foldseek database downloaded on\n30 July 2026 for paired structural searches, MGnify90 and MGnify90 2026\\_07, GEM,\nSMAG (Zenodo deposit 8223844), SPIRE, and the five desert studies listed in',
}
for a,b in changes.items():
 if s.count(a)!=1:raise ValueError((s.count(a),a))
 s=s.replace(a,b)
p.write_text(s)
print('Applied',len(changes),'source-grounded prose refinements')
