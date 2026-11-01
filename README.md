# compare-phage-genomes
Using GenBank files, compare the nucleotide and protein content of two annotated phage genomes.

Pipeline:
1. Parses CDS features from both GenBank files.
2. Pairs orthologous proteins between the genomes
3. Globally aligns each protein pair (BLOSUM62) and walks the alignment codon by codon.
4. Classifies every aligned position:
     - identical codon              
     - synonymous differences
     - nonsynonymous differences
   (Residues that fall in alignment gaps are counted separately as indel residues.)
5. Reports per-gene results, genome-wide totals, and a per-site list of differences.

Requires biopython >= 1.80.
