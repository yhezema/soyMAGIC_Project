# soyMAGIC_Project
# Reference-Guided De Novo Assembly in Soybean Using Short-Read Whole Genome Sequencing

SoyMAGIC is a novel soybean population derived from eight diverse parental lines, represents a great source for soybean breeding. Accurate and complete genome assemblies of these parents are essential to enhance SNP detection, imputation, and genomic analyses of SoyMAGIC population. Traditional reference-guided assembly is biased towards the used reference, leading to biased or incomplete variant discovery, while de novo assembly remains computationally demanding for complex plant genomes. We developed a reference-guided de novo assembly pipeline to generate a genome assembly for the soybean parents. The approach used was modified from  Lischer and Shimizu (2017) and Schneeberger et al. (2011). The pipeline included reference-guided mapping, combined with de novo assembly of unmapped reads to capture conserved and divergent genomic regions via nine steps, including:
 1-read trimming
 2-reference mapping
 3-block definition
 4-de novo assembly of superblocks and unmapped reads
 5-merge contigs
 6-supercontigs mapping
 7-de novo assembly of unmapped reads
 8-error correction
 9-gap closing

<img width="974" height="880" alt="image" src="https://github.com/user-attachments/assets/2865eafa-ba67-48f9-8495-b8181aa5a796" />
