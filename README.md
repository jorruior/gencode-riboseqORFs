## README

This readme is a guideline for any user that wants to use the main Methods to generate an unified set of Ribo-seq ORFs in this article:

While this script is designed to unify independent sets of ORFs and mapped them to a specific Ensembl version, it is not a tool to analyze ribosome-profiling data. However, GENCODE plans a Phase II to re-analyze ribosome-profiling datasets and generate consistent sets of ORFs.


DEPENDENCIES: 

-**gffread** (v0.1.10) http://ccb.jhu.edu/software/stringtie/dl/gffread-0.12.6.tar.gz

-**bedtools** (v2.27.1 or newer) https://github.com/arq5x/bedtools2

-**Python packages: Biopython**



1) **getRegions.py**: Extract collapsed conserved and non-conserved regions for every gene from a set of transcripts:
```
python3 getRegions.py --gtf <TRANSCRIPT_GTF> --fasta <TRANSCRIPT_FASTA> --db <BLAST_DB> -O <ORTHOLOGS_LIST> -m <PSEUDOGENES_GTF> --out <OUT_NAME> 
```
TRANSCRIPT_GTF: Ensembl GTF or similar. Only 'exon' lines are parsed.

TRANSCRIPT_FASTA: Transcript FASTA. IDs in GTF and FASTA should be similar.

BLAST_DB: Blast database (makeblastdb) built over a transcript FASTA on a different species to assess for conservation.

ORTHOLOGS_LIST (optional): List of genes that are known orthologs and should be considered as conserved genes regardless of BLAST.

PSEUDOGENES_GTF (optional): GTF with genes that are known pseudogenes and should be masked. Pseudogenes will not be considered if the biotype is found in the transcript_gtf as well.

OUT_NAME: Unique name for output files. They will be generated in the 'out' folder.


2) **getRNP.py**: Get a list of RNPs for a list of transcripts or regions:
```
python3 getRNP.py --input <TRANSCRIPT_PRED> --sam <SAMFILE_PLUS,SAMFILE_REV/SAMFILE> --cds <BED_CDS> --out <OUT_NAME>
```
TRANSCRIPT_PRED: Transcript coordinates in PRED format for Rfoot analyses. It is recommended to only include non-translated regions (ORFs can be substracted using bedtools).

SAMFILE_PLUS,SAMFILE_REV/SAMFILE: Mapped SAM file with Ribo-Seq reads. If stranded, separate SAM strands into two files and specify both of them separated by comma.

CDS: CDS coordinates in BED format used as training model and to mask translated regions from putative RNPs.

OUT_NAME: Unique name for output files. They will be generated in the 'out' folder.


3) **featureCov.py**: Compute the overlap in the previously computed regions for a specific feature (BED file): RNA-seq and Ribo-seq, promoter, ORF, RNP, or CLIP-seq overlap:
```
python3 featureCov.py --input <REGIONS_OUTPUT_BED> -f <BAM/BED_FEATURES> --stranded <yes/no> --out <OUT_NAME> 
```
REGIONS_OUTPUT_BED: Transcript/region coordinates in a BED file. The output generated by getRegions.py can be used here.

BAM/BED_FEATURES: BAM/BED file including features to check for coverage in the main file. e.g. RNA-seq reads, Ribo-seq reads, promoters, RNPs, or ORFs.

STRANDED: If 'yes', coverage is limited to the same strand. Otherwise, specify 'no'.

OUT_NAME: Unique name for output files. They will be generated in the 'out' folder.


**Guidelines for reproducibility of methods in Ruiz-Orera et al.:**

- The initial mouse transcript dataset corresponds to the annotated version in Ensembl v.89. Repeats were masked using RepeatMasker.

- The human transcript dataset used for building the database in BLAST was obtained from: https://figshare.com/articles/Ruiz-Orera_et_al_2017_/4702375

- The final list in Ruiz-Orera et al. was curated by eliminating lncRNA regions that had protein-coding orthologs in mouse (possible unannotated pseudogenes), putative misannotated UTR regions (located within 4kb from a sense protein-coding gene and/or with evidence of being part of the same gene using RNA-Seq data), or regions with a RNA-seq coverage lower than 56.38 reads/kb.

- The full coordinates of translated ORFs can be found in the input folder (mmu89_t_orfs.tar). getDNDS.py allows to compute dn/ds on a list of ORFs. Both species1 and species2 FASTA are needed (genomic alignments, two fasta files). For this study, the genomic alignments between mouse and human ORFs were used. The alignment of the 9 peptide candidates in lncRNAs can be reproduced by:
```
python3 getDNDS.py -1 input/candidate_peptides_sp1.fa -2 input/candidate_peptides_sp2.fa -o candidates
```

- The folder 'tables' contain all raw data to reproduce the figures in the paper:

*Table 1 contains data for all considered regions (for Figs 1A,2B,3C,4B)*

*Table 2 contains data for equally sized gene subregions (for Fig 1B)*

*Table 3 contains Ribo-Seq coverage data for regions divided by read length (for Fig 4A)*

*Table 4 contains data for all considered genes (for Figs 2A,3B)*

