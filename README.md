## README

This readme is a guideline for any user that wants to generate an unified set of Ribo-seq ORFs and reproduced the Phase I ORF list publised in this article: https://doi.org/10.1038/s41587-022-01369-0

While this script is designed to unify independent sets of ORFs and map them to a specific Ensembl version, it is not a tool to analyze ribosome-profiling data. However, GENCODE plans to develop a Phase II to re-analyze ribosome-profiling datasets and call ORFs.

**Version 1.1 (2022.09.19):**
-The updated version also includes and considers annotated CDS sequences without annotated start and/or stop codons. Phase I contains 48 Ribo-seq ORFs overlapping to incomplete proteins.

-Unitary pseudogenes are now considered as a separate category and included in the output files. GENCODE plans to include this category in next updates, but they were not considered initially in Phase I.


**DEPENDENCIES:**

This script is based on python3 and shell, requiring some additional packages to correctly run all steps:

-**gffread** (v0.1.10) http://ccb.jhu.edu/software/stringtie/dl/gffread-0.12.6.tar.gz

-**bedtools** (v2.27.1 or newer) https://github.com/arq5x/bedtools2

-**Python3 packages: Biopython**


**INPUT ANNOTATION FILES**: 

This pipeline requires a series of files in the correct format to analyze the data. First of all, the user will need to collect all annotation files in a single folder <FOLDER> that will be required as input for the main script. These files can be downloaded from Ensembl or GENCODE; we also included a bash script to automatically download and convert all files for a specific Ensembl release:
```
bash scripts/retrieve_ensembl_data.sh <ENSEMBL_RELEASE (e.g. 101)> <GENOME_ASSEMBLY (e.g. GRCh38)>
```
Files will be stored in a new <FOLDER> named 'Ens + number of release'. The files are:
       
PROTEOME_FASTA: Fastq file with all annotated proteins. The protein ID should be the only element in the header.

TRANSCRIPTOME_FASTA: Fasta file with all transcripts. It should include both mRNA and ncRNA, and the transcript ID should be the only element in the header.

SORTED_TRANSCRIPTOME_GTF: Sorted GTF file with exon and CDS sequences for the transcripts included in TRANSCRIPTOME_FASTA -Ensembl format-. Transcript IDs (field 'transcript_id') and protein IDs (field 'protein_id') should match with the IDs in the fasta files. A GTF can be sorted running this command in bash:
```
sort -k1,1 -k4,4n -k5,5n <TRANSCRIPTOME_GTF> > <SORTED_TRANSCRIPTOME_GTF>
```
TRANSCRIPT_SUPPORT: Tab-delimited file containing all transcripts IDs, Transcript Support Level (TSL) scores, and APPRIS support scores. This file will be used to prioritize transcripts that translate each ORF. **<This file has to be downloaded from Ensembl Biomart, a downloaded compatible file based on Ensembl v.101 is located in the folder Ens101>**

PSITES_BED: File including coordinates of P-sites for the annotated proteins. This file can be obtained from the fasta file by running:
```
python3 ../scripts/calculate_frame_bed.py <SORTED_TRANSCRIPTOME_GTF>
```

**INPUT ORF FILES:**

In addition, the pipeline will require input files with the protein sequences and coordinates of the ORFs that will be inspected and unified, tagged by study. The name of the ORF needs to be concatenated with the name of the study by '--' in the FASTA file (to avoid name duplication from different studies). e.g.:
```
>A1BG_58858387_46aa--6_Chen2020
MQPRAQGAVGVLRSAGDSGLAPSPPVAAQGRGLWGAGEASLIPPRN*
>A1BG_58858945_25aa--6_Chen2020
MPSCAARDPSPTSPSSCCARARRRP*
>A1BG:ENST00000598345.1.3797--3_Raj2016
MPSCAARDPSPTSPSSCCARARRRP*
```
Additionally, both the name of the ORF and the study have to be present in the fourth and fifth fields of the BED files. e.g.:
```
19      58346882        58347022        A1BG_58858387_46aa      6_Chen2020      -
19      58347503        58347580        A1BG_58858945_25aa      6_Chen2020      -
19      58347503        58347580        A1BG:ENST00000598345.1.3797     3_Raj2016       -
```

Ideally both sequences and exonic coordinates should be included as input for the main script, but ORF studies are very heterogeneous and sometimes only sequences or coordinates are available. If only the protein sequences are available, please use the **-a** option in the main script. If only the exonic coordinates are available, it is possible to run this in-house bash script to convert **1-based (Important: All BED/GTF coordinate files should be 1-based)** BED coordinates (<BED_FILE>) into protein sequences with the format needed to run the pipeline:
```
bash scripts/bed1_to_fasta.sh <BED_FILE>
```
Outputs: <BED_FILE>.nucl.fa and <BED_FILE>.prot.fa   


**MAIN SCRIPT - ORF_mapper_to_GENCODE.py**:
```
python3 ORF_mapper_to_GENCODE.py -d <FOLDER> -f <ORFS_FA_FILE> -b <ORFS_BED_FILE (1-based)>
Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -d FOLDER, --input_dir=FOLDER
                        (Required) Directory with required Ensembl files
                        (transcriptome gft and fasta, proteome fasta, tab file
                        with APPRIS and TSL support, protein psites generated
                        by 'calculate_frame_bed.py'
  -f ORFS_FA_FILE, --input_fasta=ORFS_FA_FILE
                        (Required) File with all translated candidate ORFs. A
                        FASTA can be generated from a BED file using the
                        script 'bed1_to_fasta.sh'
  -b ORFS_BED_FILE, --bed=ORFS_BED_FILE
                        (Required) File with 1-based BED coordinates of all
                        translated candidate ORFs. ORF names should match in
                        both fasta and bed files. If -a is activated, the BED
                        coordinates will be written into this file.
  -o OUT_NAME, --output=OUT_NAME
                        Tag name for the generated output files. Default: BED
                        filename
  -l LEN_CUTOFF, --min_len_cutoff=LEN_CUTOFF
                        Minimum ORF length threshold in amino acids (without
                        stop codon). default=16
  -L MAX_LEN_CUTOFF, --max_len_cutoff=MAX_LEN_CUTOFF
                        Maximum ORF length threshold in amino acids (without
                        stop codon). default=none
  -c COL_THR, --collapse_cutoff=COL_THR
                        Minimum required fraction to collapse ORFs with
                        similar stretches of overlapping amino-acid sequences.
                        default=0.9
  -m METHOD, --collapse_method=METHOD
                        Method to cluster ORF variants. 'longest_string' will
                        collapse ORFs if the longest shared string is above -c
                        threshold and the start codon and/or stop codon are 
                        shared (default). 'psite_overlap' is a slower method 
                        that collapses ORFs if the fraction of shared psites
                        is above -c threshold
  -a CALCULATE_COORDINATES, --make_annot_bed=CALCULATE_COORDINATES
                        If ORF BED file is not available, generate it from the
                        FASTA file. WARNING: BED file will be writen to the
                        filename given by -b/--bed. (ATG/NTG/XTG/no, default =
                        no)
  --multiple=MULT       If -a option is activated, include ORFs that map to
                        multiple genomic coordinates. (yes/no, default = yes)
  -g GENOMIC, --genomic=GENOMIC
                        If -a option is activated, this optional argument uses
                        a BED file with ORF genomic coordinates to help
                        mapping ORF sequences to the correct exon-intron
                        positions. Use 'none' if no file is given or -a is not
                        activated. (default = none)
  -G FGENOMIC, --force_genomic=FGENOMIC
                        If -a option is activated, this optional argument uses
                        a BED file with ORF genomic coordinates to FORCE
                        mapping ORF sequences to the correct exon-intron
                        positions. ORFs that cannot be mapped to these genomic
                        regions will be discarded. Use 'none' if no file is
                        given or -a is not activated. (default = none)
  -C CDS_CASES, --add_cds=CDS_CASES
                        If 'yes', include CDSs and non-unitary pseudogenes in
                        the output GTF. (default = 'no')

                        
```
If the BED file including ORF coordinates is not available, the option **-a** should be enabled and the BED coordinates will be written into the file specified by the option **-b**. If **-a** is set to *no* (default), **-b** BED file will be used as input for extracting ORF coordinates. **-a** has three different writing options: **ATG** will write all ORF starting with ATG codon; **NTG** will write all ORFs starting with non-ATG codons, **XTG** will write all ORFs regardless of the translation initiation codon. The first Phase I ORF dataset was generated including ATG ORFs. Additionally, there is the possibility of including an additional 'guide' BED file (**-g**) to limit ORF mapping to the ranges defined by this file. Unlike the main BED file, this annotation file should not include intron coordinates and should define the genomic loci coordinates where the ORF can be found, discarding ORFs that mapped to alternative genomic loci. This parameter is optional and included for the studies that only generated output sequences and genomic coordinates without introns.

Moreover, the script includes the parameter **-m** which offers two options for collapsing ORFs located in the same locus and sharing some degree of similarity. **'longest_string'** is the default parameter used to calculate Phase I ORFs and will collapse pairs of overlapping ORFs when the length of the longest shared amino acid string divided by the length of the short ORF is higher or equal than **-c** (default: 0.9). **'psite_overlap'** will collapse pairs of overlapping ORFs sharing a minimum fraction of P-site positions higher or equal than **-c** (default: 0.9). In both cases the longest ORF will be selected as representative and the shorter ORFs will be considered as variants, with the possibility of a variant being assigned to two or more independent ORFs.

Outputs: 

<OUT_NAME>.orfs.fa: Protein sequences for all unified ORFs  (including pseudogenes and annotated CDS)

<OUT_NAME>.orfs.bed: 1-based BED file with all unified ORFs (including pseudogenes and annotated CDS)

<OUT_NAME>.orfs.gtf: 1-based GTF file with all unified ORFs (excluding pseudogenes and annotated CDS)  

<OUT_NAME>.orfs.out: Table with features of all unified ORFs:

**orf_id:**	*Unique identifier for each Ribo-seq ORF*

**phase_id:**     *Identifer for the Ribo-seq ORF in Phase I ('unknown' otherwise)*

**chrm:**	*Chromosome*

**starts:**	*Start exon coordinates*

**ends:**	*End exon coordinates*

**strand:**	*Strand*

**trans:**	*Ensembl transcript ID (v.101) of the main host transcript assigned to the ORF. When more than one transcript could be assigned to an ORF, the one with the highest APPRIS or Transcript Support Level (TSL) support was prioritized*

**gene:**	*Ensembl gene ID (v.101) of the main host gene*

**gene_name:**	*Ensembl gene name (v.101) of the main host gene*

**orf_biotype:**	*ORF biotype based on 6 outlined categories (Box 1): Upstream ORFs (uORFs), Upstream overlapping ORFs (uoORFs), Downstream ORFs (dORFs), Downstream overlapping ORFs (doORFs), Internal out-of-frame ORFs (intORFs), Long non-coding RNA ORFs (lncRNA-ORFs)*

**gene_biotype:**	*Ensembl gene biotype (v.101) of the main host gene: protein-coding or non-coding. LncRNA-ORFs in protein-coding genes are mapped to processed transcripts*

**pep:**	*Amino acid protein sequence*

**orf_length:**	*ORF length, it includes the stop codon*

**n_datasets:**	*Number of Ribo-seq datasets with evidence of translation of the ORF (main ORF or any alternative clustered ORF)*

**X_study:**	*1: The ORF was identified in this study. 0: The ORF was not identified in this study. A column per study is added to the table.*

**pseudogene_ov:** *1: The ORF overlaps a pseudogene. 0: The ORF does not overlap any pseudogene*

**CDS_ov:**	*1: The ORF overlaps an annotated CDS in any frame. 0: The ORF does not overlap any CDS*

**CDS_as_ov:**	*1: The ORF overlaps an annotated CDS in an alternative antisense frame. 0: The ORF does not overlap any antisense CDS*

**all_trans:**	*Ensembl transcript IDs (v.101) of all the transcripts that could be assigned to the ORF*

**all_genes:**	*Ensembl gene IDs (v.101) of all the genes that could be assigned to the ORF*

**all_gene_names:**	*Ensembl gene names (v.101) of all the genes that could be assigned to the ORF*

**n_variants:**	*Total number of additional clustered ORF isoforms*

**seq_variants**:	*Amino acid protein sequences of the alternative ORF isoforms*

**all_orf_names**:	*Original names of the identified ORFs in each Ribo-seq dataset*
       
**phase_id:**     *Identifer for the Ribo-seq ORF in Phase I*
       
**phase_biotype:**     *Identifer for the Ribo-seq ORF biotype in Phase I*

If **-a** is enabled, the script generates two additional files '<OUT_NAME>.**altmapped**' and '<OUT_NAME>.**unmapped**' with additional statistics about ORFs that did not map to the transcriptome or that map to multiple genomic regions (in this case, a random region is selected if **--multiple** is set to yes). 

Example of an **altmapped** output file:
```
BET1L_180332_35aa--6_Chen2020   altmap  MAPPEGARSRPNQLPAEVPEAARSFNFLRAARGSC*    138230979       138231086       180226  180333  ENST00000446912
BRCA1:NM_007297:chr17:-:194:274--1_Ji2015       altexons        MDLSALRVEEVQNVINAMQKILECPI*     43121676;43124017       43121679;43124096       43106533;43124017       43106536;43124096       ENST00000652672
```
In this example, the first ORF maps to two alternative genomic loci (the selected X:138230979-138231086 and the alternative X:180226-180333). The second ORF maps to the same genomic loci but there are two different possible 1st exon combinations (the selected two-exon ORF X:43121676-43121679 and X:43124017-43124096 and the alternative X:43106533-43106536 and X:43124017-43124096).

Example of an **unmapped** output file:
```
ENST00000075120_71_149--4_VanHeesch2019 unmapped        MIPARPRDAVMVRLWILPEDVEKTCC*
ENST00000616540_42_66--7_Gaertner2020   nanoORF MSARSPTS*
A1BG-AS1_58865324_38aa--6_Chen2020      NTG     LVHSAWCDPPRLWSRISTQVIQLRPALPRPTRDMCSVT*
```
In this example, the first ORF was not fully map to any transcript isoform (**unmapped**), the second ORF was shorter than the minimum length cut-off (**nanoORF**), and the third ORF started with a near-cognate ORF (**NTG**, these cases can be included with the option -a NTG or -a XTG). If **--multiple=no**, altmapped ORFs are also included in this file (tag **multiple_coords**). If **--max_len_cutoff** is used, excluded ORFs are also included in this file (tags **longORF**).



**TEST: GENERATING THE PHASE I RIBO-SEQ ORF LIST -Ensembl v.101-**
```
bash scripts/retrieve_ensembl_data.sh 101 GRCh38
python3 ORF_mapper_to_GENCODE.py -d Ens101 -f phaseI/sORFs_genomic_hg38.prot.20210329.fa -b phaseI/sORFs_genomic_hg38.prot.20210329.bed -a no
#Pseudogenes and annotated CDS should be excluded from the output
```
