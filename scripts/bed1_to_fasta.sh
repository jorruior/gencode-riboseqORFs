bed=$1 #ORFs in 1-based coordinate file
genome=$2 #Genome fasta file

gffread <(awk '{print $1"\texon\texon\t"$2"\t"$3"\t.\t"$6"\t.\t gene_id \""$4"--"$5"\"; transcript_id \""$4"--"$5"\";"}' $bed) -g $genome -w $bed.nucl.fa
python3 translate_orfs.py $bed.nucl.fa > $bed.prot.fa
