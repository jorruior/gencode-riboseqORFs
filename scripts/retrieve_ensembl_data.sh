#First field: Number of Ensembl release (e.g. 103)
#Second field: Genome assembly (e.g. GRCh38)

mkdir Ens$1
cd Ens$1

wget http://ftp.ensembl.org/pub/release-$1/gtf/homo_sapiens/Homo_sapiens.$2.$1.gtf.gz
wget http://ftp.ensembl.org/pub/release-$1/fasta/homo_sapiens/pep/Homo_sapiens.$2.pep.all.fa.gz
wget http://ftp.ensembl.org/pub/release-$1/fasta/homo_sapiens/cdna/Homo_sapiens.$2.cdna.all.fa.gz
wget http://ftp.ensembl.org/pub/release-$1/fasta/homo_sapiens/ncrna/Homo_sapiens.$2.ncrna.fa.gz

gunzip -d *gz

sort -k1,1 -k4,4n -k5,5n Homo_sapiens.$2.$1.gtf > Homo_sapiens.$2.sorted.gtf
cat Homo_sapiens.$2.cdna.all.fa Homo_sapiens.$2.ncrna.fa | cut -d"." -f1,1 > tmpfile; mv tmpfile Homo_sapiens.$2.trans.fa
cut -d"." -f1,1 Homo_sapiens.$2.pep.all.fa > tmpfile; mv tmpfile Homo_sapiens.$2.pep.all.fa
cut -d"." -f1,1 Homo_sapiens.$2.trans.fa > tmpfile; mv tmpfile Homo_sapiens.$2.trans.fa
python3 ../scripts/calculate_frame_bed.py Homo_sapiens.$2.sorted.gtf

#Change biomart to the corresponding version
wget -O ENST_support.txt 'http://ensembl.org/biomart/martservice?query=<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="0" count="" datasetConfigVersion="0.6">     
<Dataset name="hsapiens_gene_ensembl" interface="default">
    <Attribute name="ensembl_transcript_id"/>
    <Attribute name="transcript_tsl"/>
    <Attribute name="transcript_appris"/>
</Dataset>
</Query>'
