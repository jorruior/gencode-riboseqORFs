#!/usr/bin/env python
import sys
import string
import subprocess
import os
from Bio import SeqIO
from Bio.Seq import Seq
from datetime import datetime
from optparse import OptionParser

import functions

__author__ = "Jorge Ruiz-Orera"
__contributor__="..."
__copyright__ = ""
__credits__ = []
__license__ = ""
__version__="1.0.0"
__maintainer__ = "Jorge Ruiz-Orera"
__email__ = "jorruior@gmail.com"


def check_arg (arg_str,s):
	'''Check if arg was written'''
	if not arg_str:   # if filename is not given
		print("Error: " + str(s) + " argument not given\n")
		exit()

def check_file (file_str):
	'''Check if input really exists'''
	try:
		open("%s" %file_str)
	except:
		print("Error: " + file_str + " input not found\n")
		exit()

def main():
	usage = "\n%prog  [options]"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-d","--input_dir",action="store",dest="folder",help="(Required) Directory with required Ensembl files (transcriptome gft and fasta, proteome fasta, tab file with APPRIS and TSL support, protein psites generated by 'calculate_frame_bed.py'")
	parser.add_option("-f","--input_fasta",action="store",dest="orfs_fa_file",help="(Required) File with all translated candidate ORFs. A FASTA can be generated from a BED file using the script 'bed1_to_fasta.sh'")
	parser.add_option("-b","--input_bed",action="store",dest="orfs_bed_file",help="File with 1-based BED coordinates of all translated candidate ORFs. ORF names should match in both fasta and bed files. If -m is activated, the BED coordinates willbe written into this file.")
	parser.add_option("-o","--output",action="store",dest="out_name",default="input_bed",help="Tag name for the generated output files. Default: BED filename")	
	parser.add_option("-l","--len_cutoff",action="store",type=int,dest="len_cutoff",default="16",help="Minimum ORF length threshold in amino acids (without stop codon). default=16")
	parser.add_option("-c","--collapse_cutoff",action="store",type=float,dest="col_thr",default="0.9",help="Minimum required fraction to collapse ORFs with similar stretches of overlapping amino-acid sequences. default=0.9")
	parser.add_option("-m","--collapse_method",action="store",dest="method",default="longest_string",help="Method to cluster ORF variants. 'longest_string' will collapse ORFs if the longest shared string is above -c threshold (default). 'psite_overlap' is a slower method that collapses ORFs if the fraction of shared psites is above -c threshold")
	parser.add_option("-a","--make_annot_bed",action="store",dest="calculate_coordinates",default="no",help="If ORF BED file is not available, generate it from the FASTA file. (ATG/NTG/XTG/no, default = no)")
	parser.add_option("--multiple",action="store",dest="mult",default="yes",help="If -a option is activated, include ORFs that map to multiple genomic coordinates. (yes/no, default = yes)")

	(opt,args)=parser.parse_args()

	check_arg(opt.folder,"--folder")
	check_arg(opt.orfs_fa_file,"--input_fasta")
	check_arg(opt.orfs_bed_file,"--input_bed")
	check_file(opt.orfs_fa_file)
	if opt.calculate_coordinates != "no":
		print("A new BED file will be generated in " + opt.orfs_bed_file)
	else:
		check_file(opt.orfs_bed_file)
	folder = opt.folder
	orfs_fa_file = opt.orfs_fa_file
	orfs_bed_file = opt.orfs_bed_file
	len_cutoff = opt.len_cutoff
	col_thr = opt.col_thr
	method = opt.method
	calculate_coordinates = opt.calculate_coordinates
	if opt.out_name != "input_bed":
		out_name = orfs_bed_file
	else:
		out_name = opt.out_name
	if not calculate_coordinates in ("ATG","NTG","XTG","no"):
		print("Error: " + method + " is not a valid -a argument\n")	
		exit()
	if not method in ("longest_string","psite_overlap"):
		print("Error: " + calculate_coordinates + " is not a valid -m argument\n")		
		exit()
	if not opt.mult in ("yes","no"):
		print("Error: " + opt.mult + " is not a valid --multiple argument\n")
		exit()
	elif (calculate_coordinates == "no") and (opt.mult != "yes"):
		print("Error: " + opt.mult + " is not a valid --multiple argument when -a is not enabled\n")
		exit()
	else:
		mult = opt.mult

	#Annotation files
	transcriptome_fa_file = folder + "/Homo_sapiens.GRCh38.trans.fa" #include both mRNA and ncRNA, only transcript ID in header, e.g. cat Ens103/Homo_sapiens.GRCh38.cdna.all.fa Ens103/Homo_sapiens.GRCh38.ncrna.fa | cut -d"." -f1,1 > prueba; mv prueba Ens103/Homo_sapiens.GRCh38.trans.fa
	transcriptome_gtf_file =  folder + "/Homo_sapiens.GRCh38.sorted.gtf" #sort -k1,1 -k4,4n -k5,5n Ens103/Homo_sapiens.GRCh38.gtf > Ens103/Homo_sapiens.GRCh38.sorted.gtf
	psites_bed_file = folder + "/Homo_sapiens.GRCh38.sorted.gtf_psites.bed" #Run calculate_frame_bed.py on the sorted gtf
	proteome_fa_file = folder + "/Homo_sapiens.GRCh38.pep.all.fa" #include both mRNA and ncRNA, only protein ID in header. e.g. cut -d"." -f1,1 Ens103/Homo_sapiens.GRCh38.pep.all.fa > prueba; mv prueba Ens103/Homo_sapiens.GRCh38.pep.all.fa
	t_support = folder + "/ENST_support.txt" #Biomart: APPRIS and TSL annotation	

	check_file(transcriptome_fa_file)
	check_file(transcriptome_gtf_file)
	check_file(psites_bed_file)
	check_file(proteome_fa_file)
	check_file(t_support)

	try:
   		os.mkdir(folder + "/tmp/")
	except:
		pass 

	(orfs_fa,transcriptome_fa,proteome_fa) = functions.load_fasta(orfs_fa_file,transcriptome_fa_file,proteome_fa_file)
	gtf = functions.parse_gtf(transcriptome_gtf_file,"exon")
	if calculate_coordinates != "no":
		functions.make_bed(orfs_fa,transcriptome_fa,gtf,len_cutoff,calculate_coordinates,orfs_bed_file,out_name,mult)
	(appris,supp) = functions.read_support(t_support)
	(overlaps,overlaps_cds,other_overlaps,total_studies) = functions.insersect_orf_gtf(orfs_bed_file,transcriptome_gtf_file,folder)
	other_overlaps = functions.pseudo_or_cds_ov(orfs_bed_file,transcriptome_gtf_file,other_overlaps,folder)
	(candidates,trans_orfs,second_names,coord_psites) = functions.orf_tags(overlaps,overlaps_cds,orfs_fa,transcriptome_fa,proteome_fa,gtf,len_cutoff,folder)
	(exc,variants,variants_names,datasets) = functions.exclude_variants(trans_orfs,col_thr,candidates,method,coord_psites,folder)
	functions.write_output(orfs_bed_file,candidates,exc,variants,variants_names,datasets,appris,supp,gtf,transcriptome_fa,second_names,len_cutoff,col_thr,total_studies,other_overlaps,psites_bed_file,folder,out_name)


if __name__ == '__main__':
	main()

exit(0)

#Tmp folder file for multiple calculations and remove at the end


