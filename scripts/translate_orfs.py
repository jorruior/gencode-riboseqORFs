#! /usr/bin/python
__author__="jruiz"
__date__ ="$Aug 27, 2013 15:18:00 PM$"
'''Script to check how many SNPs are in a set of transcripts/ORFs, and computing some density and PN/PS stats.
Argument 1: GTF file with transcripts
Argument 2: GTF file with ORFs
Argument 3: BED file of SNPs (one chromosome)
'''

import sys
import re

try:
	sys.path.append("/home/jruiz/lib/python/")
	sys.path.append("/home/jruiz/lib64/python/")
except:
	pass

try:
	sys.path.append("/data/users/jruiz/lib/python/")
	sys.path.append("/data/users/jruiz/lib64/python/")
except:
	pass

import os
import Bio
from Bio import SeqIO
from Bio.Seq import Seq

##ARGUMENTS

fasta = sys.argv[1]
seqs = SeqIO.index(fasta, "fasta")

for seq in seqs:
	prot = seqs[seq].seq.translate()
	print(">" + seq + "\n" + str(prot) + "\n", end = "")

exit(0)
