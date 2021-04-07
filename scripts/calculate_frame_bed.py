import sys
import os

bed = sys.argv[1]
cds = "../../"
out = open(bed + "_psites.bed","w+")


if bed.endswith("bed"):
	phase = {}
	for line in open(bed):
		name = line.split("\t")[3]
		if not name in phase:
			if "\t+" in line:
				phase[name] = 0
			elif "\t-" in line:
				phase[name] = 2
		for coord in range(int(line.split("\t")[1]),int(line.split("\t")[2])+1):
			if phase[name] % 3 == 2:
				out.write(line.split("\t")[0] + "\t" + str(coord) + "\t" + str(coord) + "\t" + name + "\t" + str(phase[name]) + "\t" + line.split("\t")[5].rstrip("\n") + "\n")
			phase[name] += 1
elif bed.endswith("gtf"):
	phase = {}
	status = {}
	for line in open(bed):
		if "\tstart_codon\t" in line:
			if not line.split('transcript_id "')[1].split('"')[0] in status:
				status[line.split('transcript_id "')[1].split('"')[0]] =  0
			status[line.split('transcript_id "')[1].split('"')[0]] = status[line.split('transcript_id "')[1].split('"')[0]] + 1
		elif "\tstop_codon\t" in line:
			if not line.split('transcript_id "')[1].split('"')[0] in status:
				status[line.split('transcript_id "')[1].split('"')[0]] =  0
			status[line.split('transcript_id "')[1].split('"')[0]] = status[line.split('transcript_id "')[1].split('"')[0]] + 1
	for line in open(bed):	
		if not "\tCDS\t" in line:
			continue
		name = line.split('protein_id "')[1].split('"')[0]
		trans = line.split('transcript_id "')[1].split('"')[0]
		if not trans in status:
			continue
		if status[trans] < 2:
			continue
		if not name in phase:
			if "\t+" in line:
				phase[name] = 0
			elif "\t-" in line:
				phase[name] = 2
		for coord in range(int(line.split("\t")[3]),int(line.split("\t")[4])+1):
			if phase[name] % 3 == 2:
				out.write(line.split("\t")[0] + "\t" + str(coord) + "\t" + str(coord) + "\t" + name + "\t" + str(phase[name]) + "\t" + line.split("\t")[6] + "\n")
			phase[name] += 1	

out.close()
exit(0)
