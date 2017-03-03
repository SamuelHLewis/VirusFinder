#!/usr/bin/env python3

import os
import sys
import subprocess
import random as random
import shutil
import argparse
import re

#############################
## DEFAULT ARGUMENT VALUES ##
#############################
# number of cores to use
Cores = 1
# minimum length of viral contigs
MinContigLength = 100

###########################
## USER ARGUMENT PARSING ##
###########################
parser = argparse.ArgumentParser(description='Read arguments')
parser.add_argument('-l', '--infile1', type=str, help='input RNASeq file 1 (left reads)')
parser.add_argument('-r', '--infile2', type=str, help='input RNASeq file 2 (right reads)')
parser.add_argument('-g', '--genomeDB', type=str, help='path to bowtie2 index for host genome (to map reads to initially)')
parser.add_argument('-v', '--virusDB', type=str, help='path to viral BLAST database (must be built using DIAMOND)')
parser.add_argument('-n', '--nrDB', type=str, help='path to nr BLAST database (must be built using DIAMOND)')
parser.add_argument('-m', '--minlength',  type=int, help='minimum length of viral contigs')
parser.add_argument('-c', '--cores',  type=int, help='number of cores')
args = parser.parse_args()
# input fasta file parsing
InFile1 = args.infile1
if InFile1 is not None:
	if InFile1.endswith('.fq') or InFile1.endswith('.fastq') or InFile1.endswith('.fa') or InFile1.endswith('.fas') or InFile1.endswith('.fasta') or InFile1.endswith('.fq.gz') or InFile1.endswith('.fastq.gz') or InFile1.endswith('.fa.gz') or InFile1.endswith('.fas.gz') or InFile1.endswith('.fasta.gz'):
		print('Input RNASeq file 1 (left reads) is ' + InFile1)
	else:
		print('Input RNASeq file 1 (left reads) ' + InFile1 + ' does not have normal fasta/fastq filename ending - please check file format')
		sys.exit(0)
else:
	print('ERROR: no input RNASeq file 1 (-l) specified')
	sys.exit(0)
# input fasta file parsing
InFile2 = args.infile2
if InFile2 is not None:
	if InFile2.endswith('.fq') or InFile2.endswith('.fastq') or InFile2.endswith('.fa') or InFile2.endswith('.fas') or InFile2.endswith('.fasta') or InFile2.endswith('.fq.gz') or InFile2.endswith('.fastq.gz') or InFile2.endswith('.fa.gz') or InFile2.endswith('.fas.gz') or InFile2.endswith('.fasta.gz'):
		print('Input RNASeq file 2 (right reads) is ' + InFile2)
	else:
		print('Input RNASeq file 2 (right reads) ' + InFile2 + ' does not have normal fasta/fastq filename ending - please check file format')
		sys.exit(0)
else:
	print('ERROR: no input RNASeq file 2 (-r) specified')
	sys.exit(0)
# host genome bowtie2 index parsing
GenomeDB = args.genomeDB
if GenomeDB is not None:
	print('Path to host genome bowtie2 index is ' + GenomeDB)
else:
	print('ERROR: no host genome bowtie2 index (-g) specified')
	sys.exit(0)
# viral BLAST database parsing
VirusDB = args.virusDB
if VirusDB is not None:
	print('Path to viral BLAST database is ' + VirusDB)
else:
	print('ERROR: no viral BLAST database (-v) specified')
	sys.exit(0)
# non-redundant protein BLAST database parsing
nrDB = args.nrDB
if nrDB is not None:
	print('Path to non-redundant protein BLAST database is ' + nrDB)
else:
	print('ERROR: no non-redundant protein BLAST database (-n) specified')
	sys.exit(0)
# minimum viral contig length parsing
if args.minlength is None:
	print('Using default minimum length for viral contigs (' + str(MinContigLength) + ')')
else:
	MinContigLength = args.minlength
	if MinContigLength > 0:
		print('minimum length for viral contigs set to ' + str(MinContigLength))
	else:
		print('ERROR: minimum length for viral contigs (-m) must be >0')
		sys.exit(0)
# cores parsing
if args.cores is None:
	print('Using default number of cores (' + str(Cores) + ')')
else:
	Cores = args.cores
	if Cores > 0:
		print(str(Cores) + ' cores specified')
	else:
		print('ERROR: 1 or more cores (-c) required')
		sys.exit(0)

# function to map reads using bowtie2 (requires bowtie2 mapping databases to exist already, discards unmapped reads by default)
def BowtieMapper(InFile1,InFile2,Index,MappedOut,UnmappedOut,cores=1):
	# map fastq reads
	print('Mapping ' + InFile1 + ' & ' + InFile2 + ' to index ' + Index)
	cmd = 'bowtie2 --fast -p ' + str(cores) + ' -q -x ' + Index + ' -1 ' + InFile1 + ' -2 ' + InFile2 + ' -S ' + MappedOut + '.sam --un-conc ' + UnmappedOut
	subprocess.call(cmd,shell=True)
	# rename unmapped read files
	shutil.move(UnmappedOut+'.1',UnmappedOut+'.1.fq')
	shutil.move(UnmappedOut+'.2',UnmappedOut+'.2.fq')
	print('Mapped ' + InFile1 + ' & ' + InFile2 + ' to index ' + Index + ', unmapped reads to ' + UnmappedOut + '.1.fq & ' + UnmappedOut + '.2.fq')
	return(UnmappedOut+'.1.fq',UnmappedOut+'.2.fq')

# function to assembly paired-end reads with Trinity
def TrinityAssembler(InFile1,InFile2,MinContigLength,cores=1):
	# remove TrinityAssembly dir if already exists and replace it with new (blank) dir
	if os.path.isdir('TrinityAssembly') is True:
		shutil.rmtree('TrinityAssembly')
		print('Old TrinityAssembly tempfiles removed')
	os.makedirs('TrinityAssembly')
	# assemble NonHost reads
	print('Assembling ' + InFile1 + ' & ' + InFile2)
	cmd = 'Trinity --left ' + InFile1 + ' --right ' + InFile2 + ' --seqType fq --max_memory 10G --min_contig_length ' + str(MinContigLength) + ' --output ./TrinityAssembly --CPU ' + str(cores)
	subprocess.call(cmd,shell=True)
	# extract longest ORF from each transcript, rename this file, and move other files
	cmd = 'TransDecoder.LongOrfs -t ./TrinityAssembly/Trinity.fasta -m ' + str(int(MinContigLength/3))
	subprocess.call(cmd,shell=True)
	shutil.move('./Trinity.fasta.transdecoder_dir/longest_orfs.cds','./Trinity_min50.fasta')
	shutil.move('./Trinity.fasta.transdecoder_dir/','./TrinityAssembly/')
	print('Assembled contigs to Trinity_min50.fasta')
	return('./Trinity_min50.fasta')

# function to assembly paired-end reads with Trinity
def BLASTer(InFile,BlastDB,OutFile,cores=1):
	# blastx
	cmd = 'diamond blastx -p ' + str(cores) + ' -d ' + BlastDB + ' -q ' + InFile + ' -o ' + OutFile + ' -f 6 qseqid evalue length stitle'
	subprocess.call(cmd,shell=True)
	return(OutFile)

# master function to identify and extract viral contigs
def VirusFinder(InFile1,InFile2,HostGenome,VirusDB,nrDB,MinContigLength,Threads=1):
	# map reads to host genome and keep non-host reads
	NonHost1,NonHost2 = BowtieMapper(InFile1=InFile1,InFile2=InFile2,Index=HostGenome,MappedOut='Host',UnmappedOut='NonHost',cores=Threads)
	# assemble non-host reads
	Contigs = TrinityAssembler(InFile1=NonHost1,InFile2=NonHost2,MinContigLength=MinContigLength,cores=Threads)
	# blast assembled contigs against NCBI viruses
	VirusBLASTOut = BLASTer(InFile='Trinity_min50.fasta',BlastDB=VirusDB,OutFile='ncbi_viruses.out',cores=Threads)
	# check whether potential viruses have been found
	fileinfo = os.stat(VirusBLASTOut)
	if fileinfo.st_size == 0:
		print('No viral candidates found')
		return('No viral candidates found')
	else:
		# count number of viral candidates
		count = 0
		for line in open(VirusBLASTOut):
			count += 1
		print(str(count) + ' viral candidates found')
		# grab sequences of viral candidates
		cmd = 'cat ' + VirusBLASTOut + ' | cut -f 1 | sort | uniq | grep -A 1 --no-group-separator -F -f - Trinity_min50.fasta > candidates.fasta'
		subprocess.call(cmd,shell=True)
		# blast the viral candidates against NCBI nr
		BLASTer(InFile='candidates.fasta',BlastDB=nrDB,OutFile='candidates_ncbi_nr.out',cores=Threads)
		# grab all candidate names and print to temp file
		cmd = 'grep ">" candidates.fasta | grep -o "^[^ ]*" | sed \'s/>//\' > candnames.txt'
		subprocess.call(cmd,shell=True)
		# for each candidate name, print the top blast hit to file
		for candidate in open('candnames.txt'):
			cmd = 'grep -m 1 ' + candidate.strip('\n') + ' candidates_ncbi_nr.out >> topblasthits.txt'
			subprocess.call(cmd,shell=True)
		# print to file the names of all candidates whose top hit was a virus, excluding those coming from transposons
		cmd = 'grep \'[Vv]irus\' topblasthits.txt | grep -v \'transpos*\' | cut -f 1 | sort | uniq > verified_viruses.txt'
		subprocess.call(cmd,shell=True)
		# check whether there are any verified viruses 
		fileinfo = os.stat('verified_viruses.txt')
		if fileinfo.st_size == 0:
			print('No verified viruses found')
			# remove temp files
			os.remove('candnames.txt')
			os.remove('topblasthits.txt')
			return('No verified viruses found')
		else:
			# print to file the name and sequence of all candidates whose top hit was a virus (& remove temp files)
			cmd = 'cat verified_viruses.txt | grep -A 1 --no-group-separator -F -f - Trinity_min50.fasta > verified_viruses.fasta'
			subprocess.call(cmd,shell=True)
			# grab name of each candidate and its top blast hit
			CandidatesAndHits = {}
			for candidate in open('candidates_ncbi_nr.out'):
				if not candidate.split('\t')[0] in CandidatesAndHits:
					CandidatesAndHits[candidate.split('\t')[0]] = candidate.split('\t')[3].strip('\n')
			# grab each verified viral contig and its top blast hit
			VerifiedAndHits = {}
			for contig in open('verified_viruses.txt'):
				VerifiedAndHits[contig.strip('\n')] = CandidatesAndHits[contig.strip('\n')]
			# make tab-delimited output of each verified contig and its top blast hit
			AnnotatedNames = 'Contig\tTop hit\n'
			for contig in VerifiedAndHits:
				entry = contig + '\t' + VerifiedAndHits[contig] + '\n'
				AnnotatedNames+=entry
			# output compound names to file
			CompoundOut = open('compoundnames.txt','wt')
			CompoundOut.write(AnnotatedNames)
			CompoundOut.close()
			# remove temp files
			os.remove('candnames.txt')
			os.remove('topblasthits.txt')
			# count number of verified viruses
			count = 0
			for i in open('verified_viruses.txt'):
				count+=1
			print(str(count) + ' verified viruses found')
			return(str(count) + ' verified viruses found')

VirusFinder(InFile1=InFile1,InFile2=InFile2,HostGenome=GenomeDB,VirusDB=VirusDB,nrDB=nrDB,MinContigLength=MinContigLength,Threads=Cores)

