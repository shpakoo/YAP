########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2012 Sebastian Szpakowski, J.Craig Venter Institute.
########################################################################################



#!/usr/bin/python
#################################################
## 	A new program
#################################################
import sys, shlex, glob, time
from optparse import OptionParser
from string import maketrans
from subprocess import *


_author="Sebastian Szpakowski"
_date="2011/01/01"
_version="Version 1"

#################################################
##		Classes
##
	#################################################
	### Iterator over input file.
	### every line is converted into a dictionary with variables referred to by their 
	### header name
class 	GeneralPurposeParser:
	def	__init__(self, file, skip=0, sep="\t"):
		self.filename = file
		self.fp = open(self.filename, "r")	
		self.linecounter = 0
		self.sep=sep
	
		while skip>0:
			skip = skip-1
			self.next()
			
	def	__iter__(self):
		return (self)
	
	def	next(self):
		for currline in self.fp:
			currline = currline.strip().split(self.sep)
			self.linecounter = self.linecounter + 1
			return(currline)					
		raise StopIteration
					
	def	__str__(self):
		count = 0
					
		return ("%s\n[LINE %s]" % (self.filename, self.linecounter))
					
	#################################################
	### Iterator over input fata file.
	### Only reading when requested
	### Useful for very large FASTA files
	### with many sequences
	
class	FastaParser:
	def	__init__ (self, x):
		self.filename = x
		self.fp = open(x, "r")	
		self.currline = "" 
		self.currentFastaName = ""
		self.currentFastaSequence = ""
		self.lastitem=False	
			
	def	__iter__(self):
		return(self)	
				
		##### 
	def	next(self):
		for self.currline in self.fp:
			if self.currline.startswith(">"):
				self.currline = self.currline[1:]
				if self.currentFastaName == "":
					self.currentFastaName = self.currline
				else:
					otpt = (self.currentFastaName.strip(), self.currentFastaSequence.strip())
					self.currentFastaName = self.currline
					self.currentFastaSequence = ""	
					self.previoustell = self.fp.tell()
					return (otpt)
				
			else:
				self.addSequence(self.currline)	
		
		if not self.lastitem:
			self.lastitem=True			
			return (self.currentFastaName.strip(), self.currentFastaSequence.strip())
		else:
			raise StopIteration	
			       				
       	def	addSequence(self, x):
       		self.currentFastaSequence = "%s%s" % (self.currentFastaSequence, x.strip())			
       					
	
	def	__str__():
		return ("reading file: %s" %self.filename)	

		
		

#################################################
##		Functions
##
def	qsub(num):
	batches.append(num)
	command = "qsub -P %s -N tmp.%s.in -cwd -M sszpakow@jcvi.org -m a \"/home/sszpakow/bin/python %s -i tmp.%s.in -o tmp.%s.out -s -p REFSEQ \"" % (options.projectid, num, sys.argv[0], num, num   )
	
	#command = "qsub -P %s -N tmp.%s.in -cwd -l \"fast\" -M sszpakow@jcvi.org -m a \"/home/sszpakow/bin/python %s -i tmp.%s.in -o tmp.%s.out -s -p REFSEQ \"" % (options.projectid, num, sys.argv[0], num, num   )
	# 	   qsub -P 810013 -N tmp.2.in -cwd -l "fast" -M sszpakow@jcvi.org -m a "/home/sszpakow/bin/python /home/sszpakow/scripts/summarizeAlignment.py -i tmp.2.in -o tmp.2.out -s -p REFSEQ "
	#print shlex.split(command)
	p = Popen(shlex.split(command), stdout=PIPE)	
	
	print command
#################################################
##		Arguments
##
parser = OptionParser()

parser.add_option("-p", "--pattern", dest="pattern",
                  help="the first sequence name that starts with P will be used as a reference for coordinates", metavar="P")

parser.add_option("-i", "--input", dest="infilename",
                  help="NAME of an input file (fasta format)", metavar="NAME")

parser.add_option("-o", "--output", dest="outfilename",
                  help="NAME of an output file", metavar="NAME")

parser.add_option("-t", "--tmpsize", type="int", dest="size",
                  help="N fasta sequences will be put in each temporary file [default: %default]", metavar="N", default=10000)

parser.add_option("-s", "--slave",
                  action="store_false", dest="master", default=True,
                  help="qsub worker")
                  
parser.add_option("-P", "--project", dest="projectid", default="810013",
                  help="projectid", metavar="num")                  
                  
(options, args) = parser.parse_args()




#################################################
##		Begin
##


if options.master:
	#### init
	curname = 0
	batches = list()

	#### find pattern in the large input file
	file = FastaParser(options.infilename)
	refhed = ""
	seq= ""

	print "Searching for %s..." % (options.pattern)

	for head, seq in file:
		if head.startswith("%s" % (options.pattern)) and refhed == "":
			refhed = head
			refseq = seq
			print "\t^--found !"
			break
	
	#### create temporary files with the reference found and submit qsub jobs
	file = FastaParser(options.infilename)
	
	counter=0
	otpt = open("tmp.%s.in"% (curname), "w")
	otpt.write(">REFSEQ\n%s\n" % (refseq))
	
	for head, seq in file:
		if counter == options.size:		
			otpt.close()
			qsub(curname)
			curname +=1 
			otpt = open("tmp.%s.in"% (curname), "w")
			otpt.write(">REFSEQ\n%s\n" % (refseq))
			counter=0		
		otpt.write(">%s\n%s\n" % (head, seq))
		
		counter+=1	
					
	otpt.close()		
	qsub(curname)		
	
	#### wait until the jobs finish, merge whatever is available
	
	c1 = list()
	c2 = list()
	c3 = list()
	
	previousnums=set()
	while len(batches)>0:
		for file in glob.glob("*.out.done"):
			print file
			cols = dict()
			for line in GeneralPurposeParser(file[:-5]):
				
				#line = line.strip().split("\t")
				
				for column, val in enumerate(line):
					if not cols.has_key(column):
						cols[column]=list()
												
					cols[column].append(int(val))
			
			#print cols.keys()
			
			if len(c1)==0:
				c1 = cols[0]
				c2 = cols[1]
				c3 = cols[2]
			else:
				tmp = cols[2]
				
				c3 = [c3[index]+val for index, val in enumerate(tmp)]
			
			
			batch = file.strip().split(".")[1]			
			batches.remove(int(batch))
			
			for XXX in glob.glob("tmp.%s.*" % (batch)):
				print "-%s" % (XXX)
				p = Popen(shlex.split("rm %s" % (XXX)),  stdout=PIPE)

		if len(batches) in previousnums:
			pass
		else:	
			print "waiting for %s more..." % (len(batches))
			previousnums.add(len(batches))
			
		time.sleep(5)	
			
	#### output results
	otpt = open(options.outfilename, "w")
	for index, val in enumerate(c1):
		line = "%s\t%s\t%s\n" % (c1[index], c2[index], c3[index])
		otpt.write(line)	
	otpt.close()
else:

	file = FastaParser(options.infilename)
	coords = list()	
	inttab = "acgtACGTN.-"
	outtab = "11111111000"
	transtab = maketrans(inttab, outtab)
	refhed = ""
	
	print "Searching for %s..." % (options.pattern)
	for head, seq in file:
		if head.startswith("%s" % (options.pattern)) and refhed == "":
			refhed = head
			refseq = list(seq.translate(transtab))
			refseq = [int(y) for y in refseq]
			print "\t^--found !"
		else:
			seq = list(seq.translate(transtab))
			seq = [int(y) for y in seq]
			#print seq
			if len(coords)==0:
				coords=seq
			else:			
				coords = [coords[index]+val for index, val in enumerate(seq)]
	
				
	counter=0	
	otpt = open(options.outfilename, 'w')
	for index, val in enumerate(refseq):
		if val==1:
			otpt.write("%s\t%s\t%s\n" % (index+1, counter+1, coords[index] ) )
			counter+=1		
	otpt.close()

	open("%s.done" % (options.outfilename), 'w')
	

		
#################################################
##		Finish
#################################################
