########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2012 J.Craig Venter Institute.
########################################################################################



#!/usr/bin/python
#################################################
## 	A new program
#################################################
import sys
from optparse import OptionParser

_author="Sebastian Szpakowski"
_date="2011/01/01"
_version="Version 1"

#################################################
##		Classes
##
	#################################################
	### Iterator over input fata file.
	### Only reading when requested
	### Useful for very large FASTA files
	### with many sequences
class	FastaParser:
	def	__init__ (self, x, quals=False):
		self.filename = x
		self.fp = open(x, "r")	
		self.currline = "" 
		self.currentFastaName = ""
		self.currentFastaSequence = ""
		self.lastitem=False	
		if quals:
			self.linesep=" "	
		else:
			self.linesep=""	
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
       		self.currentFastaSequence = "%s%s%s" % (self.currentFastaSequence,self.linesep, x.strip())			
       					
	
	def	__str__():
		return ("reading file: %s" %self.filename)	
		
#################################################
##		Functions
##
	################################################
	### Read in a file and return a list of lines
	###
def	loadLines(x):
	try:
		fp = open(x, "r")
		cont=fp.readlines()
		fp.close()
		#print "%s line(s) loaded."  % (len(cont))
	except:
		cont=""
		#print "%s cannot be opened, does it exist? " % ( x )	
	return cont
	
	
	
#################################################
##		Arguments
##
parser = OptionParser()

#parser.add_option("-f", "--file", dest="filename",
#                  help="write report to FILE", metavar="FILE")
#parser.add_option("-q", "--quiet",
#                  action="store_false", dest="verbose", default=True,
#                  help="don't print status messages to stdout")


(options, args) = parser.parse_args()

#################################################
##		Begin
##

### need list file
### need taxonomy file
### need group file	
### need fasta file

### these two use OTUid 
OTUs = dict()
TAXONOMY = dict()
### link seqid to OTUid
SEQid = dict()
### keep track of groups
GROUPS = dict()


counter = 0
#list
for line in loadLines(sys.argv[2]):
	line = line.strip().split("\t")
	dist = line[0]
	ids = line[1:]
	if dist == sys.argv[1] :	
			for idgroup in ids:
				OTUs[counter] = idgroup.split(",")
				for id in idgroup.split(","): 
					SEQid[id]=counter
				counter +=1
				
			
#taxonomy, skip header			
for line in loadLines(sys.argv[3])[1:]:
	OTUid, size, classification = line.strip().split("\t")
	OTUid = int(OTUid)
	
	TAXONOMY[OTUid] = classification		
	#print TAXONOMY[OTUid], int(size), len( OTUs[OTUid])


# GROUPS
for line in loadLines(sys.argv[4]):
	id, group = line.strip().split("\t")
	GROUPS[id]=group

# fasta
for head, seq in FastaParser(sys.argv[5]):
	tmp = set()
	### OTUs	
	for K in OTUs[SEQid[head]]:
		tmp.add(GROUPS[K])
	tmp = list(tmp)
	tmp.sort()
	otpt = ">%s|%s|%s|%s|%s\n%s" % ( head, TAXONOMY[SEQid[head]],  ";".join(tmp), len(OTUs[SEQid[head]]), len(tmp), seq)
	
	print otpt

#################################################
##		Finish
#################################################
