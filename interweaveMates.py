########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2013 Sebastian Szpakowski
########################################################################################


#################################################
## 	Interwaeve two mates
#################################################
import sys
from optparse import OptionParser
from collections import defaultdict

_author="Sebastian Szpakowski"
_date="2011/01/01"
_version="Version 1"

#################################################
##		Classes
##
	#################################################
	### Iterator over input fasta file.
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

class	FastaLikeParser:
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
       		self.currentFastaSequence = "%s\n%s" % (self.currentFastaSequence, x.strip())			
       					
	def	__str__():
		return ("reading file: %s" %self.filename)		
		
		
#################################################
##		Functions
##

def parseMate(infile):
	seqs = dict()
	for head, seq in FastaParser(infile):
		id = head.strip().split("/")[0].strip().split()[0]
		seqs[id]= seq

	print "loaded %s reads" % (len(seqs.keys())) 
	return seqs

def	parseCluster(infile):
	repr = defaultdict(set)
	for clusterid, contents in FastaLikeParser(infile):
		representative= ""
		descendants = set()
		for line in contents.strip().split("\n"):
			id = line.split(">")[1].split("...")[0].split()[0]
			id = id.strip().split("/")[0]
			if line.endswith("*"):
				representative = id
			
			### self is always a descendant of itself!
			descendants.add(id)
			
		for k in descendants:
			repr[k]=descendants
	
	return (repr)
	
def	mergeClusters(c1, c2):
	unique = set()
	all = set(c1.keys()) | set(c2.keys())
	for id in all:
		A1 = c1[id]
		A2 = c2[id]
		
		representative = A1 & A2
		representative = list(representative)
		representative.sort()
		if len(representative)>0:
			representative = representative[0]
			unique.add(representative)
		else:
			representative = ""
			
		nondups = A1 ^ A2
		
		unique = unique | nondups
# 
# 		if representative !="":
# 			print id, c1[id], c2[id],
# 			print representative, nondups
# 			print 
	return (unique)
#################################################
##		Arguments
##
parser = OptionParser()

parser.add_option("-f", "--fasta", dest="fn_fasta",
                  help="write report to FILE", metavar="FILE")
                  
parser.add_option("-c", "--clusters", dest="fn_cluster", default="",
                  help="write report to FILE", metavar="FILE")   
                  
#parser.add_option("-q", "--quiet",
#                  action="store_false", dest="verbose", default=True,
#                  help="don't print status messages to stdout")


(options, args) = parser.parse_args()

#################################################
##		Begin
##

file1, file2 = options.fn_fasta.strip().split(",")
seq1 = parseMate(file1)	
seq2 = parseMate(file2)

A = set(seq1.keys()).intersection(set(seq2.keys()))
print "\t%s sequences mated." % (len(A))

clusters = options.fn_cluster.strip().split(",")
if len(clusters)==2:
	c1, c2 = clusters
	c1 = parseCluster(c1)
	c2 = parseCluster(c2)

	unique = mergeClusters(c1, c2)
	print "\t%s unique mates." % (len(unique))
	
	A = A & unique

otptfilename = "%s.weave.fasta" % ".".join(file1.strip().split(".")[:-1])

otpt = open(otptfilename, "w")
for id in A:
	line =  ">%s/1\n%s\n>%s/2\n%s\n" % (id.strip().split()[0], seq1[id], id.strip().split()[0], seq2[id])
	otpt.write(line)
otpt.close()		
#################################################
##		Finish
#################################################
