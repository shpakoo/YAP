########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2013 Sebastian Szpakowski
########################################################################################

#################################################
## rename anything with "uncultured" with a species that is most abundant in a cluster
#################################################
import sys
from optparse import OptionParser
from collections import defaultdict

_author="Sebastian Szpakowski"
_date="2012/03/30"
_version="Version 1"

#################################################
##		Classes
##

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
def	findName(names):
	tmp = defaultdict(int)
	for n in names:
		if n.lower().find("uncultured")==-1 and n.lower().find("associated")==-1:
			n = "_".join(n.split("_")[1:])
			n = "%s" % (n)
			tmp[n]+=1
	
	if len(tmp.keys())>0:
		otpt = defaultdict(list)	
		for n,count in tmp.items():
			otpt[count].append(n)
			
		keys = otpt.keys()
		keys.sort()
		keys.reverse()
		newname = otpt[keys[0]]
		return newname	
	else:
		return list()
	

#################################################
##		Arguments
##
parser = OptionParser()

parser.add_option("-c", "--clusters", dest="fn_clstr",
                 help="clustering from CDHIT", metavar="FILE")

parser.add_option("-s", "--sequences", dest="fn_seqs", default ="",
                 help="names from CDHIT", metavar="FILE")
                 
parser.add_option("-o", "--output", dest="fn_out", default ="",
                 help="names from CDHIT", metavar="FILE")                 

(options, args) = parser.parse_args()

#################################################
##		Begin
##

names = defaultdict(set)
clusters = defaultdict(set)
remapping = defaultdict(list)

## parse clusters
for clusterid, contents in FastaLikeParser(options.fn_clstr):
	representative= ""
	descendants = set()
	for line in contents.strip().split("\n"):
		id = line.split(">")[1].split("...")[0]
		if line.endswith("*"):
			representative = id
		
		### self is always a descendant of itself!
		descendants.add(id)
		[descendants.add(x) for x in  names[id]]
	clusters[representative] = descendants

### find new names
for K in clusters:
	if K.lower().find("uncultured")>-1 or K.lower().find("associated")>-1:
		remapping[K] = findName(clusters[K])

## rename fasta	
otptfile = open(options.fn_out, "w")
for head, seq in FastaParser(options.fn_seqs):
	if len(remapping[head])>0:
		newname = "%s_%s" % (head.split("_")[0], remapping[head][0])
		print head, "->" , newname, remapping[head]

	else:
		newname = head
		
	otptfile.write(">%s\n%s\n" % (newname, seq))	
otptfile.close		
		



#################################################
##		Finish
#################################################
