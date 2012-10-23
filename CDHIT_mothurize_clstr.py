########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2012 J.Craig Venter Institute.
########################################################################################



#!/usr/bin/python
#################################################
## 	create or merge with existing names file an output from
##	CD-HIT
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
	### Iterator over input file.
	### every line is converted into a dictionary with variables referred to by their 
	### header name
class 	GeneralPurposeParser:
	def	__init__(self, file, skip=0, sep="\t"):
		self.filename = file
		self.fp = open(self.filename, "r")	
		self.sep = sep
		self.linecounter = 0
		self.currline=""
		
	def	__iter__(self):
		return (self)
	
	def	next(self):
		otpt = dict()
		for currline in self.fp:
			currline = currline.strip().split(self.sep)
			self.currline = currline
			self.linecounter = self.linecounter + 1
			return(currline)			
		raise StopIteration
					
	def	__str__(self):
		return "%s [%s]\n\t%s" % (self.filename, self.linecounter, self.currline)
		
		
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

#################################################
##		Arguments
##
parser = OptionParser()

parser.add_option("-c", "--cluster", dest="fn_clstr",
                 help="clustering from CDHIT", metavar="FILE")

parser.add_option("-n", "--names", dest="fn_name", default ="",
                 help="names from CDHIT", metavar="FILE")
                 
parser.add_option("-o", "--output_mode", dest="out_mode", default = "name",
                 help="output mode, i.e. either output list or names files", metavar="FILE")
                 
#parser.add_option("-q", "--quiet",
#                  action="store_false", dest="verbose", default=True,
#                  help="don't print status messages to stdout")


(options, args) = parser.parse_args()

#################################################
##		Begin
##

names = defaultdict(set)
outputnames = defaultdict(set)

if options.fn_name != "" :
	for id, desc in GeneralPurposeParser(options.fn_name, skip=0, sep="\t"):
		[names[id].add(x) for x in desc.split(",")]
				
for clusterid, contents in FastaLikeParser(options.fn_clstr):
	representative= ""
	descendants = set()
	for line in contents.strip().split("\n"):
		line = line.strip()
		id = line.split(">")[1].split("...")[0].split()[0]
		if line.endswith("*"):
			representative = id
			#print "repr"
		### self is always a descendant of itself!
		descendants.add(id)
		[descendants.add(x) for x in  names[id]]	
		
	outputnames[representative] = descendants
			
if options.out_mode =="name":
	prefix = options.fn_clstr.strip().split("/")[-1]
	otptfile = open("%s.name" % (prefix), "w")	
	
	for key, vals in outputnames.items():
		otptfile.write("%s\t%s" % (key, key))
		
		vals.discard(key)
		if len(vals)>0:
			otptfile.write(",%s\n" % (",".join(list(vals))))
		else:
			otptfile.write("\n")
		
	otptfile.close()	
	
else:
	otus = defaultdict(list)
	prefix = options.fn_clstr.strip().split("/")[-1]
	listfile = open("%s.list" % (prefix), "w" )
	listfile.write("%s\t%s" % (options.out_mode, len(outputnames.keys() )))
	
	rafile =open("%s.rabund" % (prefix), "w")
	rafile.write("%s\t%s" % (options.out_mode, len(outputnames.keys() )))
	
	safile = open("%s.sabund" % (prefix), "w")
	safile.write("%s" % (options.out_mode))
	
	mx = 0
	for key, vals in outputnames.items():
		otus[len(vals)].append(vals)
		
		### the OTUs must start with a representative, followed by rest
		#print key, vals
		vals.remove(key)
		curotu = "%s,%s" % (key, (",").join(list(vals)))
		listfile.write("\t%s" % curotu.strip(","))
		
		### we removed the representative, add one!)
		rafile.write("\t%s" % (len(vals)+1))
		mx = max (mx, len(vals)+1)
	
	safile.write("\t%s" % (mx) )
	for x in range (max(otus.keys())):
		x = x+1
		safile.write("\t%s" % len(otus[x]))
		
	listfile.write("\n")	
	rafile.write("\n")	
	safile.write("\n")	
	
	
	listfile.close()
	rafile.close()
	safile.close()	
	

#################################################
##		Finish
#################################################
