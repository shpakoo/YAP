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
#################################################
##		Functions
##

#################################################
##		Arguments
##
parser = OptionParser()

#parser.add_option("-f", "--file", dest="filename",
#                  help="write report to FILE", metavar="FILE")
#parser.add_option("-q", "--quiet",
#                  action="store_false", dest="verbose", default=True,
#                  help="don't print status messages to stdout")

parser.add_option("-f", "--fasta", dest="fn_fasta", default="",
                  help="load fasta FILE", metavar="FILE")
                  
parser.add_option("-l", "--list", dest="fn_list", default = "",
				  help="load list FILE", metavar="FILE")
				  
parser.add_option("-n", "--names", dest="fn_names", default = "",
				  help="load names FILE (optional)", metavar="FILE")				  
				  
parser.add_option("-g", "--groups", dest="fn_groups",
                  help="load group FILE", metavar="FILE")                  

parser.add_option("-o", "--output", dest="fn_output",
                  help="load and clean fasta FILE", metavar="FILE")



(options, args) = parser.parse_args()


#################################################
##		Begin
##
groups = dict()
for id, group in GeneralPurposeParser(options.fn_groups, sep="\t"):
	groups[id] = group
otpt = open(options.fn_output, "w")
scratch = open("%s.missing" % (options.fn_groups), "w" )

if options.fn_fasta!="":
	for head, seq in FastaParser(options.fn_fasta):
		if groups.has_key(head):
			otpt.write("%s\t%s\n" % (head, groups[head]))
		else:
			scratch.write("%s\n" % (head) )
	if options.fn_names!=""	
		for name, merged in GeneralPurposeParser(options.fn_names, sep="\t"):
			for m in merged.split(","):
				if groups.has_key(head):
					otpt.write("%s\t%s\n" % (head, groups[head]))
				else:
					scratch.write("%s\n" % (head) )	
	otpt.close()
	
elif options.fn_list!="":
	tooutput=set()
	for line in GeneralPurposeParser(options.fn_list, sep="\t"):
		for k in line[2:]:
			for z in k.strip("").split(","):
				tooutput.add(z)
	for o in tooutput:
		if groups.has_key(o):
			otpt.write("%s\t%s\n" % (o, groups[o]))
	
	otpt.close()	
scratch.close()

#################################################
##		Finish
#################################################
