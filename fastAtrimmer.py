########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2012 Sebastian Szpakowski, J.Craig Venter Institute.
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

#################################################
##		Arguments
##
parser = OptionParser()

parser.add_option("-e", "--end", dest="end", type="int",
                  help="output fastA or quals up to END position", metavar="END")
parser.add_option("-l", "--lucy",
                  action="store_true", dest="lucymode", default=False,
                  help="lucy files")


(options, args) = parser.parse_args()

#################################################
##		Begin
##


if options.lucymode:
	fp = FastaParser(sys.argv[-2],quals=False)
	qp = FastaParser(sys.argv[-1],quals=True)
	
	outputfasta = open("%s.trim.lucy.fasta" % ( sys.argv[-3]), "w")
	outputqual =  open("%s.trim.lucy.qual" % ( sys.argv[-3]), "w")
	
	for head, seq in fp:
		qhead, qseq = qp.next()
		name, minclone, medclone, maxclone, left, right = head.strip().split()
		if name == qhead:
		
			
			left = int(left)
			right= int(right)
#			print "orig"
#			print seq
#			print qseq
			seq = seq[left:right]
			qseq = " ".join(qseq.split()[left:right])
#			print "----"
#			print seq
#			print qseq
#			print " "
			
			outputfasta.write(">%s\n%s\n" % (name, seq))
			outputqual.write(">%s\n%s\n" % (name, qseq))
			
		else:
			print name. qhead
			
	outputfasta.close()
	outputqual.close()		

else:
	qual = sys.argv[-1].strip().split(".")[-1]=="qual"
	#print qual
	fp = FastaParser(sys.argv[-1], quals=qual)
	for head, seq in fp:
		
		### deal with qual files
		if qual:
			seq=seq.split()
			sep=" "
		else:
			sep=""	
			
		### deal with length	
		if len(seq)>options.end:
			if head.find("length=")>-1:
				head="%s-%s=%s" % (head, len(seq) - options.end, options.end)
			seq = seq[0:options.end]
			
		#if head.startswith("F0K1MYC02GCRVD"):	
		print ">%s\n%s" % (head, sep.join(seq))



#################################################
##		Finish
#################################################
