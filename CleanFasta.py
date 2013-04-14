########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2013 Sebastian Szpakowski
########################################################################################

#################################################
## 	run mothur alignments
#################################################
import sys, time, random, shlex, os, itertools, gc
from optparse import OptionParser
from Bio import SeqIO, Seq
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from threading import *
from subprocess import *

_author="Sebastian Szpakowski"
_date="2011/06/01"
_version="Version 1"

#################################################
##		Classes
##
class	CommandWorker (Thread)  :
	def __init__(self, command):
		Thread.__init__(self)
		self.command = command
		
	def	run(self):	
		pool_sema.acquire()
		gc.collect()
		p = Popen(shlex.split(self.command), stdout=PIPE)
		gc.collect()
		pool_sema.release()	

	def	__str__(self):
		otpt =  "%s\n%s" % (self.command, len(result))
		return (otpt)
	
	
	#################################################
	### Iterator over input fasta file.
	### Only reading when requested
	### Useful for very large FASTA files
	### with many sequences
	
class FastaParser:
	def __init__ (self, x, quals=False):
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
	def __iter__(self):
		return(self)    
				
		##### 
	def next(self):
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
									
	def addSequence(self, x):
				self.currentFastaSequence = "%s%s%s" % (self.currentFastaSequence,self.linesep, x.strip())            
	def __str__():
		return ("reading file: %s" %self.filename) 
				
	
		
#################################################
##		Functions
##
def	cleanup(recordlist):
	otpt=list()
	for rec in recordlist:
		seq = str(rec.seq).replace(".", "").replace("-", "").replace(" ", "").replace("\n", "")
		id = rec.id.strip().split()[0]
		
		if len(id.split("|"))==5:
			id = id.split("|")[1]
		tmp = SeqRecord(seq=Seq(seq), id=id, name=id, description=id)
		otpt.append(tmp)
	return (otpt)

#################################################
##		Arguments
##
parser = OptionParser()

#parser.add_option("-f", "--file", dest="filename",
#                  help="write report to FILE", metavar="FILE")
#parser.add_option("-q", "--quiet",
#                  action="store_false", dest="verbose", default=True,
#                  help="don't print status messages to stdout")

parser.add_option("-i", "--input", dest="fn_input",
                  help="load and clean fasta FILE", metavar="FILE")

parser.add_option("-o", "--output", dest="fn_output",
                  help="load and clean fasta FILE", metavar="FILE")
		  
(options, args) = parser.parse_args()

#################################################
##		Begin
##


if options.fn_input != None :
	print "cleaning up..."
	#C = list(SeqIO.parse(open(options.fn_input), "fasta"))
	#cleana = cleanup(C)
	output_handle = open(options.fn_output, "w")
	#SeqIO.write(cleana, output_handle, "fasta")
	
	for head, seq in FastaParser(options.fn_input):	
		
		seq = seq.replace(".", "").replace("-", "").replace(" ", "").replace("\n", "")
		head = head.strip().split()[0]
		
		output_handle.write(">%s\n%s\n" % (head, seq))
	output_handle.close()
#################################################
##		Finish
#################################################
