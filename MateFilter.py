#!/usr/bin/python
#################################################
## 	A new program
#################################################
import sys
from optparse import OptionParser
from Bio import SeqIO

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

parser.add_option("-i", "--input", dest="fn_input",
                  help="fastQ file name (or names, comma separated)", metavar="FILE")


parser.add_option("-f", "--filter", dest="fn_filter",
                  help="list of names to remove", metavar="FILE")
                  

parser.add_option("-t", "--format", dest="file_format", default="fastq",
                  help="list of names to remove", metavar="FILE")

(options, args) = parser.parse_args()

#################################################
##		Begin
##

names = [x.strip().split()[0].split(":")[0] for x in loadLines(options.fn_filter)]
print names
for file in options.fn_input.strip().split(","):
	file_in = open(file, "r")
	
	filename = file.strip().split("/")[-1]
	ext = filename.split(".")[-1]
	filename = "%s.filt.%s" % (".".join(filename.split(".")[:-1]), ext)
	
	outputs = list()
	
	for record in SeqIO.parse(file_in, options.file_format) :
		
		if record.id.split()[0] in names:
			print "filtering out:", record.id
		else:
			outputs.append(record)
					
	file_in.close()
	file_out = open(filename, "w")
	SeqIO.write(outputs, file_out, options.file_format)	
	file_out.close()
#################################################
##		Finish
#################################################
