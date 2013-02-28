########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2013 Sebastian Szpakowski
########################################################################################

#################################################
## 	filter mates to keep pairs
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


parser.add_option("-f", "--filter", dest="fn_filter",default = "",
                  help="list of names to remove", metavar="FILE")

parser.add_option("-k", "--keep", dest="fn_keep",default = "",
                  help="list of names to keep", metavar="FILE")
           
parser.add_option("-t", "--format", dest="file_format", default="fastq",
                  help="format of the input file", metavar="FILE")

parser.add_option("-s", "--suffix", dest="suffix", default="filt",
                  help="suffix to add right before output file extension", metavar="FILE")

(options, args) = parser.parse_args()

#################################################
##		Begin
##

if options.fn_filter != "":
	keep = False
	names = [x.strip().split()[0].strip("@") for x in loadLines(options.fn_filter)]
	names = set(names)
	
elif options.fn_keep != "":
	keep = True
	names = [x.strip().split()[0].strip("@") for x in loadLines(options.fn_keep)]
	names = set(names)	

else:
	print "keep or filter required"
	parser.print_help()
	sys.exit(1)

for file in options.fn_input.strip().split(","):
	
	
	filename = file.strip().split("/")[-1]
	ext = filename.split(".")[-1]
	filename = "%s.%s.%s" % (".".join(filename.split(".")[:-1]), options.suffix, ext)
	
	file_in = open(file, "r")
	file_out = open(filename, "w")
	
	for record in SeqIO.parse(file_in, options.file_format) :
		recid = record.id.split()[0].strip("@")
		found = recid in names
		if ( found and keep ) or ( not found and not keep ) :
			#outputs.append(record)	
			SeqIO.write([record], file_out, options.file_format)	
			names.discard(recid)		
		else:
			print "filtering out:", record.id
								
	file_in.close()	
	file_out.close()
#################################################
##		Finish
#################################################
