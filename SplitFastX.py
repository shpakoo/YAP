#!/usr/bin/python
#################################################
## 	A new program
#################################################
import sys
from optparse import OptionParser
from Bio import SeqIO

_author="Sebastian Szpakowski"
_date="2012/03/29"
_version="Version 1"

#################################################
##		Classes
##

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

parser.add_option("-c", "--chunk", dest="chunk_size", type=int,
                  help="put N sequences per output file", metavar="N")

parser.add_option("-f", "--format", dest="file_format", default="fasta",
                  help="list of names to remove", metavar="FILE")

(options, args) = parser.parse_args()

#################################################
##		Begin
##

counter_F = 0
outputs= list()
file_in = open(options.fn_input, "r")

for record in SeqIO.parse(file_in, options.file_format) :	
	outputs.append(record)
	if len(outputs) >= options.chunk_size:
		newfilename = "%s.%s.chunk.%s" % (".".join(options.fn_input.strip().split(".")[:-1]), counter_F, options.file_format)
		file_out = open(newfilename, "w")
		SeqIO.write(outputs, file_out, options.file_format)	
		file_out.close()
		outputs = list()
		counter_F += 1

if len(outputs)>0:
	newfilename = "%s.%s.chunk.%s" % (".".join(options.fn_input.strip().split(".")[:-1]), counter_F, options.file_format)
	file_out = open(newfilename, "w")
	SeqIO.write(outputs, file_out, options.file_format)	
	file_out.close()
		

file_in.close()
	

#################################################
##		Finish
#################################################
