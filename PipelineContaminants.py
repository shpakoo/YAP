########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2013 Sebastian Szpakowski
########################################################################################

#################################################
## 	Remove Contaminants
#################################################
import sys, os
from optparse import OptionParser, OptionGroup
from StepsLibrary import *
from StepsLibrary_EXP import *
from collections import defaultdict

_author="Sebastian Szpakowski"
_date="2011/01/01"
_version="Version 1"

#################################################
##		Classes
##

class	InfoParser:
	def __init__(self, filename):
		self.file = filename
		self.sampleinfo = defaultdict(list)
		
		id = ""
		for line in GeneralPurposeParser(self.file, sep=" "):
			if line[0] == "Lane":
				id = line[2]
			else:
				self.sampleinfo[id].append(line[0])	
				
	def	getSamples(self):
		return self.sampleinfo.keys()
		
	def	getFiles(self, id):
		return self.sampleinfo[id][0:2]
			

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

group = OptionGroup(parser, "Required", description="Will not run without these !")

group.add_option("-P", "--PROJECT", dest="project", default="",
                 help="project code", metavar="#")
group.add_option("-E", "--EMAIL", dest="email", default="",
                 help="e-mail address", metavar="@")                 
group.add_option("-r", "--reference", dest="reference", default="human",
                 help="remove reads that match X (human, hg19 or mouse, mm9)", metavar="human")
group.add_option("-i", "--inputs", dest="fn_inputs", default="",
                 help="A file with ABSOLUTE paths directories or individual files one per line. Each directory will be searched for fasta files", metavar="inputs.txt")

parser.add_option_group(group)

group = OptionGroup(parser, "Optional Configuration", description="parameters to alter if necessary")
parser.add_option_group(group)

group = OptionGroup(parser, "Technical", description="could be useful sometimes")
group.add_option("-C", "--NODESIZE", dest="nodesize", default=30,
                 help="maximum number of grid node's CPUs to use\n[%default]", metavar="#")
parser.add_option_group(group)
	
(options, args) = parser.parse_args()


#################################################
##		Begin
##

if options.fn_inputs =="" or options.email == "" or options.project =="":
	parser.print_help()
	sys.exit(1)
		
init(options.project, options.email)

filters = list()
processed = list()

for path in loadLines(options.fn_inputs):
	path = path.strip()
	files = list()
	if os.path.isdir(path):
		files = glob.glob("%s/*.seq" % (path))
		files.extend(glob.glob("%s/*.fa" % (path)))
		files.extend(glob.glob("%s/*.fasta" % (path)))
		files.extend(glob.glob("%s/*.fna" % (path)))
	else:
		files.append( path )
		
	if options.reference =="mouse":
		ref= "/usr/local/projects/KFLAB/sszpakow/annotation/GENOMES/MM9/mm9"	
	else:
		ref = "/usr/local/projects/KFLAB/sszpakow/annotation/GENOMES/HG19/hg19"
	
	for file in files:
		
		INS = {"fasta": ["%s" % (file)]}
		input = FileImport(INS)
		
		### remove spaces from header to make contig names unique and make all 50 base fragments (overlapping by 25)
		b0 = TilingFasta(dict(), [input])
		ARGS = {
			" "	:  ref,
			"-v"	: "3",
			"-k" : "1",
			"-5" : "1",
			"-3" : "1",
			"--mm" : "",
			"-f" : "",
			#"--suppress": "2,3,4,5,6,7,8",
			"-p": "%s" % options.nodesize
		}				
		
	 	b1 = Bowtie1(options.nodesize, dict(), ARGS, [b0])
	 	filters.append(b1)
	 	
	 	b2 = ContaminantRemoval(INS, dict(), [b1])	
	 	
	 	processed.append(b2)	 	
	 	
OutputStep("FILTERS", "bowtie1alignment", filters)		
OutputStep("FASTAS", "fasta", processed)	

#################################################
##		Finish
#################################################
