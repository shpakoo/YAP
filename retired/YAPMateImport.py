#!/usr/bin/python
#################################################
## 	A new program
#################################################
import sys
from optparse import *
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

parser.add_option("-i", "--info", dest="fn_info",
                 help="a list of lanes from JIRA", metavar="allinfo.csv")

parser.add_option("-r", "--reference", dest="reference", default="human",
                 help="reference genome of the host (to filter out contaminants)", metavar="human")

parser.add_option("-P", "--projectid", dest="projectid", default="810013",
                 help="projectID", metavar="#")

parser.add_option("-I", "--imported", dest="importflag", default=False, action="store_true",
                 help="By default this should be run on the archive 1 to import the data.")



group = OptionGroup(parser, "Technical", "\"behind-the-scene\" options")

group.add_option("-T", "--threads", 
					dest="maxThreads",type="int", default=4,
                  	help="maximum number of threads to spawn. [default: %default]")
                  	
group.add_option("-N", "--nodes", 
					dest="maxNodes",type="int", default=125,
                  	help="maximum number of grid-nodes to use. [default: %default]")   
     
group.add_option("-F", "--filesOpen", 
					dest="maxFiles",type="int", default=50,
                  	help="maximum number of files opened at a given time. [default: %default]")       

group.add_option("-C", "--NODESIZE", dest="nodesize", default=30,
                 help="maximum number of grid node's CPUs to use", metavar="#")

                  	
parser.add_option_group(group)

(options, args) = parser.parse_args()

#################################################
##		Begin
##

init(options.projectid, "sszpakow@jcvi.org")
info = InfoParser(options.fn_info)

ASSEMBLIES = list()
ORFS = list()
ORFSfilt = list()

if options.reference =="mouse":
	ref= "/usr/local/projects/KFLAB/sszpakow/annotation/GENOMES/MM9/mm9"
else:
	ref = "/usr/local/projects/KFLAB/sszpakow/annotation/GENOMES/HG19/hg19"

for sample in info.getSamples():
	
	M1, M2 = info.getFiles(sample)
	samplename = M1.strip().split("/")[-1].strip().split("_")[0]
	INS = {"mate1": ["%s~%s" % (M1, sample)], "mate2": ["%s~%s" % (M2, sample)]}
	x = FileImport(INS)
	



#################################################
##		Finish
#################################################
