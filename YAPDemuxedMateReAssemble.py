########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2013 Sebastian Szpakowski
########################################################################################

#################################################
## 	Pipeline for assemblying Mates
#################################################
import sys, subprocess
from optparse import *
from StepsLibrary import *
from StepsLibrary_EXP import *
from collections import defaultdict

_author="Sebastian Szpakowski"
_date="2012/12/17"
_version="Version 2"

#################################################
##		Classes
##

class	InfoParser:
	def __init__(self, dirname):
		
		self.sampleinfo=defaultdict(list)
		
		for mate in (1,2):
			for file in glob.glob("%s/*.mate%s" % (dirname, mate)):
				filename = file.strip().split("/")[-1]
				id = filename.split(".")[-2]
				self.sampleinfo[id].append("%s/%s" % (dirname, filename))
			
				
	def	getSamples(self):
		return self.sampleinfo.keys()
		
	def	getFiles(self, id):
		return self.sampleinfo[id][0:2]
			
		
#################################################
##		Functions
##

def countLines(fname):
	#http://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
    p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE, 
                                              stderr=subprocess.PIPE)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])
	



#################################################
##		Arguments
##
parser = OptionParser()

parser.add_option("-i", "--infodir", dest="dir_info",
                 help="a directory with demultiplexed mate1 and mate2", metavar="dir")

parser.add_option("-r", "--reffasta", dest="fn_ref",
                 help="a file with contigs", metavar="refdir")

parser.add_option("-P", "--PROJECT", dest="project", default="",
                 help="project code", metavar="#")

parser.add_option("-E", "--EMAIL", dest="email", default="",
                 help="e-mail address", metavar="@") 

parser.add_option("-I", "--imported", dest="importflag", default=False, action="store_true",
                 help="Unless specified, this script will only import the data. Useful for running the script first on the archive1 server.")

parser.add_option("-H", "--head", dest="head", default=0, type="int",
                 help="For dry runs, import only # of lines from the input files")

parser.add_option("-S", "--sequences", 
					dest="minSeqs",type="int", default=10,
                  	help="Make sure the input files have at least n sequences [default: %default]")
parser.add_option("-x", "--independent", dest="matesindependent", action = "store_true", default=False,
                 help="""If Specified, mates will be treated independently, otherwise complete pairs will be used. [%default]""", metavar="#") 


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


QC=list()
CLEAN= list()
ASSEMBLIES = list()
COVERAGES = list()
READS = list()

init(options.project, options.email, maxnodes = options.maxNodes)

INS = {"fasta": [options.fn_ref] }
ref = FileImport(INS)
 
for fasta in glob.glob("%s/*.fasta" % (options.dir_info)):
	if countLines(fasta)/2 >  options.minSeqs:
		sample = ".".join(fasta.strip().split("/")[-1].split(".")[1:-1])
	 	INS = {"fasta": ["%s~%s" % (fasta, sample)]}
	 	
	 	#print M1, M2
	  	S1 = FileImport(INS)
		
		
		if options.matesindependent:
			ARGS = 	{
			"-p": "no",
			"--cpus": "10", 
			"-m" : "0.5",
			"--no-progress": ""				
		}
		else:
			ARGS = 	{
			"-p": "fb ss 180 500",
			"--cpus": "10", 
			"-m" : "0.5",
			"--no-progress": ""			
			}
		
		INS = {"properties": ["/usr/local/packages/clc-ngs-cell/license.properties"]}
		S10 = CLC_Assemble_Ref(INS, ARGS, [S1, ref] )
		S11 = CLC_Assemble_Info({}, [S10])	
		S12 = CLC_Assemble_Table({}, [ref, S10])
		S13 = SingletonsFishOut({},{}, [S12])
			
		COVERAGES.append(S11)
		READS.append(S13)
	
  
PLOTS2 = R_OTUplots( {}, {}, COVERAGES)
 	
OutputStep("COVERAGES", "png,pdf,clcassemblystats", PLOTS2)
OutputStep("SINGLETONS", "fasta", READS)
 
#################################################
##		Finish
#################################################?
