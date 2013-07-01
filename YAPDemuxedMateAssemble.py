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

parser.add_option("-r", "--reference", dest="reference", default="human",
                 help="reference genome of the host (to filter out contaminants)", metavar="human")

parser.add_option("-P", "--PROJECT", dest="project", default="",
                 help="project code", metavar="#")

parser.add_option("-E", "--EMAIL", dest="email", default="",
                 help="e-mail address", metavar="@") 

parser.add_option("-I", "--imported", dest="importflag", default=False, action="store_true",
                 help="Unless specified, this script will only import the data. Useful for running the script first on the archive1 server.")

parser.add_option("-H", "--head", dest="head", default=0, type="int",
                 help="For dry runs, import only # of lines from the input files")

parser.add_option("-S", "--sequences", 
					dest="minSeqs",type="int", default=1000,
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


info = InfoParser(options.dir_info)

QC=list()
CLEAN= list()
ASSEMBLIES = list()
COVERAGES = list()
READS = list()

init(options.project, options.email, maxnodes = options.maxNodes)

if options.reference =="mouse":
	ref = "/usr/local/projects/KFLAB/sszpakow/annotation/GENOMES/MM10/mm10"
else:
	ref = "/usr/local/projects/KFLAB/sszpakow/annotation/GENOMES/HG19/hg19"
 
for sample in info.getSamples():
	#print sample
 	M1, M2 = info.getFiles(sample)
 	#print M1, M2
 	if countLines(M1)/4 >= options.minSeqs and countLines(M2)/4 >= options.minSeqs:
	 	INS = {"mate1": ["%s~%s" % (M1, sample)], "mate2": ["%s~%s" % (M2, sample)]}
	 	if options.head == 0:
	 		x = FileImport(INS)
	 	else:
	 		x = FileMiniImport(INS, {"lines": options.head})	
	 	
		if options.importflag or options.head>0:
			Q = getQ(M1)
	 		
			ARGS = {
				 "-h": "20",
				 "-m": "",
				 "-v": ""
				}
			sqa1 = SQA(ARGS, [x]) 
			QC.append(sqa1)
	 		
			# -> 200 000 sequences per file (assuming fastq)
			ARGS = {
						"types": "mate1,mate2",
						"chunk":  "800000"
			}
			chunks = FileSplit(ARGS, [x])
	 	
		 	stats1 = fastx_quality_stats(dict(), {"-Q": "%s" %(Q) } , [x])			
	 	 	
			ARGS = {
					"-x"	:  ref,
					"--nondeterministic": "",
					"--phred%s" % (Q): "" ,
					"--no-unal": "",
					"--no-hd": "",
					"--no-head": "",
					"--no-sq": "",
					"--very-fast" : "",
					"--no-discordant": "",
					"-p": "4" 
			}						
			S1 = Bowtie2(dict(), ARGS, [chunks])
	 			
			ARGS = {
						"type": "sam",
						"newtype": "filter",
						"awk" : " NF>2 {print $1}",
						"uniq" : "",
						"sort" : "",
						"postprocess" : " tr -d \"@\" | sort | uniq "
					}
			S2 = AwkCommand({}, ARGS, [S1] )	 
			S3 = ContaminantRemoval(dict(), dict(), [S2])	
	 			 		
			ARGS = {
			 "-h": "20"
			}
			S4 = SQAtrim(ARGS, [S3])
	 			 		
			ARGS = {
			 "-l": "30"
			}
			S4a = SQAlenfil(ARGS, [S4])
	 		
			ARGS = { 
					"-Q": Q
			}
			S5 = fastq2fasta(dict(), ARGS, [S4a])
	 			
	 		if options.matesindependent:
	 			S6 = MateMerge(dict(), dict(), [S5], prefix=sample)
	 		else:	
		 		S6 = mateInterweave(dict(), dict(), [S5])	
	 		
		 	S7 = FileMerger("fasta", [S6], prefix=sample)
		 	CLEAN.append(S7)
 	else:
		 print "skipping:\n\t%s\n\t%s" % (M1, M2)

if options.importflag or options.head>0:
	
	############################################################################# 	
	### assemble using all samples
	### 	
	
	SALL = FileMerger("fasta", CLEAN, prefix = "ALL")
	
	INS = {"properties": ["/usr/local/packages/clc-ngs-cell/license.properties"]}
	
	if options.matesindependent:
			ARGS = 	{
			"-p": "no",
			"--cpus": "4",
			"-v" : "", 
			"--no-progress": ""				
		}
	else:
		ARGS = 	{
				"-p": "fb ss 180 500",
				"--cpus": "4",
				"-v" : "", 
				"--no-progress": ""				
		}
	S8 = CLC_Assemble(INS, ARGS, [SALL])
	
	args = {
			"types": "fasta", 
			"old": "contig", 
			"new": "ALL"
	}
 		
 	S9 = SED_replace({}, args, [S8])	
		
	OutputStep("CLEAN", "fasta", CLEAN) 	 	
	OutputStep("ASSEMBLED", "fasta", [S9])


 
#################################################
##		Finish
#################################################?
