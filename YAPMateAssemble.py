########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2013 Sebastian Szpakowski
########################################################################################

#################################################
## 	Pipeline for assemblying Mates
#################################################
import sys
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

parser.add_option("-P", "--PROJECT", dest="project", default="",
                 help="project code", metavar="#")

parser.add_option("-E", "--EMAIL", dest="email", default="",
                 help="e-mail address", metavar="@") 

parser.add_option("-I", "--imported", dest="importflag", default=False, action="store_true",
                 help="Unless specified, this script will only import the data. Useful for running the script first on the archive1 server.")

parser.add_option("-H", "--head", dest="head", default=0, type="int",
                 help="For dry runs, import only # of lines from the input files")


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


info = InfoParser(options.fn_info)

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
	
	M1, M2 = info.getFiles(sample)
	
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
	 	
				
#		ARGS = {
#				" "	:  ref,
#				"-v"	: "3",
#				"-k" : "1",
#				"-5" : "1",
#				"-3" : "1",
#				"--phred%s-quals" % Q : "",
#				"--mm" : "",
#				"--suppress": "2,3,4,5,6,7,8",
#				"-p": "4" 
#		}						
#		b1 = Bowtie1(dict(), ARGS, [chunks])	#		


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
		
		
#		
##		ARGS = {
##				"-Q": Q ,
##				"-q": "25",
##				"-p": "95"
##		}		
##		S4 = fastq_quality_filter(dict(), ARGS, [S3])
#		
		ARGS = { 
				"-Q": Q
		}
		S5 = fastq2fasta(dict(), ARGS, [S4a])
			
	 	S6 = mateInterweave(dict(), dict(), [S5])	
		
	 	S7 = FileMerger("fasta", [S6], prefix=sample)
	 	CLEAN.append(S7)
	 	
		INS = {"properties": ["/usr/local/packages/clc-ngs-cell/license.properties"]}
		ARGS = 	{
					"-p" : "no" ,
					"--cpus": "4",
					"-v" : "", 
					"--no-progress": ""				
		}
		S8 = CLC_Assemble(INS, ARGS, [S7])
	
		args = {
					"types": "fasta", 
					"old": "contig", 
					"new": sample
				}
		
		S9 = SED_replace({}, args, [S8])	
	
		ARGS = 	{
					"-p" : "no" ,
					#### changed from 10 to 4
					"--cpus": "4", 
					"-m" : "0.5",
					"--no-progress": ""				
		}
		
		S10 = CLC_Assemble_Ref(INS, ARGS, [S7,S9])
		S11 = CLC_Assemble_Info({}, [S10])
		S11a =ContigCoverageUpdate({}, [S11])
		
		S12 = CLC_Assemble_Table({}, [S7, S10])
		S13 = SingletonsFishOut({},{}, [S12])
				
	 	ASSEMBLIES.append(S11a)
	 	COVERAGES.append(S11)
	 	READS.append(S13)

if options.importflag or options.head>0:	
	PLOTS1 = FastaSummaryRPlots({},ASSEMBLIES)
	PLOTS2 = R_OTUplots( {}, {}, COVERAGES)
	
	OutputStep("QC", "pdf,png,matrix", QC)
	OutputStep("CLEAN", "fasta", CLEAN)
	OutputStep("ASSEMBLIES", "pdf,png,fasta", PLOTS1)
	OutputStep("COVERAGES", "png,pdf,clcassemblystats", PLOTS2)
	OutputStep("SINGLETONS", "fasta", READS)



#################################################
##		Finish
#################################################?
