########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2013 Sebastian Szpakowski
########################################################################################

#################################################
## Pipeline for CLC asssembly of fasta files
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
		
#################################################
##		Functions
##

#################################################
##		Arguments
##
parser = OptionParser()

parser.add_option("-d", "--directory", dest="dn_fasta", default="",
                 help="directory that will be scanned to fasta files", metavar="path/to/dir")

parser.add_option("-f", "--fasta", dest="fn_fasta", default="",
                 help="fasta file", metavar="path/to/file")

parser.add_option("-r", "--reference", dest="reference", default="human",
                 help="reference genome of the host (to filter out contaminants)", metavar="human")

parser.add_option("-P", "--PROJECT", dest="project", default="",
                 help="project code", metavar="#")
parser.add_option("-E", "--EMAIL", dest="email", default="",
                 help="e-mail address", metavar="@") 



group = OptionGroup(parser, "Technical", "\"behind-the-scene\" options")

group.add_option("-T", "--threads", 
					dest="maxThreads",type="int", default=4,
                  	help="maximum number of threads to spawn. [default: %default]")
                  	
group.add_option("-N", "--nodes", 
					dest="maxNodes",type="int", default=125,
                  	help="maximum number of grid-nodes to use. [default: %default]")   
     
group.add_option("-F", "--filesOpen", 
					dest="maxFiles",type="int", default=5,
                  	help="maximum number of files opened at a given time. [default: %default]")       

group.add_option("-C", "--NODESIZE", dest="nodesize", default=30,
                 help="maximum number of grid node's CPUs to use", metavar="#")

                  	
parser.add_option_group(group)

(options, args) = parser.parse_args()

#################################################
##		Begin
##

if options.fn_fasta=="" and options.dn_fasta== "" or options.project=="":
	parser.print_help()
	sys.exit(1)

CLEAN= list()
ASSEMBLIES = list()
COVERAGES = list()
READS = list()


if options.reference =="mouse":
	ref= "/usr/local/projects/KFLAB/sszpakow/annotation/GENOMES/MM9/mm9"
else:
	ref = "/usr/local/projects/KFLAB/sszpakow/annotation/GENOMES/HG19/hg19"

files = list()
if options.fn_fasta != "":
	files.append(options.fn_fasta)

if options.dn_fasta != "":
	for ext in ("fa", "fasta", "fna"):
		tmp = glob.glob("%s/*.%s" % (options.dn_fasta, ext) )	
		if len(tmp)>0:
			files.extend(tmp)
			
init(options.project, options.email)

for file in files:
	#sample =  ".".join(file.strip().split(".")[1:-1])
	tmp = file.strip().split("/")[-1].split(".")
	if len(tmp)>=3:
		sample = tmp[1]
	else:
		sample = tmp[0]	
	
	INS= {"fasta": ["%s~%s" % (file, sample)]}
	
	S1 = FileImport(INS)
	S2 = FastaSort({"types": "fasta"}, [S1])
	S3 = FastaSplit({"types": "fasta", "chunk": "100000"}, [S2])
				
#	ARGS = {
#			" "	:  ref,
#			"-v"	: "3",
#			"-k" : "1",
#			"-5" : "1",
#			"-3" : "1",
#			"--mm" : "",
#			"--suppress": "2,3,4,5,6,7,8",
#			"-p": "4" 
#	}						
#	S4 = Bowtie1(dict(), ARGS, [S3])	

	
	ARGS = {
			"-x"	:  ref,
			"-k" : "1",
			"--nondeterministic": "",
			"--local" : "",
			"--no-unal": "",
			"--no-hd": "",
			"--no-sq": "",
			"-N": "1",
			"-p": "4" 
	}						
	S4a = Bowtie2(dict(), ARGS, [S3])	
	
	ARGS = {
					"type": "sam",
					"newtype": "filter",
					"awk" : " NF>2 {print $1}",
					"uniq" : "",
					"sort" : "",
					"postprocess" : " tr -d \"@\" | sort | uniq "
				}
	
	S4 = AwkCommand({}, ARGS, [S4a] )

	S5 = ContaminantRemoval(dict(), {"-t": "fasta" }, [S4])	
	S6 = FileMerger("fasta", [S5], prefix=sample)
		
	CLEAN.append(S6)
	 	
	INS = {"properties": ["/usr/local/packages/clc-ngs-cell/license.properties"]}
	ARGS = 	{
				"-p" : "no" ,
				"--cpus": "4",
				"-v" : "", 
				"--no-progress": ""				
	}
	S7 = CLC_Assemble(INS, ARGS, [S6])

	args = {
				"types": "fasta", 
				"old": "contig", 
				"new": sample
			}
	
	S8 = SED_replace({}, args, [S7])	

	ARGS = 	{
				"-p" : "no" ,
				"--cpus": "10", 
				"-m" : "0.5",
				"--no-progress": ""				
	}
	
	S9 = CLC_Assemble_Ref(INS, ARGS, [S6,S8])
	S10 = CLC_Assemble_Info({}, [S9])
	S11 = CLC_Assemble_Table({}, [S6, S9])
	S12 = SingletonsFishOut({},{}, [S11])
	
	
 	ASSEMBLIES.append(S8)
 	COVERAGES.append(S10)
 	READS.append(S12)
 	
PLOTS1 = FastaSummaryRPlots({},ASSEMBLIES)
PLOTS2 = R_OTUplots( {}, {}, COVERAGES)

OutputStep("CLEAN", "fasta", CLEAN)
OutputStep("ASSEMBLIES", "pdf,png,fasta", PLOTS1)
OutputStep("COVERAGES", "png,pdf,clcassemblystats", PLOTS2)
OutputStep("SINGLETONS", "fasta", READS)

##################################################
###		Finish
##################################################?
