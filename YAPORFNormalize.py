########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2013 Sebastian Szpakowski
########################################################################################

#################################################
## 	Pipeline for normalization of ORFs
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

parser.add_option("-S", "--dir", dest="dn_STRUCTURE", default="",
                  help="directory that will be scanned for contig,singletons,ORF and PROK files", metavar="path/to/dir")
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

if  options.dn_STRUCTURE== "" :
	parser.print_help()
	sys.exit(1)


assemblies = dict()
### ASSEMBLIES 
for f in glob.glob("%s/*ASSEMBLIES*/*.fasta" % (options.dn_STRUCTURE)):
	id = f.strip().split("/")[-1].split(".")[1].split("-")[0]
	assemblies[id] = f

reads_A = dict()
reads_S = dict()
### READS 
for f in glob.glob("%s/*SINGLETONS*/*.fasta" % (options.dn_STRUCTURE)):
	id = f.strip().split("/")[-1].split(".")[1].split("-")[0]
	if f.find("assembled")>-1:
		reads_A[id]=f
	else:
		reads_S[id]=f	

orfs = dict()
### ORFs
for f in glob.glob("%s/*ASSEMBLIES*/ORFs/*/*/frag-gene-scan/*.faa" % (options.dn_STRUCTURE)):
	id = f.strip().split("/")[-1].split(".")[1].split("-")[0]
	orfs[id]=f
	
proks = dict()
### PROKs	
for f in glob.glob("%s/*ASSEMBLIES*/PROKs/*.faa" % (options.dn_STRUCTURE)):
	
	id = f.strip().split("/")[-1].split(".")[0].split("-")[0]
	proks[id]="%s/prok-annotation/annotation_rules.combined.out" % f
						
init(options.project, options.email)


all = list()
proknorm = list()
for k in orfs.keys():
	INS= {"fasta": ["%s~%s.contigs" % (assemblies[k], k)]}
	Ia = FileImport(INS)
	
	INS= {"fasta": ["%s~%s.assembled" % (reads_A[k], k)]}
	Ib = FileImport(INS)
	
	INS= {"fasta": ["%s~%s.orfs" % (orfs[k], k)]}
	Ic = FileImport(INS)
	
	INS= {"txt": ["%s~%s.prok" % (proks[k], k)]}
	Id = FileImport(INS)
		 	
		 	
	INS = {"properties": ["/usr/local/packages/clc-ngs-cell/license.properties"]}
	ARGS = 	{
				"-p" : "no" ,
				"--cpus": "10", 
				"-m" : "0.5",
				"--no-progress": ""				
	}
	
	S1 = CLC_Assemble_Ref(INS, ARGS, [Ia,Ib])
	S2 = CLC_Assemble_Table({}, [Ib, S1])
	S3 = ORFCoverage ({},{},[Ic,S2])
	S4 = ORFCoverageNorm({},{},[S3])
	all.append(S4)

	S5 = PROKModify({},{},[Id, S4])
	proknorm.append(S5)
	
P1 = ORFCoverageNorm({}, {}, all)
	
OutputStep("WEIGHTS", "normweight", all) 	
OutputStep("PLOTS_INDIVIDUAL", "png", all) 
OutputStep("PLOTS_GLOBAL", "png", P1) 
OutputStep("NORMALIZED", "jcvinorm", proknorm) 

##################################################
###		Finish
##################################################?
