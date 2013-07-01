########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2013 Sebastian Szpakowski
########################################################################################

#################################################
##    A pipeline for cleaning up sff files
#################################################

import sys
from optparse import OptionParser, OptionGroup
from StepsLibrary import *
from StepsLibrary_EXP import *
from collections import defaultdict
from Queue import Queue

_author="Sebastian Szpakowski"
_date="2012/08/27"
_version="Version 1"

#################################################
##        Classes
##
    #################################################
    ### Iterator over input file.
    ### every line is converted into a dictionary with variables referred to by their 
    ### header name
class GeneralPurposeParser:
    def    __init__(self, file, skip=0, sep="\t"):
        self.filename = file
        self.fp = open(self.filename, "rU")    
        self.sep = sep
        
        self.linecounter = 0
        self.currline=""
        
    def __iter__(self):
        return (self)
    
    def next(self):
        otpt = dict()
        for currline in self.fp:
            currline = currline.strip().split(self.sep)
            self.currline = currline
            self.linecounter = self.linecounter + 1
            return(currline)            
        raise StopIteration
                    
    def    __str__(self):
        return "%s [%s]\n\t%s" % (self.filename, self.linecounter, self.currline)


            
#################################################
##        Functions
##
    ################################################
    ### Read in a file and return a list of lines
    ###
def loadLines(x):
    try:
        fp = open(x, "rU")
        cont=fp.readlines()
        fp.close()
        #print "%s line(s) loaded."  % (len(cont))
    except:
        cont=""
        #print "%s cannot be opened, does it exist? " % ( x )    
    return cont

def getOligoFlip(filename):
    if filename.find(".R.")>-1:
        return "T"
    else:
        return "F"
    




        
#################################################
##        Arguments
##

parser = OptionParser()

group = OptionGroup(parser, "Required", description="Will not run without these !")

group.add_option("-P", "--PROJECT", dest="project", default="",
                 help="project code", metavar="#")
group.add_option("-E", "--EMAIL", dest="email", default="",
                 help="e-mail address", metavar="@")  

group.add_option("-b", "--database", dest="db",
                 help="absolute path to a database", metavar="blast")

parser.add_option("-d", "--directory", dest="dn_fasta", default="",
                 help="directory that will be scanned to fasta files", metavar="samples.group")

parser.add_option("-f", "--fasta", dest="fn_fasta", default="",
                 help="fasta file", metavar="reads.fasta")   

            
parser.add_option_group(group)

group = OptionGroup(parser, "Optional Configuration", description="parameters to alter if necessary")

group.add_option("-m", "--mode", dest="blast_mode", default="blastn",
                 help="which blast application to use#\n[%default]", metavar="?blast?")

group.add_option("-s", "--split", dest="split_size", default=5000, type="int",
                 help="split input fasta into chunks with n sequences#\n[%default]", metavar="#")
            
parser.add_option_group(group)

group = OptionGroup(parser, "Technical", description="could be useful sometimes")
group.add_option("-T", "--threads", 
                    dest="maxThreads",type="int", default=4,
                      help="maximum number of threads to spawn. [default: %default]")
                      
group.add_option("-N", "--nodes", 
                    dest="maxNodes",type="int", default=200,
                      help="maximum number of grid-nodes to use. [default: %default]")   
     
group.add_option("-F", "--filesOpen", 
                    dest="maxFiles",type="int", default=50,
                      help="maximum number of files opened at a given time. [default: %default]")       

group.add_option("-C", "--NODESIZE", dest="nodesize", default=30,
                 help="maximum number of grid node's CPUs to use", metavar="#")

parser.add_option_group(group)
    
(options, args) = parser.parse_args()

#################################################
##        Begin
##

if (options.fn_fasta == "" and options.dn_fasta=="") or options.email == "" or options.project =="" or options.db =="":
    parser.print_help()
    sys.exit(1)
    
files = list()
if options.fn_fasta != "":
    files.extend(options.fn_fasta.split(","))

if options.dn_fasta != "":
    for ext in ("fa", "fasta", "fna"):
        tmp = glob.glob("%s/*.%s" % (options.dn_fasta, ext) )    
        if len(tmp)>0:
            files.extend(tmp)
   
blasts = list()     
            
init(options.project, options.email, maxnodes = options.maxNodes)

for file in files:
    tmp = file.strip().split("/")[-1].split(".")
    if len(tmp)>=3:
        sample = tmp[1]
    else:
        sample = tmp[0]
        
    INS= {"fasta": ["%s~%s" % (file, sample)]}
    S1 = FileImport(INS)
    #S2 = FastaSort({"types": "fasta"}, [S1])
    S3 = FastaSplit({"types": "fasta", "chunk": options.split_size}, [S1])
    args = { "mode" : options.blast_mode,  
             "-db" : options.db,
             "-max_target_seqs" : "10",
             "-html": "",
             "-num_threads": "1"
             
              }
    S4 = BLAST(args, [S3])
    S5 = FileMerger("blast6", [S4], prefix = sample)
    blasts.append(S5)
    
OutputStep("FIN", "blast6", blasts)
    
##########################################################################    
#  
#################################################
##        Finish
#################################################
