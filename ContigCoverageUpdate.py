########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2013 Sebastian Szpakowski
########################################################################################


#################################################
##    need clcstatstable
##    need analogous contig file
##    add _[reads=x]_[_sites=y]_[cov=z]  in the header of each fasta entry  
#################################################
import sys
from optparse import OptionParser

_author="Sebastian Szpakowski"
_date="Jan 10, 2013"
_version="V1.00"

#################################################
##        Classes
##

    #################################################
    ### Iterator over input fasta file.
    ### Only reading when requested
    ### Useful for very large FASTA files
    ### with many sequences
    
class FastaParser:
    def __init__ (self, x, quals=False):
        self.filename = x
        self.fp = open(x, "r")    
        self.currline = "" 
        self.currentFastaName = ""
        self.currentFastaSequence = ""
        self.lastitem=False    
        if quals:
            self.linesep=" "    
        else:
            self.linesep=""    
    def __iter__(self):
        return(self)    
                
        ##### 
    def next(self):
        for self.currline in self.fp:
            if self.currline.startswith(">"):
                self.currline = self.currline[1:]
                if self.currentFastaName == "":
                    self.currentFastaName = self.currline
                else:
                    otpt = (self.currentFastaName.strip(), self.currentFastaSequence.strip())
                    self.currentFastaName = self.currline
                    self.currentFastaSequence = ""    
                    self.previoustell = self.fp.tell()
                    return (otpt)
                
            else:
                self.addSequence(self.currline)    
        
        if not self.lastitem:
            self.lastitem=True            
            return (self.currentFastaName.strip(), self.currentFastaSequence.strip())
        else:
            raise StopIteration    
                                   
    def addSequence(self, x):
               self.currentFastaSequence = "%s%s%s" % (self.currentFastaSequence,self.linesep, x.strip())            
                           
    
    def __str__():
        return ("reading file: %s" %self.filename)    
    

#################################################
##        Functions
##

    ################################################
    ### Read in a file and return a list of lines
    ###
def loadLines(x):
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
##        Arguments
##
parser = OptionParser()

#parser.add_option("-f", "--file", dest="filename",
#                  help="write report to FILE", metavar="FILE")
#parser.add_option("-q", "--quiet",
#                  action="store_false", dest="verbose", default=True,
#                  help="don't print status messages to stdout")


(options, args) = parser.parse_args()

#################################################
##        Begin
##

stats = dict()
for line in loadLines(sys.argv[-2]):
    
    line = line.strip().split()
    if len(line)==4:
        id = line[0].strip()
        stats[id] = "{0}_{1}_[cov={2}]".format(line[1], line[2], line[3]) 
    
        
outputfilename = "{0}.cov.fasta" .format( ".".join(sys.argv[-1].split("/")[-1].split(".")[:-1]) )
otpt = open(outputfilename, "w")
for head, seq in FastaParser(sys.argv[-1]):
    contigid = head.split("_")[-1]
    newcontigid = "{0}_{1}".format(head, stats[contigid])
    otpt.write(">{0}\n{1}\n".format(newcontigid, seq))
   
otpt.close()   



#################################################
##        Finish
#################################################