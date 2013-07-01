########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2013 Sebastian Szpakowski
########################################################################################

#################################################
##     Coverage of ORF
#################################################
import sys
from optparse import OptionParser
from collections import defaultdict

_author="Sebastian Szpakowski"
_date="Jan 29, 2013"
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
def overlap(s,e,S,E):
    return not (e<S or s>E) 
    
#################################################
##        Arguments
##
parser = OptionParser()

parser.add_option("-o", "--orf", dest="fn_orf",
                  help="ORFs", metavar="FILE")
parser.add_option("-a", "--anno", dest="fn_anno",
                  help="annotation", metavar="FILE")

parser.add_option("-e", "--outputfile", dest="fn_output",
                  help="output file", metavar="FILE")

#parser.add_option("-q", "--quiet",
#                  action="store_false", dest="verbose", default=True,
#                  help="don't print status messages to stdout")


(options, args) = parser.parse_args()

#################################################
##        Begin
##

### check which ORFS don't need to have coverage recalculated.

needed = list()
lookup = defaultdict(list)
cache = defaultdict(float)
otpt = open(options.fn_output, "w")

for head, seq in FastaParser(options.fn_orf):
    head = head.strip().split("_")
    id, tid, length, reads, cov, start, end, strand = head
    length = float(length)
    start = float(start)
    end = float(end)
    
    if end - start > length-10 and end - start < length+10 :
        otpt.write( "%s\t%s\n" % ( "_".join(head), cov.strip("]").strip("[").split("=")[-1]))
        pass
    else: 
        needed.append(int(tid))  
        lookup[int(tid)].append(("_".join(head), int(start), int(end)))
        
print "easy part done..."
#### cache reads for the "needed" contigs  
#   0         SOLEXA4:32:C13WTACXX:1:1102:10260:39080/1  101    0   51 13360      6632      6683  0   1  0  51
#   1         SOLEXA4:32:C13WTACXX:1:1102:10260:39080/2  101    2  101 23163       254       353  1   1  0  99 
for x in GeneralPurposeParser(options.fn_anno, sep=None):
    tid = int(x[5])+1
    
    if (float(x[0]) % 1000000 == 1):
            print  x 
                     
    for curorf in lookup[tid]:
        if overlap(int(x[6])+1,int(x[7])+1, curorf[1], curorf[2]):
            cache[curorf[0]]+= (int(x[4]) - int(x[3]))
    
            print curorf
            print x
            print cache[curorf[0]]
            print
            
            
#           
#        if len(x)==12:
#            qnum, qid, qlen, qstart, qend, tid, tstart, tend, qrev, qmatch, paired, score = x 
#        else:
#            qnum, qid, qlen, qstart, qend, tid, tstart, tend, qrev, qmatch, score = x
#        ### tid is a number, starting at 0!
#        contigs[int(tid)+1] += 1
#        

# process cache  

for orfid in cache.keys():
    orf = orfid.strip().split("_")
    s = int(orf[-3])
    e = int(orf[-2])
    
    otpt.write( "%s\t%.2f\n" % (orfid, cache[orfid]/(e-s))  )
    
otpt.close()     
print "\n"
   

#################################################
##        Finish
#################################################