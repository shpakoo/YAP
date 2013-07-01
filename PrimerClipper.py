#################################################
##     path
#################################################
import sys, tempfile, shlex, glob, os, stat, hashlib, time, datetime, re
from subprocess import *
from optparse import OptionParser
from collections import defaultdict

_author="Sebastian Szpakowski"
_date="Apr 8, 2013"
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
        

class fuzznuc:
    def __init__(self, templates, pattern, mm):
        self.templates = templates
        self.pattern = pattern 
        self.trim = defaultdict(list)
        self.mm = mm
        
        self.runFuzznuc()  
        
          
    def runFuzznuc(self):
        command="fuzznuc -sequence %s -pattern %s -rformat tagseq -outf stdout -pmismatch %s"  % (self.templates, self.pattern, self.mm) 
        #print (command)
        p = Popen(shlex.split(command), stdout=PIPE, stderr=PIPE, close_fds=True )
        out,err = p.communicate()
        self.parse(out)
        
    def parse(self, lines):
        inside = True
        b1 = ""
        b2 = ""
        b3 = ""
        counter = 0
        
        for line in  lines.split("\n"):
            line = line.strip("\n") 
             
            inside = not line.startswith("#")
            
            if inside :
                if line == "":
                    counter = -1
                else:    
                    counter+=1
                    if counter % 3 == 0:
                        b1 = b1 + line[len("nucleotide_motif")+1:]
                    elif counter % 3 == 1:
                        b2 = b2 + line[len("nucleotide_motif")+1:]
                    elif counter % 3 == 2:
                        b3 = b3 + line[len("nucleotide_motif")+1:]
                    else:
                        print "???", line                
            else:
                if len(b2)>0:  
                    #print len(b2), b2
                    #print len(b3), b3
                    self.trim[b2] = [b3.find("+"), b3.rfind("+")]
                    b1 = ""
                    b2 = ""
                    b3 = ""

#################################################
##        Functions
##

def revComp(string):
    global transtab
    string=string.upper()
    #reverse
    string = string [::-1]
    return string.translate(transtab)

#################################################
##        Arguments
##
parser = OptionParser()

parser.add_option("-i", "--infile", dest="fn_in",
                  help="input fastA file", metavar="FILE")

parser.add_option("-f", "--forward", dest="pri_f",
                  help="input primer F", metavar="forward")

parser.add_option("-r", "--reverse", dest="pri_r",
                  help="input primer R", metavar="reverse")

parser.add_option("-m", "--mm", dest="mismatches", type = "int", default = 1,
                                    help="input primer R", metavar="reverse")

(options, args) = parser.parse_args()

#################################################
##        Begin
##

from string import maketrans
inttab=  "ACGTN"
outtab = "TGCAN"
transtab = maketrans(inttab, outtab)

forward = options.pri_f
reverse = options.pri_r


fA = fuzznuc(options.fn_in, forward, options.mismatches)
#fB = fuzznuc(options.fn_in, revComp(forward), options.mismatches)
#fC = fuzznuc(options.fn_in, reverse, options.mismatches)
fD = fuzznuc(options.fn_in, revComp(reverse), options.mismatches)

a= set()
#for key, val in fA.trim.items():
#    a.add( key[val[0] : val[1]+1])
#print a    
    
newfilename = "%s.noprimers.fasta" % (".".join(options.fn_in.split(".")[:-1])) 
otpt = open(newfilename, "w")
    
for head, seq in FastaParser(options.fn_in):
    
    A = fA.trim[seq]
    #B = fB.trim[seq]  
    #C = fC.trim[seq]
    D = fD.trim[seq]
    if len(A)>0:
        start = A[1]+1
    else:
        start = 0 
        
    if len(D)>0:
        end = D[0]
    else:
        end = len(seq) 

    newseq = seq[start:end]    
    otpt.write(">%s\n%s\n" % (head, newseq)) 
otpt.close()  
    

   

#################################################
##        Finish
#################################################