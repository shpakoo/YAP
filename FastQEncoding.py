########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2013 Sebastian Szpakowski
########################################################################################


#################################################
##        File: FastQEncoding
#################################################
import sys
import argparse

_author="Sebastian Szpakowski"
_date="Oct 12, 2012"
_version="Version 1"

#################################################
##        Classes
##
        
    #################################################
    ### Iterator over input fastq file.
    ### Only reading when requested
    ### Useful for very large files
    ### with many sequences
class    FastqParser:
    def    __init__ (self, x, quals=False):
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
    def    __iter__(self):
        return(self)    
                
        ##### 
    def    next(self):
        for self.currline in self.fp:
            if self.currline.startswith("+") or self.currline.startswith("@"):
                #self.currline = self.currline[1:]
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
    def __str__(self):
        return ("reading file: %s" % self.filename)    


#  SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................
#  ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX......................
#  ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................
#  .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ......................
#  LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL....................................................
#  !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
#  |                         |    |        |                              |                     |
# 33                        59   64       73                            104                   126
#
# S - Sanger        Phred+33,  raw reads typically (0, 40)
# X - Solexa        Solexa+64, raw reads typically (-5, 40)
# I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
# J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
#    with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold) 
#    (Note: See discussion above).
# L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)
class FastQDecoder:
    def __init__(self):  
        self.encodings = {   
                        "64"  : set(list("""@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"""        )), 
                        "33"  : set(list("""!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ"""     ))
                     }
        
    def guessEncoding(self, input):
        guesses = set()
        #print ""
        #print "\t",  self.format(input)
        for name, quals in self.encodings.items():
            tmp = input.difference(quals)
#            print name, "\t", self.format(quals)
#            print "\t", self.format(tmp)
            #print name, "\t", self.format(tmp)
            if len(tmp)==0:
                guesses.add(name)
            
        return guesses
        
    def format(self, input):
        tmp = list(input)
        tmp.sort
        return "".join(tmp)

#################################################
##        Functions
##

    
    

#################################################
##        Arguments
##

#import argparse
#parser = argparse.ArgumentParser(description='Process some integers.')
#parser.add_argument('integers', metavar='N', type=int, nargs='+',
#                  help='an integer for the accumulator')
#parser.add_argument('--sum', dest='accumulate', action='store_const',
#                   const=sum, default=max,
#                  help='sum the integers (default: find the max)')
#args = parser.parse_args()
#print args.accumulate(args.integers)

#################################################
##        Begin
##

observedQs = set()
decoder = FastQDecoder()

encodings = set(["64", "33"])
count = 0

guessing = True

for head, seq in FastqParser(sys.argv[-1]):
    if head.startswith("+"):
        observedQs = observedQs.union(list(seq))
        encodings = decoder.guessEncoding(observedQs)
        if len(encodings)==1:
            count+=1
            
        if count >100 and len(encodings)==1:
            print list(encodings)[0]
            guessing = False
            break;    
 
if guessing:
    print "33\n"       

#################################################
##        Finish
#################################################