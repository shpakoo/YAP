########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2013 Sebastian Szpakowski
########################################################################################

#################################################
##     Calculate coverage of an ORF
#################################################
import sys
from optparse import OptionParser
from collections import defaultdict

_author="Sebastian Szpakowski"
_date="Jan 30, 2013"
_version="V1.00"

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

    
#################################################
##        Arguments
##
parser = OptionParser()

parser.add_option("-p", "--prok", dest="fn_prok",
                  help="ORFs", metavar="FILE")
parser.add_option("-w", "--weights", dest="fn_weight",
                  help="annotation", metavar="FILE")

parser.add_option("-o", "--outputfile", dest="fn_output",
                  help="output file", metavar="FILE")

(options, args) = parser.parse_args()

#################################################
##        Begin
##

cache=dict()
for line in GeneralPurposeParser(options.fn_weight):
    cache[line[0]] = line[1:]
    
otpt = open(options.fn_output, "w")    
    
for line in GeneralPurposeParser(options.fn_prok):
    if cache.has_key(line[0]):
        otpt.write("%s\t%s\t%s\n" % (line[0], "\t".join(cache[line[0]]), "\t".join(line[1:]) ) ) 
    else:
       print line[0]
otpt.close()
#################################################
##        Finish
#################################################