########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2013 Sebastian Szpakowski
########################################################################################

#################################################
##    A pipeline for finding DAT files in the mascot directory
##    filtering which of these files match a pattern of IDs
##    finally using a web crawler to cache and generate text report.
#################################################
import sys
from optparse import OptionParser, OptionGroup
from StepsLibrary import *
from StepsLibrary_EXP import *
from collections import defaultdict
from Queue import Queue

_author="Sebastian Szpakowski"
_date="2012/11/29"
_version="Version 1"

#################################################
##        Classes
##
    #################################################
    ### Iterator over input file.
    ### 
class GeneralPurposeParser:
    def __init__(self, file, skip=0, sep="\t"):
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
        self.end()
        raise StopIteration
    
    def end(self):
        self.fp.close() 
        
    def getHeader(self,n=10):
        header = n
        otpt = list()
        for currline in self.fp:
            header -=1
            currline = currline.strip().split(self.sep)
            otpt.append(currline) 
            
            if header==0:
                return (otpt)
            
        self.end()        
        return otpt
                       
    def __str__(self):
        return "%s [%s]\n\t%s" % (self.filename, self.linecounter, self.currline)

class InfoParser:
    def __init__(self, filename, directory):
        self.filename = filename
        self.info = GeneralPurposeParser(filename, sep="\t")
        self.patterns = [ "".join(k).strip() for k in self.info ]
        self.patterns = [k for k in self.patterns if len(k)>5]
        
        self.dir = directory
        
        # tuples filename/serverpath/sampleid
        self.q=Queue()
        if len(self.patterns)>0:
            self.findFiles()
    
    def findFiles(self):
        print "Searching for files..."
        files = glob.glob("%s/2011*/*.dat" % (self.dir))
        files.extend(glob.glob("%s/2012*/*.dat"  % (self.dir)))
        files.extend(glob.glob("%s/2013*/*.dat"  % (self.dir)))

        for file in files:
            header = 100
            serverfilepath = file[len(self.dir):]
            filename = file.split("/")[-1]
            
            try:
                p = GeneralPurposeParser(file, sep="=")
                for line in p.getHeader(n=10):
                    if line[0].strip() == "COM":
                        info = line[1]
                        info = info[info.find(" from ")+6: info.find(" by ")]
                        info = info.replace(" " , "-")
                        
                        if filename in self.patterns:
                            print "[file]", serverfilepath, info
                            self.q.put((serverfilepath, info))
                        else:
                            for k in self.patterns:
                                #if info.find(k)>-1 and info.find("2TO20")>-1:
                                if info.find(k)>-1:   
                                    print "[patt]", serverfilepath, info
                                    self.q.put((serverfilepath, info)) 
                            
                            
                p.end()  
            except:
                pass
            
            
    def __iter__(self):   
        return (self)

    def next(self):
        if not self.q.empty():
            return self.q.get()
        else:    
            raise StopIteration
            
#################################################
##        Functions
##
    ################################################
    ### Read in a file and return a list of lines
    ###
def    loadLines(x):
    try:
        fp = open(x, "rU")
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

group = OptionGroup(parser, "Required", description="Will not run without these !")

group.add_option("-P", "--PROJECT", dest="project", default="",
                 help="project code", metavar="#")
group.add_option("-E", "--EMAIL", dest="email", default="",
                 help="e-mail address", metavar="@")                 
group.add_option("-i", "--info", dest="fn_info", default="",
                 help="a file with patterns (IDs) that the .dat files in mascot directory will be searched for. One id/pattern per line. For example ASIRIW will match any file with a substring ASIRIW as in 2DLCT1D_ASIRIW_120928VP-L-2TO20. Patterns shorter than 5 characters are not allowed.", metavar="patterns.txt")

parser.add_option_group(group)

group = OptionGroup(parser, "Optional Configuration", description="parameters to alter if necessary")
group.add_option("-m", "--mascot", dest="dir_mascot", default="/local/pfgrc_mascot/data/",
                 help="directory that stores the mascot dat files\n[%default]", metavar="mascot.dat")
parser.add_option_group(group)

group = OptionGroup(parser, "Technical", description="could be useful sometimes")
group.add_option("-C", "--NODESIZE", dest="nodesize", default=30,
                 help="maximum number of grid node's CPUs to use\n[%default]", metavar="#")
parser.add_option_group(group)
    
(options, args) = parser.parse_args()

#################################################
##        Begin
##

if options.email == "" or options.project == "":
    parser.print_help()
    sys.exit(1)
    
   
    
############################
IP =  InfoParser(options.fn_info, options.dir_mascot)
init(options.project, options.email) 


final = list()
for file, id in IP:
    args = {"file": file, "id": id}
    S1 = MascotReportLifter(args, [])
    final.append(S1)
    
OutputStep("REPORTS", "txt", final)
OutputStep("PDFscreens", "pdf", final)
OutputStep("PNGscreens", "png", final)
OutputStep("LOGS", "log", final)

   
##########################################################################    
#  
#################################################
##        Finish
#################################################
