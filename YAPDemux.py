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

class    InfoParser:
    def    __init__(self, filename):
        self.filename = filename
        self.info = GeneralPurposeParser(filename, sep=",")
        self.store = defaultdict(list)
        for line in self.info:
            if line[0].endswith("/"):
                line[0] = line[0][:-1]
            path = "%s/%s" % (line[0].strip("\""),line[1].strip("\""))
            self.store[path].append(line[2:])
            
        self.primers = set()    
    
    def    makeOligoFile(self, x):
        tmp=defaultdict(list)
        ### group barcodes and primers
        for line in self.store[x]:
            barcode = line[0].strip().strip("\"")
            forward = line[1].strip().strip("\"")
            reverse = line[2].strip().strip("\"")
            use = line[3].upper().strip().strip("\"")
            needTwoPrimers = line[4].upper().strip().strip("\"")
            sampleID =line[5].strip().strip("\"")
                                    
            if needTwoPrimers in ("Yes", "yes", "Y", "YES"):
                if forward == "" or reverse =="":
                    print "%s: two primers needed, only one supplied F:'%s'-R:'%s' %s\n check %s " % (x, forward, reverse, self.name)
                    sys.exit(2)
                bin = "B-%s-%s" % (forward, reverse)
                tmp[bin].append((barcode, sampleID))
                
            elif use in ("N", "NONE"):
                 bin ="N-none"
                 tmp[bin].append((barcode, sampleID))
                 
            else:
                if forward != "" and "F" in use:
                    bin ="F-%s" % (forward)
                    tmp[bin].append((barcode, sampleID))
                elif forward == "" and "F" in use:
                    print "%s: demultiplexing using forward primer requested but primer sequence not specified" % (x)
                    sys.exit(3)     
                
                if reverse != "" and "R" in use:
                    bin ="R-%s" % (reverse)
                    tmp[bin].append((barcode, sampleID))
                elif reverse == "" and "R" in use:
                    print "%s: demultiplexing using reverse primer requested but primer sequence not specified" % (x)
                    sys.exit(4)
                           
        oligofiles=list() 
        for bin in tmp:
            otpt=""
            fn_otpt=""
            
            if bin.startswith("B"):
                fn_otpt = "%s.B.%s.%s.oligos" % (x.split("/")[-1], F, R)
                T, F, R = bin.strip().split("-")
                otpt += "forward\t%s\n" % (F)
                otpt += "reverse\t%s\n" % (R)
                otpt += "# %s both primers\n"

            elif bin.startswith("N"):
                fn_otpt = "%s.N.oligos" % (x.split("/")[-1])            
                otpt += "# no primers\n"

            else:
                D, P = bin.strip().split("-")
                if D == "F":
                    fn_otpt = "%s.F.%s.oligos" % (x.split("/")[-1], P)
                    otpt += "forward\t%s\n" % (P)
                    otpt += "# %s this is the FORWARD primer\n" % (x)
                else:
                    fn_otpt = "%s.R.%s.oligos" % (x.split("/")[-1], P)
                    otpt += "forward\t%s\n" % (P)
                    otpt += "# %s this is the REVERSE primer\n" % (x)

            for barcode, sampleID in tmp[bin]:
                if bin.startswith("N"):
                    otpt+= "barcode\t%s\t%s\n" % (revComp(barcode), sampleID)
                otpt+= "barcode\t%s\t%s\n" % (barcode, sampleID)
                
            otptfile = open(fn_otpt, "w")
            otptfile.write(otpt)
            otptfile.close()
            oligofiles.append(fn_otpt)
         
        
        return (oligofiles)    
                       
    def removeNonFiles(self, x):
        return x != "path/file"
    
    def getFiles(self):
        return (filter(self.removeNonFiles,  self.store.keys()))
    
    def getPrimerFilename(self):
        primerfilename =  "primers.fasta"
        
        if len(self.primers)>4:
            print "The annotation file has more than 2 primers !"
            for p in self.primers:
                print "%s" % (p.strip())
            sys.exit(5)
        
        primerfile = open(primerfilename , "w")     
        
        for p in self.primers:
            primerfile.write(p) 
        primerfile.close() 

        return (primerfilename)
            
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
    
def    demultiplex():
    forprocessing = InfoParser(options.fn_info)
    DEMUX = list()
    
    if options.strictlevel==2:
        bdiffs="0"
        pdiffs="1"
    elif options.strictlevel==3:
        bdiffs="0"
        pdiffs="2"
    elif options.strictlevel==4:
        bdiffs="1"
        pdiffs="2" 
               
    else:
        bdiffs="0"
        pdiffs="0"
        
    for file in forprocessing.getFiles():
        oligos = forprocessing.makeOligoFile(file)
        sffinputs = {"sff": ["%s" % (file)]}
    
        print oligos
        # #### extract fasta and qual files
        
        if options.useflows:
            A = MothurStep("sffinfo", options.nodesize, sffinputs, {"flow":"T"}, list())
            
        else:
            D = SFFInfoStep(sffinputs, dict(), list())    
    
    
        for oligofile in oligos:
            oligoinputs = {"oligos": ["../%s" % (oligofile)]}
            
            if options.useflows:
                args = { 
                        "bdiffs": bdiffs, 
                        "pdiffs": pdiffs, 
                        "maxhomop":"7",
                        "maxhomop":"7", 
                        "minflows":"360", 
                        "maxflows":"720"
                }        
                B = MothurStep("trim.flows", options.nodesize, oligoinputs, args, [A])
                C = MothurSHHH([B])
                D = FileMerger("fasta,name,qfile", [C])
                
            #### deconvolution
            args = {    "flip":getOligoFlip(oligofile), 
                        "bdiffs": bdiffs, 
                        "pdiffs": pdiffs, 
                        "minlength": "%s" % (options.minlength),
                        "maxlength": "2500",
                        "qtrim": "F", 
                        "maxambig":"0", 
                        "maxhomop":"7"
            }    
            if options.useflows:
                args["force"] =  "name,fasta,oligos"
            
            E = MothurStep("trim.seqs", options.nodesize, oligoinputs, args, [D])    
            
            if  options.useflows:
                DEMUX.append(E)
            
            else:
                #### generate trimming coordinates using LUCY
                F = LUCYcheck(options.nodesize, [E])

                #### trim the reads using LUCY-generated coordinates
                G = LUCYtrim([F])
                        
                DEMUX.append(G)
    

    A1 = FileMerger("fasta,name,group,qfile", DEMUX)
    A2 = MatchGroupsToFasta(dict(), [A1])
    
    args = {"mingroupmembers": options.mingroupmembers, 
            "report": "failing"}
    A3 = GroupRetriever(args, A2)
    
    args = {
            "force" : "fasta,name,group",
            "find": "groups" 
            }      
    A4 = MothurStep("remove.groups", options.nodesize, dict(), args, [A3])
    
    return ([A4, forprocessing.getPrimerFilename()])



        
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
                 help="mapping: file, barcode, primer, sample information. File should be in CSV format", metavar="allinfo.csv")

parser.add_option_group(group)

group = OptionGroup(parser, "Optional Configuration", description="parameters to alter if necessary")


group.add_option("-a", "--annotations", dest="dir_anno", default="/usr/local/devel/ANNOTATION/sszpakow/ANNOTATION/",
                 help="directory that stores auxilliary files\n[%default]", metavar="annotations")
group.add_option("-m", "--minlen", dest="minlength", default=220, type="int",
                 help="what is the minimum length of reads to process\n[%default]", metavar="#")     

group.add_option("-g", "--mingroupsize", dest="mingroupmembers", default=100, type="int",
                 help="after demultiplexing, discard groups with fewer reads than #\n[%default]", metavar="#")
            
group.add_option("-x", "--strict", dest="strictlevel", default=2, type="int",
                 help="""how strict to be at demultiplexing: 
1 very strict (barcode no mismatches, primer no mismatches) 
2 less strict (barcode no mismatches, primer 1 mismatch allowed)
3 even less strict (barcode no mismatches, primer 2 mismatches allowed)
4 last resort (barcode 1 mismatch, primer 2 mismatches allowed)
[%default]""", metavar="#")                 

group.add_option("-F", "--useFlows", dest="useflows", action="store_true", default=False,
                 help="""if specified flows will be used and pyronoise will be employed. Otherwise qual values will be used with LUCY""", metavar="#") 

parser.add_option_group(group)

group = OptionGroup(parser, "Technical", description="could be useful sometimes")
group.add_option("-C", "--NODESIZE", dest="nodesize", default=30,
                 help="maximum number of grid node's CPUs to use\n[%default]", metavar="#")

parser.add_option_group(group)
    
(options, args) = parser.parse_args()

#################################################
##        Begin
##

if options.fn_info == "" or options.email == "" or options.project =="":
    parser.print_help()
    sys.exit(1)
        
init(options.project, options.email)

############################

supplementary = list()
READY, primerfile = demultiplex()
OutputStep("DEMUX", "groupstats,fasta,group,name,list,pdf,svg,tiff,taxsummary,globalsummary,localsummary", [READY])

gs = GroupSplit({}, {}, [READY])
OutputStep("GROUPSPLIT", "fasta", [gs])

    
##########################################################################    
#  
#################################################
##        Finish
#################################################
