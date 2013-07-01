########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2013 Sebastian Szpakowski
########################################################################################

#################################################
##    The next steps...
#################################################
import sys
from optparse import OptionParser, OptionGroup
from StepsLibrary import *
from StepsLibrary_EXP import *
from collections import defaultdict
from Queue import Queue

_author="Sebastian Szpakowski"
_date="2012/10/01"
_version="Version 1"

#################################################
##        Classes
##
    #################################################
    ### Iterator over input file.
    ### every line is converted into a dictionary with variables referred to by their 
    ### header name
class     GeneralPurposeParser:
    def    __init__(self, file, skip=0, sep="\t"):
        self.filename = file
        self.fp = open(self.filename, "rU")    
        self.sep = sep
        
        self.linecounter = 0
        self.currline=""
        
    def    __iter__(self):
        return (self)
    
    def    next(self):
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
            
            if reverse =="" or forward =="":
                print "%s: please provide both primers for barcode:'%s' " % (x, barcode)
                sys.exit(1) 
            else:    
                self.primers.add(">_primer_F\n%s\n" % (forward))
                self.primers.add(">_primer_F_rc\n%s\n" % (revComp(forward)))
                self.primers.add(">_primer_R\n%s\n" % (reverse))
                self.primers.add(">_primer_R_rc\n%s\n" % (revComp(reverse)))
                        
            if needTwoPrimers in ("Yes", "yes", "Y", "YES"):
                if forward == "" or reverse =="":
                    print "%s: two primers needed, only one supplied F:'%s'-R:'%s' %s\n check %s " % (x, forward, reverse, self.name)
                    sys.exit(2)
                bin = "B-%s-%s" % (forward, reverse)
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
                fn_otpt = "%s.B.oligos" % (x.split("/")[-1])
                T, F, R = bin.strip().split("-")
                otpt += "forward\t%s\n" % (F)
                otpt += "reverse\t%s\n" % (R)
                otpt += "# %s both primers\n"

            else:
                D, P = bin.strip().split("-")
                if D == "F":
                    fn_otpt = "%s.F.oligos" % (x.split("/")[-1])
                    otpt += "forward\t%s\n" % (P)
                    otpt += "# %s this is the FORWARD primer\n" % (x)
                else:
                    fn_otpt = "%s.R.oligos" % (x.split("/")[-1])
                    otpt += "forward\t%s\n" % (P)
                    otpt += "# %s this is the REVERSE primer\n" % (x)

            for barcode, sampleID in tmp[bin]:
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

parser.add_option_group(group)

group = OptionGroup(parser, "Optional Configuration", description="parameters to alter if necessary")

group.add_option("-Y", "--Yap", dest="mode", default="16S",
                 help="""Which Pipeline: 16S ITS [%default]""", metavar="#") 
group.add_option("-a", "--annotations", dest="dir_anno", default="/usr/local/devel/ANNOTATION/sszpakow/ANNOTATION/",
                 help="directory that stores auxilliary files\n[%default]", metavar="annotations")

parser.add_option_group(group)

group = OptionGroup(parser, "Technical", description="could be useful sometimes")
group.add_option("-C", "--NODESIZE", dest="nodesize", default=30,
                 help="maximum number of grid node's CPUs to use\n[%default]", metavar="#")
parser.add_option_group(group)
    
(options, args) = parser.parse_args()

#################################################
##        Begin
##

    
if options.email == "" or options.project =="":
    parser.print_help()
    sys.exit(1)
    
if not options.mode in ("16S", "ITS"):
    parser.print_help()
    sys.exit(2)    

### parameters specific to YAP incarnations

### 16S V1-V3    
if options.mode=="16S":
    ### file in the annotations directory that has reference sequences
    _referenceseq = "thermo.fasta"
    
### ITS NSI1 - NLB4 (barcoded)    
if options.mode=="ITS":
    _referenceseq = "yeastITS.fasta"
            
init(options.project, options.email)

print "We are in %s mode" % (options.mode) 

############################
######################
#### reference: 

inputs = {"fasta": ["%s/%s" % (options.dir_anno, _referenceseq)] }
REF = FileImport(inputs)

end =  list()

for f in glob.glob("*clean.fasta"):
    
    dist = f.split("otureps_")[1].split("_")[0].strip("_")
    inputs = {"fasta": ["%s/%s" % (os.getcwd(),f)] }
    S1 = FileImport(inputs)
    S2 = FileMerger("fasta", [S1, REF], prefix = "OTUDIST_%s" % (dist))
    args = {"-mode": "quickaln"}
    #args = {"-mode": "expresso",
    #       "-email=%s" % (options.email) : " "
    #        }
    S3 = TCOFFEE({}, args, [S2])
    
    args = {
            "-tree": "",
            "-bootstrap=10000": "",
            "-clustering=NJ": "",
            "-outputtree=phylip": "",
            "-QUIET": ""
            }
    
    S4 = CLUSTALW2({}, args, [S3])
    
    end.append(S4)

OutputStep("GUIDETREES", "tre,tree,dnd", end)  
OutputStep("PHYLOTREES", "ph", end)

forotutable = defaultdict(list)
for f in glob.glob("*.group"):
    forotutable["group"].append("%s/%s" % (os.getcwd(), f))
for f in glob.glob("*.list"):
    forotutable["list"].append("%s/%s" % (os.getcwd(), f))   

INS = FileImport(forotutable)
x = OtuTable({}, {}, [INS])
OutputStep("OTUTABLES", "otutable", x)

   

    
##########################################################################    
#  
#################################################
##        Finish
#################################################
