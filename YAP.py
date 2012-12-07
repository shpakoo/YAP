########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2012 J.Craig Venter Institute.
########################################################################################

#################################################
##    A pipeline for cleaning up sff files and producing 
##    OTUs (certain regions of 16S and ITS supported)
#################################################
import sys, os.path
from optparse import OptionParser, OptionGroup
from StepsLibrary import *
from collections import defaultdict
from Queue import Queue

_author="Sebastian Szpakowski"
_date="2012/12/06"
_version="Version 4"

#################################################
##        Classes
##
    #################################################
    ### Iterator over input file.
    ### every line is converted into a dictionary with variables referred to by their 
    ### header name
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
        raise StopIteration
                    
    def __str__(self):
        return "%s [%s]\n\t%s" % (self.filename, self.linecounter, self.currline)

class InfoValidator:
    def __init__(self,filename):
        self.filename = filename
        self.info = GeneralPurposeParser(filename, sep=",")
        self.URI = "http://confluence/display/~sszpakow/YAP"
        self.epilogue = "\n***\tPlease correct before continuing...\n***\t{0}\n".format(self.URI) 
        self.header = ""
        
        self.files,  self.barcodes ,self.primersF, self.primersR, self.sampleIDs = self.parse()
        print ("***\tValidation complete, no obvious errors found.\n")
       
        
        
    def parse(self):
        counter=0;
        print ("\n***\tValidating your template\n\t{0} ...\n".format(self.filename))
        files = set()
        barcodes = set()
        primersF = set()
        primersR = set()
        sampleIDs = set()
                
        for line in self.info:
            if counter == 0:
                
                self.header = line
                has =  ",".join (self.header[:8])
                needed = "path,file,barcode,forward,reverse,use,both,SampleID"
                if not has.lower() == needed.lower():
                    self.error( "Your template's header is incorrect or missing:\nhas :\t{0}\nneed:\t{1}".format(has, needed), 101)
                if not self.header[7] == "SampleID":    
                    self.error( "Your template has\n\t'{0}' instead of \n\t'SampleID' as 8th column's header.".format(self.header[7]), 102)
                    
            else:
                files.add("{0}/{1}".format(line[0], line[1]))   
                barcodes.add(line[2])   
                primersF.add(line[3])    
                primersR.add(line[4])
                sampleIDs.add(line[7])
            counter+=1 
            
        ##### files
        for f in files:
            if not os.path.isfile(f):
                self.error("file doesn't exist\n\t{0}".format(f), 103)

        ##### F primers
        if len(primersF)>1:
            self.error("Multiple forward primers specified:\n\t{0}\n\tnot supported in the current version of YAP".format("\n\t".join(primersF)), 104)
        
        if list(primersF)[0].strip() =="" :
            self.error("Forward primer should not be empty", 104)
        
        
        ##### R primers
        if len(primersF)>1:
            self.error("Multiple reverse primers specified:\n\t{0}\n\tnot supported in the current version of YAP".format("\n\t".join(primersR)), 105)
        
        if list(primersR)[0].strip() =="" :
            self.error("Reverse primer should not be empty", 105)
        
        ##### sampleIDs
        spaces = set()
        ill = ("\\","/", "~", "-", "+", "#")
        illegalchars = set()
        digitstart = set()
        for s in sampleIDs:
            if s.count(" ")>0:
                spaces.add(s)
            for k in ill:
                if s.count(k)>0:
                    illegalchars.add(s)
            if s[0].isdigit():
                digitstart.add(s)
         
        hint = "*You could create two columns: \n\tSampleID, compliant with YAP (excel function: SUBSTITUTE()) and\n\tOriginalIDs, where any character is allowed."    
        if len(spaces)>0:
            M = "The following samplesID(s) have spaces in them:\n\t"
            for s in spaces:
                M = "{0}'{1}',".format(M, s) 
            M = "{0}\n\n\t{1}".format(M, hint)    
            self.error(M, 106)    
            
        if len(illegalchars)>0:
            M = "The following samplesID(s) have illegal chars in them {0}:\n\t".format(", ".join(ill))
            for s in illegalchars:
                M = "{0}'{1}',".format(M, s) 
            
            M = "{0}\n\n\t{1}".format(M, hint)    
            self.error(M, 107)   
            
        if len(digitstart)>0:
            M = "The following samplesID(s) start with numbers:\n\t".format(", ".join(ill))
            for s in digitstart:
                M = "{0}'{1}',".format(M, s) 
                 
            M = "{0}\n\n\t{1}".format(M, hint)    
            self.error(M, 108)  
            
            
        return (files, barcodes, primersF, primersR, sampleIDs)    
                       
                  
    def error(self, message, code):
        print "!!!\t{0}\n{1}".format(message, self.epilogue)
        sys.exit(code)
        
    def getTrimpoints(self):
        primers = self.primersF.union(self.primersR)
        if "AGAGTTTGATYMTGGCTCAG" in primers and "ATTACCGCGGCTGCTGG" in primers:
            return "1044", "13127", "1044-13127"
        else:
            return "0", "0", "unknown"

       
class InfoParser:
    def __init__(self, filename):
        self.filename = filename
        self.info = GeneralPurposeParser(filename, sep=",")
        self.store = defaultdict(list)
        for line in self.info:
            if line[0].endswith("/"):
                line[0] = line[0][:-1]
            path = "%s/%s" % (line[0].strip("\""),line[1].strip("\""))
            self.store[path].append(line[2:])
            
        self.primers = set()    
    
    def makeOligoFile(self, x):
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
                sys.exit(11) 
            else:    
                self.primers.add(">_primer_F\n%s\n" % (forward))
                self.primers.add(">_primer_F_rc\n%s\n" % (revComp(forward)))
                self.primers.add(">_primer_R\n%s\n" % (reverse))
                self.primers.add(">_primer_R_rc\n%s\n" % (revComp(reverse)))
                        
            if needTwoPrimers in ("Yes", "yes", "Y", "YES"):
                if forward == "" or reverse =="":
                    print "%s: two primers needed, only one supplied F:'%s'-R:'%s' %s\n check %s " % (x, forward, reverse, self.name)
                    sys.exit(12)
                bin = "B-%s-%s" % (forward, reverse)
                tmp[bin].append((barcode, sampleID))
            else:

                if forward != "" and "F" in use:
                    bin ="F-%s" % (forward)
                    tmp[bin].append((barcode, sampleID))
                elif forward == "" and "F" in use:
                    print "%s: demultiplexing using forward primer requested but primer sequence not specified" % (x)
                    sys.exit(13)     
                
                if reverse != "" and "R" in use:
                    bin ="R-%s" % (reverse)
                    tmp[bin].append((barcode, sampleID))
                elif reverse == "" and "R" in use:
                    print "%s: demultiplexing using reverse primer requested but primer sequence not specified" % (x)
                    sys.exit(14)
                           
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
            sys.exit(15)
        
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
    
def demultiplex():
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

def finalize(input):
    clean = CleanFasta(dict(), [input])
    
    ####### remove sequences that are too short, and with ambiguous bases 
    args = { "minlength" : "%s" % ( options.minlength ),
             "maxambig" : "0",
             "force": "fasta,name,group"}
    clean2 = MothurStep("screen.seqs", options.nodesize, dict(), args, [clean])

    args = {"mingroupmembers": 0, 
            "report": "passing"}
    clean2a = GroupRetriever(args, [clean2])
    OutputStep("2-PRE454", "groupstats,fasta,group,name,list,svg,pdf,tiff,taxsummary,globalsummary,localsummary", clean2a)

    ###################### CDHIT-454
    #### unique and de-noise
#    args=         { 
#                    "c" : "1.0",
#                    "b" : "8",
#                    "aS": "1.0",
#                    "g" : "1",
#                    "M"    : "50000",
#                    "T" : "%s" % (options.nodesize)
#                }
    
    ### aggressive denoising:
    args=         { 
                    "c" : "0.98",
                    "b" : "10",
                    "aS": "0.0",
                    "g" : "1",
                    "M"    : "0",
                    "T" : "%s" % (options.nodesize)
                }
    #### de-noise/unique collapse            
    CD_1 = CDHIT_454(options.nodesize, args, [clean2])
    CD_2 = CDHIT_Mothurize(dict(), CD_1)
    
    
    args = {"mingroupmembers": 0, 
            "report": "passing"}
    CD_2a = GroupRetriever(args, [CD_2])
    OutputStep("3-UNIQUE", "groupstats,tre,fasta,group,name,list,svg,pdf,tiff,taxsummary,globalsummary,localsummary", CD_2a)  
    
    #### add reference sequences to the merged experiments' file
    CD_3 = FileMerger("fasta,name,group,qfile", [CD_2, REF_1, REF_2, REF_3, PRI_1, PRI_2, PRI_3])
    
    #### align to reference database
    inputs = {"reference": ["%s/%s" % (options.dir_anno, _alignment)] }
    args = {    "flip":"t", 
                "ksize": "8"
            }   
     
    CD_4 = MothurStep("align.seqs", options.nodesize, inputs, args, [CD_3])
    
    #### AlignmentSummary determining alignment trimming options 
    #### sets trimstart and trimend variables that can be used by in subsequent steps.
    #### threshold means to keep the center part of the alignment with at least 
    #### the fraction of maximum coverage  
    args = {"ref": _referenceseqname, "thresh": options.dynthresh}
    CD_5 = AlignmentSummary(args,[CD_4])
    
    #### alignment plots  
    if _trimstart != _trimend:
        args = {"ref": _referenceseqname, 
                "trimstart" : _trimstart,  
                "trimend" : _trimend
                }
    else:  
        args = {"ref": _referenceseqname, 
                "trimstart" : "find",  
                "trimend" : "find"
                }        
    CD_6 = AlignmentPlot(args,[CD_5])
    

    #supplementary.append(CD_5)
    supplementary.append(CD_6)
    ###########################
    
    args = {"mingroupmembers": 0, 
            "report": "passing"}
    CD_4a = GroupRetriever(args, [CD_4])
    OutputStep("4-ALIGNED", "groupstats,tre,fasta,group,name,list,svg,pdf,tiff,taxsummary,globalsummary,localsummary", CD_4a)    
       
    cleanCD = cleanup(CD_5)
    args = {"mingroupmembers": 0, 
            "report": "passing"}
    cleanCDa = GroupRetriever(args, [cleanCD])
    OutputStep("5-CLEAN", "groupstats,fasta,group,name,list,svg,pdf,tiff,taxsummary,globalsummary,localsummary", cleanCDa)
    
    clusterCD = CDHITCluster(cleanCD)
    
    x = plotsAndStats(clusterCD)
    INS = {"annotation" : [options.fn_info]}
    ARGS = {"dist": "0.03"}
    output1 = R_defaultplots(INS, ARGS, x)
    output2 = AnnotateClusters(dict(), dict(), output1)
        
    return (output2)

def cleanup(input):

    ### remove the "ref" group
    args = {
            "force" : "fasta,name,group",
            "groups": "ref" 
            }
            
    s15 = MothurStep("remove.groups", options.nodesize, dict(), args, [input])
    
    ####### remove sequences that are too short (bad alignment?)  
    args = {
                "minlength" : "%s" % (options.minlength), 
                "maxambig" : "0",
                "force" : "fasta,name,group" ,
            }
    s16 = MothurStep("screen.seqs", options.nodesize, dict(), args, [s15])
    
    ####### find chimeric sequences  
    toremove = list()
    for ch in [ "uchime" ]:
        ### chimeras against reference
        args = {"force" : "fasta,reference"}
        inputs = {"reference": ["%s/%s" % (options.dir_anno, _alignment)] }
        
        A = MothurStep("chimera.%s" % (ch),options.nodesize, inputs, args, [s16])    
        toremove.append(A)
        
        if not options.quickmode:
            ### chimeras against self
            args ={"force": "name,group,fasta"}
            inputs = {}
            
            A = MothurStep("chimera.%s" % (ch),options.nodesize, inputs, args, [s16])    
            toremove.append(A)
        
    ### merge all accnos files and remove ALL chimeras    
    allchimeras = FileMerger("accnos", toremove)
    s17 = MothurStep("remove.seqs",options.nodesize, dict(), dict(), allchimeras)
    
    #### if primer trimming points are not unknown
    if _trimstart!=_trimend:
        ### primer cut
        args = {
                    "s" : _trimstart, 
                    "e": _trimend,
                }
    else:
        args = {
                "s" : "find:trimstart",  
                "e" : "find:trimend"
        }
     
        
    s18a = AlignmentTrim(dict(), args, [s17])
            
    ####### remove sequence fragments, bad alignments (?) 
    args = { "minlength" : "%s" % (options.minlength),
             "force": "fasta,name,group"}
    s18b = MothurStep("screen.seqs", options.nodesize, dict(), args, [s18a])
    
    ### build a tree
    #s18b_tree = ClearcutTree({}, s18b)
    
    ####### remove empty columns
    args = {"vertical" : "T"}
    s19 = MothurStep("filter.seqs",options.nodesize, dict(), args, [s18b]) 
    
    ####### taxonomy
    inputs = {    "reference": ["%s/%s" % (options.dir_anno,_trainset)],
                "taxonomy": ["%s/%s" % (options.dir_anno, _taxonomy )]
            }
            
    args = {    "iters" : "100",
            "cutoff":  "60"
            }
    s20 = MothurStep("classify.seqs", options.nodesize, inputs, args, [s19])
    
    ### remove - and . for subsequent clustering efforts 
    s21 = CleanFasta(dict(), [s20])
    
    return (s21)

def CDHITCluster(input):
    cdhits = list()
    for arg in ["0.99", "0.97", "0.95", "0.90"]:
        args = {"c": arg,
                "n": "8",
                "g": "1",
                "M": "10000",
                "T": "%s" % (options.nodesize) 
                }
        
        CD_1 = CDHIT_EST(options.nodesize, args, [input])
        
        ### make sth. analogous to mothur's labels
        arg = 1.0 - float(arg)
        if arg == 0:
            arg = "unique"
        else:
            arg = "%s" % (arg)
        
        args = {"mode": arg    
                }
        CD_2 = CDHIT_Mothurize(args, CD_1)
        CD_2a = CDHIT_Perls({}, CD_2)            
        cdhits.append(CD_2)
                
    READY = FileMerger("list,rabund,sabund", cdhits)    
    SORTED = FileSort("list,rabund,sabund", READY)
    return (SORTED)

def plotsAndStats(input):
    
    ### all groups!
    args = {"mingroupmembers": 0, 
            "report": "passing"}
    s23 = GroupRetriever(args, [input])
    
    ######## make a shared file 
    args = {"label" : "0.01-0.03-0.05-0.1", "find": "groups"}
    s24 = MothurStep("make.shared", options.nodesize, dict(), args, [s23])


    args = {
            "label" : "0.01-0.03-0.05-0.1",
            "basis" : "otu"
            }
            
    s25a= MothurStep("classify.otu",  options.nodesize, dict(), args, [s24])
    args = {
                "taxonomy": "otu.taxonomy",
                "taxsummary": "otu.taxsummary"                
            }
    s25aa = FileType(args, [s25a])
    
    args = {
            "label" : "0.01-0.03-0.05-0.1",
            "basis" : "sequence"
            }
            
    s25b = MothurStep("classify.otu",  options.nodesize, dict(), args, [s24])
    args = {
                "taxonomy": "seq.taxonomy",
                "taxsummary": "seq.taxsummary"                
            }
    s25bb = FileType(args, [s25b])
    
    args = {"force" : "list", "calc": "nseqs-sobs-simpson-invsimpson-chao-shannon-shannoneven-coverage"}
    s26 = MothurStep("summary.single",options.nodesize, dict(), args, [s25bb])
    
    args = {"summary": "globalsummary"}
    s26a = FileType(args, [s26])
    
    args = {"force" : "shared", "calc": "nseqs-sobs-simpson-invsimpson-chao-shannon-shannoneven-coverage"}
    s27 = MothurStep("summary.single", options.nodesize, dict(), args, [s25bb])
    
    args = {"summary": "localsummary"}
    s27a = FileType(args, [s27])
    
    args = {"force" : "shared", "calc": "thetayc-jclass-braycurtis"}
    s28 = MothurStep("tree.shared", options.nodesize, dict(), args, [s24]) 
    
    supplementary.append(s28)
    
    args = {"force" : "list", "calc": "nseqs-sobs-simpson-invsimpson-chao-shannon-shannoneven-coverage", "freq": "0.01"}
    s29 = MothurStep("rarefaction.single", options.nodesize, dict(), args, [s24])
    #return ([s23, s24, s25aa, s25bb, s26a, s27a, s28, s29])
    
    args = {"force" : "shared", "calc": "nseqs-sobs-simpson-invsimpson-chao-shannon-shannoneven-coverage", "freq": "0.05"}
    s30 = MothurStep("rarefaction.single",options.nodesize, dict(), args, [s24]) 
    return ([s23, s24, s25aa, s25bb, s26a, s27a, s28, s29, s30])
    
    
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

group.add_option("-Y", "--Yap", dest="mode", default="16S",
                 help="""Which Pipeline: 16S ITS [%default]""", metavar="#") 

group.add_option("-D", "--dynamic", dest="dynamic", action = "store_true", default=False,
                 help="""If specified, alignment will be scanned for primer locations and trimmed accordingly. Otherwise a database of known primers and trimming points will be used. [%default]""", metavar="#") 

group.add_option("-d", "--thresh", dest="dynthresh", default=0.1, type="float",
                 help="""in conjunction with -D, otherwise this is ignored. This allows to specify how much of the alignment to keep using the per-base coverage. The [%default] value indicates that ends of the alignment are trimmed until a base has a coverage of [%default] * peak coverage.""", metavar="#") 

group.add_option("-a", "--annotations", dest="dir_anno", default="/usr/local/devel/ANNOTATION/sszpakow/ANNOTATION/",
                 help="directory that stores auxilliary files\n[%default]", metavar="annotations")
group.add_option("-S", "--SAMPLE", dest="sampletimes", default=0, type="int",
                 help="perform sub.sampling of all reads based on the number of reads in smallest group. if 0 - all reads are used. if 1 - the sampling will be performed once, if 2 or more, then 2 or more independent samplings are going to be performed.\n[%default]", metavar="#")                 
group.add_option("-m", "--minlen", dest="minlength", default=220, type="int",
                 help="what is the minimum length of reads to process\n[%default]", metavar="#")     

group.add_option("-g", "--mingroupsize", dest="mingroupmembers", default=100, type="int",
                 help="after demultiplexing, discard groups with fewer reads than #\n[%default]", metavar="#")
    
group.add_option("-q", "--quick", dest="quickmode", action = "store_true", default=False,
                 help="""If specified, only single, reference DB based chimera checking will be used. [%default]""", metavar="#") 
              
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
  
   
    
if not options.mode in ("16S", "ITS"):
    parser.print_help()
    sys.exit(2)    

### parameters specific to YAP incarnations

### 16S V1-V3    
if options.mode=="16S":
    ### file in the annotations directory that has reference sequences
    _referenceseq = "ecolis.fasta"
    ### which fasta ID use as the reference (if file has more than one)
    _referenceseqname = "e_coli2_genbank"
    ### mothur's compendium of ALIGNED 16S sequences
    _alignment = "silva.bacteria.fasta"
    ### mothur's curated version of RDP's curated train set and corresponding taxonomy
    _trainset = "trainset9_032012.pds.fasta"
    _taxonomy = "trainset9_032012.pds.tax"
    ### until automatic primer detection is implemented, these are coordinates of primers
    ### when aligned to the silva.bacteria.fasta (for in-silico PCR and subsequent primer trimming)
    #_trimstart = "1044"
    #_trimend = "13127"
    
### ITS NSI1 - NLB4 (barcoded)   
elif options.mode=="ITS":
    _referenceseq = "yeastITS.fasta"
    _referenceseqname = "AF293_reference"
    _alignment = "FungalITSseed.092012.1.aln.fasta"
    _trainset = "FungalITSdb.092012.1.fasta"
    _taxonomy = "FungalITSdb.092012.1.tax"
    #_trimstart = "1716"
    #_trimend = "2795"    

else:
    parser.print_help()
    sys.exit(2)
                   
validator = InfoValidator(options.fn_info)  
_trimstart , _trimend, _region = validator.getTrimpoints()    
                                                              
BOH = init(options.project, options.email)
BOH.toPrint("-----", "GLOBAL",  "We are in %s mode" % (options.mode)) 

if options.dynamic or _region == "unknown":
    BOH.toPrint("-----", "GLOBAL",  "experimental dynamic alignment trimming enabled")
    BOH.toPrint("-----", "GLOBAL",  "alignment will be trimmed using %s * peak coverage threshold" % (options.dynthresh))
    _trimstart = "0"
    _trimend = "0"
else:
    BOH.toPrint("-----", "GLOBAL",  "alignment trimming predefined: %s - %s" % (_trimstart, _trimend))


#############################
#######################
##### reference: 
inputs = {"fasta": ["%s/%s" % (options.dir_anno, _referenceseq)] }
REF = FileImport(inputs)
REF_1 = MakeNamesFile([REF])
REF_2 = MakeGroupsFile([REF], "ref")
REF_3 = MakeQualFile  ([REF], "40" )
##############################
#
supplementary = list()
READY, primerfile = demultiplex()

#######################
##### primers: 
inputs = {"fasta": ["%s/%s" % (os.getcwd(), primerfile)] }
PRI = FileImport(inputs)
PRI_1 = MakeNamesFile([PRI])
PRI_2 = MakeGroupsFile([PRI], "ref")
PRI_3 = MakeQualFile([PRI], "40")

OutputStep("1-DEMUX", "groupstats,fasta,group,name,list,pdf,svg,tiff,taxsummary,globalsummary,localsummary", [READY])

if options.sampletimes==0:
    tmp = finalize(READY)    
    y = R_rarefactions(dict(), dict(), tmp)
    z = R_OTUplots(dict(), dict(), tmp)
    supplementary.append(y)
    supplementary.append(z)
    OutputStep("6-ENTIRE", "groupstats,fasta,group,name,list,pdf,svg,tiff,taxsummary,globalsummary,localsummary,phylotax", [tmp])
    OutputStep("8-TBC", "phylotax,group,list,fasta", [tmp])
    
else:
    thefinalset = list()
    for k in xrange(0, options.sampletimes):
        args =     {
                    "force" : "fasta,name,group",
                    "persample": "T",
                    "iter": "%s" % (k)
                }            
        sampled = MothurStep("sub.sample", options.nodesize, dict(), args, [READY])
        tmp = finalize(sampled)    
        y = R_rarefactions(dict(), dict(), tmp)
        z = R_OTUplots(dict(), dict(), tmp)
        supplementary.append(y)
        supplementary.append(z)
        OutputStep("SAMPLED_%s" % (k), "groupstats,fasta,group,name,list,pdf,svg,tiff,taxsummary,globalsummary,localsummary", [tmp])
        thefinalset.append(tmp)
    
OutputStep("7-SUPP_PLOTS", "tre,pdf,svg,tiff,r_nseqs,rarefaction,r_simpson,r_invsimpson,r_chao,r_shannon,r_shannoneven,r_coverage", supplementary)
    
    
###########################################################################    
##  
##################################################
###        Finish
##################################################
