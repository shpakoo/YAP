########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2013 Sebastian Szpakowski
########################################################################################

#################################################
##    A pipeline for miseq data 
##    OTUs (certain regions of 16S and ITS supported)
##    This is for demultiplexed MiSeq data
#################################################
import sys, os.path
from optparse import OptionParser, OptionGroup
from StepsLibrary import *
from StepsLibrary_EXP import *
from collections import defaultdict
from Queue import Queue

_author="Sebastian Szpakowski"
_date="2013/04/01"
_version="Version 5"

#################################################
##        Classes
##
 
class InfoValidator:
    def __init__(self,filename):
        self.filename = filename
        self.info = GeneralPurposeParser(filename, sep=",")
        self.URI = "http://confluence/display/~sszpakow/YAP"
        self.epilogue = "\n***\tPlease correct before continuing...\n***\t{0}\n".format(self.URI) 
        self.header = ""
        
        self.tech = ""
        
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
                has =  ",".join (self.header)
                needed454 = "path,file,barcode,forward,reverse,use,both,SampleID"
                neededMiSeq = "path,file1,file2,forward,reverse,SampleID"
                
                if has.lower().startswith( needed454.lower()) :
                    self.tech = "454"
                elif  has.lower().startswith( neededMiSeq.lower()) :
                    self.tech = "MiSeq"
                else:
                    self.error( "Your template's header is incorrect or missing:\nhas :\t{0}\nneed (454):\t{1}\n\t(illumina)\t{2}".format(has, needed454, neededMiSeq), 101)
                
                if not ("SampleID" in self.header):    
                    self.error( "Your template has\n\t'{0}' instead of \n\t'SampleID' in the column's header.".format(self.header[7]), 102)
                    
            else:
                files.add("{0}/{1}".format(line[0], line[1].strip()))   
                if self.tech == "454":
                    barcodes.add(line[2])   
                    primersF.add(line[3])    
                    primersR.add(line[4])
                    sampleIDs.add(line[7])
                elif self.tech == "MiSeq":
                    if line[2].strip() != "":
                        files.add("{0}/{1}".format(line[0], line[2].strip())) 
                    primersF.add(line[3])    
                    primersR.add(line[4])
                    sampleIDs.add(line[5])
                        
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
    def getTech(self):
        return self.tech
     
class InfoParserMiSeq:
    def __init__(self, filename):
        self.filename = filename
        self.info = GeneralPurposeParser(filename, sep=",", skip=1)
        self.store = list()
        self.IDs = defaultdict(str)
        self.primers = set()
        self.forward = ""
        self.reverse = ""
        
        for line in self.info:
            path = line[0]
            file1 = line[1]
            file2 = line[2]
            forward = line[3]
            reverse = line[4]
            
            if path.endswith("/"):
                path = path[:-1]
            
            path1 = "%s/%s" % (path, file1)
            path2 = "%s/%s" % (path, file2)
                       
            if file2=="":
                self.store.append([path1])
                self.IDs[path1] = line[5]
            else:
                self.store.append([path1, path2])
                self.IDs[path1] = line[5]
                self.IDs[path2] = line[5]

            if reverse =="" or forward =="":
                print "%s: please provide both primers for file(s):'%s' " % (x, ",".join(file1, file2))
                sys.exit(11) 
            else:    
                self.primers.add(">_primer_F\n%s\n" % (forward))
                self.primers.add(">_primer_F_rc\n%s\n" % (revComp(forward)))
                self.primers.add(">_primer_R\n%s\n" % (reverse))
                self.primers.add(">_primer_R_rc\n%s\n" % (revComp(reverse)))
                self.forward = forward
                self.reverse = reverse

    def getFiles(self):
        return (self.store)       
    def getSampleID(self, file):
        return self.IDs[file]
       
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
    
def preprocess():
    forprocessing = InfoParserMiSeq(options.fn_info)
    PREPROCESS = list()

    for files in forprocessing.getFiles():
        INS = {}
        
        
        if len(files) == 2:
            M1 = files[0]
            M2 = files[1]
            sampleid = forprocessing.getSampleID(M1)
            INS = {"mate1": ["%s~%s" % (M1, sampleid)], "mate2": ["%s~%s" % (M2, sampleid)]}
        else:
            M1 = files[0]
            sampleid = forprocessing.getSampleID(M1)
            INS = {"fastq": ["%s~%s" % (M1, sampleid)]}
          
        #### import files    
        if options.head == 0:
            x = FileImport(INS)
        else:
            x = FileMiniImport(INS, {"lines": options.head})    

        #### determine the encoding of fastQ
        Q = getQ(M1)
        
       
        if Q == "":
            print (Q)
            print "Q issues"
            print files
            sys.exit(1)
        ### generate quality information:
        ARGS = {
             "-h": options.minqual,
             "-m": "",
             "-v": ""
            }
        qc = SQA(ARGS, [x])
        supplementary.append(qc) 
        
        ### split into smaller files for parallelization
        ### 100,000 sequences (x4 since fastq)
        ARGS = {
                    "types": "mate1,mate2,fastq",
                    "chunk":  "400000"
        }
        P0 = FileSplit(ARGS, [x]) 
            
        #### trim fastQ files
        ARGS = {
         "-h": options.minqual,
        }
        P1 = SQAtrim(ARGS, [P0])
        
        #### overlap mates if available
        if len(files)==2:
            ARGS = {
             "-M": "200",
             "-p": Q,
             "-r": "250"
            }
            P2 = Flash({}, ARGS, [P1])
        else:    
            P2 = P1
               
        #### convert fastq to fasta
        ARGS = { 
                "-Q": Q
        }
        P3 = fastq2fasta(dict(), ARGS, [P2])
        
        #### use fuzznuc to find cut primer sequences
        ARGS = {
                "-f": forprocessing.forward,
                "-r": forprocessing.reverse,
                "-m": "1"
        }
        P4 = PrimerClipper ( {}, ARGS, [P3])

        ### make fastA headers less problematic
        P5 = FastaHeadHash({}, {}, [P4])
        P6 = FileMerger("fasta", [P5])
        P7 = MakeGroupsFile([P6], sampleid)
        P8 = MakeNamesFile([P6])
        
        PREPROCESS.extend([P6,P7,P8])
       
    
    A1 = FileMerger("fasta,group,name", PREPROCESS)
    
    
    args = {"mingroupmembers": options.mingroupmembers, 
            "report": "failing"}
    A2 = GroupRetriever(args, [A1])
    
    args = {
            "force" : "fasta,name,group",
            "find": "groups" 
            }   
    A3 = MothurStep("remove.groups", options.nodesize, dict(), args, [A2])   

    return (A3)

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
    OutputStep("2-NOISY", "groupstats,fasta,group,name,list,svg,pdf,tiff,taxsummary,globalsummary,localsummary", clean2a)

    ###################### CDHIT-454
    #### unique and de-noise
    args = {}
    
    ### strictly unique collapsing
    if options.strictlevel==1:
        args=         { 
                        "c" : "1.0",
                        "b" : "8",
                        "aS": "1.0",
                        "g" : "1",
                        "M"    : "50000",
                        "T" : "%s" % (options.nodesize)
                       
                    }
    
    ### aggressive de-noising:
    elif options.strictlevel==2:
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
    CD_3 = FileMerger("fasta,name,group,qfile", [CD_2, REF_1, REF_2, REF_3])
    
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
    args = {}
    if options.dynamic:
        args = { "minlength" : "50" ,
                 "force": "fasta,name,group"}
    else:
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
                "d" : "0",
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
    
    if options.quickmode:
        return ([s23, s24, s25aa, s25bb, s26a, s27a, s28, s29])
    else:
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

group.add_option("-d", "--thresh", dest="dynthresh", default=0.75, type="float",
                 help="""in conjunction with -D, otherwise this is ignored. This allows to specify how much of the alignment to keep using the per-base coverage. The [%default] value indicates that ends of the alignment are trimmed until a base has a coverage of [%default] * peak coverage.""", metavar="#") 

group.add_option("-a", "--annotations", dest="dir_anno", default="/usr/local/devel/ANNOTATION/sszpakow/ANNOTATION/",
                 help="directory that stores auxilliary files\n[%default]", metavar="annotations")
group.add_option("-S", "--SAMPLE", dest="sampletimes", default=0, type="int",
                 help="perform sub.sampling of all reads based on the number of reads in smallest group. if 0 - all reads are used. if 1 - the sampling will be performed once, if 2 or more, then 2 or more independent samplings are going to be performed.\n[%default]", metavar="#")                 
group.add_option("-m", "--minlen", dest="minlength", default=200, type="int",
                 help="what is the minimum length of reads to process\n[%default]", metavar="#")     

group.add_option("-g", "--mingroupsize", dest="mingroupmembers", default=100, type="int",
                 help="after demultiplexing, discard groups with fewer reads than #\n[%default]", metavar="#")

group.add_option("-Q", "--minqual", dest="minqual", default=30, type="int",
                 help="Keep stretches of reads this good or better #\n[%default]", metavar="#")

group.add_option("-q", "--quick", dest="quickmode", action = "store_true", default=False,
                 help="""If specified, only single, reference DB based chimera checking will be used. [%default]""", metavar="#") 
 
parser.add_option("-H", "--head", dest="head", default=0, type="int",
                 help="For dry runs, import only # of lines from the input files")

              
group.add_option("-x", "--strict", dest="strictlevel", default=2, type="int",
                 help="""how strict to be at pre-clustering: 
1 very strict, conservative denoising (precluster identical sequences) 
2 less strict, aggresive denoising (precluster using 98% similarity)
[%default]""", metavar="#")                 

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
_tech = validator.getTech()
                                                              
BOH = init(options.project, options.email)
BOH.toPrint("-----", "GLOBAL",  "We are in %s mode" % (options.mode)) 
BOH.toPrint("-----", "GLOBAL",  "We will be processing %s data" % (_tech)) 

if options.dynamic or _region == "unknown":
    BOH.toPrint("-----", "GLOBAL",  "Dynamic alignment trimming enabled")
    BOH.toPrint("-----", "GLOBAL",  "Alignment will be trimmed using %s * peak coverage threshold" % (options.dynthresh))
    _trimstart = "0"
    _trimend = "0"
else:
    BOH.toPrint("-----", "GLOBAL",  "Alignment trimming predefined: %s - %s" % (_trimstart, _trimend))

#############################
#######################
##### reference: 
inputs = {"fasta": ["%s/%s" % (options.dir_anno, _referenceseq)] }
REF = FileImport(inputs)
REF_1 = MakeNamesFile([REF])
REF_2 = MakeGroupsFile([REF], "ref")
REF_3 = MakeQualFile  ([REF], "40" )
##############################

supplementary = list()
READY = preprocess()
OutputStep("1-PREPROCESS", "groupstats,fasta,group,name,list,pdf,svg,tiff,taxsummary,globalsummary,localsummary", READY)

if options.sampletimes==0:
    tmp = finalize(READY)    
    y = R_rarefactions(dict(), dict(), tmp)
    z = R_OTUplots(dict(), dict(), tmp)
    supplementary.append(y)
    supplementary.append(z)
    OutputStep("6-ENTIRE", "groupstats,fasta,group,name,list,pdf,svg,tiff,taxsummary,globalsummary,localsummary,phylotax", [tmp])
    OutputStep("8-TBC", "phylotax,group,list,fasta", [tmp])
    
#else:
#    thefinalset = list()
#    for k in xrange(0, options.sampletimes):
#        args =     {
#                    "force" : "fasta,name,group",
#                    "persample": "T",
#                    "iter": "%s" % (k)
#                }            
#        sampled = MothurStep("sub.sample", options.nodesize, dict(), args, [READY])
#        tmp = finalize(sampled)    
#        y = R_rarefactions(dict(), dict(), tmp)
#        z = R_OTUplots(dict(), dict(), tmp)
#        supplementary.append(y)
#        supplementary.append(z)
#        OutputStep("SAMPLED_%s" % (k), "groupstats,fasta,group,name,list,pdf,svg,tiff,taxsummary,globalsummary,localsummary", [tmp])
#        thefinalset.append(tmp)
#    
OutputStep("7-SUPP_PLOTS", "tre,pdf,png,svg,tiff,r_nseqs,rarefaction,r_simpson,r_invsimpson,r_chao,r_shannon,r_shannoneven,r_coverage", supplementary)
    
    
###########################################################################    
##  
##################################################
###        Finish
##################################################
