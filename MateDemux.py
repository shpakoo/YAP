#################################################
##     /Volumes/JCVI/home/sszpakow/DEVEL/YAP/scripts/Sources/MateDemux.py
#################################################
import sys, argparse
from collections import defaultdict

_author="Sebastian Szpakowski"
_date="May 23, 2013"
_version="V1.00"

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

class FastQParser:
    def __init__(self, file, skip=0, sep="\t"):
        self.filename = file
        self.fp = open(self.filename, "rU") 
        
        self.linecounter = 0
        self.currline=""

    def __iter__(self):
        return (self)
    
    def next(self):
        buffer = list()
        id, seq, qual = "","", ""
        
        for currline in self.fp:
            buffer.append(currline.strip())
            self.currline = currline
            self.linecounter = self.linecounter + 1
            if len(buffer)==4:
                return(buffer[0][1:], buffer[1], buffer[3])            
        raise StopIteration
                    
    def __str__(self):
        return "%s [%s]\n\t%s" % (self.filename, self.linecounter, self.currline)
              
#################################################
##        Functions
##

#################################################
##        Arguments
##

_arguments = argparse.ArgumentParser(description='remove barcode, create a groups file and scrap file', epilog='if all else fails: shpakoo@gmail.com')
_arguments.add_argument('-f', '--fasta', metavar='fasta', type=str, dest = "fn_file")
_arguments.add_argument('--format', metavar='fastq', type=str, dest = "format", default="fastq")
_arguments.add_argument('-m', '--mapping', metavar='mapping', type=str, dest = "fn_mapping")

args = _arguments.parse_args()


#################################################
##        Begin
##
barcodes = dict()
files = dict()
barlength = set()
bardups = defaultdict(list)
bartotal = defaultdict(int)

for sample, barcode in GeneralPurposeParser(args.fn_mapping, sep=","):
    if sample != "SampleID":
        barcodes[barcode]  = sample
        barlength.add(len(barcode))
        newfilename = "%s.%s.%s"  % (".".join(args.fn_file.split("/")[-1].split(".")[:-1]), sample, args.format)
        tmp = open(newfilename,  "w")
        files [barcode] = tmp
    
scrapfile = "%s.scrap.%s.%s"  % (".".join(args.fn_file.split("/")[-1].split(".")[:-1]), args.format, "fasta")
groupsfile = "%s.%s"  % (".".join(args.fn_file.split("/")[-1].split(".")[1:]), "groups")    
    
files["scrap"] = open(scrapfile, "w")
files["groups"] = open(groupsfile, "w")
  
count = 0
 
if len(barlength)>1:
    print "barcodes of different lengths"
    sys.exit(1) 
else:
    barlength = list(barlength)[0]
   
for id, seq, qual in FastQParser(args.fn_file):
    bar = seq[0:barlength]
    countbar = seq.count(bar)
    if files.has_key(bar):
        sample = barcodes[bar]
        bartotal[bar]+=1
        if countbar >1:
            files["scrap"].write(">%s|%s|b#%s\n%s\n" % (id, sample, countbar, seq))
            bardups[bar].append(countbar)
        else:
            files[bar].write("@%s\n%s\n+\n%s\n" % (id, seq[barlength:], qual[barlength:]))   
            
    else:
        files["scrap"].write(">%s|b\n%s\n" % (id, seq))         
    
    
for key in barcodes.keys():
    if len(bardups[key]) >0:
        print "%s\t%s\t%s\t%s\t%s\t%s" % ( key,   barcodes[key], bartotal[key], len(bardups[key]),  float(sum(bardups[key]))/len(bardups[key]), ",".join(map(str, bardups[key])))
    else:
        print "%s\t%s\t%s\t%s\t%s\t%s" % ( key,   barcodes[key], bartotal[key], len(bardups[key]),  0, ",".join(map(str, bardups[key])))
  
    #print key, barcodes[key], bartotal[key], len(bardups[key]),  float(sum(bardups[key]))/len(bardups[key])

    
#################################################
##        Finish
#################################################