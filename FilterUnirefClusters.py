#################################################
##     /Volumes/JCVI/home/sszpakow/DEVEL/YAP/scripts/Sources/FilterUnirefClusters.py
#################################################
import sys, argparse
from collections import defaultdict

_author="Sebastian Szpakowski"
_date="Jul 15, 2013"
_version="V1.00"
_description= """given a list of Uniref clusters, 
        find all members of the clusters and output ones that are of origin of interest"""
_epilog ="if all else fails: shpakoo@gmail.com"

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
    def __str__(self):
        return ("reading file: %s" %self.filename) 
        
    #################################################
    ### Iterator over input file.
    ### every line is converted into a dictionary with variables referred to by their 
    ### header name
class GeneralPurposeParser:
    def    __init__(self, filename, skip=0, sep="\t"):
        self.filename = filename
        self.fp = open(self.filename, "rU")    
        self.sep = sep
        
        self.skip=skip
        
        self.linecounter = 0
        self.currline=""
        
    def __iter__(self):
        return (self)
    
    def next(self):
        for currline in self.fp:
            currline = currline.strip().split(self.sep)
            self.currline = currline
            self.linecounter = self.linecounter + 1
            return(currline)            
        raise StopIteration
                    
    def    __str__(self):
        return "%s [%s]\n\t%s" % (self.filename, self.linecounter, self.currline)

class UnirefCluster:
    def __init__(self, clusterid):
        self.clusterid = clusterid
        self.proteinids = list()
        self.descriptors = dict()
        self.sequences = dict()
        self.lengths = dict()
        
    def add(self, id, desc, seq, term ):
        #print self.clusterid, id, desc, term, seq
        if desc.lower().find(term.lower())>-1:
            self.proteinids.append(id)
            self.descriptors[id] = desc
            self.sequences[id] = seq
            self.lengths[len(seq)] = id
            
    def getRepresentative(self):
        if len(self.proteinids)>0:
            
            tmp = self.lengths.keys()
            print "\t", tmp,
            tmp.sort()
            tmp.reverse()
            print "->", tmp
            id = self.lengths[tmp[0]]
            tmp = ">%s %s [%s]\n%s\n" % (self.clusterid, self.descriptors[id].split("|")[2],id, self.sequences[id] )
            return (tmp)
        else:
            return None
     
#################################################
##        Functions
##
    
#################################################
##        Arguments
##

_arguments = argparse.ArgumentParser(description=_description, epilog=_epilog)

_arguments.add_argument('-c', '--clusters', 
                            metavar="uniref50", 
                            type=str, 
                            dest = "fn_clusterlist")

_arguments.add_argument('-f', '--fasta', 
                            metavar='N', 
                            type=str, 
                            dest = "fn_fasta",
                            help="complete UNIPROT database to output sequences later")

_arguments.add_argument('-m', '--mapping', 
                            metavar='N', 
                            type=str, 
                            dest = "fn_mapping", 
                            help="mapping of protein IDs to clusters")

_arguments.add_argument('-k', '--keep', 
                            metavar='Homo sapiens', 
                            type=str, 
                            dest = "keep", 
                            help="Which taxon to keep from the clusters")


args = _arguments.parse_args()


#################################################
##        Begin
##

tokeep = defaultdict(str)
clusters = defaultdict(set)
sequences = defaultdict(str)
descriptors = defaultdict(str)

unirefversion = args.fn_clusterlist.split("/")[-1].split(".")[0]
if unirefversion.lower()=="uniref100":
    uni = "Uniref100"
    offset = 8
elif unirefversion.lower()=="uniref90":
    uni = "Uniref90"
    offset = 9
elif unirefversion.lower()=="uniref50":
    uni = "Uniref50"
    offset = 10
else:
    print "Not sure which Uniref this [%s] is..." % (args.fn_clusterlist) 
    sys.exit(1)
    
print "%s detected" % (uni)

otptfiltered = open("%s.filtered.%s.fasta" % (uni, args.keep.replace(" ", "_")), "w")
otptmissing  = open("%s.missing.%s.list" % (uni, args.keep.replace(" ", "_")) , "w")
    
print ("Caching useful clusters")
for head, seq in FastaParser(args.fn_clusterlist):
    head = head.strip()
    desc = " ".join(head.split(" ")[1:]).lower()
    if desc.find("[fragment]") >-1 or desc.find("(fragment)") >-1 or desc.find("[fragments]")>-1 or desc.find("(fragments)")>-1:
       pass
    else:
        tokeep[head.split(" ")[0]]=desc
        #print (desc)

    #if len(tokeep)==100:
    #    break    
          
print ("Found: %s" % len(tokeep))

print ("Caching cluster mappings...")
# read in cluster mapping
for line in GeneralPurposeParser(args.fn_mapping):
    id, cluster = line[0], line[offset]
    #if cluster in tokeep.keys():
    #print "keeping %s - %s" % (id, cluster)
    clusters[cluster].add(id)

        #if len(clusters.keys()) == len(tokeep.keys()):
        #    break

print (len(clusters.keys()))

print ("Caching sequences...")
for head, seq in FastaParser(args.fn_fasta):
    id = head.strip().split("|")[1]
    sequences[id] = seq
    descriptors[id] = head

print ("Retrieving Clusters...")
for cluster in clusters.keys():
    
    if tokeep[cluster]!="":
        print "_____"
        print cluster, tokeep[cluster]
        tmp = UnirefCluster(cluster)
        
        for proteinid in clusters[cluster]:
            tmp.add(proteinid, descriptors[proteinid], sequences[proteinid], args.keep)
        
        x = tmp.getRepresentative()
        if x is None:
            print "\tNo %s found in %s" % (args.keep, cluster) 
            otptmissing.write("%s\t%s\n" % (cluster, tokeep[cluster]))
        else:
            otptfiltered.write("%s" % x )
            print "\t%s" % (x.strip().split("\n")[0])
    else:
        print "Skipping fragment %s" % (cluster)
        
#################################################
##        Finish
#################################################