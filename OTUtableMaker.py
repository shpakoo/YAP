#################################################
##     OTUtableMaker.py
#################################################
import sys, argparse
from collections import Counter

_author="Sebastian Szpakowski"
_date="Mar 1, 2013"
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

_arguments = argparse.ArgumentParser(description='A new program', epilog='if all else fails: shpakoo@gmail.com')
_arguments.add_argument("-l", "--fn_list", type=str, help="mothur list file")
_arguments.add_argument("-g", "--fn_group", type=str, help="mothur groups file")

args = _arguments.parse_args()



#################################################
##        Begin
##

seqid2group = dict()

for seqid, group in GeneralPurposeParser(args.fn_group):
        seqid2group[seqid] = group

groups=list(set(seqid2group.values()))

dists = dict()
for line in GeneralPurposeParser(args.fn_list):
    if  int(line[1]) == len(line[2:]):
        dists[line[0]] = line[2:]
    else:
        len(line [2:])
    
for dist in dists.keys():
    print (dist)
    otpt = open("dist_{0}.otutable".format(dist), "w")
    otpt.write("otuid\t{0}\n".format("\t".join(groups)))
    for OTUs in dists[dist]:

        seqids = OTUs.split(",")
        line = "{0}".format(seqids[0])
        cnt = Counter([seqid2group[x] for x in seqids])
        
        for g in groups:
            line = "{0}\t{1}".format(line, cnt[g])
        otpt.write("{0}\n".format(line))
    otpt.close()  
        


 
#################################################
##        Finish
#################################################