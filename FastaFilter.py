#################################################
##     /Volumes/JCVI/home/sszpakow/DEVEL/YAP/scripts/Sources/FastaFilter.py
#################################################
import sys, argparse

_author="Sebastian Szpakowski"
_date="May 20, 2013"
_version="V1.00"

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
#################################################
##        Arguments
##

_arguments = argparse.ArgumentParser(description='Given a list of IDs, Filter a fastA file.', epilog='if all else fails: shpakoo@gmail.com')

_arguments.add_argument('-f', '--fasta', metavar='N', type=str, dest = "fn_fasta")
_arguments.add_argument('-l', '--list', metavar='N', type=str, dest = "fn_list")
_arguments.add_argument('-m', '--mode', metavar='N', type=str, dest = "mode", help="either keep or filter to remove the contents of the list ")

args = _arguments.parse_args()


#################################################
##        Begin
##

if args.mode=="filter":
    keep=True
else:
    keep=False    

newfilename = "%s.filter.fasta" % ("".join(args.fn_fasta.split("/")[-1].split(".")[:-1]))
print newfilename
otpt = open(newfilename, "w")

filters = list()
counter =0 

for line in loadLines(args.fn_list):
    filters.append(line.strip())
   
for head, seq in FastaParser(args.fn_fasta):
    tmp = head.split("|")[0]

    if keep:
        if tmp in filters and keep:
            otpt.write(">%s\n%s\n" % (head, seq))
            counter+= 1
    else:        
         if not (tmp in filters) :
            otpt.write(">%s\n%s\n" % (head, seq))
            counter+= 1
            
        
 
otpt.close()
print "%s sequences written" % (counter)        
        
        


#################################################
##        Finish
#################################################