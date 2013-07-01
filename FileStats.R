########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2013 Sebastian Szpakowski
########################################################################################

library(reshape2)

# list directories
files = dir("./", pattern=glob2rx("*.fasta") , recursive=TRUE) 
files = append(files, dir("./", pattern=glob2rx("*.mate1") , recursive=TRUE))
files = append(files, dir("./", pattern=glob2rx("*.mate2") , recursive=TRUE))

x = data.frame(file = files, stringsAsFactors=FALSE)
x$step = unlist(lapply(x$file, function(x)
					{ 
						
						tmp = unlist(strsplit(x, "/")) 
						tmp = unlist(strsplit(tmp[1],"_"))
						tmp = tmp[-1]
						tmp = tmp[-length(tmp)]
						tmp = paste(tmp, sep="_", collapse="_")
						tmp						
					} ))

x$filename = unlist(lapply(x$file, function(x)
					{ 	
						
							tmp = unlist(strsplit(x, "/")) 
							tmp = tmp[length(tmp)]
							return(tmp)
							
					} ))
	
x$filetype = unlist(lapply(x$filename, function(x)
					{ 		
						
						if (regexpr( "fasta.mate1",x)>-1)
						{
							return("fasta (mate1)")
						}
						else if (regexpr( "fasta.mate2",x)>-1)
						{
							return("fasta (mate2)")
						}
						else if (regexpr( "contig",x)>-1)
						{
							return("fasta (contigs)")
						}
						else if (regexpr("singletons",x )>-1)
						{
							return("fasta (singletons)")
						}
						else if (regexpr("assembled",x)>-1)
						{
							return("fasta (assembled)")
						}
						else if (regexpr("fasta",x)>-1)
						{
							return("fasta")
						}
						
						else if (regexpr("mate1",x)>-1)
						{
							return("fastq (mate1)")
						}
						else if (regexpr("mate2",x)>-1)
						{
							return("fastq (mate2)")
						}
						
						else
						{
							tmp = unlist(strsplit(x, "\\.")) 
							tmp = tmp[length(tmp)]
							tmp 
						}	
						
					} ))	
	
x$id = unlist(lapply(x$filename, function(x)
					{ 		
						tmp = unlist(strsplit(x, "\\.")) 
						tmp = tmp[-length(tmp)]
						tmp = tmp[-1]
						tmp = tmp[1]
						tmp					
						
					} ))

x$seqs = 	unlist(lapply(x$file, function(x)
				{ 		
					
					if (regexpr("fasta", x)>-1)
					{
						tmp = system(paste("grep -c \"^>\" ", x, sep="" ), intern=TRUE, ignore.stderr=TRUE)				
						tmp = as.numeric(tmp)	
					}
					else 
					{
						tmp = system(paste("wc -l ", x, sep="" ), intern=TRUE, ignore.stderr=TRUE)				
						tmp = unlist(strsplit(tmp, " "))
						
						tmp = as.numeric(tmp[1])
						if (regexpr("mate", x)>-1)
						{
							tmp = tmp/4
						}
					}
					return(tmp)
					
				} ))

#x$wc = ifelse(x$filetype=="mate1", x$wc/4, x$wc)
#x$wc = ifelse(x$filetype=="mate2", x$wc/4, x$wc)


x$idtype = 	apply(x[,c("id", "filetype")], 1, function(x)
		{
			x = paste(x, sep="\t", collapse="\t")	
			x
		})


y = aggregate(x$seqs, x[, c("idtype", "step")], sum)

names(y)[length(names(y))]=c("sum")
z = dcast(y, idtype ~ step)

names(z)[1] = "id\ttype"
z[is.na(z)]=""
write.table(z, "summary.dir.tab.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )



















