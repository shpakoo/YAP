########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2013 Sebastian Szpakowski
########################################################################################

makePlots = function(tmp, dataset, f)
{
	for ( stat in c("length", "meanq", "medq", "R50.meanq", "cov")) 
		{
			print (stat)
			if ( max(dataset[,stat])>10 )
			{
				br = seq(0-1.5, max(dataset[,stat])+1.5, by=1)
			}
			else
			{
				br = 100
			}
			#print (br)
			hist (tmp[,stat], br = br, main=paste(f, stat, sep="\n"), xlim=c(0, max(dataset[,stat]) ) ) 
		}
		
		x=sample(1:nrow(tmp), 0.1* nrow(tmp) )
		plot (tmp$length[x], tmp$meanq[x], pch=19, col="gray", cex=0.1, xlim=c(0, max(dataset[,"length"]) ), ylim=c(0, max(dataset[,"meanq"]) ))
		points(lowess(tmp$length, tmp$meanq, f=0.05) , type="l", lwd=6, col="gray50")
		abline (v= mean(tmp[,"length"]), col="green", lwd=3)
		abline (h= mean(tmp[,"meanq"]), col="blue", lwd=3)
		legend ("bottom", legend=c(round(mean(tmp[,"length"]),0), round(mean(tmp[,"meanq"]),0), nrow(tmp)), fill=c("green", "blue", "black"), cex=1.5 , ncol=3)
}


processFiles = function(filelist, prefix, parentchild=list(), raster=TRUE)
{
	dataset=data.frame()
	for (filename in filelist)
	{
		cat("Loading ", filename,"...", sep="")
		tmp = read.table(filename, header=TRUE, as.is=TRUE)

		if (nrow(dataset)==0)
		{
			dataset=tmp
		}
		else
		{
			dataset=rbind(dataset, tmp)
		}
		cat("\n")
	}

	dataset$cov = dataset$sdq / dataset$meanq
	files = unique(dataset$file)
	panels = length(files)+length(parentchild)
	if (raster)
	{
		png(paste(prefix, ".stats.png", sep=""), width=26, height=5*(panels), unit="in", res=400)
	}
	else
	{
		pdf(paste(prefix, ".stats.pdf", sep=""), width=24, height=5*(panels), paper="special")	
	}
	par(mfrow= c(panels,6))
	for (f in files)
	{
		tmp= dataset[dataset$file==f, ]
		makePlots(tmp, dataset, f)
		
		for (k in names(parentchild))
		{
			if (k == f)
			{
				cat("Loading ", parentchild[[k]], " for ", k, "...", sep="")
				subsets = read.table(parentchild[[k]], header=TRUE, as.is=TRUE)
				tmp = tmp[tmp$read %in% subsets$QueryName,]
				makePlots(tmp, dataset, parentchild[[k]])
				cat("\n")
			}
		}
		
	}
	dev.off()
	invisible(dataset)
}

# lista=list()
# lista["T1.trim.qual"] = "T1.trim.unique.fragclust.align.report"
# filelist=c("T1.trim.qual.stats.tab.txt")
# processFiles(filelist, "test", lista)

