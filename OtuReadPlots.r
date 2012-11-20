########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2012 J.Craig Venter Institute.
########################################################################################

library(ggplot2)
library(grid)
library(scales)

multiplot <- function(..., plotlist=NULL, cols=1, layout=NULL)
{
	require(grid)
	
	# Make a list from the ... arguments and plotlist
	plots = c(list(...), plotlist)
	
	numPlots = length(plots)
	print (numPlots)
	
	# If layout is NULL, then use 'cols' to determine layout
	if (is.null(layout)) 
	{
		layout = grid.layout(ceiling(numPlots/cols), cols)
	}
	
	# create a panel that will be used to map ids later
	
	panel = matrix(1:(layout$nrow * layout$ncol), nrow = layout$nrow, ncol = layout$ncol, byrow=TRUE )
	
	if (numPlots==1) 
	{
		print(plots[[1]])
	} 
	else 
	{
		# Set up the page
		grid.newpage()
		pushViewport(viewport(layout = layout))
		
		for (i in 1:numPlots) 
		{
			matching <- as.data.frame(which(panel == i, arr.ind = TRUE))
			
			print ( plots[[i]], vp = viewport(	layout.pos.row = matching$row, 
							layout.pos.col = matching$col)
			)
		}
	}
}


makePlot = function(filename = "read.groups.tab.txt", minreads = 0, mingroups = 0, nlevels = 500, mode="OTU", xlim=c(0,0), ylim=c(0,0))
{
	if (mode=="OTU")
	{
		x = read.table(filename, as.is=T, header=F, sep="\t")
		x=x[-1,]
		names(x) = c("reads", "groups")
		
		ylab = "number of samples showing an OTU"
		xlab = "OTU size (i.e. number of reads in OTU, log10 scale)"
		
		histlab="OTU size"
		
		groups = unique(x$groups)
		labs = groups
		
		alldata = x
		subsetdata =  x[x$reads>minreads & x$groups>mingroups,]
		
		v = NA
		h = NA
		
		forhist = subsetdata$reads
		forhist[forhist>100]=100
		
	}
	else if ( mode=="coverage")
	{
		x = read.table(filename, as.is=T, header=F, sep="", skip=34)
		#print (dim(x))
		#x=x[-1,]
		
		###  Contig     Sites     Reads   Coverage
		names(x) = c("x1", "reads", "x2", "groups")
		
		ylab = "coverage"
		xlab = "Contig Size (log10 scale)"
		
		histlab="Contig Coverage"
				
		x$groups = log10(x$groups)
		
		#groups = round(seq(0,ceiling(max(x$groups)), length.out=20),2)
		groups = log10(c(1,10,100,1000,10000))
		labs = round(10^groups,0)
		
		alldata = x
		subsetdata =  x[x$reads>minreads & x$groups>mingroups,]
		
		forhist = 10^subsetdata$groups
		
		h = mean(subsetdata$reads)
		v = sum(subsetdata$x2) / sum(subsetdata$reads) * 100
		
		xlim = log10(c(100,100000))
		ylim = log10 (c(0.1,10000))
		
		forhist[forhist>100]=100
		print (v)
		
	}
	else
	{
		print ("unknown plotting mode.")
		return (0)
	}
	
	alldata = alldata[is.finite(alldata$groups),]
	subsetdata = subsetdata[is.finite(subsetdata$groups),]
	
	### setup  canvas using all data

	nf <- layout(matrix(c(1,2,3,0), 2, 2, byrow=TRUE), respect=FALSE, width=c(0.8,0.2))
	
	if (xlim[1]==xlim[2])
	{
		ylim = c(min(alldata$groups)-0.5,max(alldata$groups)+0.5)
		plot (log10(alldata$reads), alldata$groups, type= "n", ylab = ylab, xlab = xlab, ylim=ylim, axes=FALSE, main=filename)

	}
	else
	{
		plot (log10(alldata$reads), alldata$groups, type= "n", ylab = ylab, xlab = xlab, ylim=ylim, xlim=xlim,  axes=FALSE, main=filename)
	}
	
	
	
	axis(2, at=groups, lab=labs, las=2)
	A = c(1:9)
	B =seq(0, log10(max(alldata$reads)), by=1)
	C =expand.grid(A, B)
	C$at = C[,1] * 10^(C[,2])
	C$lab = ifelse (C[,1]%%10==1, 10^C[,2], "")
	
	axis(1, at=log10(C$at), lab=rep("", nrow(C)) , col.ticks="gray", lwd.ticks=2 )
	C2 = C[C$lab!="", ]
	axis(1, at=log10(C2$at), lab=C2$lab , col.ticks="black", lwd.ticks= 4 )

	### use only subset data from now on:
	subsetdata$reads = log10(subsetdata$reads)
		
	# contour
	require(MASS)
	xykde = MASS::kde2d(subsetdata$reads, subsetdata$groups, lims = c( min(subsetdata$reads)-5, max(subsetdata$reads)+5, min(subsetdata$groups)-5, max(subsetdata$groups)+5) )
	
	if (sum(is.nan(xykde$z))>0)
	{
		cat("have to nudge groups...\n")
		b = subsetdata$groups + sample(c(-0.1, 0.1), length(subsetdata$groups), replace=T)
		xykde = MASS::kde2d(subsetdata$reads, b, lims = c( min(subsetdata$reads)-5, max(subsetdata$reads)+5, min(subsetdata$groups)-5, max(subsetdata$groups)+5) )
	}
	
	zlim = range(xykde$z, finite = TRUE)
	lev = seq(zlim[1], zlim[2], le = nlevels)
	col = rev(heat.colors(length(lev)))
	contour(xykde, add = TRUE, levels = lev,
		drawlabels = FALSE, col=col  )        
	
			
	# trendline
	trendline = lowess((subsetdata$reads), subsetdata$groups, f=0.1)
	points(trendline , type="l", lwd=5, col="gray30")
	
	# points
	points(subsetdata$reads, subsetdata$group, pch=19, col="gray60", cex=0.4) 
	
	# cutoff lines	if (mingroups>0)
	{
		abline(h=mingroups, col="red", lty=2, lwd=0.5)
	}
	if (minreads>0)
	{
		abline(v=log10(minreads), col="red", lty=2, lwd=0.5)
	}
	box(lwd=5)
	
	# legend
	plot (0,100, xlim=c(0,1), ylim = c(min(lev), max(lev)), xlab="", ylab="", axes=F, main="")
	rect(0, lev[-length(lev)], 1, lev[-1L], col = col, border=NA)  

	#print (lev)
	
	at = unique(round(seq(0, max(lev), length.out=10),2))
	lab = paste(round((1-(at/ceiling(max(at)))) * 100,0) , "%", sep="")
	axis(4, las=1, cex.axis=0.5)
	axis(2, las=1, at = at, lab=lab)
	box()
		
	
	# histogram
	hist ((forhist), col="gray", br=seq(0.5,100.5, by=1), axes =FALSE, xlab = histlab, main="")
	axis(2)
	at = c(1:10, seq(20,100, by=10) )
	labs = paste(at)
	labs [2:9]=""
	
	labs[length(labs)] = paste(labs[length(labs)], "+", sep="")
	axis(1, at = at, lab=labs, las=3)
	
	if (! is.na(v))
	{	
		#abline(v= v, lwd=5, col="green")
		#abline(h= h, lwd=5, col="blue")
		legend("top", c(
							paste( "Average coverage: ", round(v,2), sep=" "),
							paste( "Average length: ", round(h,2), sep=" ")
							), 
							
		fill=c("green", "blue"))
	}

	#print ( nrow(x[x$groups<=1 & x$reads<=0.1,]) )
	#print ( nrow(x))

	 
}


getData = function(files, mode="clccoverage")
{
	global = data.frame()
	for (f in files)
	{
		
		print (f)
		if (mode =="clccoverage")
		{
			tmp = read.table(f, as.is=T, header=F, sep="", skip=34)
			names(tmp) <- c("Contig", "Sites", "Reads", "Coverage")
			
		}
		else if (mode == "otu")
		{
			tmp = read.table(f, as.is=T, header=F, sep="\t")
			names(tmp) <- c("Reads", "Groups")
			
		}
		else
		{
			return (global)
		}
		
		a = strsplit( f, "\\.")[[1]][2]
		b = strsplit( a,"-")[[1]][1]
		cat(f, "->", a, "->", b, "\n")
		tmp$File = rep(b, nrow(tmp))
		
		if (nrow(global)==0)
		{
			global=tmp
		}
		else
		{
			global=rbind(tmp,global) 
		}
			
	}
	return (global)
}

makeGGplotOTU = function(global)
{
	print ("not there yet")
}
	
makeGGplotCOVERAGE = function(global)
{
	
	### how do files compare w/r to contig sizes
	plots = list()
	
	for (S in c("Reads", "Sites", "Coverage"))
	{
		x1 = (ggplot(global,  aes_string(y=S, x="File", color="File"))
					+ geom_boxplot() 
					+ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
							labels = trans_format("log10", math_format(10^.x)))
					+ coord_flip()
					+ theme(axis.ticks = element_blank(), axis.text.y = element_blank())
					
					)
		
		x2 = (ggplot(global,  aes_string(x=S, y="..count..",  fill="File", group="File" )) 
					+ stat_bin(aes(y=..count..), geom="area", position="stack", weight=2)
					+ scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
							labels = trans_format("log10", math_format(10^.x)))
					)
		
		x3 = (ggplot(global,  aes_string(x=S, y="..count..",  fill="File", group="File" )) 
					+ stat_bin(aes(y=..count..), geom="area", position="dodge", alpha=0.5, weight=2)
					+ scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
							labels = trans_format("log10", math_format(10^.x)))
					)
		
		x4 = (ggplot(global,  aes_string(x=S, y="..count..", fill="File", group="File" )) 
					+ geom_bar(position="fill")
					+ scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
							labels = trans_format("log10", math_format(10^.x)))
					)
	
		
	 	plots = append(list(x1, x2, x3, x4), plots)
		
	}
	y1 = (ggplot(global,  aes(y=Reads / Sites, x=File, color=File))
				+ geom_boxplot() 
				+ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
						labels = trans_format("log10", math_format(10^.x)))
				+ coord_flip()
				+ theme(axis.ticks = element_blank(), axis.text.y = element_blank())
				
				)
	y2 = (ggplot(global,  aes(x= Reads/Sites , y=..count..,  fill=File, group=File )) 
				+ stat_bin(aes(y=..count..), geom="area", position="stack", weight=2)
				+ scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
						labels = trans_format("log10", math_format(10^.x)))
				)
	
	y3 = (ggplot(global,  aes(x= Reads/Sites , y=..count..,  fill=File, group=File ))  
				+ stat_bin(aes(y=..count..), geom="area", position="dodge", alpha=0.5, weight=2)
				+ scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
						labels = trans_format("log10", math_format(10^.x)))
				)
	
	y4 = (ggplot(global,  aes(x= Reads/Sites , y=..count..,  fill=File, group=File ))  
				+ geom_bar(position="fill")
				+ scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
						labels = trans_format("log10", math_format(10^.x)))
				)
		
	plots = append(plots, list(y1, y2, y3, y4), )
	
	
	
	#return (plots)
	
#	B = (ggplot(global,  aes(File, Sites, color=File)) 
#				+ geom_boxplot() 
#				+ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
#						labels = trans_format("log10", math_format(10^.x)))
#				+ coord_flip()
#				+ theme(axis.ticks = element_blank(), axis.text.y = element_blank())
#				
#				)
#	C = (ggplot(global,  aes(File, Coverage, color=File)) 
#				+ geom_boxplot() 
#				+ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
#						labels = trans_format("log10", math_format(10^.x)))
#				+ coord_flip()
#				+ theme(axis.ticks = element_blank(), axis.text.y = element_blank())
#				
#				)
	
	
	png("AllAssemblyStats_GG.png", width=20, height = 12, units="in", res=250)
	multiplot(plotlist=plots, cols=4)
	dev.off()
	invisible(plots)

}

makeBatchOTU = function(filename) 
{
	pdf(paste(filename,"otusizes.pdf", sep="."), paper="special", width=11, height=10)
	makePlot(filename = filename, minreads=0, mingroups=0)
	makePlot(filename = filename, minreads=1, mingroups=0)
	makePlot(filename = filename, minreads=1, mingroups=1)
	makePlot(filename = filename, minreads=5, mingroups=1)
	dev.off()
}

makeBatchCoverage = function(filename)
{
	pdf(paste(filename,"coveragestats.pdf", sep="."), paper="special", width=11, height=10)
	makePlot(filename = filename, minreads=0, mingroups=0, mode="coverage", nlevels=50)
	makePlot(filename = filename, minreads=0, mingroups=log10(10), mode="coverage", nlevels=50)
	makePlot(filename = filename, minreads=1000, mingroups=log10(10), mode="coverage", nlevels=50)
	dev.off()
}


###################################################################
##### OTU stats
files = dir( pattern=glob2rx("*.otustats"))
for (f in files)
{
	print (f)
	makeBatchOTU (f)
}

if (length(files)>0)
{
	inputdata = getData(files, mode="otu") 
	makeGGplotOTU (inputdata)
}

##### COVERAGE stats
files = dir( pattern=glob2rx("*.clcassemblystats"))
for (f in files)
{
	print (f)
	makeBatchCoverage (f)
}

if (length(files)>0)
{
	inputdata = getData(files, mode="clccoverage") 
	makeGGplotCOVERAGE (inputdata)
}
