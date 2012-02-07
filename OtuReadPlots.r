########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2012 Sebastian Szpakowski, J.Craig Venter Institute.
########################################################################################


makePlot = function(filename = "read.groups.tab.txt", minreads = 0, mingroups = 0, nlevels = 500)
{
	x = read.table(filename, as.is=T, header=F, sep="\t")
	x=x[-1,]
	names(x) = c("reads", "groups")
	
	#print(x)
	### setup
	nf <- layout(matrix(c(1,2,3,0), 2, 2, byrow=TRUE), respect=FALSE, width=c(0.8,0.2))
	plot (log10(x$reads), x$groups, type= "n", ylab = "number of samples showing an OTU", xlab = "OTU size (i.e. number of reads in OTU, log10 scale)", ylim = c(min(x$groups)-0.5, max(x$groups)+0.5), axes=FALSE)
	groups = unique(x$groups)
	axis(2, at=groups, lab=groups, las=2)
	
	A = c(1:9)
	B =seq(0, log10(max(x$reads)), by=1)
	C =expand.grid(A, B)
	C$at = C[,1] * 10^(C[,2])
	C$lab = ifelse (C[,1]%%10==1, 10^C[,2], "")
	
	axis(1, at=log10(C$at), lab=rep("", nrow(C)) , col.ticks="gray", lwd.ticks=2 )
	C2 = C[C$lab!="", ]
	axis(1, at=log10(C2$at), lab=C2$lab , col.ticks="black", lwd.ticks= 4 )
	
	
	#######DEBUG CODE
	#b = rnorm(length(x$groups), mean=0, sd=0.001)
	#x$groups = x$groups+b
	#x$groups = sample(c(1:2), length(x$groups), replace=T)
	########
	
	x = x[x$reads>minreads & x$groups>mingroups,]
	
	#return (x)
	origreads = x$reads
	origreads[origreads>100]=100
	
	x$reads = log10(x$reads)
	
	
	# contour
	require(MASS)
	xykde = MASS::kde2d(x$reads, x$groups, lims = c( min(x$reads)-5, max(x$reads)+5, min(x$groups)-5, max(x$groups)+5) )
	
	if (sum(is.nan(xykde$z))>0)
	{
		cat("have to nudge groups...\n")
		b = x$groups + sample(c(-0.1, 0.1), length(x$groups), replace=T)
		xykde = MASS::kde2d(x$reads, b, lims = c( min(x$reads)-5, max(x$reads)+5, min(x$groups)-5, max(x$groups)+5) )
	}
	
	zlim = range(xykde$z, finite = TRUE)
	lev = seq(zlim[1], zlim[2], le = nlevels)
	col = rev(heat.colors(length(lev)))
	contour(xykde, add = TRUE, levels = lev,
		drawlabels = FALSE, col=col  )        
	
			
	# trendline
	trendline = lowess((x$reads), x$groups, f=0.1)
	points(trendline , type="l", lwd=5, col="gray30")
	
	# points
	points(x$reads, x$group, pch=19, col="gray70", cex=0.5) 
	
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

	
	at = unique(round(seq(0, max(lev), length.out=10),2))
	lab = paste(round((1-at) * 100,0) , "%", sep="")
	axis(4, las=1, cex.axis=0.5)
	axis(2, las=1, at = at, lab=lab)
	box()
		
	
	# histogram
	hist ((origreads), col="gray", br=seq(0.5,100.5, by=1), axes =FALSE, xlab="OTU size")
	axis(2)
	at = c(1:10, seq(20,100, by=10) )
	labs = paste(at)
	labs [2:9]=""
	
	labs[length(labs)] = paste(labs[length(labs)], "+", sep="")
	axis(1, at = at, lab=labs, las=3)

	#print ( nrow(x[x$groups<=1 & x$reads<=0.1,]) )
	#print ( nrow(x))

	 
}

makeBatch = function(filename) 
{
	pdf(paste(filename,"otusizes.pdf", sep="."), paper="special", width=11, height=10)
	makePlot(filename = filename, minreads=0, mingroups=0)
	makePlot(filename = filename, minreads=1, mingroups=0)
	makePlot(filename = filename, minreads=1, mingroups=1)
	makePlot(filename = filename, minreads=5, mingroups=1)
	dev.off()
}
