library(Biostrings)
library(ggplot2)
library(grid)
library(scales)

multiplot <- function(..., plotlist=NULL, cols=1, layout=NULL)
{
	require(grid)
	
	# Make a list from the ... arguments and plotlist
	plots <- c(list(...), plotlist)
	numPlots = length(plots)
	
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


histopanel = function(stat, filename, xlab, xlim, color, extension)
{
	
	par(mar =c(0,4,3,2))
	boxplot(stat, horizontal=TRUE, ylim=xlim, ylim = c(0,1), lwd=1, main="", ylab="",  col=color, axes=FALSE,  xlim=c(0.8, 1.2), outline=F)
	labs = c(median(stat), quantile(stat, c(0.25, 0.75)))
	axis(1, at=labs, lab = round(labs,1) ,  las=2)
	
	par(mar =c(5,4,4,2))
	hist(stat,main="", xlab=xlab, col=color, xlim=xlim, br=seq(xlim[1]-extension, xlim[2]+extension, length.out=150), axes=FALSE, ylab=paste("frequency"))
	labs = seq(min(stat), max(stat), length.out=10)
	axis(1, at = labs, lab=round(labs,1), las=0)
	axis(2)
	box()
	
}

makePlot = function(fastafile, fastas=list(),  main="Fasta", prefix="",...)
{
	layout( matrix(c(1,2,3,4), 2, 2, byrow = FALSE), heights=c(0.1, 0.9), widths=c(0.15, 0.85))
	
	cols = c("blue", "lightblue", "red", "green", "seagreen", "lightgreen", "tomato1", "violetred")
	random= round(runif(1, min=1, max=length(cols)),0)
	
	if (length(fastas)==0)
	{
		#fastas= readFASTA(fastafile, strip.descs=TRUE)
		fastas = readDNAStringSet(fastafile, use.names=FALSE)
		lengths = width(fastas)
		freqs = alphabetFrequency(fastas, baseOnly=TRUE)
		actualchars = apply(alphabetFrequency(fastas, baseOnly=TRUE, collapse=FALSE)[,1:4],1,sum)
	}
	else
	{
		lengths = fastas$lengths
		actualchars = fastas$actualchars
	}
		
	stat= log10(lengths)
	desc = "sequence lengths"
	xlim=c(min(stat),max(stat))
	histopanel(stat, fastafile, desc, xlim, cols[random], 20)
	
	#~ stat= log10(lengths)
	#~ desc = "log 10 of sequence length"
	#~ xlim=c(min(stat),max(stat))
	#~ histopanel(stat, fastafile, desc, xlim, cols[random], log10(20))
	
	stat = log10(actualchars)
	desc = "lengths excluding excluding ,.-"
	xlim=c(min(stat),max(stat))
	histopanel(stat, fastafile, desc, xlim, cols[random], 20)
	mtext(paste("Distribution of fragment sizes in\n", main, "\nn=", length(stat), sep=""), outer=TRUE, side=3, line=-5)
	
	return (list(lengths=lengths, actualchars=actualchars))
	
			
}
makeGGplot=function(fastafiles, fastas=list(), main=fastafile, prefix="",...)
{
	global = data.frame()
	for (f in fastafiles)
	{
		fastas = readDNAStringSet(f, use.names=FALSE)
		lengths = width(fastas)
		freqs = as.data.frame(alphabetFrequency(fastas, baseOnly=TRUE))
		actualchars = apply(alphabetFrequency(fastas, baseOnly=TRUE, collapse=FALSE)[,1:4],1,sum)
		alldata = cbind(freqs, actualchars, lengths)
		alldata$gcrelative = (alldata$C + alldata$G) / (alldata$A + alldata$C + alldata$G + alldata$T) 
		alldata$sizecat=as.character(floor(log10(actualchars)))
		
		a = strsplit( f, "\\.")[[1]]
		if (length(a)>2)
		{
			a = a[2]
			b = strsplit( a,"-")[[1]]
			if (length(b)>2)
			{
				f=b[1]
			}
			else
			{
				b=f
			}
		}	
		else
		{
			b=f
		}
		cat(f, "->", a, "->", b, "\n")
	
		
		alldata$file=rep(b, nrow(alldata))
		if (nrow(global)==0)
		{
			global=alldata
		}
		else
		{
			global=rbind(alldata,global) 
		}
	}
	
	#maintitle = paste("Distribution of sequence sizes in\n", main, "\n(n=", length(lengths),")", sep="")
	
	A = (ggplot(global,  aes(file, lengths, color=file)) 
			+ geom_boxplot() 
			+ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
							labels = trans_format("log10", math_format(10^.x)))
			+ coord_flip()
			+ theme(axis.ticks = element_blank(), axis.text.y = element_blank())
			
			)		
			
	B = (ggplot(global,  aes(sizecat , (C + G) / (A + C + G + T), color=file )) 
			+ geom_jitter(alpha = I(1/100)) 
			+ geom_boxplot()
			#+ scale_x_discrete(labels=math_format( "[",  10^.x , 10^.x+1,  ")") )
			+ scale_x_discrete(labels=math_format( expr=paste("[",10^.x, ",", 10^(.x+1), ")"  )))
			)
	


	C1 = (ggplot(global,  aes(lengths, ..count..,  fill=file, group=file )) 
			+ stat_bin(aes(y=..count..), geom="area", position="stack", weight=2)
			+ scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
							labels = trans_format("log10", math_format(10^.x)))
			)

	C2 = (ggplot(global,  aes(lengths, ..count..,  fill=file, group=file )) 
			+ stat_bin(aes(y=..count..), geom="area", position="dodge", alpha=0.5, weight=2)
			+ scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
					labels = trans_format("log10", math_format(10^.x)))
			)

	
	
#	C1 = (ggplot(global,  aes(lengths, ..count..,  fill=file, group=file )) 
#			+ geom_bar(position="stack")
#			+ scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
#							labels = trans_format("log10", math_format(10^.x)))
#			)
#
#	C2 = (ggplot(global,  aes(lengths, ..density..,  fill=file, group=file )) 
#			+ geom_bar(position="stack")
#			+ scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
#							labels = trans_format("log10", math_format(10^.x)))
#			)
#
#	D1 = (ggplot(global,  aes(lengths, ..count..,  fill=file )) 
#			+ geom_bar(position = "dodge")
#			+ scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
#							labels = trans_format("log10", math_format(10^.x)))
#			)
#
#	D2 = (ggplot(global,  aes(lengths, ..density..,  fill=file )) 
#			+ geom_bar( position= "dodge")
#			+ scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
#							labels = trans_format("log10", math_format(10^.x)))
#			)
#
#	E1 =  (ggplot(global,  aes(lengths, fill=file )) 
#			+ stat_bin(aes (ymax = ..count.. ))	
#				)



#	D1 = (ggplot(global,  aes(lengths ,  fill=file )) 
#				+ geom_bar(position="dodge")
#				+ scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
#						labels = trans_format("log10", math_format(10^.x)))
#				)
#	
#	D2 = (ggplot(global,  aes(lengths, ..density.. ,  fill=file )) 
#			+ geom_bar(position="dodge")
#			+ geom_density(alpha=0.5, weight=2)
#			+ scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
#					labels = trans_format("log10", math_format(10^.x)))
#			)
#
#
	F = (ggplot(global,  aes(lengths, fill=file, group=file )) 
				+ geom_bar(position="fill")
				+ scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
						labels = trans_format("log10", math_format(10^.x)))
				)
	
	

	
	
	
	#pdf("test.pdf", width=10, height = 15)
	png("FastaFileStats_GG.png", width=10, height = 20, units="in", res=250)
	multiplot(A, B,  C1, C2,  F, ncol=1)
	dev.off()
	
	invisible(alldata)
}



files = dir(pattern=glob2rx("*.muscle"))
files = append(files, dir(pattern=glob2rx("*.fa")))
files = append(files, dir(pattern=glob2rx("*.align")))
files = append(files, dir(pattern=glob2rx("*.fasta")))


pdf(paste("FastaFileStats_indi.pdf", sep=""), width = 20, height = 8)	

L = c()
A = c()

for (f in files)
{
	print (f)
	x = makePlot(f, main=f)
	L = append(L, x$lengths)
	A = append(A, x$actualchars)
	
}
dev.off()

pdf(paste("FastaFileStats_comb.pdf", sep=""), width = 20, height = 8)	
makePlot("", fastas=list(lengths=L, actualchars=A), main=paste("combined ", length(files)," files", sep=""))
dev.off()


makeGGplot(files)
