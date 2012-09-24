library(Biostrings)
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
		fastas= readFASTA(fastafile, strip.descs=TRUE)
		lengths = unlist(lapply(fastas, function(x){nchar(x$seq) }))
		actualchars = unlist(lapply(fastas, function(x){ x=unlist(strsplit(x$seq, "[,.-]")) ; sum(nchar(x)) }))
	}
	else
	{
		lengths = fastas$lengths
		actualchars = fastas$actualchars
	}
		
	stat= lengths
	desc = "sequence lengths"
	xlim=c(min(stat),max(stat))
	histopanel(stat, fastafile, desc, xlim, cols[random], 20)
	
	#~ stat= log10(lengths)
	#~ desc = "log 10 of sequence length"
	#~ xlim=c(min(stat),max(stat))
	#~ histopanel(stat, fastafile, desc, xlim, cols[random], log10(20))
	
	stat = actualchars
	desc = "lengths excluding excluding ,.-"
	xlim=c(min(stat),max(stat))
	histopanel(stat, fastafile, desc, xlim, cols[random], 20)
	mtext(paste("Distribution of fragment sizes in\n", main, "\nn=", length(stat), sep=""), outer=TRUE, side=3, line=-5)
	
	return (list(lengths=lengths, actualchars=actualchars))
	
			
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
