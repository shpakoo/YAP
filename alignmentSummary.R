########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2012 J.Craig Venter Institute.
########################################################################################

reassemble=function(L)
{	
	otpt = data.frame()
	for (k in L)
	{
		if (nrow(otpt) ==0)
		{
			otpt = k
		}
		else
		{
			otpt = rbind(otpt, k)
		}
	}
	return (otpt)
}

paralellize=function(input, func, cpus=8, ...)
{
	library(multicore)	
	input$xxxid  = 1:nrow(input) 
	input$xxcpu = input$xxxid %% cpus
	
	#print (input[1:10, c("xxxid", "xxcpu")])
		
	allmyprocesses = list()
	counter=0
	
	for (core in unique(input$xxcpu))
	{
		gc()
		cat(core)
		tmp = input[input$xxcpu==core,]		
		o <- parallel( func(tmp,  ... ) 	)
		counter = counter + 1
		allmyprocesses[[counter]] = o
		cat(".", sep="")	
	}
	
	cat("waiting...", sep="")
	tmp = collect(allmyprocesses, wait = TRUE)
	cat("reassembling...", sep="")
	tmp = reassemble(tmp)			
	cat("\n")
	return (tmp)
}

convertAlignment=function(input)
{	
	otpt = data.frame()
	
	
	for (k in 1:nrow(input))
	{
		
		#print (names(input))
		#print (input$desc[1:10])
		acgt = input$seq[k]
		head = input$desc[k]	
		
		#cat(name, string, "\n", sep=" ")
		
		mask = chartr("acgtACGT.-", "1111111120", acgt)
		mask = as.numeric(unlist(strsplit(mask, "")))
		mask[mask==2]=-1
		
		#cat(head, acgt,  "\n",sep=" ")
		
		tmp = matrix( mask, nrow=1, ncol=nchar(acgt), byrow=TRUE )
		tmp = data.frame(tmp, stringsAsFactors=FALSE)
		
		row.names(tmp)<-c(head)
		names(tmp)<-(1:nchar(acgt))
		#print (tmp)
		#print (dim(tmp))
		
		if (nrow(otpt)==0)
		{
			otpt = tmp
		}
		else
		{
			otpt = rbind(otpt, tmp)
		}
		
	}		
	return (otpt)	
}
makeAlignmentHistogram=function(input, ref)
{
	reference = input[regexpr(ref, row.names(input))>-1,] 
	allelse = input[! (row.names(input) %in% row.names(reference)), ]
	
	print (dim(reference))
	print (dim(allelse))
	
	reference = reference[1,]
	
	locations = names(reference)[reference[1,]==1]
	#print (locations)
	
	xs = c()
	labs=c()
	
	ys = c()
	
	### mothur decided to cut the first 19 bases off
	counter=20
	for (L in locations)
	{
		tmp = allelse[,L]
		
		xs = c(xs, counter)
		ys = c(ys, length(tmp[tmp==1]))
		labs=c(labs, L)
		
		counter=counter+1
	}
	
	plot(0, 0 , type="n",  axes=FALSE, xlab="position on e-coli's 16S consensus", xlim =range(xs), ylim = range(ys), ylab="count of aligned bases, (\"reads\" / species) at a position ")
	abline()
	
	pattern= (1:10)%/%10 == 1
	abline(v=xs[pattern], col="gray90")
	
	### V8 
	coords = c(1243, 1294)
	info = "V8"
	rect(coords[1], min(ys), coords[2], max(ys), col="gray80", border="gray90")
	text(coords[2], min(ys), paste(coords[1], "-", coords[2],"\n" , info, sep=""), cex=0.75, pos=4)
	
	### V7 
	coords = c(1117, 1173)
	info = "V7"
	rect(coords[1], min(ys), coords[2], max(ys), col="gray80", border="gray90")
	text(coords[2], min(ys), paste(coords[1], "-", coords[2],"\n" , info, sep=""), cex=0.75, pos=4)
	
	### V6 
	coords = c(986, 1043)
	info = "V6"
	rect(coords[1], min(ys), coords[2], max(ys), col="gray80", border="gray90")
	text(coords[2], min(ys), paste(coords[1], "-", coords[2],"\n" , info, sep=""), cex=0.75, pos=4)

	### V5 
	coords = c(822, 879)
	info = "V5"
	rect(coords[1], min(ys), coords[2], max(ys), col="gray80", border="gray90")
	text(coords[2], min(ys), paste(coords[1], "-", coords[2],"\n" , info, sep=""), cex=0.75, pos=4)

	### V4 
	coords = c(576, 682)
	info = "V4"
	rect(coords[1], min(ys), coords[2], max(ys), col="gray80", border="gray90")
	text(coords[2], min(ys), paste(coords[1], "-", coords[2],"\n" , info, sep=""), cex=0.75, pos=4)
		
	### V3 
	coords = c(433, 497)
	info = "V3"
	rect(coords[1], min(ys), coords[2], max(ys), col="gray80", border="gray90")
	text(coords[2], min(ys), paste(coords[1], "-", coords[2],"\n" , info, sep=""), cex=0.75, pos=4)

	### V2 
	coords = c(137, 242)
	info = "V2"
	rect(coords[1], min(ys), coords[2], max(ys), col="gray80", border="gray90")
	text(coords[2], min(ys), paste(coords[1], "-", coords[2],"\n" , info, sep=""), cex=0.75, pos=4)
	
	### V1 
	coords = c(69,99)
	info = "V1"
	rect(coords[1], min(ys), coords[2], max(ys), col="gray80", border="gray90")
	text(coords[2], min(ys), paste(coords[1], "-", coords[2],"\n" , info, sep=""), cex=0.75, pos=4)

	### 27F
	coords = c(8, 27)
	info = "27F"
	lines(coords, c(median(ys), median(ys)), lwd=5, col="red" )
	text(coords[2], median(ys), paste(coords[1], "-", coords[2],"\n" , info, sep=""), cex=0.75, pos=4)
	
	### 534R
	coords = c(518, 534)
	info = "534R"
	lines(coords, c(median(ys), median(ys)), lwd=5, col="red" )
	text(coords[2], median(ys), paste(coords[1], "-", coords[2],"\n" , info, sep=""), cex=0.75, pos=4)
	
	### 1492R
	coords = c(1498, 1498+30)
	info = "1492R"
	lines(coords, c(median(ys), median(ys)), lwd=5, col="red" )
	text(coords[2], median(ys), paste(coords[1], "-", coords[2],"\n" , info, sep=""), cex=0.75, pos=2)
	
	
	# read counts
	points(xs, ys, col="yellowgreen",lwd=2, type="h")
	
	# trendline
	points(xs, lowess(ys, f=0.1)$y , type="l", lwd=6, col="gray50")
	
	axis(1, at=xs[pattern], lab=xs[pattern], las=2, cex.axis=0.5, col="blue")
	axis(3, at=xs[pattern], lab=labs[pattern], las=2, cex.axis=0.5, col="lightblue")
	axis(2)

	box()
	
	legend ("bottomright", 
		c("read counts", "trendline", "primer locations", "\"V\" regions", "E.coli coordinates", "SILVA coordinates"),
		fill=c("yellowgreen", "gray50", "red", "gray80", "blue", "lightblue"), 	
		ncol=1, cex=0.55, bg="ivory")
	
}

makeAlignmentCountPlot=function(input, ref)
{
	reference = input[regexpr("e_coli", row.names(input))>-1,] 
	allelse = input[! (row.names(input) %in% row.names(reference)), ]
	
	print (dim(reference))
	print (dim(allelse))
	
	
	locations = names(reference)[reference[1,]==1]
	print (locations)
	
	xs = c()
	labs=c()
	
	ys = c()
	
	### mothur decided to cut the first 19 bases off
	counter=20
	for (L in locations)
	{
		tmp = allelse[,L]
		
		xs = c(xs, counter)
		ys = c(ys, length(tmp[tmp==1]))
		labs=c(labs, L)
		
		counter=counter+1
	}
}

batch=function(inputfile, description, pattern="e_coli", cpus = 16)
{


	#~ ### with gaps
	library(Biostrings)
	fastas = readFASTA(inputfile, strip.descs=TRUE)
	descs = unlist(lapply (fastas, function(x){x$desc}))
	seqs  = unlist(lapply (fastas, function(x){x$seq}))

	seqlookup = data.frame(lookupid = 1:length(fastas), desc = descs, seq = seqs, stringsAsFactors=FALSE )
	refid = seqlookup$lookupid[regexpr(pattern, seqlookup$desc)>-1]
	refid = refid[1]

	out = paralellize(seqlookup, cpus=cpus, convertAlignment )
	save(list=c("out"), file=paste("ecoli.alignment.",description,".RData", sep=""))

	pdf(paste("ecoli.alignment.",description,".pdf", sep=""), width=15, height=6)
	makeAlignmentHistogram(out, pattern)
	dev.off()

}


makeAlignmentHistogram2=function(input, ref, trimstart=0, trimend=0)
{

	labs = input$global
	xs = input$local
	### mothur decided to cut the first 19 bases off
	if (ref=="e_coli2_genbank")
	{
		xs= xs + 27
	}

	ys = input$val
	
	print (range(xs))
	print (range(ys))
		
	###############################################################################################	
	### for labeling purposes make sure that the Y axis goes 10 percent higher and lower than data	
	Ry = abs(diff(range(ys))) * 0.1
	ylim =  c( min(ys) - Ry, max(ys) + Ry)
	
	Rx = abs(diff(range(xs))) * 0.06
	xlim =  c( min(xs), max(xs) + Rx )
	###############################################################################################
	plot(0, 0 , type="n",  axes=FALSE, xlab="position on REFERENCE consensus", xlim =xlim, ylim = ylim, ylab="count of aligned bases, (\"reads\" / species) at a position ")
	abline()
	
	#pattern= (1:10)%%10 == 1
	pattern = seq(xlim[1], xlim[2], length.out=50 )
	print (pattern)
	
	abline(v=xs[pattern], col="gray90")
	
	if (ref=="e_coli2_genbank")
	{
		### V8 
		coords = c(1243, 1294)
		info = "V8"
		rect(coords[1], ylim[1], coords[2], ylim[2], col="gray80", border="gray90")
		text(coords[1], ylim[1], paste(coords[1], "-", coords[2],"\n" , info, sep=""), cex=0.75, pos=4)
		
		### V7 
		coords = c(1117, 1173)
		info = "V7"
		rect(coords[1], ylim[1], coords[2], ylim[2], col="gray80", border="gray90")
		text(coords[1], ylim[1], paste(coords[1], "-", coords[2],"\n" , info, sep=""), cex=0.75, pos=4)
		
		### V6 
		coords = c(986, 1043)
		info = "V6"
		rect(coords[1], ylim[1], coords[2], ylim[2], col="gray80", border="gray90")
		text(coords[1], ylim[1], paste(coords[1], "-", coords[2],"\n" , info, sep=""), cex=0.75, pos=4)
	
		### V5 
		coords = c(822, 879)
		info = "V5"
		rect(coords[1], ylim[1], coords[2], ylim[2], col="gray80", border="gray90")
		text(coords[1], ylim[1], paste(coords[1], "-", coords[2],"\n" , info, sep=""), cex=0.75, pos=4)
	
		### V4 
		coords = c(576, 682)
		info = "V4"
		rect(coords[1], ylim[1], coords[2], ylim[2], col="gray80", border="gray90")
		text(coords[1], ylim[1], paste(coords[1], "-", coords[2],"\n" , info, sep=""), cex=0.75, pos=4)
			
		### V3 
		coords = c(433, 497)
		info = "V3"
		rect(coords[1], ylim[1], coords[2], ylim[2], col="gray80", border="gray90")
		text(coords[1], ylim[1], paste(coords[1], "-", coords[2],"\n" , info, sep=""), cex=0.75, pos=4)
	
		### V2 
		coords = c(137, 242)
		info = "V2"
		rect(coords[1], ylim[1], coords[2], ylim[2], col="gray80", border="gray90")
		text(coords[1], ylim[1], paste(coords[1], "-", coords[2],"\n" , info, sep=""), cex=0.75, pos=4)
		
		### V1 
		coords = c(69,99)
		info = "V1"
		rect(coords[1], ylim[1], coords[2], ylim[2], col="gray80", border="gray90")
		text(coords[1], ylim[1], paste(coords[1], "-", coords[2],"\n" , info, sep=""), cex=0.75, pos=4)
	
		### 27F
		coords = c(8, 27)
		info = "27F"
		lines(coords, rep(ylim[2],2), lwd=5, col="red" )
		text(coords[2], ylim[2], paste(coords[1], "-", coords[2],"\n" , info, sep=""), cex=0.5, pos=4)
		abline(v=coords[2], lty=3, col="red", lwd=1)
		
		### 534R
		coords = c(518, 534)
		info = "534R"
		lines(coords, rep(ylim[2],2), lwd=5, col="red" )
		text(coords[1], ylim[2], paste(coords[1], "-", coords[2],"\n" , info, sep=""), cex=0.5, pos=2)
		abline(v=coords[1], lty=3, col="red", lwd=1)
		
		### 1492R
		coords = c(1498, 1498+30)
		info = "1492R"
		lines(coords, rep(ylim[2],2), lwd=5, col="red" )
		text(coords[1], ylim[2], paste(coords[1], "-", coords[2],"\n" , info, sep=""), cex=0.5, pos=2)
		abline(v=coords[1], lty=3, col="red", lwd=1)
	}	
	
	# read counts
	points(xs, ys, col="yellowgreen",lwd=1, type="h")
		
	# trendline
	points(xs, lowess(ys, f=0.05)$y , type="l", lwd=5, col="gray40", lty=1)

	
	axis(1, at=xs[pattern], lab=xs[pattern], las=2, cex.axis=0.5, col="blue")
	axis(3, at=xs[pattern], lab=labs[pattern], las=2, cex.axis=0.5, col="lightblue")
	axis(2)

	if (trimstart!=0)
	{
		print (xs[labs==trimstart])
		abline(v=xs[labs==trimstart], col="gray80", lwd=2, lty=3)
		
	}
	if (trimend !=0)
	{
		print (xs[labs==trimend])
		abline(v= xs[labs==trimend], col="gray80", lwd=2, lty=3)
	}
	
	box()
	
	if (ref=="16S")
	{
		legend ("right", 
			c("read counts", "trendline", "primer locations", "\"V\" regions", "REFERENCE coordinates", "SILVA coordinates"),
			fill=c("yellowgreen", "gray50", "red", "gray80", "blue", "lightblue"), 	
			ncol=1, cex=0.75, bg="ivory")
	}
	else
	{
		legend ("right", 
				c("read counts", "trendline", "REFERENCE coordinates", "ALIGNMENT coordinates"),
				fill=c("yellowgreen", "gray50", "blue", "lightblue"), 	
				ncol=1, cex=0.75, bg="ivory")
	}
	
}

batch2 = function(inputfile, ref="e_coli2_genbank", trimstart=0, trimend=0)
{
	input = read.table(inputfile, as.is=TRUE, header=FALSE)
	names(input) <- c("global", "local", "val")
	pdf(paste(inputfile,".pdf", sep=""), width=15, height=6)
	makeAlignmentHistogram2(input, ref, trimstart, trimend)
	dev.off()
}





