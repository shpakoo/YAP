########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2013 Sebastian Szpakowski
########################################################################################

getRandomColors=function(num)
{
	
	
	if (num==1)
	{
		if (exists ( ".colorcounter" , env = .GlobalEnv))
		{	
			.colorcounter	= get(".colorcounter", envir = .GlobalEnv)
		}
		else
		{	
			.colorcounter	= 0
		}
	
		library(RColorBrewer)
		mypalette<-((brewer.pal(12, "Paired")))
		.colorcounter	= (.colorcounter %% 12) +1
		K = mypalette[.colorcounter]		
		assign(".colorcounter", .colorcounter, envir = .GlobalEnv)
		return (K)
		
	}
	else if (num>12)
	{
		return (sample(rainbow(num)))
	}
	else
	{
		library(RColorBrewer)
		mypalette<-(brewer.pal(12, "Paired"))
		
		return (sample(mypalette, num))
	}

	return (sample(rainbow(num)))

}

getColors=function(x)
{
	if (exists ( "colorlookup" , env = .GlobalEnv))
	{	
		colorlookup	= get("colorlookup", envir = .GlobalEnv)
	}
	else
	{	
		colorlookup	= list()
	}

	otpt = list()
	for (taxon in x)
	{
	
		if (regexpr("\\.\\.\\.", taxon)>-1)
		{
			colorlookup[taxon] = "black"
		}
		
		else if (taxon %in% names(colorlookup))
		{
			#print (colorlookup[taxon])
		}
		else
		{
			colorlookup[taxon] = getRandomColors(1)
		}
			
		otpt = append(otpt, colorlookup[taxon])
		#print (otpt)
	}
	
	
	assign("colorlookup", colorlookup, envir = .GlobalEnv)
	return (unlist(otpt))
}


makeCurveSet=function(input, main, xlim=c(0,0), ylim=c(0,0), title="Distance to form OTU", ylab="# of OTUs")
{

		#input = input[order(input$numsampled),]
		if (xlim[1] == 0 && xlim[2] ==0)
		{
			xlim = range(c(input$numsampled), na.rm=TRUE)
		}
		if (ylim[1] == 0 && ylim[2] ==0)
		{
			
			ylim = range(input[,seq(4, ncol(input), by=3)], na.rm=TRUE)
		}

		plot (0,0, type="n", main = main, ylab=ylab, xlab="# sampled", xlim=xlim, ylim=ylim)
		dists = list()
		cols =  list()
		
		
		# plot confidence intervals first (start from last columns - furthers distance) 
		for (idx in ((ncol(input)-1)/3):1)
		{
			K = 2+(idx-1)*3	
			#column = names(input)[K]
			#color= getColors(column)
			
			#dists = append(dists,column)
			#cols = append(cols, color) 
	
			low = aggregate(input, by=list(input$numsampled), min)
			high = aggregate(input, by=list(input$numsampled), max)
			
			segments (low$numsampled, low[,K+2], high$numsampled, high[,K+3], col="gray90", lwd=0.5)
					
			#lines (input$numsampled, input[,column],  cex=0.2, col=color)
			#points (input$numsampled, input[,column], pch=19, cex=0.2, col=color)
				
		}
		
		# plot curves on the top 
		for (idx in ((ncol(input)-1)/3):1)
		{
			K = 2+(idx-1)*3	
			column = names(input)[K]
			color= getColors(column)
			
			dists = append(dists,column)
			cols = append(cols, color) 
			
#			low = aggregate(input, by=list(input$numsampled), min)
#			high = aggregate(input, by=list(input$numsampled), max)
			
			#segments (low$numsampled, low[,K+2], high$numsampled, high[,K+3], col="gray80")
			
			lines (input$numsampled, input[,column],  lwd=2, col=color)
			#points (input$numsampled, input[,column], pch=19, cex=0.2, col=color)
			
		}
		
		legend ("bottom", fill=unlist(cols), unlist(dists), title=title, cex=0.75, ncol=3)
}

#### rarefaction, r_coverage, r_chao

stats = c("r_nseqs","rarefaction","r_simpson","r_invsimpson","r_chao","r_shannon","r_shannoneven","r_coverage")


for (stat in stats)
{
	files = dir(path=".", pattern=glob2rx(paste("*sorted.*",stat, sep="")) )
	
	if (stat == "rarefaction")
	{
		ylab = "# of OTUs"
	}
	else
	{
		ylab = unlist( strsplit(stat, "_")[[1]][2])
	}
	
	if (length(files)>0)
	{
		input = data.frame()
		for (f in files)
		{
				
			input = read.table(f, sep="\t", header=T, as.is=T, fill=T, skip=0, check.names=F)

			filelab = unlist( strsplit(f, "\\.")[[1]][1])
			
			### remove missing rows
			input = input[! is.na(input[,1]),]
			
			### remove missing columns
			input = input[, ! names(input) %in% c( "X", "" )]	
			
			####### rarefactions, global
			if (regexpr("groups", f)==-1)
			{
				cat("global: ", f,"\n", sep="")
				pdf(paste(f, "curves.pdf", sep="."), paper="special", width=6, height=5)
				par(mfrow=c(1,1))	
				makeCurveSet(input, paste(stat, " (",filelab, ")\nall samples", sep="" ), title="Distance to form OTU", ylab = ylab)
				dev.off()
			}	
			####### rarefactions local, make plots per sample
			else
			{
				cat("local: ", f,"\n", sep="")
				samples = names(input)[seq(2, length(names(input)), by=3)]
				samples = unique(unlist(lapply(samples, function(x) { strsplit(x, "-")[[1]][2]} )))
				
				### first, establish the same scale for all plots
				
				xlim = range(c(input$numsampled), na.rm=TRUE)
				ylim = range(input[,seq(4, ncol(input), by=3)], na.rm=TRUE)			
					
				print (xlim)
				print (ylim)
				
				#pdf(paste(f, "curves.pdf", sep="."), paper="special", width=6, height=5*length(samples))
				#par(mfrow=c(length(samples), 1))
				pdf(paste(f, "curves.pdf", sep="."), paper="special", width=6, height=5)
				par(mfrow=c(1,1))
				for (s in samples)
				{	
					tmp = input[,names(input)=="numsampled"| regexpr(s, names(input))>-1 ]
					#print (names(tmp))
					names(tmp)[2:length(names(tmp))] = unlist(lapply(names(tmp)[2:length(names(tmp))], function(x) { strsplit(x, "-")[[1]][1]} ))
					makeCurveSet(tmp, paste(stat, " (", filelab, ")\n", s, "'s reads", sep="" ), title="Distance to form OTU" , ylab = ylab,  xlim=xlim, ylim=ylim)
				}
				dev.off()
			}
				
		}	
		
	}

}	

