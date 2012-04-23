########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2012 J.Craig Venter Institute.
########################################################################################



# getName = function(x)
# {
# 	if (regexpr("/", x)>-1)
# 	{
# 		return (substr(x, 6,regexpr("/", x)-1 ))
# 	}
# 	else
# 	{
# 		return (x)
# 	}
# 	
# }
# 
# 
# dataset=data.frame()
# for (f in dir(pattern=glob2rx("*trim*"), recursive=TRUE))
# {
# 	if (regexpr("phylip.fn.rarefaction", f)>-1)
# 	{
# 		#print (f)
# 		x=getName(f)
# 		tmp = read.table(f, as.is=TRUE, header=TRUE)
# 		if (nrow(dataset)==0)
# 		{
# 			dataset = tmp[,c("numsampled", "X0.03", "lci.3", "hci.3")]
# 			Ns= paste(x, c("0.03", "lci", "hci"),sep="_")
# 			names(dataset) <- c("numsampled", Ns) 
# 		}			
# 		else
# 		{
# 			tmp = tmp[,c("numsampled", "X0.03", "lci.3", "hci.3")]
# 			Ns= paste(x, c("0.03", "lci", "hci"),sep="_")
# 			names(tmp) <- c("numsampled", Ns) 
# 			dataset = merge(dataset, tmp, by = "numsampled", all.x=TRUE, all.y=TRUE)
# 		}
# 	}
# }
# 
# xlim = c(0,50000)
# 
# labs = c(names(dataset)[regexpr("0.03",names(dataset))>-1])
# cols = c("magenta", "lightskyblue4", "salmon", "gray85", "yellow3")
# #cols = colors()[runif(length(labs), min=1, max=length(colors()))]
# pdf( "rarefaction.pdf", paper="special", width=10, height=8)
# par (mfrow = c(2,1))
# 
# plot (0,0, xlim=xlim,  ylim = range(dataset[,2:ncol(dataset)], na.rm=TRUE), type="n", xlab = "# sequences selected", ylab="# of OTUs", main = "rarefaction" )
# for (N in labs)
# {
# 	points (dataset[,1], dataset[, N], type="l", col=cols[labs==N], lwd=2)
# }
# legend("top", legend = labs, fill=cols, ncol=6, cex=0.8)
# 
# plot (0,0, xlim=xlim,  ylim = c(0,1), type="n", xlab = "# sequences selected", ylab="% coverage ", main="convergence to method's total OTU count" )
# for (N in labs)
# {
# 	points (dataset[,1], dataset[, N]/max(dataset[,N], na.rm=TRUE), type="l", col=cols[labs==N], lwd=2)
# }
# legend("bottom", legend = labs, fill=cols, ncol=6, cex=0.8)
# 
# 
# 
# dev.off()


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
		mypalette<-((brewer.pal(9, "Set1")))
		.colorcounter	= (.colorcounter %% 9) +1
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

getSampleName=function(x, statistic)
{
	S = regexpr("sorted\\.", x)
	E = regexpr(paste("\\.",statistic, sep=""), x)
	
	return (substr(x, S+7, E))	
}
makeCurveSet=function(input, main, xlim=c(0,0), ylim=c(0,0))
{

		#input = input[order(input$numsampled),]
		if (xlim[1] == 0 && xlim[2] ==0)
		{
			xlim = range(c(input$numsampled))
		}
		if (ylim[1] == 0 && ylim[2] ==0)
		{
			ylim = range(input[,seq(3, ncol(input), by=3)])
		}
		
		plot (0,0, type="n", main = main, ylab="# of OTUs", xlab="# sampled", xlim=xlim, ylim=ylim)
		dists = list()
		cols =  list()
		
		for (idx in 1:((ncol(input)-1)/3))
		{
			K = 2+(idx-1)*3	
			column = names(input)[K]
			color= getColors(column)
			
			dists = append(dists,column)
			cols = append(cols, color) 
	
			low = aggregate(input, by=list(input$numsampled), min)
			high = aggregate(input, by=list(input$numsampled), max)
			
			segments (low$numsampled, low[,K+2], high$numsampled, high[,K+3], col="gray80")
					
			lines (input$numsampled, input[,column],  cex=0.2, col=color)
			#points (input$numsampled, input[,column], pch=19, cex=0.2, col=color)
				
		}
		
		legend ("bottom", fill=unlist(cols), unlist(dists), title="Distance to form OTU", cex=0.75, ncol=3)
}
makeCurveSetDEBUG=function(input, main)
{

		input = input[order(input$numsampled),]
		xlim = range(c(input$numsampled))
		ylim = range(input[,seq(3, ncol(input), by=3)])
		
		plot (0,0, type="n", main = main, ylab="# of OTUs", xlab="# sampled", xlim=c(0,50000), ylim=c(0,10000))
		dists = list()
		cols =  list()
		
		for (idx in 1:((ncol(input)-1)/3))
		{
			K = 2+(idx-1)*3	
			column = names(input)[K]
			color= getColors(column)
			
			dists = append(dists,column)
			cols = append(cols, color) 
	
			low = aggregate(input, by=list(input$numsampled), min)
			high = aggregate(input, by=list(input$numsampled), max)
			
			#segments (low$numsampled, low[,K+2], high$numsampled, high[,K+3], col="gray80")
					
			lines (input$numsampled, input[,column],  cex=0.2, col=color)
			#points (input$numsampled, input[,column], pch=19, cex=0.2, col=color)
				
		}
		
		legend ("bottom", fill=unlist(cols), unlist(dists), title="Distance to form OTU", cex=0.5, ncol=4)
}
#### rarefaction, r_coverage, r_chao

stats = c("r_nseqs","rarefaction","r_simpson","r_invsimpson","r_chao","r_shannon","r_shannoneven","r_coverage")

#pdf("curves.pdf", paper="special", width=5, height=2.5*length(stats))
#par(mfrow=c(length(stats),1))
#for (stat in stats)
#{
#	####### rarefactions, global
#	files = dir(path=".", pattern=glob2rx(paste("*sorted.*",stat, sep="")) )
#	if (length(files)>0)
#	{
#	
#		input = data.frame()
#		curves = 0
#		for (f in files)
#		{
#			if (regexpr("groups", f)==-1)
#			{
#				print (f)
#				curves = curves+1
#				tmp = read.table(f, sep="\t", header=F, as.is=T, fill=T, skip=1)
#				headers = read.table(f, sep="\t", header=T, nrow=1, check.names=F)
#				headers = names(headers)
#				if (ncol(input)==0)
#				{
#					input = tmp	
#				}
#				else
#				{	
#					input = rbind(input, tmp)
#				}
#			}	/
#		}
#		
#		### remove missing rows
#		input = input[! is.na(input[,1]),]
#		
#		### remove missing columns
#		toremove = list()
#		for (k in 1:ncol(input))
#		{
#			if (sum( input[,k]== "") == nrow(input))
#			{
#				toremove = append(toremove, k)
#			}
#		}
#		if (length(toremove)>0)
#		{
#			toremove = -unlist(toremove)
#			input = input[,toremove]
#		}
#		### name the columns
#		names(input)=headers
#		
#		#print ((input[,2]))
#		
#		if (! is.numeric(input[,2]))
#		{
#			for (g in unique(input[,2]))
#			{
#				tmp = input[input[,2]==g,-2]
#				makeCurveSet(tmp, paste(g, ", ",stat, "\n", curves, " curves per distance", sep="" ))
#			}
#			
#		}
#		else
#		{
#			makeCurveSet(input, paste(stat, "\n", curves, " curves per distance", sep="" ))
#		}
#		
#	
#		
#	}
#
#
#	
#}	
#dev.off()

for (stat in stats)
{
	####### rarefactions, global
	files = dir(path=".", pattern=glob2rx(paste("*sorted.*",stat, sep="")) )
	if (length(files)>0)
	{
	
		input = data.frame()
		for (f in files)
		{
				curves = 1
				tmp = read.table(f, sep="\t", header=F, as.is=T, fill=T, skip=1)
				headers = read.table(f, sep="\t", header=T, nrow=1, check.names=F)
				headers = names(headers)
				
				input = tmp	
				
		
				### remove missing rows
				input = input[! is.na(input[,1]),]
				
				### remove missing columns
				toremove = list()
				for (k in 1:ncol(input))
				{
					if (sum( input[,k]== "") == nrow(input))
					{
						toremove = append(toremove, k)
					}
				}
				if (length(toremove)>0)
				{
					toremove = -unlist(toremove)
					input = input[,toremove]
				}
				### name the columns
				names(input)=headers
				
				#pdf(paste(f, "curves.pdf", sep="."), paper="special", width=5, height=2.5*length(stats))
				#par(mfrow=c(length(stats),1))
				
				pdf(paste(f, "curves.pdf", sep="."), paper="special", width=5, height=5)
				par(mfrow=c(1,1))
				
				if (! is.numeric(input[,2]))
				{
					maxx = 1000000
					maxy = 0
					for (g in unique(input[,2]))
					{
						tmp = input[input[,2]==g,-2]
						m = max(tmp$numsampled)
						maxx = min(m, maxx)

					}
					
					for (g in unique(input[,2]))
					{
						tmp = input[input[,2]==g,-2]	
						tmp = tmp[tmp$numsampled<maxx,]
						maxy = max(max(tmp[,seq(3, ncol(input), by=3)]), maxy)
					}
					
					for (g in unique(input[,2]))
					{
						tmp = input[input[,2]==g,-2]
						makeCurveSet(tmp, paste(g, ", ",stat, "\n", curves, " curves per distance", sep="" ), xlim = c(0, maxx), ylim = c(0, maxy))
					}
					
				}
				else
				{
					makeCurveSet(input, paste(stat, " (",f, ")\n", curves, " curves per distance", sep="" ))
				}
				
				dev.off()
				
		}	
	
		
	}


}	




