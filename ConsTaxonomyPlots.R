########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2012 J.Craig Venter Institute.
########################################################################################

library(fdrtool)
library(RColorBrewer)
library(ade4)
library(vegan)
library(randomForest)
library(gplots)

collookup=function(x, annotation, column="SampleID")
{
	val = annotation[annotation$SampleID==x, column] 	
	library(RColorBrewer)
	allvals = sort(unique(annotation[, column]))
	#mypalette <- colorRampPalette(brewer.pal(12, "Paired"), space="rgb")
	mypalette = rainbow(length(allvals))
	return (mypalette[allvals == val])
}

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
		#mypalette<-((brewer.pal(9, "Set1")))
		#.colorcounter	= (.colorcounter %% 9) +1
		
		mypalette <- colorRampPalette(brewer.pal(9, "Set1"), space="rgb")
		mypalette = mypalette(20)
		
		.colorcounter	= (.colorcounter %% 20) +1
		
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
	if (exists ( "taxoncolors" , env = .GlobalEnv))
	{	
		taxoncolors	= get("taxoncolors", envir = .GlobalEnv)
	}
	else
	{	
		taxoncolors	= list()
	}

	otpt = list()
	for (taxon in x)
	{
	
		if (regexpr("\\.\\.\\.", taxon)>-1)
		{
			taxoncolors[taxon] = "black"
		}
		
		else if (taxon %in% names(taxoncolors))
		{
			#print (taxoncolors[taxon])
		}
		else
		{
			taxoncolors[taxon] = getRandomColors(1)
		}
		
		
		otpt = append(otpt, taxoncolors[taxon])
		#print (otpt)
	}
	
	
	assign("taxoncolors", taxoncolors, envir = .GlobalEnv)
	return (unlist(otpt))
}

dist.binary = function(x, ... )
{
	return (dist(x, method="binary", ...))
}

fundist= function(x, ... )
{
	return (vegdist(x, method="bray", ...))
}

funhclust = function(x, ... )
{
	return (hclust(x, method="average", ...))
}


makeLevelDataset=function(dataset, level=4, cumulative = FALSE)
{
	if (! cumulative)
	{
		orig = dataset
	 	toreturn = dataset[dataset$depth==level,]
	 	
	 	OTHER = data.frame(rbind (rep(0, ncol(toreturn))))
	 	names(OTHER) <- names(dataset)
	 	
	 	toreturn=as.data.frame(rbind(toreturn, OTHER))	
	 	toreturn[nrow(toreturn), "label"] = "Unclassified"
	 	
	 	
	 	for (taxonid in unique(dataset[dataset$depth==level,"taxonid"]))
	 	{
	 		tmp = orig[orig$depth>level, ]
	 		tmp = tmp[regexpr(paste(taxonid, "\\.", sep=""), tmp$taxonid)==1,] 		
	
	 		for (samp in names(toreturn)[4:ncol(toreturn)])
	 		{
	 			x = sum(tmp[,samp])
	 			toreturn [toreturn$taxonid==taxonid, samp] = as.numeric(toreturn [toreturn$taxonid==taxonid, samp]) + x			
	 		} 			
	 	} 	
	 	
	 	for (samp in names(toreturn)[4:ncol(toreturn)])
	 	{
	 			
	 			x = sum(as.numeric(dataset [, samp])) -  sum(as.numeric(toreturn [, samp]))
	 			toreturn[toreturn$label=="Unclassified", samp] = x
	 	}
	 	return (toreturn)
	}
	else
	{
		
		tmp = dataset[dataset$depth==level,]
		return (tmp)
		
	}	
 		
}


ordering = function(...)
{
	return (reorder(..., agglo.FUN = "max"))
}
makeHeatMap=function(dataset, annotation, highlight="SampleID", main="Heatmap", perc=TRUE)
{

	tmp = prepareDataset(dataset)
	
	mypalette <- colorRampPalette(c("blue","yellow", "red"), space="rgb")
	mypalette = mypalette(20)

	#### normalize	
	if (perc)
	{
		for (k in 1:ncol(tmp))
		{
			### the dataset is per taxon level!
			### i.e. NOT cummulative  
			tmp[,k] = tmp[,k]/sum(tmp[,k])
		}
		lab = paste(round(seq(0,max(tmp*100), length.out=11),0), "%")
		at = seq(0,max(tmp), length.out=11)	
		
		sumsfilter = (apply(tmp, 1, sd)>0.01)	
	}
	else
	{
		sumsfilter = (apply(tmp, 1, sd)>0)
	}
	
	
	if (sum(sumsfilter)<3)
	{
		plot(0,0, axes=FALSE, type="n", xlab="", ylab="", main=paste(main, "nothing to plot", sep="\n") )
	}
	
	else
	{
	
		tmp = tmp[sumsfilter,]		
		tmp = data.matrix(tmp)
		tmp = t(tmp)
		tmp = tmp[row.names(tmp) %in% annotation$SampleID,]
			
		colors = unlist( lapply(dimnames(tmp)[[1]], collookup, annotation, column=highlight) )
		
		legendstuff = data.frame( fill = colors, SampleID = dimnames(tmp)[[1]], stringsAsFactors=F )
		minianno = annotation[,c("SampleID", highlight)]
		legendstuff = merge(legendstuff, minianno, by="SampleID", all.x=T, all.y=F)
		legendstuff = aggregate(legendstuff, list(legendstuff[,highlight]), unique)
	
		description = paste(main, ifelse(perc, "OTUs normalized per sample", "Raw counts"), "Bray-Curtis, Average Linkage", sep="\n")
		
		#heatmap(tmp, col = mypalette, margins=c(20,20), keep.dendro=T, distfun = fundist, hclustfun = funhclust, scale=NULL, RowSideColors=colors, xlab="taxons", ylab="samples", main=paste(main, sep="\n"))
		#legend("left", legend = legendstuff[,highlight], fill = legendstuff$fill, col = legendstuff$fill, cex=0.5)
	
		#heatmap.2( xlab="taxons", ylab="samples", main=paste(main, sep="\n"))
		heatmap.2(tmp, margins=c(25,5), distfun = fundist, hclustfun = funhclust, col=mypalette, RowSideColors=colors, trace="none", scale="none", tracecol="white", main = description)
		
		legend("bottomleft", legend = legendstuff[,highlight], fill = legendstuff$fill, col = legendstuff$fill, cex=0.8, horiz=T )
	
	}
}

mykde2d = function (dfxy, xax = 1, yax = 2, pch = 20, cpoint = 1, neig = NULL, 
    cneig = 2, xlim = NULL, ylim = NULL, grid = TRUE, addaxes = TRUE, 
    cgrid = 1, include.origin = TRUE, origin = c(0, 0), sub = "", 
    csub = 1.25, possub = "bottomleft", pixmap = NULL, contour = NULL, 
    area = NULL, add.plot = FALSE, ... ) 
{
    if (!require(MASS)) 
        stop("library MASS required for kde2d")
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    s.label(dfxy, xax = xax, yax = yax, clab = 0, pch = pch, 
        cpoint = cpoint, neig = neig, cneig = cneig, xlim = xlim, 
        ylim = ylim, grid = grid, addaxes = addaxes, cgrid = cgrid, 
        include.origin = include.origin, origin = origin, sub = sub, 
        csub = csub, possub = possub, pixmap = pixmap, contour = contour, 
        area = area, add.plot = add.plot)
    x <- as.numeric(dfxy[, xax])
    y <- as.numeric(dfxy[, yax])
    xykde = MASS::kde2d(x, y, lims = par("usr"))
    zlim = range(xykde$z, finite = TRUE)
    lev = seq(zlim[1], zlim[2], le = 8)
    lev = lev[2:7]
    contour(xykde, add = TRUE, levels = lev, 
        drawlabels = FALSE, ...)
    invisible(match.call())
}

PCAPlot2=function(dataset, annotation, subs = "Sex", main="")
{
	library(ade4)
	
	tmp = prepareDataset(dataset)
	
	sumsfilter = (apply(tmp, 1, sum)>0)
	tmp=tmp[sumsfilter,]
	
	tmp = data.matrix(tmp)
	tmp = t(tmp)
	tmp = tmp[row.names(tmp) %in% annotation$SampleID,]
		
	ids = c()
	for (k in row.names(tmp))
	{
		ids = append(ids, annotation[annotation$SampleID==k ,subs])
	}	
	ids = as.factor(ids)
	
	###anomalous, need GlobalEnv)
	assign("ids", ids, envir=.GlobalEnv)	
	assign("tmp", tmp, envir=.GlobalEnv)
	
	pca1  <- dudi.pca(tmp,  center = TRUE,  scale = FALSE, scan = FALSE)
	pca2  <- dudi.pca(tmp,  center = TRUE,  scale = TRUE,  scan = FALSE)
	
	colors = unique(unlist( lapply(dimnames(tmp)[[1]], collookup, annotation, column=subs) ))
		
	par(mfrow = c(2,2))
	mykde2d(pca1$li, col="gray90", lwd=0.8, lty=1)

	s.class(pca1$li, ids, cpoint = 1, sub = paste(main, "original", sep="\n"), csub=1, col=colors, add.p=TRUE, label="")
	if (subs=="SampleID")
	{
		### no legend, dots will have text.
	}
	else if (length(levels(ids))>10)
	{
		legend("bottomright", levels(ids), fill=colors, ncol=2, cex=0.5 )
	}
	else
	{
		legend("bottomright", levels(ids), fill=colors, ncol=1, cex=0.7 )
	}
	
	text(pca1$li, rownames(pca1$li), cex=0.5, col="gray30", pos=3)

	#s.corcircle(pca1$co, lab = dimnames(tmp)[[2]], full = TRUE, box = TRUE, sub = paste(main, "Contribution of variables to axes.",sep="\n"), csub=1)
	s.arrow(    pca1$c1, lab = dimnames(tmp)[[2]], sub = paste(main, "Importance of variable.", sep="\n") , csub=1, possub = "bottom")
	
	mykde2d(pca2$li, col="gray90", lwd=0.8, lty=1)
	s.class(pca2$li, ids, cpoint = 1, sub = paste(main,"scaled", sep="\n"), csub=1, col=colors, add.p=TRUE, label="")
	if (subs=="SampleID")
	{
		### no legend, dots will have text.
	}
	else if (length(levels(ids))>10)
	{
		legend("bottomright", levels(ids), fill=colors, ncol=2, cex=0.5 )
	}
	else
	{
		legend("bottomright", levels(ids), fill=colors, ncol=1, cex=0.7 )
	}
	text(pca2$li, rownames(pca2$li), cex=0.5, col="gray30", pos=3)
	
	s.corcircle(pca2$co, lab = dimnames(tmp)[[2]], full = TRUE, box = TRUE, sub = paste(main, "Contribution of variables to axes.",sep="\n"), csub=1, possub="bottom")
	#s.arrow(    pca2$c1, lab = dimnames(tmp)[[2]], sub = paste(main, "Importance of variable.", sep="\n") , csub=1)

	par(mfrow = c(1,1))
	
	return (pca1)
}


fixlabels=function(dataset)
{
	for (r in 1:nrow(dataset))
	{
		tmp = dataset[r,]
		
		label = tmp$label
		id = unlist(strsplit(tmp$taxonid, "\\."))
		parentlabel = ""

		while (label == "unclassified")
		{
			parentid = paste(id[1:(length(id)-1)], collapse=".", sep=".")	
			#cat(label, id, "[",parentid, "] -> ")
			parentlabel = dataset$label[dataset$taxonid==parentid]
			#cat(parentlabel)			
			id   = unlist(strsplit(parentid, "\\."))
			label = parentlabel
			#cat ("\n")
		}
#		print (tmp$label)
#		print (tmp$depth)
#		print (parentlabel)
#		
		if (parentlabel=="")
		{
			label = paste(tmp$label, " [", tmp$depth, "]", sep = "")
		}
		else
		{
			label = paste("Unclassified", " ", parentlabel, " [", tmp$depth, "]", sep = "")
		}
		
		dataset[r, "newlabel"] = label
		
	}
	
	dataset$label = dataset$newlabel
	return (dataset[,-ncol(dataset)])
}

prepareDataset = function(dataset)
{
	tmp = dataset[,4:ncol(dataset)]
	x = aggregate (tmp, list(dataset$label), sum)
	labs = x[,1]
	x = x[,-1]
	row.names(x) <- labs
	return (x)
}
ordering =function(x)
{

	return (sum(abs(diff(x))))
}

getPvals = function(x, mapping)
{

	tmp = x
	allheaders = unlist(mapping)
	map = names(tmp)

	for (k in 1:length(mapping))
	{
		if (length(mapping[[k]])==1)
		{
			for (id in mapping[[k]])
			{
				map = map[-which(map==id)]
				tmp = tmp[-which(names(tmp)==id)]
			}
		}
		else
		{
			for (id in mapping[[k]])
			{
				map[map==id]= k-1
			}
		}	
	}
	
	pval = NA
	
#	print ("...")
#	print (mapping)
#	print (tmp)
#	print (as.numeric(map))
#	print (sum(as.numeric(map)))
#	print (sd(as.numeric(map)))
	
	
	if (length(tmp)>0 && !is.na(sum(as.numeric(map))) && sd(as.numeric(map))>0  )
	{

		pval = kruskal.test(tmp, as.numeric(map) )$p.value 
		if (is.nan(pval))
		{
			return (NA)
		}
	}
	
	return (pval)
}

profilePlot = function(dataset, annotation, desc, column, topmost=30, perc=TRUE, orderclasses=c(), legendtitle = paste("OTUs per ", "taxon", sep=""), statsfilename = "")
{
	if (length(unique(annotation[,column]))==1)
	{
		return(data.frame())
	}
	
	tmp = prepareDataset(dataset)


	if (length(orderclasses)>0)
	{
		
	}
	else
	{
		orderclasses = unique(annotation[,column])
	}
	
	#print ("")
	#print (tmp)
	#print (orderclasses)
	
	######## normalize
	if (perc)
	{
		for (k in 1:ncol(tmp))
		{
			#print (".")
			### the dataset is per taxon level!
			### i.e. NOT cummulative  
			if (sum(tmp[,k])==0)
			{
			
			}
			else
			{
				tmp[,k] = tmp[,k]/sum(tmp[,k])
			}
				
		}
		lab = paste(round(seq(0,max(tmp*100), length.out=11),0), "%")
		at = seq(0,max(tmp), length.out=11)

	}
	
	######### subset SampleIDs per group
	x = aggregate(annotation[, c("SampleID")], list(annotation[,column]), function(x){ return(levels(x)[x]) } )
	
	tmpsd=tmp
	
	################
	#print ("...")
	#print (x)
	pvals = apply(tmp[,names(tmp)%in% unlist(x[,2])], 1, getPvals, x[,2])
	
	
	###############
	
	
	######## aggregate statistics
	Ns = list()
	
	for (G in orderclasses)
	{
		# print (G)
		group = x[x[,1]==G,]
		members=unlist(group[,2])
		Ns$new = length(members)

		if (length(members)>1)
		{
			tmp$new = apply( tmp[, members], 1, mean)
			tmpsd$new = apply(tmpsd[, members], 1, sd)
			tmpsd$new = 1.96 * tmpsd$new /sqrt(length(members))		
		}
		else
		{
			tmp$new=tmp[,members]
			tmpsd$new = 0
		}
		names(tmp)[names(tmp)=="new"]<-c(group[,1])	
		names(tmpsd)[names(tmpsd)=="new"]<-c(group[,1])
		names(Ns)[names(Ns)=="new"]<-c(group[,1])
		
	}	
	
	
	######## keep only the aggregated statistics
	tmp=tmp[,(ncol(tmp)-nrow(x)+1):ncol(tmp)]
	tmpsd=tmpsd[,(ncol(tmpsd)-nrow(x)+1):ncol(tmpsd)]
	
	
	###### table	
	otpt = data.frame(taxon = row.names(tmp) )
	#otpt$taxon = row.names(tmp)
	for (G in orderclasses)
	{
	
		#print (names(tmp))
		#print (names(tmpsd))
		#print (names(Ns))
		#print (dim(otpt))
		
		otpt.tmp = data.frame(tmp[,G], tmpsd[,G], Ns[[G]])
		names(otpt.tmp) = c( paste(G, c("mean", "SD", "n"), sep="_" ) )
		if (nrow(otpt)==0)
		{
			otpt = otpt.tmp
		}
		else
		{
			otpt = cbind(otpt, otpt.tmp)
		}
	}
	otpt$pval = pvals	
	pvals = otpt$pval[!is.na(otpt$pval)]
	if (length(pvals)>100)
	{
		U = fdrtool(pvals, statistic="pvalue", plot=FALSE, color.figure=FALSE, verbose=FALSE)
		otpt$pval2[!is.na(otpt$pval)] = U$pval
		otpt$qval[!is.na(otpt$pval)] = U$qval
		otpt$lfdr[!is.na(otpt$pval)] = U$lfdr
	}
	else
	{
		otpt$pval2= NA
		otpt$qval = NA
		otpt$lfdr = NA
	}
	
	if (statsfilename != "")
	{
		write.table(otpt, paste(statsfilename, column, "stats.tab.txt", sep="." ), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE )
	}
	
	######## order experiments (based on most differences between classes)
	
	fororder = 	apply(tmp, 1, ordering)
	tmp = tmp[order(fororder, decreasing=TRUE),]
	tmpsd = tmpsd[order(fororder, decreasing=TRUE),]
	
	#### plot only topX and accumulate the rest to form "Other"
	
	
	if (topmost < nrow(tmp))
	{
		tmpA = tmp[1:(topmost+1)    ,]
		tmpB = tmp[(topmost+1):nrow(tmp),]	
		
		tmpA[topmost+1,] = apply(tmpB, 2, sum)
		row.names(tmpA)[topmost+1] = paste("... [remaining" , nrow(tmpB),"]", sep=" ")
	
		tmpsdA = tmpsd[1:(topmost+1)    ,]
		
		tmpsdA[topmost+1,] = rep(0, ncol(tmpsdA))
		row.names(tmpsdA)[topmost+1] = paste("... [remaining" , nrow(tmpB),"]", sep=" ")
	
	}
	
	else
	{
		tmpA=tmp
		tmpsdA = tmpsd
	}

	
	### reverse order so that the MOST changed are at the top
	tmpA = tmpA[rev(1:nrow(tmpA)),]	
	tmpsdA = tmpsdA[rev(1:nrow(tmpsdA)),]	
	
	#print (dim(tmpA))
	#print (dim(tmpsdA))
	

	### flip so that data is grouped by taxa
	tmpA = t(data.matrix(tmpA))
	tmpsdA = t(data.matrix(tmpsdA))
	
	### reverse the order of classes so that it matches the ordercolumn from top to bottom 
	tmpA = tmpA[rev(1:nrow(tmpA)),]	
	tmpsdA = tmpsdA[rev(1:nrow(tmpsdA)),]	
	
	#cols = getColors(row.names(tmpA))
	
	barx = barplot(tmpA, 
			main=paste("Microbial profile on the taxonomic level of ",desc, sep=""),
			legend = row.names(tmpA), 
			axes=FALSE,
			cex.axis=0.8,
			beside=TRUE,
			horiz=TRUE, 
			las=1,
			xlim = c(0,2.5*max(at)),
			args.legend=list(cex=0.8, bg="ivory", title=legendtitle))
	
	
	error.bar(barx, tmpA,tmpsdA, horiz=TRUE, col="black", lwd=3)
	
	axis(1, at=at, lab=lab, las=3)
	ticks = (1:ncol(tmpA))*(nrow(tmpA)+1) -  nrow(tmpA)/2
	axis(2, pos=0, at=ticks, lab=rep("",length(ticks)), lwd=5, col="gray80")
	
	
	####repeat)
	par(new=TRUE)
	
	barx = barplot(tmpA, 
			main=paste("Microbial profile on the taxonomic level of ",desc, sep=""),
			legend = row.names(tmpA), 
			axes=FALSE,
			cex.axis=0.8,
			beside=TRUE,
			horiz=TRUE, 
			las=1,
			xlim = c(0,2.5*max(at)),
			args.legend=list(cex=0.8, bg="ivory", title=legendtitle))	
	
	
			
	invisible (list(means = tmpA, sds=tmpsdA))		
	
	#### pvalues and q values
	
	forlab = dimnames(tmpA)[[2]]
	forlab = data.frame(x = 1:length(forlab), taxon=forlab)
	pvals = otpt[otpt$taxon %in% forlab$taxon, c("taxon", "pval", "qval")]
	forlab = merge(forlab, pvals, by="taxon", all.x=FALSE, all.y=FALSE, sort=FALSE)
	
	if (length(dimnames(tmpA)[[2]]) == nrow(forlab))
	{
		labs = c(round(forlab$pval,4))
	}
	else
	{
		labs = c(1, round(forlab$pval,4))
	}
	labs = unlist(lapply(labs, function(x) { ifelse( !is.na(x) && x<0.05, paste("p = ", format(x, nsmall=4, scientific=FALSE, zero.print=TRUE)), "" ) } ))	
	axis(2, pos=max(at)+0.15, at=ticks, lab=labs, lwd=0.5, col="gray90", las=2, cex.axis=1)
	
	if (length(dimnames(tmpA)[[2]]) == nrow(forlab))
	{
		labs = c(round(forlab$qval,4))
	}
	else
	{
		labs = c(1, round(forlab$qval,4))
	}	
	labs = unlist(lapply(labs, function(x) { ifelse( !is.na(x) && x<0.05, paste("q = ", format(x, nsmall=4, scientific=FALSE, zero.print=TRUE)), "" ) } ))
	axis(4, pos=max(at)+0.15, at=ticks, lab=labs, lwd=0.5, col="gray90", las=2, cex.axis=1)


	
}

comparisonPlot = function(dataset, annotation, desc, column, columnb, funcname="diff", topmost=30, perc=TRUE, orderclasses=c(), legendtitle = paste(column), ordervars=c())
{
	if (ncol(dataset)<5)
	{
		return (data.frame())
	}
	
	#print (dim(dataset))
	
	tmp = prepareDataset(dataset)	
	func = get(funcname)
	
	if (length(orderclasses)>0)
	{
		
	}
	else
	{
		orderclasses = unique(annotation[,column])		
	}
	######## normalize
	if (perc)
	{
		for (k in 1:ncol(tmp))
		{
			### the dataset is per taxon level!
			### i.e. NOT cummulative  
			if (sum(tmp[,k])==0)
			{
				
			}
			else
			{
				tmp[,k] = tmp[,k]/sum(tmp[,k])
			}
		}	
	}
		
	######### subset SampleIDs per group
	x = aggregate(annotation[, c("SampleID")], list(annotation[,column ]), function(x){ return(levels(x)[x]) } )
	y = aggregate(annotation[, c("SampleID")], list(annotation[,columnb]), function(x){ return(levels(x)[x]) } )
	
	######## aggregate statistics
	for (G in orderclasses)
	{
		# print (G)
		groupA = x[x[,1]==G,]
		membersA=unlist(groupA[,2])
		#print (membersA)
		
		# for func:
		orderclassesB = sort(unique(annotation[,columnb]))
	
		if (length(orderclassesB)==1)
		{
			return (data.frame())
		}
		for (H in orderclassesB)
		{
			groupB = y[y[,1]==H,]
			membersB=unlist(groupB[,2])
			members = membersA[membersA %in% membersB]
			
			if (length(members)==0)
			{
				return (data.frame())
			}
			else if (length(members)>1)
			{
				tmp[,H] = apply( tmp[, members], 1, mean)
			}
			else
			{
				tmp[,H]=tmp[,members]
		
			}
			
		}

		### !!!!!!!!!!!!!!!!!!!!!!!!!		
		### in case of diff, default, the order needs to be reversed, i.e. diff (1,2) = 1

		tmp$new = apply(tmp[,rev(orderclassesB)],1, func)
		names(tmp)[ncol(tmp)]<-c(G)
		tmp = tmp[,!(names(tmp) %in% orderclassesB)]
					
	}	
	
	######## keep only the aggregated statistics
	#print (ncol(tmp))
	tokeep = (ncol(tmp)-nrow(x)+1):ncol(tmp)
	if (length(tokeep)>1)
	{
		tmp=tmp[,tokeep]	
	}
	else
	{
		OO = data.frame(tmp[,tokeep], tmp[,tokeep])
		names(OO) <- c(names(tmp)[tokeep], "removeme")
		row.names(OO) <- row.names(tmp)
		tmp=data.frame(OO)

	}
	

	######## order experiments (based on most differences between classes)
	
	#print (tmp)
	if (length(ordervars)==0)
	{
		fororder = 	apply(tmp, 1, ordering)
		tmp = tmp[order(fororder, decreasing=TRUE),]
	}
	else
	{
		tmp = tmp[ordervars,]
	}

	#print (tmp)

	
	#### plot only topX and accumulate the rest to form "Other"
	if (topmost < nrow(tmp))
	{
		tmpA = tmp[1:(topmost+1)    ,]
		tmpB = tmp[(topmost+1):nrow(tmp),]	
		
		tmpA[topmost+1,] = apply(tmpB, 2, sum)
		row.names(tmpA)[topmost+1] = paste("... [remaining" , nrow(tmpB),"]", sep=" ")	
	}
	
	else
	{
		tmpA=tmp
	}

	
	### reverse order so that the MOST changed are at the top
	tmpA = tmpA[rev(1:nrow(tmpA)),]	
	

	### flip so that data is grouped by taxa
	tmpA = t(data.matrix(tmpA))
	
	### reverse the order of classes so that it matches the ordercolumn from top to bottom 
	tmpA = tmpA[rev(1:nrow(tmpA)),]	
	
	#print (tmpA)
	
	#cols = getColors(row.names(tmpA))
	
	if (perc)
	{
	
		lab = paste(round(seq(min(tmp*100, -50),max(tmp*100, 50), length.out=11),0), "%")
		at = seq(min(tmp, -0.5),max(tmp, 0.5), length.out=11)
	}
	else
	{
		lab = paste(round(seq(min(tmp, -150),max(tmp, 150), length.out=11),0))
		at = seq(min(tmp, -150),max(tmp, 150), length.out=11)
	}
	
	
	#print(row.names(tmpA))
#	print (dim(tmpA))
#	tmpA=data.frame(tmpA[row.names(tmpA)!="removeme",])
#	print (tmpA)
#	print (dim(tmpA))
	barx = barplot(tmpA, 
			main=paste(funcname, " ( ", paste(orderclassesB, collapse=" - "), " )\n", desc, sep=""),
			legend = row.names(tmpA), 
			axes=FALSE,
			cex.axis=0.8,
			beside=TRUE,
			horiz=TRUE, 
			las=1,
			xlim = c(1.5*min(at),1.5*max(at)),
			args.legend=list(x="right", cex=0.8, bg="ivory", title=legendtitle))
	
	axis(1, at=at, lab=lab, las=3)
	ticks = (1:ncol(tmpA))*(nrow(tmpA)+1) -  nrow(tmpA)/2
	axis(2, pos=1.5*min(at), at=ticks, lab=rep("",length(ticks)), lwd=5, col="gray80")
	axis(2, pos=0 , at=ticks, lab=rep("",length(ticks)), lwd=5, col="gray80")
	
	print (names(tmpA))
	
	####repeat)
	par(new=TRUE)
	
	barx = barplot(tmpA, 
			main=paste(funcname, " ( ", paste(orderclassesB, collapse=" - "), " )\n", desc, sep=""),
			legend = row.names(tmpA), 
			axes=FALSE,
			cex.axis=0.8,
			beside=TRUE,
			horiz=TRUE, 
			las=1,
			xlim = c(1.5*min(at),1.5*max(at)),
			args.legend=list(x="right", cex=0.8, bg="ivory", title=legendtitle))	
			
	invisible (tmp)		
	
}



error.bar <- function(x, y, upper, lower=upper, arlength=0.025, horiz=TRUE, ... ){
	if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
	{
		stop("vectors must be same length")
	}
		
	if (horiz)
	{
		arrows(y, x, y+upper,x,  angle=90, code=2, length=arlength, ...)
	}
	else
	{	
		arrows(x,y+upper, x, y-lower, angle=90, code=3, length=arlength, ...)
	}
}

stackedBarPlot = function(dataset, annotation, desc, column, topmost=20, perc=TRUE, orderclasses=c())
{
	tmp = prepareDataset(dataset)

	if (length(orderclasses)>0)
	{
	}
	else
	{
		orderclasses = unique(annotation[,column])
		orderclasses = orderclasses[!is.na(orderclasses)]
		
		#print (orderclasses)
	}	

	#print (annotation[,column])
	#print (orderclasses)

	x = aggregate(annotation[, c("SampleID")], list(annotation[,column]), function(x){ return(levels(x)[x]) } )

	tmp$sum = apply(tmp, 1, sum)
	
	for (G in orderclasses)
	{
		
		group = x[x[,1]==G,]		
		members=unlist(group[,2])
		
		if (length(members)>1)
		{
			tmp$new = apply( tmp[, names(tmp) %in% members], 1, mean)
		}
		else
		{			
			tmp$new=tmp[, names(tmp) %in% members]
			tmp = tmp[, ! (names(tmp) %in% members) ]
		}
		names(tmp)[ncol(tmp)]<-c(group[,1])
	}

	tmp = tmp[order(tmp$sum, decreasing=TRUE),]
	
	#### plot only topX and accumulate the rest to form "Other"
	if (topmost < nrow(tmp))
	{

		tmpA = tmp[1:(topmost+1)    , (ncol(tmp)-nrow(x)):ncol(tmp)]
		tmpB = tmp[(topmost+1):nrow(tmp), (ncol(tmp)-nrow(x)):ncol(tmp)]				
		
		tmpA[topmost+1,] = apply(tmpB, 2, sum)
		row.names(tmpA)[topmost+1] = paste("... [remaining" , nrow(tmpB),"]", sep=" ")

	}
	
	else
	{
		tmpA=tmp[,(ncol(tmp)-nrow(x)):ncol(tmp)]
	}

	if (perc)
	{
		for (k in 2:ncol(tmpA))
		{
			tmpA[,k] = tmpA[,k]/sum(tmpA[,k])
		}
		lab = paste(seq(0,100, by=10), "%")
		at = seq(0,1, by=0.1)
		
	}

	tmpA = tmpA[rev(1:nrow(tmpA)),]	
	cols = getColors(row.names(tmpA))
	
	#print (max(tmpA[,-1]))
	
	barplot(data.matrix(tmpA[,-1]), 
				legend = row.names(tmpA), 
				main=paste("Relative OTU distribution\n", desc, sep=""),
				col=cols, axes=FALSE, 
				ylim = c(0, 1.5 * sum(tmpA[,2])), 
				args.legend=list(cex=0.8, bg="ivory", title=paste("Taxon", sep="")), las=2)
	axis(2, at=at, lab=lab, las=1)
	
}

getLevel = function(desc)
{
	tab = list()
	tab$Kingdom = 1
	tab$Domain  = 1
	tab$Phylum	= 2
	tab$Class	= 3
	tab$Order	= 4
	tab$Family	= 5
	tab$Genus	= 6
	tab$Species	= 7
	tab$Strain	= 8
	
	tab$kingdom = 1
	tab$domain  = 1
	tab$phylum	= 2
	tab$class	= 3
	tab$order	= 4
	tab$family	= 5
	tab$genus	= 6
	tab$species	= 7
	tab$strain	= 8
	
	tab[["level 1"]]=  20
	tab[["level 2"]]=  21
	tab[["level 3"]] = 22
	tab[["function"]]= 23
	
	return (tab[[desc]])
}

prepareTaxonomy=function(infile, annotation, tag = "")
{
	
	origdataset = read.table(infile, sep="\t", as.is=T, header = T)
	#origdataset = origdataset[, c(-4,-5, -(ncol(origdataset)))]
	names(origdataset)[1:3] <-c("depth","taxonid","label")
	
	if (tag != "")
	{
		origdataset$label=paste(origdataset$label, tag, sep="_")
	}
	origdataset = origdataset[, names(origdataset) %in% c("depth","taxonid","label", annotation$SampleID)]
	origdataset = fixlabels(origdataset)
	return(origdataset)
}

findRanges=function(input)
{
	for (k in names(input))
	{
		off = regexpr("ranges", k) 
		if (off >-1)
		{
			#print (k)
			#print (input[,k])
			input[,k] = as.numeric(input[,k])
			input[,k] = levels(cut(input[,k], 5))[cut(input[,k],5)]
			newname = substr(k, 1, off-3)
			names(input)[names(input)==k]=c(paste(newname, " [intervals]", sep=""))
			#print (input[,c(paste(newname, " [intervals]", sep=""))])
		}
	}
	return (input)
}

normalize = function(input)
{
	#print (unlist(input))
	if (is.numeric(input))
	{
		return (input/sum(input))
	}
	else
	{
		return (input)
	}
	
}
makeDefaultBatchOfPlots=function(annotationfilename, constaxonomyfilename, fileprefix="AllFeat", norm=FALSE, taxa = c("Phylum", "Class", "Order", "Family", "Genus"), origdataset=data.frame(), origannotation = data.frame())
{
		
# 	origannotation = annotationfile
#  	origdataset = constaxonomyfile

	if (nrow(origannotation)==0)
	{
		origannotation = read.csv(annotationfilename, as.is=TRUE, header=TRUE)
		origannotation$All = rep("All", nrow(origannotation))
		origannotation = origannotation[, which(names(origannotation)=="SampleID"):ncol(origannotation)]
		### this should allow for any sample name...
		origannotation$SampleID_input = origannotation$SampleID
		origannotation$SampleID = make.names(origannotation$SampleID)
		
		FIXORDER = data.frame(orig = unique(origannotation$SampleID) )		
		origannotation = aggregate(origannotation, by=list(origannotation$SampleID), function(x) { paste(unique(x), collapse=",", sep=",") } )[,-1]
	
		FIXORDER$new = origannotation$SampleID 
		FIXORDER$undo = order(FIXORDER$orig)
		FIXORDER$undone = FIXORDER$new[order(FIXORDER$undo)]
		
		origannotation = origannotation[order(FIXORDER$undo),]
		origannotation = findRanges(origannotation)
		
 	} 
 	
 	if (nrow(origdataset)==0)
 	{
		origdataset = prepareTaxonomy( constaxonomyfilename, origannotation )
 	}
  	
 	if (norm)
 	{
 		fileprefix = paste(fileprefix, "norm", sep="")	
 	}
 
	for (taxon in taxa)
	{
		cat (taxon, "\n", sep="")
		x = makeLevelDataset(origdataset, getLevel(taxon), cumulative=TRUE)	
		
		
		if (norm)
		{
			for (k in 4:ncol(x))
			{
				x[,k] = x[,k]/sum(x[,k])
			}

		}
		
		annotation = origannotation[origannotation$SampleID %in% names(x),]
		categories = names(annotation)
		
		Orel=c()
		Oabs=c()
		
		if (nrow(x)>2)
		{
			cat ("HEAT:")
			filename = paste(fileprefix, "TAX", taxon, "HEAT.pdf", sep=".")
			pdf (filename, paper="special", width= 10 + 0.1 * nrow(x), height= 10 + 0.1 * ncol(x))
			#pdf (filename, paper="special", width=max(20), height=max(20))
			
			for (category in categories)
			{
				cat ("\t", category, sep="")
				
				curanno = annotation[annotation[,category]!= "NA" & annotation[,category]!= "" & ! is.na(annotation[,category]), ]
			
				curdata = x[,names(x) %in% c("label", "depth", "taxonid", curanno$SampleID),]
				if (ncol(curdata)>3)
				{
					
				
					mainlabel = paste("Taxon: ", taxon, ", samples grouped by: ", category, sep="")
					
					if (norm)
					{
						makeHeatMap(curdata, curanno, category, mainlabel, perc=FALSE)
					}
					else
					{
						makeHeatMap(curdata, curanno, category, mainlabel, perc=TRUE)
						makeHeatMap(curdata, curanno, category, mainlabel, perc=FALSE)
					}
					
					
				}
	
			} 	
			dev.off()
			cat("\n")
		
		
			cat ("PCA:")
			filename = paste(fileprefix, "TAX", taxon, "PCA.pdf", sep=".")
			pdf (filename, paper="special", width=12, height=9)
			for (category in categories)
			{
				cat ("\t", category, sep="")
				
				curanno = annotation[annotation[,category]!= "NA" & annotation[,category]!= "" & ! is.na(annotation[,category]), ]
				curdata = x[,names(x) %in% c("label", "depth", "taxonid", curanno$SampleID),]
				if (ncol(curdata)>5)
				{
					#print (curdata)
					#print (dim(curdata))
					#print (names(curdata))
					mainlabel = paste("Taxon: ", taxon, ", samples grouped by: ", category, sep="")
					PCAPlot2(curdata, curanno, subs = category, main=mainlabel)
				}
			} 	
			dev.off()
			cat("\n")
		
			cat ("Stack:")
			filename = paste(fileprefix,"TAX", taxon, "STACK.pdf", sep=".")
			pdf (filename, paper="special", width=12, height=15)
			#layout( matrix(c(1,2), 1, 2, byrow = FALSE), heights=c(1), widths=c(0.75, 0.25))
			for (category in categories)
			{
				cat ("\t", category, sep="")
				
				#print (annotation)
				curanno = annotation[annotation[,category]!= "NA" & annotation[,category]!= "" & ! is.na(annotation[,category]), ]
				#print (curanno)
				curdata = x[,names(x) %in% c("label", "depth", "taxonid", curanno$SampleID),]
				#print (names(x))
				#print (names(curdata))
				
				if (ncol(curdata)>3)
				{
					mainlabel = paste("Taxon: ", taxon, ", samples grouped by: ", category, sep="")
					stackedBarPlot(curdata, curanno, mainlabel, category)
		
				} 	
			}
			dev.off()
			cat("\n")
		
		
		
			cat ("Profl.:")
			filename = paste(fileprefix,"TAX", taxon, "PROFILE.pdf", sep=".")
			filename2 = paste(fileprefix,"TAX", taxon, sep=".")
			pdf (filename, paper="special", width=12, height=20)
			par(mar=c(5,14,4,2))
			for (category in categories)
			{
				cat ("\t", category, sep="")
				
				curanno = annotation[annotation[,category]!= "NA" & annotation[,category]!= "" & ! is.na(annotation[,category]), ]
				curdata = x[,names(x) %in% c("label", "depth", "taxonid", curanno$SampleID),]
				
				
				if (ncol(curdata)>3)
				{
					mainlabel = paste("Taxon: ", taxon, ", samples grouped by: ", category, sep="")
					profilePlot(curdata, curanno, mainlabel, category, topmost=20, legendtitle=category, statsfilename = filename2)		
				}
			} 	
			par(mar=c(5, 4, 4, 2))
			dev.off()
			cat("\n")		
	
			cat ("Comparisons (Absolute):\n")
			pdf (paste(fileprefix,"TAX", taxon, "COMPabs.pdf", sep="."), paper="special", width=8, height=10)
			par(mar=c(5,14,4,2))
			
			for (g in categories)
			{
				tocompare = sort(unique(annotation[,g]))
				tocompare = tocompare[tocompare != ""]
				
				if (length(tocompare)==2)
				{
					for (h in categories)
					{
						if (h!=g)
						{
							per = sort(unique(annotation[,h]))
							per = per[per!=""]
							
							cat ("\tdiff ", g, " [", paste(tocompare, collapse="-", sep="") ,"] per ", h, " [", paste(per[1:min(5, length(per))], collapse=", ", sep="") , "]", sep="")
							
							curanno = annotation[annotation[,g]!= "NA" & annotation[,g]!= "" & ! is.na(annotation[,g]) & annotation[,h]!= "NA" & annotation[,h]!= "" & ! is.na(annotation[,h]), ]
							#curanno = annotation[annotation[,g] != "" & annotation[,h] != "", ]
							curdata = x[,names(x) %in% c("label", "depth", "taxonid", curanno$SampleID),]
							if (ncol(curdata)>3)
							{				
								dataset = comparisonPlot(curdata, curanno, paste("taxon:", taxon), h, g, funcname="diff", topmost=15,  perc=FALSE, legendtitle = paste(h), ordervars = c())
								
								if (nrow(dataset)==0)
								{
									cat("\t...skip\n")
								}
								else
								{
									cat("\t...ok\n")
									if (length(Orel)==0)
									{
										Orel = row.names(dataset)
									}
								}
							}	
							
						}
					}
				}
			}
			par(mar=c(5, 4, 4, 2))
			dev.off()
			cat("\n")
			
			cat ("Comparisons Relative:\n")
			pdf (paste(fileprefix,"TAX", taxon, "COMPrel.pdf", sep="."), paper="special", width=8, height=10)
			par(mar=c(5,14,4,2))	
			for (g in categories)
			{
				tocompare = sort(unique(annotation[,g]))
				tocompare = tocompare[tocompare != ""]
				
				if (length(tocompare)==2)
				{
					for (h in categories)
					{
						if (h!=g)
						{
							per = sort(unique(annotation[,h]))
							per = per[per!=""]
							
							curanno = annotation[annotation[,g]!= "NA" & annotation[,g]!= "" & ! is.na(annotation[,g]) & annotation[,h]!= "NA" & annotation[,h]!= "" & ! is.na(annotation[,h]), ]
							#curanno = annotation[annotation[,g] != "" & annotation[,h] != "", ]
							curdata = x[,names(x) %in% c("label", "depth", "taxonid", curanno$SampleID),]
							if (ncol(curdata)>3)
							{
								cat ("\tdiff ", g, " [", paste(tocompare, collapse="-", sep="") ,"] per ", h, " [", paste(per[1:min(5, length(per))], collapse=", ", sep="") , "]", sep="")
								dataset = comparisonPlot(curdata, curanno, paste("taxon:", taxon), h, g, funcname="diff", topmost=15,  perc=TRUE , legendtitle = paste(h), ordervars = c())
								
								if (nrow(dataset)==0)
								{
									cat("\t...skip\n")
								}
								else
								{
									cat("\t...ok\n")
									if (length(Orel)==0)
									{
										Orel = row.names(dataset)
									}
								}
							}
							
						}
					}
				}
			}
			par(mar=c(5, 4, 4, 2))
			dev.off()
			cat("\n")
		}	
		
	}
}	

makeDebugBatchOfPlots=function(annotationfilename, constaxonomyfilename, fileprefix="AllFeat", norm=FALSE, taxa = c("Phylum", "Class", "Order", "Family", "Genus"), origdataset=data.frame(), origannotation = data.frame())
{

	if (nrow(origannotation)==0)
	{
		origannotation = read.csv(annotationfilename, as.is=TRUE, header=TRUE)
		origannotation = origannotation[, which(names(origannotation)=="SampleID"):ncol(origannotation)]
		origannotation = aggregate(origannotation, by=list(origannotation$SampleID), function(x) { paste(unique(x), collapse=",", sep=",") } )[,-1]
		origannotation = findRanges(origannotation)
 	} 
 	
 	if (nrow(origdataset)==0)
 	{
		origdataset = prepareTaxonomy( constaxonomyfilename, origannotation )
 	}
 
 	if (norm)
 	{
 		fileprefix = paste(fileprefix, "norm", sep="")	
 	}
 
	for (taxon in taxa)
	{
		cat (taxon, "\n", sep="")
		x = makeLevelDataset(origdataset, getLevel(taxon), cumulative=TRUE)	
		
		if (norm)
		{
			for (k in 4:ncol(x))
			{
				x[,k] = x[,k]/sum(x[,k])
			}

		}

		annotation = origannotation[origannotation$SampleID %in% names(x),]
		
		
						
	}
}	





