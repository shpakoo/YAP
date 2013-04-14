########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2013 Sebastian Szpakowski
########################################################################################

library(ggplot2)
library(reshape)
library(phyloseq)
library(fdrtool)
library(RColorBrewer)
library(ade4)
library(randomForest)
library(vegan)


########################################################################################

FFP.mothur=function(infile, tag="")
{

	origdataset = read.table(infile, sep="\t", as.is=T, header=T)
	names(origdataset)[1:3] <-c("depth","taxonid","label")
	if (tag != "")
	{
		origdataset$label=paste(origdataset$label, tag, sep="_")
	}
	origdataset = origdataset[, names(origdataset) %in% c("depth","taxonid","label", annotation$SampleID)]
	origdataset = fixlabels(origdataset)
	return(origdataset)
	
}

FFP.mgrast = function(infile, column = "abundance", mapping=data.frame())
{
	dataset = read.table(infile, as.is=T, sep="\t", header=T, check.names=F, comment.char="", strip.white=T)
	
	#### detect usable columns
	tmp = range(which(names(dataset) %in% names(getAllLevels())) )
	start = tmp[1]
	end = tmp[2]
	labs = dataset[start:end]
	otpt = data.frame(stringsAsFactors = FALSE )
	
	### make a table for merging containing depth, taxonid, and label	
	for (taxon in names(labs))
	{
		level = getLevel(taxon)
		tmp = unique(labs[,taxon])
		
		tmpfix = fixMGRastTaxonomy(tmp)
		tmpfix =  paste("[", level,"] ", tmpfix,  sep="")
		
		tmp = paste("[", level,"] ", tmp,  sep="")
		tmp = data.frame(depth = rep(level, length(tmp)), taxonid = tmpfix, label= tmp, stringsAsFactors = FALSE  )
		if (nrow(otpt)==0)
		{
			otpt = tmp  
		}
		else
		{
			otpt = rbind(otpt, tmp)
		}
		
		
	}	
	
	for (m in unique(dataset$metagenome))
	{
		tmp = dataset[dataset$metagenome == m,start:ncol(dataset)]
		
		
		cat(m , sep="")
		experiment = m
		if (nrow(mapping)>0)
		{
			experiment = mapping$mappedid[mapping$origid==m]
			cat("->", experiment, sep="")
		}
		cat(", ")
		### skip things not in the mapping if supplied
		
		onesample = data.frame(stringsAsFactors = FALSE )
		
		if (length(experiment)>0)
		{
			for (taxon in names(tmp))
			{
				
				level = getLevel(taxon)
				if (! is.null(level))
				{
					cat(taxon, sep="")	
					cat("[", level, "], ", sep="")
					
					tmptmp = aggregate(tmp[,column], list(tmp[,taxon]), sum)
					names(tmptmp)=c("label", experiment)
					tmptmp$label = paste("[", level,"] ",tmptmp$label,  sep="")
					
					if (nrow(onesample)==0)
					{
						onesample = tmptmp
					}
					else
					{
						onesample = rbind(onesample, tmptmp)
					}
				}
			}			
			otpt = merge(otpt, onesample, by="label", all.x=T, all.y=T, sort=F)
		}
		else
		{
			cat(" skip ")	
		}
		cat("\n")
	}
	
	for (k in 4:ncol(otpt))
	{
		otpt[is.na(otpt[,k]),k] = 0
	}	
	otpt$label = otpt$taxonid
	otpt$taxonid = rep(0, nrow(otpt))
	return (otpt[,c("depth", "taxonid", "label", names(otpt)[4:ncol(otpt)])])
}

FFP.metarep = function(filename, skip=56)
{
	datain = read.table(filename, as.is=TRUE, sep="\t", header = TRUE, skip=skip)
}


###http://wiki.stdout.org/rcookbook/Graphs/Multiple%20graphs%20on%20one%20page%20(ggplot2)/
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
	require(grid)
	
	# Make a list from the ... arguments and plotlist
	plots <- c(list(...), plotlist)
	
	numPlots = length(plots)
	
	# If layout is NULL, then use 'cols' to determine layout
	if (is.null(layout)) {
		# Make the panel
		# ncol: Number of columns of plots
		# nrow: Number of rows needed, calculated from # of cols
		layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
				ncol = cols, nrow = ceiling(numPlots/cols))
	}
	
	if (numPlots==1) {
		print(plots[[1]])
		
	} else {
		# Set up the page
		grid.newpage()
		pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
		
		# Make each plot, in the correct location
		for (i in 1:numPlots) {
			# Get the i,j matrix positions of the regions that contain this subplot
			matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
			
			print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
							layout.pos.col = matchidx$col))
		}
	}
}

assemblePhyloseq=function(path, annotation, distance)
{
	cat(".")
	
	files = dir(path = path, pattern=glob2rx("*.*"))
	OT = files[regexpr("otutable$", files )>-1 & regexpr(paste(distance), files )>-1]
	TR = files[regexpr("ph$", files )>-1 & regexpr(paste(distance), files )>-1] 
	TA= files[regexpr("phylotax$", files )>-1 & regexpr(paste(distance), files )>-1]
	
	cat(".")
	### create otu table
	
	### VERY slow! tapply,
	### instead use python ~/DEVEL/YAP/scripts/OTUtableMaker.py --help 
	### x = import_mothur(LL, mothur_group_file = GG, mothur_tree_file=NULL, cutoff = dist)
		
	otutab = read.table(OT, sep="\t", header = TRUE)
	otumat = as.matrix(otutab[,-1])
	colnames(otumat) = names(otutab)[-1]
	rownames(otumat) = otutab[,1]
	
	
	myotus = otu_table(otumat, taxa_are_rows=TRUE)
		
	### established taxonomy
	cat(".")
	taxa = read.table(TA, as.is=T, sep="\t", header=TRUE)
	OTUIDS=taxa[,1]
	TAXIDS=names(taxa)
	taxa$Species[is.na(taxa$Species)] = OTUIDS[is.na(taxa$Species)]
	taxa = as.matrix(taxa[,-1])
	row.names(taxa)=OTUIDS
		
	tokeep = apply(taxa, 2, function(x) { sum(is.na(x))!=length(x) } )  
	taxa = taxa[,tokeep]
	taxa = tax_table(taxa)
		
	### sample information
	cat(".")
	anno = read.csv(annotation, as.is=TRUE, header=TRUE)
	anno = anno[, which(names(anno)=="SampleID"):ncol(anno)]
	anno = aggregate(anno, by=list(anno$SampleID), function(x) { paste(unique(x), collapse=",", sep=",") } )[,-1]
	
	samplenames = anno[,"SampleID"]
	row.names(anno) <- samplenames
	anno = sample_data(anno)
	
	### attach the tree
	cat(".")
	TR = read.tree(TR)
	TR = root (TR, TR$tip.label[regexpr("e_coli", TR$tip.label)>-1], resolve.root=TRUE)

	x = phyloseq(myotus, anno, taxa, TR)
	
	cat("\n")
	return (x)
	
}	
getSafeTaxa=function(PHSQ)
{
	levlist=c()
	cursum=0
	for (taxon in rank_names(x))
	{
		tryCatch(	
				{
					
					tmp = get_taxa_unique(x, taxon)
					if (length(tmp) >1 )
					{ 
						levlist = append(levlist, taxon)
					}	
					
				}, 
				error = function(err)
				{
					cat(paste("No data for ",taxon, "\n", sep=""))
				}
		)
	}
	return (levlist)
}

makeTREES=function( datain, taxa = c("Phylum", "Class") )
{
	pdf("TREES.pdf", height=10, width=10)
	for (k in names(sample_data(datain)))
	{
		
		cat(k, ":\t", sep="")
		for ( taxon in taxa)
		{
			cat (taxon, " ", sep="")	
			tokeep= sample_data(datain)[,k]
			tokeep = tokeep[tokeep != ""]
			tmp = prune_samples(sample_names(tokeep), x)
			tmp_lev = tax_glom(tmp, taxon)
			title = paste("Level ", taxon, "\nSubset: \"", k, "\" with n=",nsamples(tmp_lev), sep="")
			p = plot_tree(tmp_lev, size="abundance", color = k, label.tips=taxon, title=title) 
			print (p)
		}
		cat("\n")
		
	}
	dev.off()
	
}
makeHEATMAPS=function(datain, taxa=c("Phylum", "Class"), dists=c("none", "bray", "jaccard", "unifrac"), methods=c("DCA", "RDA", "NMDS", "PCoA", "MDS"))
{
	pdf("HEATMAPS.pdf", height=10, width=10)
	for (k in names(sample_data(datain)))
	{
		cat(k, ":\t", sep="")
		for ( taxon in taxa)
		{
			cat (taxon, " ", sep="")	
			tokeep= sample_data(datain)[,k]
			tokeep = tokeep[tokeep != ""]
			tmp = prune_samples(sample_names(tokeep), datain)
			tmp_lev = tax_glom(tmp, taxon)
			
			for (method in methods)	
			{
				for (dist in dists)
				{
					if (method %in% c("DCA", "RDA") & dist=="none" | method %in% c("NMDS", "PCoA", "MDS") & dist!="none" )
					{
						desc = paste(method, ifelse(dist=="none", "-", dist), sep=" / ")	
						title = paste("Level ", taxon, "\nSubset: \"", k, "\" with n=",nsamples(tmp_lev), "\n", desc, sep="")
						#print (desc)
						p = plot_heatmap(tmp_lev, method, dist, k, taxon, title=title, low="#FFFFCC", high="#000033", na.value="white") 
						print (p)
					}
				}
			}
			
		}
		cat("\n")
	}
	dev.off()
}

################################################################################

### ln -s ../*OUTPUT_6*/*.list ./
### ln -s ../*OUTPUT_6*/*.group ./
### ln -s ../*OUTPUT_6*/*.phylotax ./
### ln -s ../*OUTPUT_8*/*PHYLOTREES*/*.ph ./
### ~/DEVEL/YAP/scripts/OTUtableMaker.py --help 

annotation = "/usr/local/projects/T1D-CNMC/sszpakow/BATCH_01_16S/annotation3.csv"
distance = 0.03
x = assemblePhyloseq("./", annotation, distance)

################################################################################


makeTREES(x)
makeHEATMAPS(x)


########################################################################################
