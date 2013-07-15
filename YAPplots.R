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

FFP.mothur=function(infile, annotation, tag="")
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


openAnno = function(annotation, as.is=TRUE)
{
	anno = read.csv(annotation, as.is=as.is, header=TRUE)
	anno = anno[, which(names(anno)=="SampleID"):ncol(anno)]
	anno = aggregate(anno, by=list(anno$SampleID), function(x) { paste(unique(x), collapse=",", sep=",") } )[,-1]
	samplenames = anno[,"SampleID"]
	row.names(anno) <- samplenames
	return (anno)
}

taxoLogic=function(x)
{
	initlength = length(x)	
	last = x[initlength]
	
	#print (x)
	x = x[regexpr("classified",x) == - 1 ]
	x = x[regexpr("sub_",x) == - 1 & regexpr("_of_",x) == - 1]
	#print (x)
	
	### if the last item is "unclassified" change it to more meaningful label
	### otherwise keep it
	if (length(x)!= initlength & x[length(x)] != last)
	{
		last = x[length(x)]
		lastname = names(last)
		label = paste("[", lastname," ", last,"]",  sep="")
	}
	else
	{
		label = x[length(x)]
	}
	
#	if (regexpr("sub-", label)>-1)
#	{
#		label = substr(label, regexpr("_of_", label)+4, nchar(label))
#	}
	
	return(label)
}

fixLabels=function(datain)
{
	x = names(datain)
	
	x = x[-1]
	x = x[-length(x)]
	x = x[-length(x)]
	
	x = rev(x)

	for (rankid in x)
	{
		tmp = datain[,2:which (names(datain)==rankid)]
		#print (ncol(tmp))
		if (!is.null(ncol(tmp)))
		{
			old = datain[,rankid]
			datain[,rankid] = apply(tmp, 1, taxoLogic)	
			new = datain[,rankid]
			#print (cbind(old, new))
		}

	}
	return (datain)
	
}

assemblePhyloseq=function(path, annotation, distance, relabelunclassified=TRUE)
{
	cat(".")
	
	files = dir(path = path, pattern=glob2rx("*.*"))
	OT = files[regexpr("otutable$", files )>-1 & regexpr(paste(distance), files )>-1]
	TR = files[regexpr("ph$", files )>-1 & regexpr(paste(distance), files )>-1] 
	TA= files[regexpr("phylotax$", files )>-1 & regexpr(paste(distance), files )>-1]
	
	print (OT)
	print (TR)
	print (TA)
	
	cat(".")
	### create otu table
	### VERY slow! tapply,
	### instead use python ~/DEVEL/YAP/scripts/OTUtableMaker.py --help 
	### x = import_mothur(LL, mothur_group_file = GG, mothur_tree_file=NULL, cutoff = dist)
		
	otutab = read.table(paste(path, OT, sep="/"), sep="\t", header = TRUE)
	otumat = as.matrix(otutab[,-1])
	colnames(otumat) = names(otutab)[-1]
	rownames(otumat) = otutab[,1]
	
	myotus = otu_table(otumat, taxa_are_rows=TRUE)
		
	### established taxonomy
	cat(".")
	taxa = read.table(paste(path, TA, sep="/"), as.is=T, sep="\t", header=TRUE)
	OTUIDS=taxa[,1]
	TAXIDS=names(taxa)
	taxa$Species[is.na(taxa$Species)] = OTUIDS[is.na(taxa$Species)]
	
	if (relabelunclassified)
	{
		taxa = fixLabels(taxa)
	}
	
	taxa = as.matrix(taxa[,-1])
	row.names(taxa)=OTUIDS
		
	tokeep = apply(taxa, 2, function(x) { sum(is.na(x))!=length(x) } )  
	taxa = taxa[,tokeep]
	#print (taxa[1:10,])
	taxa = tax_table(taxa)
		
	### sample information
	cat(".")
	anno = openAnno(annotation)
	anno = sample_data(anno)
	
	
	### attach the tree
	if ( length(TR) >0 )
	{
		cat(".")
		TR = read.tree(paste(path, TR, sep="/"))
		TR = root (TR, TR$tip.label[regexpr("thermococcus_SILVA_ABSV01001619", TR$tip.label)>-1], resolve.root=TRUE)
		x = phyloseq(myotus, anno, taxa, TR)
	}
	else
	{	
		cat("x")
		x = phyloseq(myotus, anno, taxa)
	}
	cat ("\n")
	print (x)
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

makeTREES=function( datain, taxa = c("Phylum", "Class"), categories = names(sample_data(datain)) )
{
	pdf("TREES.pdf", height=10, width=10)
	for (k in categories )
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
makeHEATMAPS=function(datain, taxa=c("Phylum", "Class"), dists=c("none", "bray", "jaccard", "unifrac"), methods=c("DCA", "RDA", "NMDS", "PCoA", "MDS"), categories = names(sample_data(datain)))
{
	pdf("HEATMAPS.pdf", height=10, width=10)
	for (k in categories)
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

makeCCA=function(datain, shape="status", color="status", minotusize = 2, taxa=c("Phylum", "Class"),  subsets = names(sample_data(datain)), constraints = names(sample_data(datain)), height=20, width=10 )
{
	
	
	availablelevels = getSafeTaxa(datain)
	labS = ifelse (length(subsets)<=2, paste(subsets, collapse="-", sep="-"), length(subsets))
	labC = ifelse (length(constraints)<=2, paste(constraints, collapse="-", sep="-"), length(constraints))
	labT = ifelse (length(taxa)<=2, paste(taxa, collapse="-", sep="-"), length(taxa))
	
	filename = paste("ORDconst", labS, labT, labC, minotusize, color, shape, "pdf", sep="." )	
	pdf(filename, height=height, width=width)
	
	for (k in subsets)
	{
		cat ("__________________________\n")
		cat(k, ":\t", sep="")
		for ( taxon in taxa)
		{
			cat (taxon, "\n", sep="")
			
			tokeep= sample_data(datain)[,k]
			tokeep = tokeep[tokeep != ""]
			tmp.samples = prune_samples(sample_names(tokeep), x)
			#print (tmp.samples)
			
			### remove empty and singletons taxa for given subset samples
			tokeep = taxa_sums(tmp.samples) 
			tokeep = tokeep[tokeep>=minotusize]
			tmp.samples.nonempty = prune_taxa(names(tokeep), tmp.samples)
			print (tmp.samples.nonempty)
			
			
			
			# for some reason the object needs to be accessible globally...
			tmp = tmp.samples.nonempty
			assign( "tmp", tmp, envir = .GlobalEnv)
			print (tmp)
			
			### constrained correspondence analysis
			for (const in constraints)
			{
				if (k != const)
				{
					
					myform = paste("tmp ~ ", const, sep="" )
					cat(myform)
					tmp.ORD <- ordinate(as.formula(myform), "CCA")
					
					desc = paste("Correspondence Analysis, constrained by ", const, sep="" )
					title = paste("All taxa\nSubset: \"", k, "\" with n=",nsamples(tmp), "\n", desc, sep="")
					
					p1 = plot_ordination(tmp, tmp.ORD, type="samples", shape=shape, color=color, label=k, title = title )	
					#p1 = p1  + geom_point(size=4) + stat_density2d(aes(alpha=..density..), geom="raster", contour=FALSE) +  geom_line() + geom_point(size=3)
					
					p2 = plot_ordination(tmp, tmp.ORD, type="scree", color = taxon, title = title )	
					p2 = p2 + theme(legend.position="none")
					
					multiplot(p1, p2, cols=1)
					
					
					p = plot_ordination(tmp, tmp.ORD, type="taxa", color=taxon)
					
					wrappingform = as.formula(paste("~", availablelevels[which(availablelevels==taxon)-1]))
					print (p + stat_density2d(aes(alpha=..density..), geom="raster", contour=FALSE, drop=TRUE) +  facet_wrap(wrappingform, ncol=3))
					
					cat("\t")	
					
				}	
				
			}
			
		}	
		cat("\n")
	}
	dev.off()
}

makeORDINATION=function(datain, shape="status", color="status", minotusize = 2, taxa=c("Phylum", "Class"),  subsets = names(sample_data(datain)), methods = c("CCA", "DCA", "DPCoA" , "RDA", "NMDS", "PCoA"), height=20, width=10 )
{
	availablelevels = getSafeTaxa(datain)
	descs = list()
	descs["CCA"] = "(unconstrained) Correspondence Analysis"
	descs["DCA"] = "Detrended Correspondence Analysis"
	descs["DPCoA"] = "Double Principal Coordinate Analysis"
	descs["RDA"] = "Redundancy Analysis / Principal Component Analysis"
	descs["NMDS"] = "Non-Metric Multidimensional Scaling "
	descs["PCoA"] = "Principal Coordinate Analysis"
	
	
	labS = ifelse (length(subsets)<=2, paste(subsets, collapse="-", sep="-"), length(subsets))
	labT = ifelse (length(taxa)<=2, paste(taxa, collapse="-", sep="-"), length(taxa))
	labM = ifelse (length(methods)<=2, paste(methods, collapse="-", sep="-"), length(methods))
	
	filename = paste("ORDuncon.", labS, labT, labM, minotusize, color, shape, "pdf", sep="." )	
	pdf(filename, height=height, width=width)
	
	
	for (k in subsets	)
	{		
		cat ("__________________________\n")
		cat(k, ":\t", sep="")
		for ( taxon in taxa)
		{
			cat (taxon, "\n", sep="")
			#print (x)
			
			tokeep= sample_data(datain)[,k]
			tokeep = tokeep[tokeep != ""]
			tmp.samples = prune_samples(sample_names(tokeep), x)
			#print (tmp.samples)
			
			### remove empty and singletons taxa for given subset samples
			tokeep = taxa_sums(tmp.samples) 
			tokeep = tokeep[tokeep>=minotusize]
			tmp.samples.nonempty = prune_taxa(names(tokeep), tmp.samples)
			#print (tmp.samples.nonempty)
			
			# for some reason the object needs to be accessible globally...
			tmp = tmp.samples.nonempty
			assign( "tmp", tmp, envir = .GlobalEnv)
			print (tmp)
			
			### unconstrained  analyses	
			
			for (m in methods)
			{
				tmp.ORD <- ordinate(tmp, m)
				cat(m)
				
				desc=descs[names(descs)==m]
				title = paste("All taxa\nSubset: \"", k, "\" with n=",nsamples(tmp), "\n", desc, sep="")
				
				p1 = plot_ordination(tmp, tmp.ORD, type="samples", shape=shape, color=color, label="SampleID", title = title )	
				p1 = p1 + geom_line() + geom_point(size=5) 
				
				p2 = plot_ordination(tmp, tmp.ORD, type="taxa", color = taxon, title = title )	
				
				multiplot(p1, p2, cols=1)
				
				p = plot_ordination(tmp, tmp.ORD, type="taxa", color=taxon) 
				wrappingform = as.formula(paste("~", availablelevels[which(availablelevels==taxon)-1]))
				
				
				print (p + facet_wrap(wrappingform, ncol=3))
				
				cat("\t")
				
				
#				p <- plot_ordination(tmp, tmp.ORD, type="split", color=k, label=k, title = title )	
#				print (p)
#				
#				p <- plot_ordination(tmp, tmp.ORD, type="samples", color=k, label=k, title = title ) + geom_line() + geom_point(size=5)		
#				print (p)
				
				# # p <- plot_ordination(tmp, tmp.ORD, type="species", color=taxon, label=NULL, title = title ) 
				# # print (p)
				
				# wrappingform = as.formula(paste("~", availablelevels[which(availablelevels==taxon)-1]))
				
				# p <- plot_ordination(tmp, tmp.ORD, type="taxa", color=taxon) 
				# print (p + facet_wrap(wrappingform, nrow=3))
				
				#print (ggplot(p$data, aes(x=CA1, y=CA2, color=Phylum)) + stat_density() )
				
				###for the purpose of plotting the contour remove singleton classes
				# xlim_ = range(p$data$CA1)
				# ylim_ = range(p$data$CA2)
				# p$data = p$data[p$data[, taxon] %in% names(table(p$data[, taxon])[table(p$data[,taxon])>1]), ]
				
				
				# print (ggplot(p$data, aes(x=CA1, y=CA2, color=Phylum)) + xlim(xlim_) + ylim(ylim_) + geom_density2d() + facet_wrap(~Phylum, nrow=3) )
				# # print (ggplot(p$data, aes(x=CA1, y=CA2, color=Phylum)) + stat_binhex(bins=100, na.rm=T) )
				#print (ggplot(p$data, aes(x=CA1, y=CA2, color=Phylum)) + stat_binhex(bins=100, na.rm=T) + facet_wrap(~Phylum, nrow=3))
				
				
				# tmp.ORD <- ordinate(tmp, "DCA")
				# desc="Detrended Correspondence Analysis"
				# title = paste("All taxa\nSubset: \"", k, "\" with n=",nsamples(tmp), "\n", desc, sep="")
				
				# p <- plot_ordination(tmp, tmp.ORD, type="split", color=k, label=k, title = title )	
				# print (p)
				
				# p <- plot_ordination(tmp, tmp.ORD, type="taxa", color=taxon) 
				# print (p + facet_wrap(wrappingform, nrow=3))
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

#annotation = "/usr/local/projects/T1D-CNMC/sszpakow/BATCH_01_16S/annotation3.csv"
#distance = 0.03
#x = assemblePhyloseq("./", annotation, distance)

################################################################################


#makeTREES(x)
#makeHEATMAPS(x)


########################################################################################
