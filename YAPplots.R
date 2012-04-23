########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2012 J.Craig Venter Institute.
########################################################################################

library(fdrtool)
library(RColorBrewer)
library(ade4)
library(randomForest)
library(vegan)

########################################################################################

parseMOTHURConsensusTax=function(filename, tag="")
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

parseMGRASTdataset = function(filename)
{
	datain = read.table(filename, as.is=T, sep="\t", header=T, check.names=F, comment.char="", strip.white=T)
}

parseMETAREPdataset = function(filename, skip=56)
{
	datain = read.table(filename, as.is=TRUE, sep="\t", header = TRUE, skip=skip)
}

########################################################################################
