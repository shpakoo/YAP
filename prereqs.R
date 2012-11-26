########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2012 J.Craig Venter Institute.
########################################################################################

init = function()
{
	### verify that all required packages are installed
	m <- getCRANmirrors(all = FALSE, local.only = FALSE)
	m$counter = 1:nrow(m)
	res = min(m$counter[m$CountryCode=="us"])+1
	URL <- m[res, "URL"]
	repos <- getOption("repos")
	cat(URL, "\n")
	repos["CRAN"] <- gsub("/$", "", URL[1L])
	options(repos = repos)
	
	source("http://bioconductor.org/biocLite.R")
	
	tryCatch(biocLite("BiocUpgrade"), error=function(e){})
}


init()
cran   = c("ade4", "RColorBrewer", "fdrtool", "randomForest", "vegan", "gplots", "multicore", "MASS")
bioc = c("Biostrings", "ggplot2", "grid", "scales")

for (p in cran)
{	
	pkg<-p
	tryCatch(library(pkg, character.only = TRUE, quietly=TRUE), error = function(e){install.packages(p, dependencies=T)})
}

for (p in bioc)
{
	pkg<-p
	tryCatch(library(pkg, character.only = TRUE, quietly=TRUE), error = function(e){biocLite(p)})
}

for (p in c(cran, bioc))
{
	pkg<-p
	tryCatch(library(pkg, character.only = TRUE), error = function(e){print (e)})
}

cat("unless you see errors above, this should have worked.\n")
