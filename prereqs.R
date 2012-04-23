########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2012 J.Craig Venter Institute.
########################################################################################

installme=function()
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
	install.packages("ade4", dependencies=T)
	install.packages("RColorBrewer", dependencies=T)
	install.packages("fdrtool", dependencies=T)
	install.packages("randomForest", dependencies=T)
	install.packages("vegan", dependencies=T)
	install.packages("gplots", dependencies=T)
}
done=FALSE
try(
	{	
		library(fdrtool)
		library(RColorBrewer)
		library(ade4)
		library(vegan)
		library(randomForest)
		library(gplots)
		done=TRUE
	},
silent=TRUE)
if (! done)
{
	installme()
}
