########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2012 Sebastian Szpakowski, J.Craig Venter Institute.
########################################################################################



installme=function()
{
	### verify that all required packages are installed
	m <- getCRANmirrors(all = FALSE, local.only = FALSE)
	m$counter = 1:nrow(m)
	res = min(m$counter[m$CountryCode=="us"])
	URL <- m[res, "URL"]
	repos <- getOption("repos")
	repos["CRAN"] <- gsub("/$", "", URL[1L])
	options(repos = repos)
	install.packages("ade4", dependencies=FALSE)
	install.packages("RColorBrewer", dependencies=FALSE)
	install.packages("fdrtool", dependencies=FALSE)
}
done=FALSE
try(
	{	
		library(ade4)
		library(RColorBrewer)
		library(fdrtool)
		done=TRUE
	},
silent=TRUE)
if (! done)
{
	installme()
}
