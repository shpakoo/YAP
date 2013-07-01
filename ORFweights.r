########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2013 Sebastian Szpakowski
########################################################################################

library(ggplot2)
library(grid)
library(scales)

multiplot <- function(..., plotlist=NULL, cols=1, layout=NULL)
{
	require(grid)
	
	# Make a list from the ... arguments and plotlist
	plots <- c(list(...), plotlist)
	numPlots = length(plots)
	
	# If layout is NULL, then use 'cols' to determine layout
	if (is.null(layout)) 
	{
		layout = grid.layout(ceiling(numPlots/cols), cols)
	}
	
	# create a panel that will be used to map ids later
	
	panel = matrix(1:(layout$nrow * layout$ncol), nrow = layout$nrow, ncol = layout$ncol, byrow=TRUE )
	
	if (numPlots==1) 
	{
		print(plots[[1]])
	} 
	else 
	{
		# Set up the page
		grid.newpage()
		pushViewport(viewport(layout = layout))
			
		for (i in 1:numPlots) 
		{
			matching <- as.data.frame(which(panel == i, arr.ind = TRUE))
			
			print ( plots[[i]], vp = viewport(	layout.pos.row = matching$row, 
												layout.pos.col = matching$col)
			)
		}
	}
}

loadWeights <- function ()
{
	x = data.frame()
	for (file in dir (pattern=  glob2rx("*.weight")))
	{
		newfile = unlist(strsplit(file, "\\."))
		newfile = newfile[-length(newfile)]
		newfile = paste(paste(newfile, collapse=".", sep=""), "normweight", sep= ".")
		print (newfile)
		tmp = read.table(file, sep="\t", header = FALSE)
		names(tmp)<-c("orfid", "coverage")
		tmp$scaled = tmp$coverage / median(tmp$coverage)
		tmp$scaledlog2 = log2(tmp$coverage / median(tmp$coverage))
		###
		write.table(tmp, newfile, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
		###
		
		tmp$file = rep(file, nrow(tmp))
		cat(file, median(tmp$coverage), mean(tmp$coverage), "\n", sep= "\t")
		if (nrow(x)==0)
		{
			x = tmp
		}
		else
		{
			x  = rbind(tmp, x)
		}
	}
	return (x)
}





################################################################################

 global = loadWeights()

 files = unique(global$file)

####### transforms

 y_trans_log10 = scale_y_log10(	breaks = trans_breaks("log10", function(x) 10^x),
 								labels = trans_format("log10", math_format(10^.x)))	

 x_trans_log10 = scale_x_log10(	breaks = trans_breaks("log10", function(x) 10^x),
 								labels = trans_format("log10", math_format(10^.x)))
	
	
 y_trans_log2 = scale_y_log10(	breaks = trans_breaks("log2", function(x) 2^x),
 								labels = trans_format("log2", math_format(2^.x)))	

 x_trans_log2 = scale_x_log10(	breaks = trans_breaks("log2", function(x) 2^x),
 								labels = trans_format("log2", math_format(2^.x)))

####### plots

 A1 = (ggplot(global,  aes(file, coverage, color=file)) 
 		+ geom_boxplot() 
 		+ coord_flip()
 		+ theme(axis.ticks = element_blank(), axis.text.y = element_blank())	
 		)	
 A2 = (ggplot(global,  aes(file, scaled, color=file)) 
 		+ geom_boxplot() 
 		+ coord_flip()
 		+ theme(axis.ticks = element_blank(), axis.text.y = element_blank())
 		)	
 A3 = (ggplot(global,  aes(file, scaledlog2, color=file)) 
 		+ geom_boxplot() 
 		+ coord_flip()
 		+ theme(axis.ticks = element_blank(), axis.text.y = element_blank())	
 		)							

 B1 = (ggplot(global,  aes(coverage, ..count..,  fill=file, group=file )) 
 		+ stat_bin(aes(y=..count..), geom="area", position="stack", weight=2)
 		)
 B2 = (ggplot(global,  aes(scaled, ..count..,  fill=file, group=file )) 
 		+ stat_bin(aes(y=..count..), geom="area", position="stack", weight=2)
 		)
 B3 = (ggplot(global,  aes(scaledlog2, ..count..,  fill=file, group=file )) 
 		+ stat_bin(aes(y=..count..), geom="area", position="stack", weight=2)
 		)

 C1 = (ggplot(global,  aes(coverage, ..count..,  fill=file, group=file )) 
 		+ stat_bin(aes(y=..count..), geom="area", position="dodge", alpha=0.5, weight=2)
 		)
 C2 = (ggplot(global,  aes(scaled, ..count..,  fill=file, group=file )) 
 		+ stat_bin(aes(y=..count..), geom="area", position="dodge", alpha=0.5, weight=2)
 		)
 C3 = (ggplot(global,  aes(scaledlog2, ..count..,  fill=file, group=file )) 
 		+ stat_bin(aes(y=..count..), geom="area", position="dodge", alpha=0.5, weight=2)
 		)


 D1 = (ggplot(global,  aes(coverage, fill=file, group=file )) 
 		+ geom_bar(position="fill")
 		)
 D2 = (ggplot(global,  aes(scaled, fill=file, group=file )) 
 		+ geom_bar(position="fill")
 		)
 D3 = (ggplot(global,  aes(scaledlog2, fill=file, group=file )) 
 		+ geom_bar(position="fill")
 		)

################ rendered output files

 png("COVERAGE_GG.png", width=30, height = 2 + 4 * length(files), units="in", res=250)
 multiplot(	A1, A2, A3, 
 			B1, B2, B3,
 			C1, C2, C3, 
 			D1, D2, D3, 
 			cols=3 
 			)
 dev.off()

 png("COVERAGE_log10counts_GG.png", width=30, height = 2 + 4 * length(files), units="in", res=250)
 multiplot(	A1 + y_trans_log10, A2 + y_trans_log10, A2 + y_trans_log10,
 			B1 + y_trans_log10, B2 + y_trans_log10, B2 + y_trans_log10, 
 			C1 + y_trans_log10, C2 + y_trans_log10, C2 + y_trans_log10, 
 			D1 + y_trans_log10, D2 + y_trans_log10, D2 + y_trans_log10,
 			cols=3 
 			)
 dev.off()

 png("COVERAGE_facet_GG.png", width=10, height = 2 + 2 * (4 * length(files)), units="in", res=250)
 multiplot(	
 			B1 + facet_grid( file ~ .) , B2 + facet_grid( file ~ .), B3 + facet_grid( file ~ .) ,
 			cols=1
 			)
 dev.off()

 png("COVERAGE_log10counts_facet_GG.png", width=10, height = 2 + 2 * (4 * length(files)), units="in", res=250)
 multiplot(	
 			B1 + y_trans_log10 + facet_grid( file ~ .) , B2 + y_trans_log10 + facet_grid( file ~ .), B3 + y_trans_log10 + facet_grid( file ~ .), 
 			cols=1 
 			)
 dev.off()


