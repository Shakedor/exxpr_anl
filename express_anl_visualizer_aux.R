#!/usr/bin/env Rscript

<<<<<<< Updated upstream
=======
#usage: Rscript expression_main.R -graph_name name -xaxis xtitle -yaxis ytitle -size 1.5 -color black -reg y -xfile xfile -yfile yfile 
#		-selected 1 orange y selected_file1 2 blue n selected_file2 -legend y -names general selected_name1 selected_name2 -legend_pos topleft -font "Arial Black"
#graph_name: Name of the graph
#xaxis: Title of horizontal axis
#yaxis: Title of vertical axis
#size: Size of points in the graph (numeric)
#color: Color of points in the graph
#reg: Draw regression line. Options: y/n (yes/no)
#xfile: Data file for horizontal axis
#yfile: Data file for vertical axis
#selected: Quadrupltes of: size color reg file (e.g. 1.5 blue n file_path) (e.g. 2 orange y file_path2). Must be in quadruplets!!
#legend: Put a legend in the graph. Options: y/n (yes/no)
#names: If legend exists - the names in each line of it (the colors of the dots are according to the colors in selected).
#legend_pos: If legend exists - position of the legend inside the graph. Options: bottomright, bottom, bottomleft, left, topleft, top, topright, right, center
#font: Font to be used (whole graph, including legend if exists)

#Number of names = Number of selected quadruplets


>>>>>>> Stashed changes
library(argparse)
library(ggplot2)

getargs <- function(c) {
	parser <- ArgumentParser(description = 'expression')
	parser$add_argument('-graph_name', required = TRUE)
	parser$add_argument('-xaxis', required = TRUE)
	parser$add_argument('-yaxis', required = TRUE)
	parser$add_argument('-size', required = TRUE)
	parser$add_argument('-color', required = TRUE)
	parser$add_argument('-reg', type="character", choices=c("y","n"), default="n")
	parser$add_argument('-xfile', required = TRUE)
	parser$add_argument('-yfile', required = TRUE)
	parser$add_argument('-selected', nargs="*")
	parser$add_argument('-names', nargs="*")
	parser$add_argument('-legend', type="character", choices=c("y","n"), default="n")	
	parser$add_argument('-legend_pos', type="character", choices=c("bottomright","bottom","bottomleft","left","topleft","top","topright","right","center"), default="")
	parser$add_argument('-font', default="")
	args <- parser$parse_args()

	if (length(args$selected)%%4 != 0) {
	    print("selected lists must come in size color reg(y/n) file quadruplets")
		return(NA)
	}
	
    #TODO add color format checker
    #TODO add file recogniser checker

	selected_sizes=c()
	selected_colors=c()
	selected_reg_bools=c()
	selected_files=c()
	
	pairNum=(length(args$selected))/4
	if (pairNum > 0){
		for(j in 1:pairNum){
			if(args$selected[4*j-1]!="y"&args$selected[4*j-1]!="n"){
				print("reg argument must contain either n or y")
				print(paste("error in arg selected",4*j-1))
				return(NA)
			}
			selected_sizes=c(selected_sizes,args$selected[4*j-3])
			selected_colors=c(selected_colors,args$selected[4*j-2])
			selected_reg_bools=c(selected_reg_bools,args$selected[4*j-1])
			selected_files=c(selected_files,args$selected[4*j])

		}	
	}
	if (args$legend == "y" & is.null(args$names)){
		args$names[1] = "All Genes"
		if (pairNum > 0){
			for (i in 1:pairNum){
				args$names[i+1] = tail(unlist(strsplit(selected_files[i], "/")), n=1) # take only the filename, not whole path
			}
		}
	}
	else if (args$legend == "y" & length(args$names) != pairNum + 1){
		print("If you want a legend, provide as many names as needed = number(selected tuples)+1")
		return(NA)	
	}
	
	res <- list("graphName"=args$"graph_name","xlab"=args$"xaxis","ylab"=args$"yaxis","size"=as.numeric(args$"size"),"color"=args$"color",
		      "reg"=args$"reg","xfile"=args$"xfile","yfile"=args$"yfile", "name"=args$"name", "legend"=args$"legend", "legend_pos"=args$"legend_pos",
			  "selected_sizes"=selected_sizes, "selected_colors"=selected_colors,"selected_reg_bools"=selected_reg_bools,
		      "selected_files"=selected_files,"selected_names"=args$names,"selectedNum"=pairNum, "font"=args$"font")
	return(res)
}


generateTable <- function(baseFileList,mutantFileList){
  
  
  geneList=(read.delim(file=as.character(baseFileList[1]),sep=",",header=FALSE))[,1]
  
  geneList=as.character(unlist(geneList))

  countDframe=data.frame(geneList)
  
  colnames(countDframe) <- "gene"
  
  
  for(name in baseFileList){
    colName=paste("base",match(name,baseFileList))
    colData=read.delim(file=name,sep=",",header=FALSE)[,2]
    colData=as.double(unlist(colData))
    countDframe <- data.frame(countDframe,colData)
    colnames(countDframe) <- c(names(countDframe)[-length(names(countDframe))],colName)
  }
  
  for(name in mutantFileList){
    colName=paste("mutant",match(name,mutantFileList))
    colData=read.delim(file=name,sep=",",header=FALSE)[,2]
    colData=as.double(unlist(colData))
    countDframe <- data.frame(countDframe,colData)
    colnames(countDframe) <- c(names(countDframe)[-length(names(countDframe))],colName)
  }
  
  return(countDframe)
  
  
}

cleanTableNoZeros <- function(table){
  cleanTable <- table[table[,2]!=0 & table[,3]!=0,]
  return(cleanTable)
}

generatePlots <- function(fullTable,args){

	if (args$"font" != "") {
		library("extrafont")
		#font_import()	
		loadfonts(quiet=TRUE)
	}
	
	#plot main data
	cleanTable <- cleanTableNoZeros(fullTable)
	
	tiff(filename = paste(args$graphName,"tiff",sep = "."), width=800,height=800,units="px",res=85)
	
	plot(cleanTable[,2],cleanTable[,3], log='xy', type = "n", main=args$graphName, xlab = args$xlab,ylab = args$ylab, family = args$"font")
	points(cleanTable[,2],cleanTable[,3],cex = args$size , col = args$color,pch=16)
<<<<<<< Updated upstream

=======
>>>>>>> Stashed changes
	if(args$reg == "y"){
		lm.r <- lm(log10(cleanTable[,3])~log10(cleanTable[,2]))
		abline(lm.r , col="red" , untf=F,lwd=6)    
	}
	op <- par(family = args$"font")
	if (args$legend == "y"){
		if (args$"legend_pos" == ""){
			legend(1,max(cleanTable[,3]), legend=as.character(args$selected_names), col=as.character(c(args$color,args$selected_colors)), pch=16)
		}
		else {
			legend(args$"legend_pos", legend=as.character(args$selected_names), col=as.character(c(args$color,args$selected_colors)), pch=16)
		}
	}
	par(op)
  
## ggplot past version  
#   baseGrid <- ggplot(fullTable,environment=environment())
#   points <- geom_point(aes(base.1+0.5,`mutant 1`+0.5))
#   logScaleY <- scale_y_log10()
#   logScaleX <- scale_x_log10()
#   xlab <- xlab(args$xlab)
#   ylab <- ylab(args$ylab)
#   title <- ggtitle(args$graphName)
#   
#   coefs <- coef(lm(log(fullTable$`mutant 1`+0.5) ~ log(fullTable$base.1+0.5) ))
#   print(coefs)
#   
#   abline <- geom_abline(intercept = coefs[1], slope = coefs[2])
#   plot <- baseGrid+points+logScaleX+logScaleY+xlab+ylab+title+abline
  
	if (args$selectedNum == 0){
		dev.off()
	}
	else {
		#plot selected data
		for(i in 1:args$selectedNum){
			selRegBool=as.character(args$selected_reg_bools[i])
			selSize=as.numeric(args$selected_sizes[i])
			selFile=as.character(args$selected_files[i])
			selColor=as.character(args$selected_colors[i])


			filterList <- (read.delim(file=selFile,sep=",",header=FALSE))
			filterList <-  as.character(unlist(filterList))


			truthVec=fullTable[,1] %in% filterList

			subTable=fullTable[truthVec,,]
			subTable <- subTable[order(subTable[,1]),]
			subTable <- cleanTableNoZeros(subTable)

			points(subTable[,2],subTable[,3],col=selColor,pch=16,cex=selSize)
			if(selRegBool=="y"){
			  lm.rs <- lm(log10(subTable[,2])~log10(subTable[,3]))
			  abline(lm.rs , col=selColor , untf=T,lwd=3)
			}
		}
		dev.off()
	} 
}


# ggplot past version
#     subTable <- data.frame(subTable,selColor)
#     colnames(subTable) <- c(names(subTable)[-length(names(subTable))],"colorCol")
# 
#     
#     Selpoints <- geom_point(data = subTable,aes(base.1+0.5,`mutant.1`+0.5,color=colorCol))
#     
#     
#     Selcoefs <- coef(lm(log(subTable$`mutant.1`+0.5) ~ log(subTable$base.1+0.5) ))
#  
#     SelAbline <- geom_abline(intercept = Selcoefs[1], slope = Selcoefs[2],color=selColor)
#     
#     plot <- plot+Selpoints+SelAbline
  
#   plot
#   ggsave(paste(args$graphName,"tiff",sep = "."),width = 50,height = 50 , units="cm")





