#!/usr/bin/env Rscript


getargs <- function(c) {
  #parse arguments
  suppressPackageStartupMessages(library("argparse"))
  
	parser <- ArgumentParser(description = 'regulation analysis')
	parser$add_argument('-base', type="character" , nargs = '*')
	parser$add_argument('-mutant', type="character", nargs='*')
	parser$add_argument('-selected', type="character", nargs='*')
	parser$add_argument('-disp', type="character", nargs='*')
	parser$add_argument('-pval', type="number", nargs='*')
	parser$add_argument('-out', type="character", required = TRUE)
	res <- parser$parse_args(c)

}

#getargs <- function(c) {
#  i <- 1
#
#  if (is.na(c[i]) | c[i] != "-base"){
#    print("must start with base")
#    return(NA)
#  } 
#  i <- i+1  
#
#  baseArgs=c()
#  while (c[i] != "-mutant") {baseArgs=c(baseArgs,c[i]);  i = i+1}
#  i <- i+1
#  
#  mutatedArgs=c()
#  while (c[i] != "-selected") {mutatedArgs=c(mutatedArgs,c[i]);  i <- i+1}
#  i <- i+1
#  
#  selected=c()
#  while (c[i] != "-disp") {selected=c(selected,c[i]);  i <- i+1}
#  i <- i+1
#  dispersion=c[i]
#  i=i+2
#  outDir=c[i]
#  
#  res <- list("base"=baseArgs,"mutant"=mutatedArgs,"selected"=selected,"disp"=dispersion,"out"=outDir)
  
#}

generateTable <- function(baseFileList,mutantFileList){
  #load all countfiles of base and mutant into a table
  
  
  geneList=(read.delim(file=baseFileList[1],sep=",",header=FALSE))[,]
  
  countDframe=data.frame(row.names = geneList[,1])
  

  for(name in baseFileList){
    colName=paste("base",match(name,baseFileList))           # column names will be the base gene names
    colData=read.delim(file=name,sep=",",header=FALSE)[,2]   # read the base file and its genecount 
    countDframe <- data.frame(countDframe,colData)           # create table linking gene to its count
    colnames(countDframe) <- c(names(countDframe)[-length(names(countDframe))],colName)
  }
  
  for(name in mutantFileList){
    colName=paste("mutant",match(name,mutantFileList))       # column names will be the mutant gene names
    colData=read.delim(file=name,sep=",",header=FALSE)[,2]   # read the mutant file and its genecount
    countDframe <- data.frame(countDframe,colData)           # create table linking gene to its count
    colnames(countDframe) <- c(names(countDframe)[-length(names(countDframe))],colName)
  }
  
  return(countDframe)
  
  
}

generateCondition <- function(baseFileList,mutantFileList) {
  # creates a condition vector, which consists of repeating "base" and "mutant" strings according to their list length
  condition=factor(c( rep.int("base",length(baseFileList))
              ,rep.int("mutant",length(mutantFileList))))
  return(condition)
}


generateBinomResults <- function(args){
  #normalize counts and calulate fault change and p value for each gene
  
  #generate condition vectors and countTables for base and mutant
  condition <- generateCondition(args$base,args$mutant)
  countTable <- generateTable(args$base,args$mutant)
  
  #load DESeq
  suppressMessages(library(DESeq,quietly = TRUE))
  
  #create countDATASET and execute binom test
  countDS <- newCountDataSet(countTable,condition)
  countDS <- estimateSizeFactors(countDS)
  countDS <- estimateDispersions(countDS, method = args$disp , sharingMode = "fit-only",fitType="local" )
  

  binomResTable <- nbinomTest(countDS,"base","mutant")
  
  
  # concatenate counttable with binomTest table
  ids=as.character(unlist(binomResTable[,1]))
 
  fullResTable=data.frame(ids,countTable,binomResTable[,2:length(colnames(binomResTable))])
  rownames(fullResTable) <- NULL
  colnames(fullResTable) <- c("id",colnames(fullResTable[2:length(colnames(fullResTable))]))
  
  return(list("fullResTable"=fullResTable,"countDS"=countDS))

}

filterTable <- function(table,idFile) {
  # get sub tables that correspond to the genes in the selected files
  filterList <- (read.delim(file=idFile,sep=",",header=FALSE))
  filterList <-  as.character(unlist(filterList))
  
  truthVec=table[,1] %in% filterList
  subTable=table[truthVec,,]
}

cleanTableFromNA_INF <- function(table){
  # remove NA and INF symbols
  cleanTable <- data.frame(na.omit(table))
  
  cleanTable <- cleanTable[is.finite(cleanTable$foldChange),]
  return(cleanTable)
}


cleanTableFromP_Value <- function(table,pval){
  
  cleanTable <- data.frame(na.omit(table))
  
  cleanTable <- cleanTable[cleanTable$padj < pval,]
  return(cleanTable)
}



analysis_main <- function(c){
  
  #parse arguments
  args <- getargs(c)	# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< the parentheses used to only contain "c"
  
  ######################################
  
  
  # get binom result tables and countDS
  deseqInfo=generateBinomResults(args)
  fullResTable <- deseqInfo$fullResTable
  countDS <- deseqInfo$countDS
  
  fullResTable <- fullResTable[with(fullResTable,order(fullResTable[,1])),]
  
  cleanTable=cleanTableFromNA_INF(fullResTable)
  cleanTable <- cleanTableFromP_Value(cleanTable,args$pval)
  
  filteredTables=list()
  for(S in args$selected){
    #for each selected file get the filtered table
    selTable <- filterTable(fullResTable,S)
    #add it to the list of filtered tables
    filteredTables[[length(filteredTables)+1]]=selTable
    
  }
  
  
  tResults <- list()
  for(table in filteredTables){
    cleanSubTable <- cleanTableFromNA_INF(table)
    cleanSubTable <- cleanTableFromP_Value(cleanSubTable,args$pval)
    
    tTestRes <- t.test(cleanTable$foldChange,cleanSubTable$foldChange)  #run t test on subtable
    tResults[[length(tResults)+1]]=tTestRes
    
  }
  
  
  ######################################
  
  
  #make results directory and enter it
  mainDir=getwd()
  subDir=args$out
  
  if (file.exists(subDir)){ #if name taken delete the old directory
    unlink(subDir,recursive = TRUE)
  }
  
  newpath <- file.path(mainDir, subDir)
  dir.create(newpath)
  
  
  
  
  ######################################
  
  #write results to files
  
  #write sizefactors to file
  sink(file.path(newpath,"sizeFactors"),append=TRUE)
  print(sizeFactors(countDS))  
  
  sink()
  
  #write ttest to file
  sink(file.path(newpath,"t_test_data"),append=TRUE)
  for(i in 1:length(args$selected)){
    selectedFile=args$selected[[i]]
    currTRes <- tResults[[i]]
    
    print(paste("t_Test Results between main data set and",selectedFile))
    print(currTRes)
    
  }
  sink()
  
  #write full table
  write.csv(format(fullResTable,digits=3),file=file.path(newpath,"raw_results_table.csv"))
  
  #write filtered tables
  for(i in 1:length(args$selected)){
    selectedFile=args$selected[[i]]
    table=filteredTables[[i]]
    
    write.csv(format(table,digits=3),file=file.path(newpath,paste(basename(selectedFile),"_filterred_results.csv")))
    
  }
  
  #write variance results
  options(scipen = -20)
  if(args$disp=="per-condition"){
    #base dispersion
    tiff(filename = file.path(newpath,paste("baseVariance","tiff",sep = ".")), width=800,height=800,units="px",res=85)
    plotDispEsts(countDS,name="base",main="base_variance",xlab="mean of normalized counts",ylab="dispersion", log = "xy") # <<<< added "log=..."
    #axis(1,labels = format(scientific = TRUE ))
    #axis(2,labels = format(scientific = TRUE ))
    dev.off()
    
    fitBase <- fitInfo(countDS,name="base")
    fitBaseDat <- data.frame(fitBase$fittedDispEsts)
    ordfitBaseDat<-cbind(rownames(fitBaseDat)[order(rownames(fitBaseDat))], fitBaseDat[order(rownames(fitBaseDat)),])
    write.csv(ordfitBaseDat,file=file.path(newpath,"base_Variance_Data.csv"))
    
    #mutant dispersion
    tiff(filename = file.path(newpath,paste("mutantVariance","tiff",sep = ".")), width=800,height=800,units="px",res=85)
    plotDispEsts(countDS,name="mutant",main="mutant_variance",xlab="mean of normalized counts",ylab="dispersion", log = "xy")	# <<<< added "log=..."
    #axis(1,labels = format(scientific = TRUE ))
    #axis(2,labels = format(scientific = TRUE ))
    dev.off()
    
    fitMutant <- fitInfo(countDS,name="mutant")
    fitMutantDat <- data.frame(fitMutant$fittedDispEsts)
    ordfitMutantDat<-cbind(rownames(fitMutantDat)[order(rownames(fitMutantDat))], fitMutantDat[order(rownames(fitMutantDat)),])
    write.csv(ordfitMutantDat,file=file.path(newpath,"mutant_Variance_Data.csv"))
  }
  else{ #no per condition only one variance table
    tiff(filename = file.path(newpath,paste("totalVariance","tiff",sep = ".")), width=800,height=800,units="px",res=85)
    plotDispEsts(countDS,main="variance",xlab="mean of normalized counts",ylab="dispersion", log = "xy")		# <<<< added "log=..."
    dev.off()
    
    fitTotal <- fitInfo(countDS)
    fitTotalDat <- data.frame(fitTotal$fittedDispEsts)
    ordfitTotalDat<-cbind(rownames(fitTotalDat)[order(rownames(fitTotalDat))], fitTotalDat[order(rownames(fitTotalDat)),])
    write.csv(ordfitTotalDat,file=file.path(newpath,"mutant_Variance_Data.csv"))
  }
 
}



