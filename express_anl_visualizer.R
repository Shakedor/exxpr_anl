source("express_anl_visualizer_aux.R")

#get arguments
c <- commandArgs(trailingOnly = TRUE)

#parse arguments
args <- getargs(c)
if (is.na(args))
	q() #quit

table <- generateTable(args["xfile"],args["yfile"])

generatePlots(table,args)

