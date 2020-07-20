#!/usr/bin/Rscript
# Module for Galaxy
# Adds sample information to a flowclustered file and binds 
#
#                  
# Version 1
# Emily Combe
#
#



#-------------#
# sample name #
#-------------#



sample_name <- function(input_file, sample_names){
	flowtext <- list()
	for(i in input_file){
		flowtext[i] <- read.table(input_file[i], header = TRUE, sep = "\t")
	}
	
	sample_names <- as.character(sample_names)
	
	for(i in length(flowtext)){
   		 flowtext[[i]]$sample <- rep(sample_names[i], nrow(flowtext[[i]]))
    	}

	if(length(flowtext) > 1){
		flowtext <- rbindlist(flowtext)
	} else {
		flowtext <- flowtext[[1]]
	}

	return(flowtext)
}

args <- commandArgs(trailingOnly = TRUE)

sample_name(args[1], args[3])

