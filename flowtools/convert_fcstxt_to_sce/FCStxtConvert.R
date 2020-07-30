#!/usr/bin/env Rscript
# GECO flow text conversion tool
# Authors: Emily Combe and Pablo Moreno
#
# This tool converts a flowtext file (or tabular file) into a SingleCellExperiment object
# The tool was written by Emily Combe and edited by Pablo Moreno
#
# There are the options to choose: the columns/markers to include in the assay, the columns to include in the meta data, descriptions of the markers and a metadata file.
# 
#
#
# Version 1
# July 2020 (Emily Combe / Pablo Moreno)


library(SingleCellExperiment)
#library(dplyr) - no longer needed

suppressPackageStartupMessages(library("optparse"))

sce <- function(input, fl_cols = list(), mtd_cols = list(), marker_type = list(), meta_data = NULL) {
    
    
    #---------------------#
    # reading in flowtext #
    #---------------------#
    
    flowtext <- read.table(input, sep = "\t", header=T) 
    
    #----------------------------------#
    # extract-marker-fluorescence data # 
    #----------------------------------#
    
    fl_cols_assay <- colnames(flowtext)
    
    if (length(fl_cols) > 0){
       
        if(length(fl_cols) > ncol(fl_cols_assay)){
            quit(save = "no", status = 13, runLast = FALSE)
        }
        fl_cols_assay <- fl_cols_assay[fl_cols]
    } else {
    
        channels_to_exclude <- c(grep(fl_cols_assay, pattern="FSC"),
                                 grep(fl_cols_assay, pattern="SSC"),
                                 grep(fl_cols_assay, pattern="FSC-A"),
                                 grep(fl_cols_assay, pattern="SSC-A"),
                                 grep(fl_cols_assay, pattern="FSC-W"),
                                 grep(fl_cols_assay, pattern="SSC-W"),
                                 grep(fl_cols_assay, pattern="FSC-H"),
                                 grep(fl_cols_assay, pattern="SSC-H"),
                                 grep(fl_cols_assay, pattern="Time", ignore.case = T),
                                 grep(fl_cols_assay, pattern="Population|flowSOM|cluster|SOM|pop|cluster", ignore.case = T),
                                 grep(fl_cols_assay, pattern="Live_Dead|live|dead", ignore.case = T))
        
        fl_cols_assay <- fl_cols_assay[-channels_to_exclude]
        
    }
    
    counts <- flowtext[, fl_cols_assay, drop = FALSE] 
    
    counts <- as.matrix(counts)

    # transpose data into assay as columns=cells and rows=features.
    counts <- base::t(counts)
    
    colnames(counts) <- 1:ncol(counts)  
    
    
    #-----------------#
    #coldata/meta data#
    #-----------------#
    
    # by default any columns with sample names or cluster results will be extracted - to over ride this user must provide a comma separated list of column name (mtd_cols)
    
    mtd_cols_assay <- colnames(flowtext)
    
    if (length(mtd_cols) > 0){
        
        if(length(mtd_cols) > ncol(mtd_cols_assay)){
            quit(save = "no", status = 14, runLast = FALSE)
        }
        
        mtd_cols_assay <- mtd_cols_assay[mtd_cols]
 
    } else {
        
        #print("Meta data columns from flowtext files not specified")
        #create warning here to the user - but without failing
        
        mtd_columns <- c(grep(markers, pattern="sample", ignore.case=T),
                         grep(markers, pattern="population|flowsom|cluster|pop|som", ignore.case=T))

        mtd_cols_assay <- mtd_cols_assay[mtd_columns]
    }
    
    md <- flowtext[, mtd_cols_assay, drop = FALSE]
    
    # if metadata available will be merged with meta data from flow text
    if(!is.null(meta_data)){
        
        #match column names so case insensitive
        md_col <- tolower(colnames(md))
        mtd_col <- tolower(colnames(metadata))
        
        #quit if < 1 or > 1 column names match
        if(length(intersect(md_col, mtd_col)) == 0){
            quit(save = "no", status = 15, runLast = FALSE)
        }
        if(length(intersect(md_col, mtd_col)) > 1){
            quit(save = "no", status = 16, runLast = FALSE)
        }
        
        md_intersect <- colnames(md)[intersect(md_col, mtd_col)]
        mtd_intersect <- colnames(meta_data)[intersect(mtd_col, md_col)]    
        
        #merge by matched column
        metadata <- merge(x = md, y = meta_data, by.x = md_intersect, by.y = mtd_intersect, all=T)
        
    } 

    #create Single Cell experiment object. SCOPE requires both counts and logcounts assays - for FLOW both assays contain the same data
    sce <- SingleCellExperiment(assays = list(counts=counts, logcounts=counts), colData=metadata)
    
    
    #-----------------#
    # row/marker data #
    #-----------------#
    
    if(length(marker_type) > 0){
        
	if(length(marker_type) != nrow(rowData(sce))){
	    quit(save = "no", status = 17, runLast = FALSE)
	}

        marker_type[marker_type == "l"] <- "lineage"
        marker_type[marker_type == "f"] <- "functional"
        
        rowData(sce)$marker_type <- marker_type
    }
    
    return(sce)
    
}  



args <- commandArgs(trailingOnly = TRUE)


fl_channels <- list()

# fluorescence markers to include in the assay
if (args[3]=="None") {
    flag_default <- TRUE
} else {
    if (args[3] == "i.e.:CD8, CD4, CD8"){
        flag_default <- TRUE
    } else {
        fl_channels <- as.character(strsplit(args[3], ",")[[1]])
        for (channel in fl_channels){
            if (is.na(channel)){
                quit(save = "no", status = 10, runLast = FALSE)
            }
        }
    }
}

# meta data columns to go into colDaa in SCE
mt_channels <- list()

if (args[4]=="None") {
    flag_default <- TRUE
} else {
    if (args[4] == "i.e.:Sample, Population"){
        flag_default <- TRUE
    } else {
        mt_channels <- as.character(strsplit(args[4], ",")[[1]])
        for (channel in mt_channels){
            if (is.na(channel)){
                quit(save = "no", status = 11, runLast = FALSE)
            }
        }
    }
}


#metadata file to add to the coldata in SCE. Must have column matching the sample column in the flowtext file
md <- NULL

if (args[5]=="None") {
    flag_default <- TRUE
} else {
    md <- read.table(args[5], header = TRUE, sep = "\t", check.names = FALSE, as.is = FALSE)
}

#comma separated list of values to define the markers included in the assay
mark_type <- list()

if (args[6]=="None") {
    flag_default <- TRUE
} else {
    if (args[6] == "i.e: lineage, lineage, functional"){
        flag_default <- TRUE
    } else {
        mark_type <- as.character(strsplit(args[6], ",")[[1]])
        for (mt in mark_type){
            if (is.na(mt)){
                quit(save = "no", status = 12, runLast = FALSE)
            } 
        } 
    }
}


sce <- sce(input = args[1], fl_cols = fl_channels, mtd_cols = mt_channels, meta_data = md, marker_type = mark_type)

saveRDS(sce, file = args[2])
