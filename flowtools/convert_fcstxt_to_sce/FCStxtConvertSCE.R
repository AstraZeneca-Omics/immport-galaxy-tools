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


suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(optparse))

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

        if(length(fl_cols) > ncol(flowtext)){
            quit(save = "no", status = 13, runLast = FALSE)
        }
        fl_cols_assay <- fl_cols_assay[fl_cols_assay %in% fl_cols]
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
        if(length(mtd_cols) > ncol(flowtext)){
            quit(save = "no", status = 14, runLast = FALSE)
        }
        mtd_cols_assay <- mtd_cols_assay[mtd_cols_assay %in% mtd_cols]
    } else {
        #print("Meta data columns from flowtext files not specified")
        #create warning here to the user - but without failing
        mtd_columns <- c(grep(marker_type, pattern="sample", ignore.case=T),
                         grep(marker_type, pattern="population|flowsom|cluster|pop|som", ignore.case=T))

        mtd_cols_assay <- mtd_cols_assay[mtd_columns]
    }

    md <- flowtext[, mtd_cols_assay, drop = FALSE]

    # if metadata available will be merged with meta data from flow text
    if(!is.null(meta_data)){

        #match column names so case insensitive
        md_col <- tolower(colnames(md))
        mtd_col <- tolower(colnames(meta_data))

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
        meta_data <- merge(x = md, y = meta_data, by.x = md_intersect, by.y = mtd_intersect, all=T)

    }

    #create Single Cell experiment object. SCOPE requires both counts and logcounts assays - for FLOW both assays contain the same data
    sce <- SingleCellExperiment(assays = list(counts=counts, logcounts=counts))
    if(!is.null(meta_data)) {
      colLabels(sce)<-metadata
    }


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

option_list = list(
  make_option(
    c("-i", "--input"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name for FCS txt file with sample information."
  ),
  make_option(
    c("-o", "--output"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name for output SCE R RDS Object."
  ),
  make_option(
    c("-f", "--fl_cols"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Comma separated list of Columns with markers to be included in the Single Cell Experiment assay"
  ),
  make_option(
    c("-m", "--metadata_columns"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Columns to be included in the metadata of the Single Cell Experiment."
  ),
  make_option(
    c("--metadata_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Optional meta data txt file to include in Single Cell Experiment.'
  ),
  make_option(
    c("--marker_type"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Marker type'
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

# fluorescence markers to include in the assay
fl_channels <- list()
if (is.na(opt$fl_cols)) {
    flag_default <- TRUE
} else {
    fl_channels <- as.character(strsplit(opt$fl_cols, ",")[[1]])
    for (channel in fl_channels){
        if (is.na(channel)){
            quit(save = "no", status = 10, runLast = FALSE)
        }
    }
}

# meta data columns to go into colDaa in SCE
mt_channels <- list()
if (is.na(opt$metadata_columns)) {
    flag_default <- TRUE
} else {
    mt_channels <- as.character(strsplit(opt$metadata_columns, ",")[[1]])
    for (channel in mt_channels){
        if (is.na(channel)){
            quit(save = "no", status = 11, runLast = FALSE)
        }
    }
}


#metadata file to add to the coldata in SCE. Must have column matching the sample column in the flowtext file
md <- NULL
if (is.na(opt$metadata_file)) {
    flag_default <- TRUE
} else {
    md <- read.table(opt$metadata_file, header = TRUE, sep = "\t", check.names = FALSE, as.is = FALSE)
}

#comma separated list of values to define the markers included in the assay
mark_type <- list()
if (is.na(opt$marker_type)) {
    flag_default <- TRUE
} else {
    mark_type <- as.character(strsplit(opt$marker_type, ",")[[1]])
    for (mt in mark_type){
        if (is.na(mt)){
            quit(save = "no", status = 12, runLast = FALSE)
        }
    }
}


sce <- sce(input = opt$input, fl_cols = fl_channels, mtd_cols = mt_channels, meta_data = md, marker_type = mark_type)

saveRDS(sce, file = opt$output)
