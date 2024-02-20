# Contains directories specifying data locations

#Load frequently used core R packages
library (data.table)
library (ggplot2)
library (ggpmisc)
library (readxl)

#dir.base must be specified in main script before calling this loader
#dir.base <- "~/Metastatic-breast-immune-profiling/"

dir.data   <- paste0(dir.base,"data/")
dir.metadata  <- paste0(dir.base,"metadata/")
dir.resources <- paste0(dir.base,"resources/")
dir.r         <- paste0(dir.base,"R/")
dir.out.plots <- paste0(dir.base,"output/plots/")
dir.out.data  <- paste0(dir.base,"output/data/")

dir.mrdarcy.bcr.mets   <- paste0(dir.data,"/processed/MrDarcy/metastatic/BCR/")
dir.mrdarcy.tcr.mets   <- paste0(dir.data,"/processed/MrDarcy/metastatic/TCR/")
dir.mrdarcy.bcr.early  <- paste0(dir.data,"/processed/MrDarcy/early/BCR/")

#load file structure
file <- list()

#BCR sequencing data list
file[["bcr.counts.mets"]]  <- paste0(dir.data, "/processed/bcr-data/BCRseq-data.metastatic.RData")
file[["bcr.counts.early"]] <- paste0(dir.data, "/processed/bcr-data/BCRseq-data.early.RData")

#DNA files
file[["dna.mhc1Nag.mets"]]   <- paste0(dir.data,"/processed/dna-data/dna.wes.mhc-i.RData")
file[["dna.mhc2Nag.mets"]]   <- paste0(dir.data,"/processed/dna-data/dna.wes.mhc-ii.RData")
file[["dna.mutations.mets"]] <- paste0(dir.data,"/processed/dna-data/dna.wes.mutations.RData")
file[["dna.onconem.mets"]]   <- paste0(dir.data,"/processed/dna-data/onconem/")

#RNA files
file[["rna.danaher.mets"]]      <- paste0(dir.data,"/processed/rna-data/rna.tme.danaher.RData")
file[["rna.expr.counts.mets"]]  <- paste0(dir.data,"/processed/rna-data/rna.counts.metastatic.Rdata")
file[["rna.expr.counts.early"]] <- paste0(dir.data,"/processed/rna-data/rna.counts.early.Rdata")
file[["rna.expr.tpm.mets"]]     <- paste0(dir.data,"/processed/rna-data/rna.bulk.tpm.RData")
file[["rna.mcp.mets"]]          <- paste0(dir.data,"/processed/rna-data/rna.tme.mcpcounter.RData")

#Resource files
file[["res.danaher"]]       <-paste0(dir.resources,"danaher-immune-cells.RData")
file[["res.ensgToHugo"]]    <- paste0(dir.resources,"EnsemblID.to.Hugo.v87.tsv.gz")
file[["res.ensgToGeneLen"]] <- paste0(dir.resources,"EnsemblID.to.genelength.v87.txt.gz")
file[["res.msigdb.c2"]]     <- paste0(dir.resources,"msigdb.7.3.c2.cgp.Rdata")
file[["res.msigdb.c5"]]     <- paste0(dir.resources,"msigdb.7.3.c5.go.bp.RData")
file[["res.msigdb.h"]]      <- paste0(dir.resources,"msigdb.7.3.hallmarks.RData")



#Load ggplot theme
source (paste0(dir.r,"theme.R"))

#Load tumour site metadata
metadata.mets <- data.frame(read_xlsx(paste0(dir.metadata,"mets-sample_metadata.xlsx")),stringsAsFactors = F)[,c(2,3)]
metadata.mets.noAbbrev <- metadata.mets
metadata.mets$Disease.site[metadata.mets$Disease.site %in% c("LN")] <- "Lymph node"
metadata.mets$Disease.site[metadata.mets$Disease.site %in% c("Brain","Meninges")] <- "Brain/meninges"

metadata.early <- data.frame(read_xlsx(paste0(dir.metadata, "early-sample_metadata.xlsx"),sheet = 1))

# returns BCR sequence encoded by cdr1->fwr4 + c_call
getBCRSequence <- function(x){
  return (gsub("-","",gsub("\\.","",paste0(x$cdr1,x$fwr2,x$cdr2,x$fwr3,x$cdr3,x$fwr4,"_",x$c_call))))
}

# returns VDJ sequence encoded by cdr1->fwr4
getVDJSequence <- function(x){
  return (gsub("-","",gsub("\\.","",paste0(x$cdr1,x$fwr2,x$cdr2,x$fwr3,x$cdr3,x$fwr4))))
}
