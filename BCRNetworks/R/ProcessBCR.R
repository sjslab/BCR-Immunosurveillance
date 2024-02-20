
#function to convert a dataframe to a fasta file
writeFasta<-function(X, fileName){
  D <- do.call(rbind, lapply(seq(nrow(X)), function(i) t(X[i, ])))
  write.table(D, fileName,row.names = FALSE, col.names = FALSE, quote = FALSE)
  if (file.exists(paste0(fileName,".gz"))) file.remove(paste0(fileName,".gz"))
  # R.utils::gzip(fileName,destname=paste0(fileName,".gz"))
}


initialiseDirs <- function(dir){
  # create working directory if does not exist
  if (!file.exists(dir)) dir.create(dir, recursive = T)

  #MRDARCY file structure
  if (!grepl("\\/$",dir)) dir <- paste0(dir,"/")
  dir.imgt.in  <- paste0(dir, "IMGT/input/")
  dir.imgt.out <- paste0(dir, "IMGT/output/")
  dir.mrdarcy.in.fasta  <- paste0(dir, "input/fasta/")
  dir.mrdarcy.in.imgt <- paste0(dir, "input/imgt/")
  dir.mrdarcy.out     <- paste0(dir, "output/")

  dir.create(dir.imgt.in,recursive = T, showWarnings = F)
  dir.create(dir.imgt.out,recursive = T, showWarnings = F)
  dir.create(dir.mrdarcy.in.fasta,recursive = T, showWarnings = F)
  dir.create(dir.mrdarcy.in.imgt,recursive = T, showWarnings = F)
  dir.create(dir.mrdarcy.out,recursive = T, showWarnings = F)

  return (list(dir.base = dir,
               dir.imgt.in          = dir.imgt.in,
               dir.imgt.out         = dir.imgt.out,
               dir.mrdarcy.in.fasta = dir.mrdarcy.in.fasta,
               dir.mrdarcy.in.imgt  = dir.mrdarcy.in.imgt,
               dir.mrdarcy.out      = dir.mrdarcy.out))
}


createImgtInput <- function(x, dir){

  if (!all(colnames(x) == c("sequence","c_call","count","sample_id","donor_id"))){
    stop("Incorrect dataframe structure. Column names expected: sequence c_call count sample_id donor_id")
  }

  saveRDS(x, paste0(dir$dir.imgt.in,"o.RData"))

  sequence <- x$sequence[!duplicated(x$sequence)]

  # annotate sequences with unique IDs
  id <- paste0(">s",c(1:length(sequence)))

  # create Master Fasta file to send to IMGT
  fasta <- data.table::data.table( name=id, seq=sequence,stringsAsFactors = F)
  fastaFileName <- paste0(dir$dir.imgt.in, "AllSequences.fasta")
  writeFasta(fasta,fastaFileName)

  # merge nucleotide sequences with IDs
  annotatedSeq      <- merge(x,fasta,by.x="sequence",by.y="seq")
  annotatedSeq$name <- gsub(">","",annotatedSeq$name)
  annotatedSeqFile  <- paste0(dir$dir.imgt.in, "AnnotatedSequences.RData")
  saveRDS(annotatedSeq, annotatedSeqFile)

  # create sample specific fastas

  x <- annotatedSeq
  x<-x[x$c_call!="",]
  class <- sort(unique(x$c_call))

  sampleFasta <- list()

  for (sample in sort(unique(x$sample_id))){
    cat("Processing: ",sample,"\n")
    y <- x[x$sample_id==sample,]
    z <- reshape2::dcast(y,sequence+name~c_call,value.var = "count")
    z[is.na(z)] <- 0

    # add IG isotypes that have not been detected
    classNotInCol <- class[!class %in% colnames(z)]
    for (i in classNotInCol) z[[i]]<-0

    # reorder column names
    ord <- c(1,2,order(colnames(z)[c(3:ncol(z))])+2)
    z <- z[, ord]
    classColNames <- colnames(z)[3:ncol(z)]

    z$id <- paste0(paste0(">",z$name,"__"),
                   apply(z[, classColNames,drop=F], 1, function(x) paste(x,collapse="_")),
                   "|",
                   paste(classColNames,collapse="_"))

    fasta <- data.frame(name=z$id, seq=z$sequence, stringsAsFactors = F)
    sampleFasta[[sample]] <- fasta
  }

  saveRDS(sampleFasta,paste0(dir$dir.imgt.in, "sampleFasta.RData"))
}



processIMGToutput <- function(dir){

  # load sample-specific fasta files
  fa <- readRDS(paste0(dir$dir.imgt.in, "sampleFasta.RData"))

  # load IMGT output
  file.imgtOutput  <- grep("3_Nt-sequences.txt", list.files(dir$dir.imgt.out,full.names = T),value = T)

  # load previous dataframe
  x <- readRDS(paste0(dir$dir.imgt.in,"o.RData"))

  if(length(file.imgtOutput)==0) {
    stop(paste0("A 3_Nt-sequences.txt file was not found in: ",dir$dir.imgt.out))
  }

  imgt <- data.table::fread(file.imgtOutput,fill=TRUE,sep="\t")
  imgt <- imgt[imgt$`V-DOMAIN Functionality`!="No results",]

  #create MRDarcy input for each sample
  toExclude <- character()
  for (sample in names(fa)){
    cat("Processing: ",sample,"\n")

    f <- fa[[sample]]
    f$id <- gsub(">","",sapply(strsplit(f$name,"__"),"[",1))
    imgtSubset   <- imgt[imgt$`Sequence ID` %in% f$id,]

    f <- f[match(imgtSubset$`Sequence ID`,f$id),]
    stopifnot(sum(f$id!=imgtSubset$`Sequence ID`)==0)
    imgtSubset$`Sequence ID`<- gsub(">","",sapply(strsplit(f$name,"\\|"),"[",1))

    if (nrow(f)>0){
      f$id <- NULL

      # write fasta and imgt files
      writeFasta(f,paste0(dir$dir.mrdarcy.in.fasta,sample,".fasta"))
      write.table(imgtSubset,paste0(dir$dir.mrdarcy.in.imgt,sample,"_Nt-sequences.txt"),
                  sep="\t",quote = F,row.names = F)
    } else {
      warning(paste0(sample," has no BCR sequences and has been discarded"))
      toExclude <- append(toExclude,sample)
    }
  }
  createMRDARCYMetadata(x, dir, toExclude)
}

selectFilterThresholds <- function(dir){

  # select thresholds: we want to maximise the number of sequences available
  # for analysis, while maximising the fragment length
  # the lower the thresholds, the more sequences rescued, but these sequences will
  # be very short

  imgtOutputFile <- paste0(dir$dir.imgt.out,"3_Nt-sequences.txt.gz")
  imgt <- fread(imgtOutputFile,fill=TRUE,sep="\t")
  imgt <- imgt[imgt$`V-DOMAIN Functionality`=="productive",]

  i <- imgt

  par(mfrow=c(2,2), mar=c(4.1, 4.1, 3.1, 0.5))
  # Assess V-region
  range <- seq(0, max(nchar(i$`V-REGION`)), 5)
  e<-numeric()
  for (r in range){
    e=append(e,sum(nchar(i$`V-REGION`)>=r))
  }
  iqr.v <- summary(nchar(i$`V-REGION`))[2]
  hist(nchar(i$`V-REGION`),main="V-region",xlab="V length")
  abline(v=iqr.v,col="red")
  plot(range,e, main=paste("Q1:",iqr.v),xlab="Cut off",ylab="Count")
  abline(v=iqr.v,col="red")
  threshold.V <- iqr.v

  # Assess with J-region
  range <- seq(0, max(nchar(i$`J-REGION`)), 5)
  e<-numeric()
  for (r in range){
    e=append(e,sum(nchar(i$`J-REGION`)>=r))
  }
  iqr.j <- summary(nchar(i$`J-REGION`))[2]
  hist(nchar(i$`J-REGION`),main="J-region",xlab="J length")
  abline(v=iqr.j,col="red")
  plot(range,e,main=paste("Q1:",iqr.j),xlab="Cut off",ylab="Count")
  abline(v=iqr.j,col="red")
  threshold.J <- iqr.j

  retained <-nrow(i[nchar(i$`J-REGION`)>=threshold.J,])

  cat("Retained: ",retained,"\n")
  cat("Lost: ",nrow(imgt)-retained,"\n")

  return (list(threshold.V=as.numeric(threshold.V),threshold.J=as.numeric(threshold.J),lost=nrow(imgt)-retained,retained=retained))

}

filterIMGToutput <- function(dir, threshold.V, threshold.J){

  # load sample-specific fasta files
  fa <- readRDS(paste0(dir$dir.imgt.in, "sampleFasta.RData"))

  # load IMGT output
  file.imgtOutput  <- grep("3_Nt-sequences.txt", list.files(dir$dir.imgt.out,full.names = T),value = T)

  # load previous dataframe
  x <- readRDS(paste0(dir$dir.imgt.in,"o.RData"))

  if(length(file.imgtOutput)==0) {
    stop(paste0("A 3_Nt-sequences.txt file was not found in: ",dir$dir.imgt.out))
  }

  imgt <- data.table::fread(file.imgtOutput,fill=TRUE,sep="\t")
  imgt <- imgt[imgt$`V-DOMAIN Functionality`!="No results",]


  #filter based on threshold lengths
  imgt <- imgt[nchar(imgt$`V-REGION`)>=threshold.V & nchar(imgt$`J-REGION`)>=threshold.J,]

  #trim sequences
  imgt$V_trim <- stringr::str_sub(imgt$`V-REGION`, start= -threshold.V)
  imgt$J_trim <- substring(imgt$`J-REGION`,1,threshold.J)
  imgt$sequence <- apply(imgt, 1, function(x)
    substring(x["V-D-J-REGION"],
              as.numeric(stringr::str_locate(x["V-D-J-REGION"],x["V_trim"])[1,1]),
              as.numeric(stringr::str_locate(x["V-D-J-REGION"],x["J_trim"])[1,2])))

  #create MRDarcy input for each sample
  for (sample in names(fa)){
    cat("Processing: ",sample,"\n")

    f <- fa[[sample]]
    f$id <- gsub(">","",sapply(strsplit(f$name,"__"),"[",1))
    imgtSubset   <- imgt[imgt$`Sequence ID` %in% f$id,]

    f <- f[match(imgtSubset$`Sequence ID`,f$id),]
    stopifnot(sum(f$id!=imgtSubset$`Sequence ID`)==0)
    imgtSubset$`Sequence ID`<- gsub(">","",sapply(strsplit(f$name,"\\|"),"[",1))

    if (nrow(f)>0){
      f$id <- NULL

      #if the sequence has been amended then this needs to be reflected in the
      #mrdarcy input
      stopifnot(sum(f$id!=imgtSubset$`Sequence ID`)==0)
      f$seq <- toupper(imgtSubset$sequence)
      imgtSubset$sequence<-NULL

      # write fasta and imgt files
      writeFasta(f,paste0(dir$dir.mrdarcy.in.fasta,sample,".fasta"))
      write.table(imgtSubset,paste0(dir$dir.mrdarcy.in.imgt,sample,"_Nt-sequences.txt"),
                  sep="\t",quote = F,row.names = F)
    }
  }

  #ok now regenerate input
  #work on trimmed files
  fastas <- list.files(dir$dir.mrdarcy.in.fasta,full.names = T)
  data <- list()

  for (f in fastas){
    sample <- gsub("Fully_reduced_","",gsub(".fasta","",basename(f)))
    donor_id <- unique(x[x$sample_id==sample,]$donor_id)

    masterFasta <- Biostrings::readDNAStringSet(f)
    masterFasta <- data.table(names=names(masterFasta),seq=as.character(masterFasta))

    id1 <- sapply(strsplit(masterFasta$names,"\\|"),"[",1)
    id2 <- sapply(strsplit(masterFasta$names,"\\|"),"[",2)

    seqID  <- sapply(strsplit(id1,"__"),"[",1)
    counts <- strsplit(sapply(strsplit(id1,"__"),"[",2),"_")
    counts <- data.frame(do.call(rbind,counts),stringsAsFactors = F)
    colnames(counts) <- unlist(strsplit(unique(id2),"_"))
    if (nrow(counts)>1){
      counts <- data.frame(apply(counts, 2, function(x) as.numeric(as.character(x))))
      counts$sequence<-masterFasta$seq
      counts$sample_id<-sample
      counts$donor_id<-donor_id

      counts <- data.table(reshape2::melt(counts,id.vars=c("donor_id","sample_id","sequence")))
      counts <- counts[,sum(value),by=.(donor_id,sample_id,sequence,variable)]
      colnames(counts)[4] <- "c_call"
      colnames(counts)[5] <- "count"
      data[[sample]] <- data.frame(counts[,c(3,4,5,2,1)],stringsAsFactors = F)
    }
  }

  data <- do.call(rbind,data)
  data <- data[data$count!=0,]

  # modify file system
  file.rename(paste0(dir$dir.base,"IMGT"),paste0(dir$dir.base,"IMGT-premerge"))
  file.rename(paste0(dir$dir.base,"input"),paste0(dir$dir.base,"input-premerge"))
  dir <- initialiseDirs(workingDirectory)
  #recreate IMGT input
  createImgtInput(data, dir)

}


createMRDARCYMetadata <- function(x, dir, toExclude){

  #Column 1: Group name
  #Column 2: Sample name
  #Column 3: Annotation file location
  #Column 4: Fasta file location
  #Column 5: Output directory
  #Column 6: Annotation file type <- use IMGT

  s <- x[,c("donor_id","sample_id")]
  s <- s[!duplicated(s),]
  s <- s[order(s$donor_id,s$sample_id),]
  col1 <- s$donor_id
  col2 <- s$sample_id
  col3 <- paste0(dir$dir.mrdarcy.in.imgt,s$sample_id,"_Nt-sequences.txt")
  col4 <- paste0(dir$dir.mrdarcy.in.fasta,s$sample_id,".fasta")
  col5 <- paste0(dir$dir.mrdarcy.out,col1,"/")
  col6 <- "IMGT"
  mdata <- data.frame(col1,col2,col3,col4,col5,col6)
  if (length(toExclude)>0){
    mdata <- mdata[!grepl(paste(toExclude,collapse="|"),mdata$col4),]
  }

  for (d in unique(mdata$col5)) dir.create(d,recursive = T,showWarnings = F)

  moreThanOneSample <- names(table(mdata$col5)[table(mdata$col5)>1])
  for (d in moreThanOneSample) dir.create(paste0(d,"CLONES"),recursive = T,showWarnings = F)
  write.table(mdata,paste0(dir$dir.base,"/metadata-mrdarcy.txt"),sep="\t",col.names = F,row.names = F,quote = F)
}


runMrDarcy <- function(dir,condaPath, MrDarcyPath){
  metadataFile <- paste0(dir$dir.base,"metadata-mrdarcy.txt")
  metadata <- read.table(metadataFile)

  donors <- unique(metadata$V1)
  numDonors<-length(donors)
  for (d in c(1:numDonors)) {

    #first merge files
    system(paste(condaPath, "run -n mrdarcy python3 ", MrDarcyPath, metadataFile, d, "1"))

    #then align sequences from largest clones
    system(paste(condaPath, "run -n mrdarcy python3 ", MrDarcyPath, metadataFile, d, "2"))

    #compress files
    dir.out <- unique(metadata[metadata$V1==donors[d],"V5"])
    mrd.files <- grep(".gz$|^CLONES$", list.files(dir.out,recursive = T,full.names = T),invert = T,value = T)
    mrd.files <- c(mrd.files, metadata[metadata$V1==donors[d],"V3"],metadata[metadata$V1==donors[d],"V4"])
    for (f in mrd.files){
      R.utils::gzip(f,destname=paste0(f,".gz"))
    }
    file.remove(paste0(grep("Tmp",mrd.files,value = T),".gz"))

  }
}





