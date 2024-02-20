
#' Check all MRDARCY files available for processing
#'
#' This function checks for the presence of Count_per_sample_clusters_merged*
#' and Clone_alignment* files in input directory which are required to perform
#' a clone centrality analysis.
#'
#' @param directory Path to directory containing MRDARCY files
#' @return A list containing the locations of a Count_per_sample_clusters_merged
#' Clone_alignment files
checkAllMrDarcyFilesExist <- function(directory){
  file.countPerSampleClusters <- grep("Count_per_sample_clusters_merged",
                                      list.files(directory,full.names = T),value = T)

  file.cloneAlignment         <- grep("_passed",
                                      grep("_summary.txt",
                                           grep("Clone_alignment_",
                                                list.files(paste0(directory,"/CLONES/"),full.names = T),value = T),
                                           value=T),invert = T, value=T)

  if (length(file.countPerSampleClusters)==0 | length (file.cloneAlignment)==0){
    warning(paste0(directory, ": A Count_per_sample_clusters_merged file or Clone_alignment file was not found in this directory."))
  } else if (length(file.countPerSampleClusters)>1 | length (file.cloneAlignment)>1){
    warning(paste0(directory, ": More than one Count_per_sample_clusters_merged file or Clone_alignment file was found in this directory."))
  } else {
    return (c(file.countPerSampleClusters, file.cloneAlignment))
  }
}

#' Filter immune clones prior to network generation
#'
#' This function retains immune clones  for the presence of Count_per_sample_clusters_merged*
#' and Clone_alignment* files in input directory which are required to perform
#' a clone centrality analysis.
#'
#' @param files Path to Count_per_sample_clusters_merged and Clone_alignment files
#' @param filter.numSamples Retain clones present in >= samples
#' @param filter.numSequences Retain clones with >= sequences
#' @return A list containing the cloneID and sizes to process
filterClones <- function(files, filter.numSamples, filter.numSequences){
  file.countPerSampleClusters <- files[1]
  file.cloneAlignment         <- files[2]

  # load cluster frequency table
  p <- data.table::fread(file.countPerSampleClusters, head=TRUE, sep="\t")
  p.cols <- colnames(p)[c(2:ncol(p))]
  p <- data.frame(p,row.names = 1,stringsAsFactors = F)
  colnames(p) <- p.cols

  # Total number of sequences per sample
  sample_sum   <- apply(p, 2, function(x) sum(as.numeric(x)) )

  # retain clones with at least x unique BCR sequences
  # and present in at least x samples
  p <- data.frame(data.table::fread(file.cloneAlignment, head=TRUE, sep="\t"),stringsAsFactors = F)
  p <- p[which(p$n_samples >= filter.numSamples),]
  p <- p[which(p$size >= filter.numSequences),]

  cloneIDs   <- as.character(p[,"clone"])
  cloneSizes <- as.numeric(p[,"size"])
  return (list(sample_sum=sample_sum, cloneIDs=cloneIDs, cloneSizes=cloneSizes))
}


checkMetadata <- function(metadata,samples){
  if (!(all(colnames(metadata) ==c("sample_id","label","colour"))==T)){
    stop("Metadata data frame should have three columns, with column names: sample_id, label, colour")
  }
  if (sum(sort(metadata$sample)!=sort(names(samples)))!=0){
    stop("Metadata sample names should match sample names in MRDARCY files")
  }
}



# align clones
alignClones <- function(inputDirectory, cloneID, cloneSize, ind, samples, plotGraph, metadata){

  bcr_dictionary <- list()
  # Find alignment file, file structure is:
  # inputDirectory,"/CLONES/Clone_alignment_",g,"__A_align_", clonesToAnalyse$cloneIDs[ind],"_", clonesToAnalyse$cloneSizes[ind],"_A.fasta.gz"
  # and filter for outliers and positions

  aln_file <- grep(paste0("__A_align_", cloneID,"_", cloneSize,"_A.fasta"),
                   grep("Clone_alignment_",
                        list.files(paste0(inputDirectory,"/CLONES/"),full.names = T),
                        value = T), value=T)

  if (length(aln_file)<1){
    cat("An alignment file was not found for clone:",cloneID, "size:",cloneSize,"\n")
  } else {

    # Alignment file found, load.
    seqs     <- ape::read.dna(aln_file, format = "fasta", as.matrix = T,as.character = T)

    # trim sequence edges
    # 95% of sequences should have an aligned base
    n_bases    <- apply(seqs, 2, function(x) length(which(x!='-')) )
    seqs       <- seqs[,which(n_bases!=0)]
    n_basesCol <- apply(seqs, 2, function(x) length(which(x!='-')) )
    min        <- min(which(n_basesCol >= nrow(seqs)*0.95))
    max        <- max(which(n_basesCol >= nrow(seqs)*0.95))
    seqs_sub1  <- seqs[,c(min:max)]
    n_basesRow <- apply(seqs_sub1, 1, function(x) length(which(x!='-')) )

    # proceed only if aligned sequence > 80 bases
    if (max-min > 80){

      # proceed only if there are three or more unique BCRs in the clone
      if(length(unique(apply( seqs_sub1 ,1, paste, collapse = ""))) >=3){
        ids  <- rownames(seqs_sub1)
        type <- sapply(strsplit(ids,"\\|"),"[",3)
        names(type) <- ids
        # type contains the sampleID for the BCR seq

        # create a distance matrix: how many bases are different?
        dm <- matrix(0, nrow=length(ids), ncol=length(ids), dimnames=list(ids,ids))
        for (i1 in c(1:length(ids))){
          for (i2 in c(i1:length(ids))){
            if(i1<i2){
              w <- intersect(which(seqs_sub1[i1,]!='-'), which(seqs_sub1[i2,]!='-'))
              d <- length(which(seqs_sub1[i1,w]!=seqs_sub1[i2,w]))
              dm[i1,i2] <- d
              dm[i2,i1] <- d
            }
          }
        }

        #group identical sequences together
        done       <- NULL
        groups     <- NULL
        ids_unique <- NULL
        for (i in c(1:nrow(dm))){
          d <- dm[i,]
          # which sequences are identical and have not been processed
          w <- setdiff(names(which(d==0)), done)
          if (length(w)>0) {
            groups <- c(groups, list(w))
            done   <- c(done, w)
            rep    <- names(which(n_basesRow[w]==max(n_basesRow[w]))[1])
            ids_unique <- c(ids_unique, rep)
          }
        }
        rm(dm)

        # generate a matrix with unique sequences only
        seq_unique           <- seqs_sub1[ids_unique,,drop=F]
        seq_unique_collapse  <- apply(seq_unique, 1, paste, collapse = "")
        rownames(seq_unique) <- seq_unique_collapse

        mat_counts <- matrix(data = 0, nrow = length(seq_unique_collapse),
                             ncol = length(samples),
                             dimnames = list(seq_unique_collapse, sort(names(samples))))
        all_counts <- mat_counts

        bcr_dictionary<-data.frame()

        #all_counts will contain a matrix of unique BCRs/sample
        for (i in c(1:length(seq_unique_collapse))) {
          w <- groups[[i]]
          t <- table(type[w])
          all_counts[seq_unique_collapse[i], names(t)] <- t
          bcr_dictionary <- rbind(bcr_dictionary,
                                  cbind(bcrSeq=as.character(seq_unique_collapse[i]),
                                        bcrID=list(unique(lapply(strsplit(groups[[i]],"__"),"[",1)))))

          q <- type[w]
          q <- strsplit(names(q),"\\|")
          s <- strsplit(sapply(strsplit(sapply(q,"[",1),"__"),"[",2),"_")
          s <- sapply(s, function(x) sum(as.numeric(x)))
          q <- data.frame(sample=sapply(q,"[",3),sum=s)
          q <- aggregate(q$sum, by=list(sample=q$sample), FUN=sum)
          h <- q$x
          names(h)<-q$sample
          mat_counts[seq_unique_collapse[i], names(h)] <- h
          mat_counts[seq_unique_collapse[i],] <- mat_counts[seq_unique_collapse[i],] / samples *100

        }
        seqs_sub1 <- ape::as.alignment(seq_unique)
        seqs_sub1 <- phangorn:::as.phyDat(seqs_sub1, type = "DNA")

        try({
          dm        <- phangorn::dist.hamming(seqs_sub1, ratio = FALSE)
          treeNJ    <- phangorn::NJ(dm)
          mst <- ape::mst(dm)


          if (plotGraph){
            fit    <- phangorn::pml(treeNJ, data = seqs_sub1)
            fitJC  <- phangorn::optim.pml(fit, TRUE)
            tree   <- fitJC$tree
            ids    <- tree$tip.label
            bsplot <- plot(tree, type = "unrooted", use.edge.length = TRUE,
                           node.pos = NULL, show.tip.label = TRUE,
                           show.node.label = FALSE, edge.color = "black",
                           edge.width = 1, edge.lty = 1, font = 3, cex = 0.001,
                           adj = NULL, srt = 0, no.margin = FALSE,
                           label.offset = 0, plot =TRUE, main=cloneID)

            size <- (rowSums(mat_counts)/100)^0.4
            size <- size*3

            cl <- metadata[metadata$sample_id %in% colnames(mat_counts),]
            stopifnot(sum(cl$sample_id!=colnames(mat_counts))==0)
            ape::tiplabels(pie = mat_counts, cex = size, piecol = cl$colour)
            legend("topleft", cl$label, pch = 21,cex = 1, bty = "n",
                   pt.bg = cl$colour, pt.lwd = 0.2)
          }

          #generate igraph representation and calculate degree per BCR
          mst <- ape::mst(dm)
          gr  <- igraph::graph_from_adjacency_matrix(mst)
          gr  <- igraph::as.undirected(gr,mode = c("collapse"))

          gr  <- data.frame(igraph::degree(gr,loops = F))
          gr$bcr <- rownames(gr)

          a          <- data.frame(all_counts)
          a$numSites <- apply(a,1, function(x) sum(x>0))
          a$bcr      <- rownames(a)

          gr <- merge(a,gr,by="bcr")
          gr$clusterIndex <- ind
          gr <- merge(gr,bcr_dictionary,by.x="bcr",by.y="bcrSeq")
          return (gr) #degree statistics


        }, silent = T)



      }
    }
  }
}


#' @export
runCentralityAnalysis <- function(dirs,
                                  filter.numSamples, filter.numSequences,
                                  plotGraph, metadata){


  directories <- list.dirs(dirs$dir.mrdarcy.out,recursive = F)

  for (mrd in directories){
    donorID <- basename(mrd)
    cat("Analysing donor: ",donorID,"\n")
    files <- checkAllMrDarcyFilesExist(mrd)
    if (length(files)>1){
      clonesToAnalyse <- filterClones(files, filter.numSamples, filter.numSequences)

      if (length(clonesToAnalyse$cloneIDs)!=0){
        # Initiate plotting device if visualisation required
        if (plotGraph){
          pdfHeight <- 20
          file.network.plots <- paste0(dirs$dir.base, donorID, "_network_align.pdf")
          pdf(file = file.network.plots, height = pdfHeight, width = pdfHeight*1.25)
          par(mfrow= c(1,1), mar = c(0,0.5,1,1))
        }

        degreeStats<-list()

        if (length(clonesToAnalyse$cloneIDs)>0){
          for (ind in c(1:length(clonesToAnalyse$cloneIDs))){
            #cat(ind,"\n")
            checkMetadata(metadata,names(clonesToAnalyse$sample_sum))
            ds <- alignClones(mrd,
                              clonesToAnalyse$cloneIDs[ind], clonesToAnalyse$cloneSizes[ind], ind,
                              clonesToAnalyse$sample_sum,plotGraph, metadata)
            if (!is.null(ds)) {
              if (class(ds)=="data.frame"){
                degreeStats[[clonesToAnalyse$cloneIDs[ind]]] <- ds
              }
            }
          }
        }
        if (plotGraph) dev.off()
        if (length(degreeStats)>0){
          degreeStats <- do.call(rbind,degreeStats)
          colnames(degreeStats)[colnames(degreeStats)=="igraph..degree.gr..loops...F."]<-"degree"
          saveRDS(degreeStats,paste0(dirs$dir.base,donorID,"-network-degree.RData"))
        }
      }
    }
  }
  # when all done, merge
  degreeList <- list()
  for (mrd in directories){
    donorID <- basename(mrd)
    file <- paste0(dirs$dir.base,donorID,"-network-degree.RData")
    if (file.exists(file)){
      f <- readRDS(file)
      degreeList[[donorID]] <- f
      file.remove(file)
    }
  }

  #present across >=2 samples and have >=10 unique BCRs in that clone
  saveRDS(degreeList, paste0(dirs$dir.base,"Network-degree.RData"))

}


