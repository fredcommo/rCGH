################################################
## HELPER FUNCTIONS (ALL INTERNAL TO THE PACKAGE)
################################################

.createGeneDB <- function(genome){

    if(genome == "hg18")
        DB <- TxDb.Hsapiens.UCSC.hg18.knownGene
    else if(genome == "hg38")
        DB <- TxDb.Hsapiens.UCSC.hg38.knownGene
    else if(genome == "hg19")
        DB <- TxDb.Hsapiens.UCSC.hg19.knownGene
    else
        stop(sprintf("'%s' is not a supported genome", genome))

    geneDB  <- genes(DB, columns=c("gene_id"))

    if(genome == "hg19"){
        geneDB <- append(geneDB,
            GRanges(
                c("chr19", "chr6"),
                IRanges(c(15270444, 32162620), c(15311792, 32191844)),
                strand = c("-", "-"),
                gene_id = c(4854, 4855)
                )
            )
        }

    if(genome == "hg38"){
        geneDB <- append(geneDB,
            GRanges(
                c("chr19", "chr6"),
                IRanges(c(15159633, 32194843), c(15200981, 32224067)),
                strand = c("-", "-"),
                gene_id = c(4854, 4855)
                )
            )
        }

    geneDB[order(as.vector(seqnames(geneDB)), start(geneDB))]

}

################################
## VALID METHODS
################################
.validAgilent <- function(filePath){
    L <- readLines(filePath, n = 25)
    v <- any(grepl("FEATURES", L))
    if(!v)
        stop("This may be not a valid Agilent FE file.")
    TRUE
}
.validSNP6 <- function(filePath){
    L <- readLines(filePath, n = 1000)
    v <- any(grepl("ProbeSet", L))
    if(!v)
        stop("This may be not a valid Affymetrix SNP6 file.")
    TRUE
}
.validCytoScan <- function(filePath){
    L <- readLines(filePath, n = 1000)
    v <- any(grepl("ProbeSetName", L))
    if(!v)
        stop("This may be not a valid Affymetrix cytoScanHD file.")
    TRUE
}
.validrCGHObject <- function(object) {
    c0 <- inherits(object, "rCGH")
    c1 <- inherits(object, "rCGH-Agilent")
    c2 <- inherits(object, "rCGH-SNP6")
    c3 <- inherits(object, "rCGH-cytoScan")
    c4 <- inherits(object, "rCGH-Illumina")
    c5 <- inherits(object, "rCGH-generic")
    if(!c0 && !c1 && !c2 && !c3 && !c4 && !c5)
        stop("Not a valid rCGH object.\n")
    TRUE
}

###############################################
## helpers called in constructors
###############################################
.readAgilentInfo <- function(filePath, verbose){

    .getAnnot <- function(aInfo, aNames, item){
        v <- aInfo[2, which(aNames == item)]
        return(as.character(v))
        }

    if(verbose)
        message('Reading information...')

    aInfo <- read.delim(filePath, header = FALSE, fill = TRUE, skip = 1, 
        nrows = 8, stringsAsFactors = FALSE, sep = "\t")
    aNames <- as.vector(aInfo[1,])
    barCode <- .getAnnot(aInfo, aNames, "FeatureExtractor_Barcode")
    gridName <- .getAnnot(aInfo, aNames, "Grid_Name")
    scanDate <- .getAnnot(aInfo, aNames, "Scan_Date")
    scanDate <- gsub("(.*)-(.*)-(.*) (.*)+", "\\3-\\1-\\2", scanDate)
    programVersion <- .getAnnot(aInfo, aNames, "Protocol_Name")
    gridGenomicBuild <- .getAnnot(aInfo, aNames, "Grid_GenomicBuild")
    ref <- 'Dual color hybridization'

    return(
        c(barCode = barCode, gridName = gridName, 
        scanDate = as.character(scanDate), programVersion = programVersion, 
        gridGenomicBuild = gridGenomicBuild, reference = ref)
        )
}

.readAgilentMatrix <- function(filePath, verbose){
    if(verbose)
        message('Reading values...')

    arrayInfo <- readLines(filePath, n = 25)
    startAt <- grep("FEATURES", arrayInfo)
    cnSet <- read.delim(filePath, header = TRUE, skip = startAt-1, sep = "\t", 
        stringsAsFactors = FALSE,
        na.strings = c("NA", "NaN", "null", "---", ""))

    cnSet <- .curateAgilentCnSet(cnSet, verbose)
    return(cnSet)
}

.getRFlags <- function(cnSet, verbose){
    flags <- which(cnSet$rIsSaturated == 1 | 
        cnSet$rIsFeatNonUnifOL == 1 | 
        cnSet$rIsWellAboveBG == 0)
    if(verbose)
        message(length(flags), ' flagged probes on chromosome ',
            unique(cnSet$ChrNum))
    return(flags)
}

.getGFlags <- function(cnSet, verbose){
    flags <- which(cnSet$gIsSaturated == 1 | 
        cnSet$gIsFeatNonUnifOL == 1 | 
        cnSet$gIsWellAboveBG == 0)
    if(verbose)
        message(length(flags), ' flagged probes on chromosome ',
            unique(cnSet$ChrNum))
    return(flags)
}

.medFlag <- function(values, flagged, minpos, maxpos){
    mf <- sapply(flagged, function(f){
        ii <- max(minpos, f-8)
        jj <- min(maxpos, f+8)
        median(values[ii:jj], na.rm=TRUE)
        })
    return(mf)
}

.replaceFlags <- function(cnSet, verbose){
    S <- split(cnSet, cnSet$ChrNum)

    if(verbose)
        message("Red channel:")
    rflags <- sapply(S, function(subset) .getRFlags(subset, verbose))

    if(verbose)
        message("\nGreen channel:")
    gflags <- sapply(S, function(subset) .getGFlags(subset, verbose))

    newR <- lapply(names(rflags), function(chr){
        chr <- as.numeric(chr)
        flagged <- rflags[[chr]]
        tmp <- S[[chr]]
        tmp$rMedianSignal[flagged] <- .medFlag(tmp$rMedianSignal,
            flagged, 1, nrow(tmp))
        as.numeric(tmp$rMedianSignal)
        })

    newG <- lapply(names(gflags), function(chr){
        chr <- as.numeric(chr)
        flagged <- gflags[[chr]]
        tmp <- S[[chr]]
        tmp$gMedianSignal[flagged] <- .medFlag(tmp$gMedianSignal, 
            flagged, 1, nrow(tmp))
        as.numeric(tmp$gMedianSignal)
        })

    cnSet$rMedianSignal <- do.call(c, newR)
    cnSet$gMedianSignal <- do.call(c, newG)

    return(cnSet)
}

.suppressFlags <- function(object, verbose){
    if(inherits(object, "rCGH-Agilent")){
        if(verbose) message('Suppressing flagged probes...')
        cnSet <- getCNset(object)
        cnSet <- .replaceFlags(cnSet, verbose)
        flagCols <- c('gIsSaturated', 'rIsSaturated', 'gIsFeatNonUnifOL', 
            'rIsFeatNonUnifOL', 'gIsWellAboveBG', 'rIsWellAboveBG')
        cnSet <- cnSet[,-which(colnames(cnSet) %in% flagCols)]
        object@cnSet <- cnSet
    }
    if(verbose) message()

    return(object)
}

.suppressDuplicProbes <- function(object, verbose){
    ## Set to NULL for CRAN
    ProbeName <- SystematicName <- NULL
    ChrNum <- ChrStart <- ChrEnd <- NULL
    rMedianSignal <- gMedianSignal <- NULL
    cnSet <- getCNset(object)
    cnSet <- cnSet[order(cnSet$ProbeName),]
    
    if (!any(colnames(cnSet) == 'ProbeName')){
        stop('None of the columns can be identifed as ProbeNames')
    }
    
    # Removing duplicated probe ids
    dup <- duplicated(cnSet$ProbeName)
    if(any(dup)){
        if(verbose) message('Suppressing duplicated probes...')
        duplicProbes <- as.character(unique(cnSet$ProbeName[dup]))
        duplicSet <- subset(cnSet, cnSet$ProbeName %in% duplicProbes)
        medianSet <- ddply(
            .data = duplicSet,
            .variables=.(ProbeName, SystematicName, ChrNum, ChrStart, ChrEnd), 
            summarize, 
            rMedianSignal = median(rMedianSignal, na.rm=TRUE), 
            gMedianSignal = median(gMedianSignal, na.rm=TRUE))
        cnSet <- rbind.data.frame(
            cnSet[!cnSet$ProbeName %in% duplicProbes,], medianSet
        )
    }
    object@cnSet <- cnSet[order(cnSet$ChrNum, cnSet$ChrStart),]
    
    return(object)
}
.suppressDuplicLocs <- function(object, verbose){
    cnSet <- getCNset(object)
    cnSet <- cnSet[order(cnSet$ProbeName),]
    
    # Removing duplicated locs
    sp <- split(cnSet, cnSet$ChrNum)
    dup <- sapply(sp, function(tmp) any(duplicated(tmp$ChrStart)) )
    any(dup)
    if(any(dup)){
        for(d in which(dup)){
            if(verbose)
                message("Suppresing dulicated locs on chr ", d, " ...")
            tmp <- sp[[d]]
            test <- duplicated(tmp$ChrStart)
            if(any(test)){
                locs <- tmp$ChrStart[which(test)]
                for(l in locs){
                    ii <- which(tmp$ChrStart == l)
                    tmp$gMedianSignal[ii] <- mean(tmp$gMedianSignal[ii],
                        na.rm = TRUE)
                    tmp$rMedianSignal[ii] <- mean(tmp$rMedianSignal[ii],
                        na.rm = TRUE)
                }
                tmp <- tmp[-which(test),]
            }
            sp[[d]] <- tmp
        }
    }
    cnSet <- as.data.frame(do.call(rbind, sp))
    cnSet <- cnSet[order(cnSet$ChrNum, cnSet$ChrStart),]
    rownames(cnSet) <- seq(1, nrow(cnSet))
    
    object@cnSet <- cnSet[order(cnSet$ChrNum, cnSet$ChrStart),]
    
    return(object)
}

.suppressDuplic <- function(object, verbose){
    object <- .suppressDuplicProbes(object, verbose)
    object <- .suppressDuplicLocs(object, verbose)

    return(object)
}

.preset <- function(object){
    object@param <- list(ksmooth=NA, Kmax=20, Nmin=160, Mwidth=2,
        UndoSD=NULL, Alpha=1e-6)
    return(object)
}

.curateAgilentCnSet <- function(cnSet, verbose){
    keepItems <- c( "ProbeName", "SystematicName", "gMedianSignal", 
        "rMedianSignal", "gIsSaturated", "rIsSaturated", 
        "gIsFeatNonUnifOL", "rIsFeatNonUnifOL", "gIsWellAboveBG", 
        "rIsWellAboveBG")
    keepCol <- which(as.character(colnames(cnSet)) %in% keepItems)

    if(verbose) message('Filtering control probes...')
    isChr = grep('^chr[^Mrandom]*$', cnSet$SystematicName)
    cnSet <- cnSet[isChr, keepCol]

    if(verbose) message('Checking chr nums...')
    systNames <- cnSet$SystematicName
    chr <- gsub(":(.*)", "", systNames)
    chrNum <- gsub("(chr)(.*):(\\d+)-(\\d+)", "\\2", systNames)
    chrNum[chrNum=="X"] <- 23
    chrNum[chrNum=="Y"] <- 24
    chrNum <- as.numeric(chrNum)
    chrStart <- as.numeric(gsub("(chr)(.*):(\\d+)-(\\d+)", "\\3", systNames))
    chrEnd <- as.numeric(gsub("(chr)(.*):(\\d+)-(\\d+)", "\\4", systNames))
    cnSet <- cbind.data.frame(ProbeName = cnSet$ProbeName,
    SystematicName = cnSet$SystematicName,
    ChrNum=chrNum, ChrStart=chrStart, ChrEnd=chrEnd,cnSet[,-c(1:2)])
    cnSet <- cnSet[order(cnSet$ChrNum, cnSet$ChrStart), ]

    return(cnSet)
}

.readSNP6 <- function(filePath, useProbes, verbose){
    if(verbose) message('Reading information...')
    aInfo <- readLines(filePath, n = 750)
    preamble <- any(grepl("GenomeWideSNP_6", aInfo))

    if(!preamble){
        arrayType <- barCode <- gridName <- scanDate <- programVersion <- NA
        ucsc <- ensembl <- gridGenomicBuild <- ref <- NA
    } else{
        arrayType <- .getTagValue(aInfo, "ArraySet")
        barCode <- .getTagValue(aInfo, "array-barcode")
        gridName <- .getTagValue(aInfo, "state-annotation-file")
        Date <- .getTagValue(aInfo, "state-time-start")
        Date <- unlist(strsplit(Date, ' '))
        scanDate <- paste(Date[5], Date[2], Date[3])
        programVersion <- .getTagValue(aInfo, "option-program-version")
        ucsc <- .getTagValue(aInfo, "genome-version-ucsc")
        ensembl <- .getTagValue(aInfo, "genome-version-ncbi")
        gridGenomicBuild <- sprintf("%s/GRCh%s", ucsc, ensembl)
        #paste(ucsc, ensembl, sep = '/')
        ref <- .getTagValue(aInfo, "state-reference-file")
    }

    infos <- c(platform=arrayType, barCode=barCode, gridName=gridName, 
        scanDate=scanDate, programVersion=programVersion, 
        gridGenomicBuild=gridGenomicBuild, reference=ref)

    startAt <- grep("ProbeSet", aInfo)
    cnSet <- .readSNP6Matrix(filePath, startAt, useProbes, verbose)

    return(list(infos=infos, cnSet=cnSet))
}

.readSNP6Matrix <- function(filePath, startAt, useProbes, verbose){
    if(verbose) message('Reading values...')
    fullSet <- read.delim(filePath, header=TRUE, skip=startAt-1, sep="\t", 
        stringsAsFactors=FALSE, na.strings = c("NA", "NaN", "null", "---", ""))
    colnames(fullSet)[1:3] <- c("ProbeName", "ChrNum", "ChrStart")

    if(!any(grepl("^SNP_A-\\d+|^CN_\\d+", fullSet$ProbeName)))
        stop("This file doesn't look like a SNP6 file.\n")

    idx <- switch(useProbes,
        snp = grep("^SNP_A-\\d+", fullSet$ProbeName),
        cn = grep("^CN_\\d+", fullSet$ProbeName),
        all = grep("^CN_\\d+|^SNP_A-\\d+", fullSet$ProbeName))

    if(length(idx)==0){
        msg <- sprintf("No %s probes in this file.", toupper(useProbes))
        stop(msg)
    }
    cnSet <- fullSet[idx,]

    idx <- which(is.na(cnSet$ChrNum) | is.na(cnSet$ChrStart) | 
    cnSet$WeightedLog2Ratio==0 | is.na(cnSet$Log2Ratio))
    if(length(idx)>0){
        cnSet <- cnSet[-idx,]
    }

    cnSet$ChrNum <- .renameChr(cnSet$ChrNum)
    cnSet <- cnSet[order(cnSet$ChrNum, cnSet$ChrStart), ]
    rownames(cnSet) <- seq(1, nrow(cnSet))

    return(cnSet)
}

.getTagValue <- function(aInfo, tag){
    x <- aInfo[grep(tag, aInfo)]
    return(unlist(strsplit(x, '='))[2])
}

.readCytoScan <- function(filePath, useProbes, verbose){
    if(verbose) message('Reading information...')
    fileName <- gsub("(.*)/", "", filePath)
    aInfo <- readLines(filePath, n = 1000)
    preamble <- any(grepl("#%affymetrix-array-type", aInfo))

    if(!preamble){
        arrayType <- barCode <- gridName <- scanDate <- programVersion <- NA
        ucsc <- ensembl <- gridGenomicBuild <- ref <- NA
    } else{
        arrayType <- .getTagValue(aInfo, "#%affymetrix-array-type")
        barCode <- .getTagValue(aInfo, "#%affymetrix-array-barcode")
        gridName <- .getTagValue(aInfo, 
            "#%affymetrix-algorithm-param-state-annotation-file")
        scanDate <- .getTagValue(aInfo, "#%affymetrix-scan-date")
        programVersion <- .getTagValue(aInfo, "#%affymetrix-algorithm-version")
        ucsc <- .getTagValue(aInfo, "genome-version-ucsc")
        ensembl <- .getTagValue(aInfo, "genome-version-ensembl")
        gridGenomicBuild <- paste(ucsc, ensembl, sep = '/')
        ref <- .getTagValue(aInfo, 
            "#%affymetrix-algorithm-param-state-reference-file")
    }

    infos <- c(platform=arrayType, barCode=barCode, gridName=gridName, 
        scanDate=format(as.Date(scanDate), "%Y-%m-%d"), 
        programVersion=programVersion, gridGenomicBuild=gridGenomicBuild, 
        reference=ref)

    startAt <- grep("ProbeSetName", aInfo)
    cnSet <- .readCytoScanMatrix(filePath, startAt, useProbes, verbose)
    
    return(list(infos=infos, cnSet=cnSet))
}

.readCytoScanMatrix <- function(filePath, startAt, useProbes, verbose){
    if(verbose) message('Reading values...')
    fullSet <- read.delim(filePath, header=TRUE, skip=startAt-1, sep="\t", 
        stringsAsFactors=FALSE, na.strings = c("NA", "NaN", "null", "---", ""))
    colnames(fullSet) <- gsub("\\.{2}(.*)+", "", colnames(fullSet))
    colnames(fullSet)[1:3] <- c("ProbeName", "ChrNum", "ChrStart")

    if(!any(grepl("S-\\d|C-\\d", fullSet$ProbeName)))
        stop("This file doesn't look like a cytoScanHD file.\n")

    adCol <- grep("AllelicDifference", colnames(fullSet)) 
    if(length(adCol)>0)
        colnames(fullSet)[adCol] <- "Allele.Difference"

    idx <- switch(useProbes,
        snp = grep("S-\\d", fullSet$ProbeName),
        cn = grep("C-\\d", fullSet$ProbeName),
        all = grep("C-\\d|S-\\d", fullSet$ProbeName))

    if(length(idx)==0){
        msg <- sprintf("No %s probes in this file.", toupper(useProbes))
        stop(msg)
    }
    cnSet <- fullSet[idx,]


    cnSet$ChrNum <- .renameChr(cnSet$ChrNum)
    
    for(ii in 3:ncol(cnSet))
        cnSet[,ii] <- as.numeric(as.character(cnSet[,ii]))

    idx <- which(is.na(cnSet$ChrNum) | is.na(cnSet$ChrStart) | 
        cnSet$WeightedLog2Ratio==0 | is.na(cnSet$Log2Ratio))
    if(length(idx)>0){
        cnSet <- cnSet[-idx,]
    }

    cnSet <- cnSet[order(cnSet$ChrNum, cnSet$ChrStart), ]
    rownames(cnSet) <- seq(1, nrow(cnSet))

    return(cnSet)
}

.renameChr <- function(ChrNum){
    if(any(ChrNum == "X")){
        ChrNum[ChrNum == "X"] <- 23
    }
    if(any(ChrNum == "Y")){
        ChrNum[ChrNum == "Y"] <- 24
    }
    ## On new ChAS version, chr23 and 24 are coded 24 and 25, resp.
    if(any(ChrNum == 25)){
        ChrNum[ChrNum == 24] <- 23
        ChrNum[ChrNum == 25] <- 24
    }
    return(as.numeric(ChrNum))
}

.readGeneric <- function(filePath){

    message("Reading data...")
    raw <- read.delim(filePath, stringsAsFactors = FALSE)
    expectedNames <- c("ProbeName", "ChrNum", "ChrStart", "Log2Ratio")

    if(!all(expectedNames %in% colnames(raw))){
        msg1 <- "Expected colnames are: \n"
        msg2 <- sprintf("%s\t", expectedNames)
        msg3 <- "\nIn your data: \n"
        msg4 <- sprintf("%s\t", colnames(raw))
        stop(c(msg1, msg2, msg3, msg4))
    }

    ChrNum <- raw$ChrNum
    if(any(is.na(ChrNum)))
        stop("Some Chr do not have any annotation. Please, fix it or remove
            the corresponding entries.")
    if(any(ChrNum == "X"))
        ChrNum[which(ChrNum == "X")] <- 23
    if(any(ChrNum == "Y"))
        ChrNum[which(ChrNum == "Y")] <- 24

    raw$ChrNum <- as.numeric(ChrNum)
    raw <- raw[order(raw$ChrNum, raw$ChrStart), expectedNames]

    return(raw)    
}

###############################################
## helpers called in adjustSignal.R
###############################################
.CyAdjust <- function(cnSet, Ref){
    ## Recall: default ref is cy3
    if(Ref=="cy3"){
        ref <- log2(cnSet$gMedianSignal)
        test <- log2(cnSet$rMedianSignal)
    } else{
        ref <- log2(cnSet$rMedianSignal)
        test <- log2(cnSet$gMedianSignal)
    }
    M <- test - ref
    A <- (test + ref)/2
    Loess <- loessFit(M, A)$fitted
    LR <- M - Loess
    cnSet$Log2Ratio <- LR
    
    return (cnSet)
}

.GCadjust <- function(cnSet){
    agilentDB <- agilentDB
    cnSet <- cnSet[order(cnSet$ProbeName),]
    idx <- match(cnSet$ProbeName, agilentDB$ProbeID)
    cnSet <- cnSet[which(!is.na(idx)), ]
    tmpDB <- agilentDB[idx[!is.na(idx)],]

    if(!all(as.character(cnSet$ProbeName) == as.character(tmpDB$ProbeID))){
        stop("An error occured in GCadjust: probeIds do not match with Agilent 
            grid.\n")
    }

    lr = cnSet$Log2Ratio
    GC <- tmpDB$GC
    adjLr <- lr - loessFit(lr, GC)$fitted
    cnSet$Log2Ratio <- adjLr
    cnSet <- cnSet[order(cnSet$ChrNum, cnSet$ChrStart), ]

    return(cnSet)
}

.dlrs <- function(x){

    if (length(x) < 3)
        stop("Vector length>2 needed for computation")

    diffs <- diff(x)
    Q <- quantile(diffs, probs = c(.025, .975), na.rm = TRUE)
    diffs <- diffs[which(diffs>=Q[1] & diffs<=Q[2])]
    dlrs <- IQR(diffs, na.rm = TRUE)/(sqrt(2)*1.34)

    return(dlrs)
}

.MAD <- function(LR){
    tmp <- abs(LR - median(LR, na.rm = TRUE))
    return(median(tmp, na.rm = TRUE))
}

.setCores <- function(nCores, verbose){
    if(!is.null(nCores) && !is.numeric(nCores))
        stop("'nCores' must be numeric.")

    maxCores <- detectCores()
    
    if(is.null(nCores)){
        nCores <- max(1, maxCores/2)
    } else if(nCores < 1){
        if(verbose) message("The number of cores must be at least 1")
        nCores <- 1
    } else if(nCores > maxCores){
        if(verbose) message("The maximum number of cores is: ", maxCores)
        nCores <- maxCores
    }
    nCores
}

# .optimalG <- function(x, G, B = 10){
#     minx <- min(x, na.rm = TRUE)
#     maxx <- max(x, na.rm = TRUE)

#     nG <- sapply(seq_len(B), function(b){
#         model <- Mclust(x,
#             G = G,
#             prior = priorControl(scale = rep(1, max(G)),
#                 mean = seq(minx, maxx, len = max(G))),
#             control = emControl(tol = 1e-6, itmax = 1000)
#             )
#         model$G
#         })
#     if(is.list(nG))
#         nG <- do.call(c, nG)
#     optimalG <- median(nG, na.rm = TRUE)
#     message("optimal G: ", optimalG)
#     return(optimalG)
# }

.modelLOH <- function(x, G, S, verbose){
    
    if(length(x)<10)
        return(x)
    
    if(any(is.na(x))){
        NAs <- is.na(x)
        x[NAs] <- rnorm(sum(NAs), 0, 0.0033)
    }
    ii <- seq(1, length(x), by=2)
    jj <- seq(2, length(x), by=2)
    xprim <- x[ii]
    
    # optimG <- .optimalG(xprim, G)
    model <- Mclust(xprim, G = G)
    K <- model$classification
    N <- as.numeric(table(K))
    pars <- model$parameters
    m <- pars$mean
    s2 <- pars$variance$sigmasq
    if(length(s2) < length(m)) s2 <- rep(s2, length(m))
    for(k in unique(K)){
        n <- N[k]
        if(!is.na(n) && n>0){
            mu <- ifelse(abs(m[k])>1.5, 1.5*sign(m[k]), m[k])
            xprim[K==k] <- rnorm(n, mu, S)
        }
    }
    x[ii] <- xprim
    x[jj] <- xprim*(-1)
    
    return(x)
}

.modelLR <- function(x, G, S, verbose){
    
    if(length(x)<10)
        return(x)
    
    if(any(is.na(x))){
        NAs <- is.na(x)
        x[NAs] <- rnorm(sum(NAs), 0, 0.0033)
    }

    minx <- min(x)
    maxx <- max(x)
    model <- Mclust(x, G = G,
        prior = priorControl(
            scale = rep(sd(x, na.rm = TRUE), max(G)),
            mean = seq(minx, maxx, len = max(G)))
        )

    K <- model$classification
    pars <- model$parameters
    m <- pars$mean
    s2 <- pars$variance$sigmasq
    if(length(s2) < length(m)) s2 <- rep(s2, length(m))

    # newx <- rnorm(length(x), 0, sd(x)*.9)
    # for(k in unique(K)){
    #     n <- sum(K == k)
    #     if(!is.na(n) && n > 0){
    #         newx[K==k] <- newx[K==k] + m[k]
    #     }
    # }

    # return(newx)

    for(k in unique(K)){
        n <- sum(K == k)
        if(!is.na(n) && n > 0){
            set.seed(111)
            x[K==k] <- rnorm(n, m[k], sqrt(s2[k])*.95)
        }
    }

    return(x)
    
}

.modelSignal <- function(signal, chr, G,
    method=c("lr", "loh"), alpha, S, nCores, verbose){
    
    options(warn = -1)
    
    method <- match.arg(method)
    switch( method,
            lr={.model <- .modelLR},
            loh={.model <- .modelLOH}
        )
    
    # alpha <- 2e3
    # if(method == "loh"){
    #     signal <- scale(signal, scale=FALSE)
    #     alpha <- 2.5e3
    # }
    
    ss <- split(signal, chr)
    
    # Force nCores to be 1 on windows
    if(.Platform$OS.type == "windows" && nCores > 1){
        warning("nCores > 1 is not supported on windows.", immediate. = TRUE)
        nCores <- 1
    }

    newSignal <- mclapply(ss, function(sss, G, alpha){
        n <- length(sss)
        l <- max(10, ceiling(n/alpha))
        idx <- round(seq(0, n, len=l))
        S <- lapply(2:length(idx), function(jj){
            tmp <- sss[(idx[jj-1]+1):idx[jj]]
            .model(tmp, G, S, verbose)
        } )
        return(do.call(c, S))
    }, G=G, alpha=alpha, mc.cores = nCores)
    
    options(warn = 0)
    
    return( do.call(c, newSignal) )
}


###############################################
## helpers called in EMnormalize.R
###############################################

.mergePeaks <- function(nG, n, m, s, p, mergeVal, verbose){

    Raw <- c(1, 1)

    while(length(Raw)!=0){

        Mdist <- matrix(0, nG, nG)
        for(i in seq_len(nG)){
            for(j in seq_len(nG)){
                Mdist[i, j] <- abs(m[i] - m[j])
            }
        }

        diag(Mdist) <- NA
        Raw <- apply(Mdist, 1, function(x) any(x<mergeVal, na.rm=TRUE) )
        Raw <- Raw*seq(1, length(Raw))
        Raw <- Raw[Raw!=0]

        if(length(Raw)!=0){
            C1 <- Raw[1]; C2 <- Raw[2]
            mu1 <- m[C1]; mu2 <- m[C2]
            s1 <- s[C1]; s2 <- s[C2]
            p1 <- p[C1]; p2 <- p[C2]
            newmu <- (p1*mu1 + p2*mu2)/(p1 + p2)
            news <- (p1*(s1+mu1^2) + p2*(s2+mu2^2))/(p1 + p2) - newmu^2
            ##news <- ((p1*n - 1)*s1 + (p2*n - 1)*s2)/(p1*n + p2*n - 2)
            newp <- p1 + p2
            m[C1] <- newmu; s[C1] <- news; p[C1] <- newp
            m <- m[-C2]; s <- s[-C2]; p <- p[-C2]
            nG <- length(m)
        }
    }

    return(list(nG = nG, m = m, s = s, p = p))
}

###############################################
## helpers called in segmentCGH.R
###############################################
.getMAD <- function(object){
    pars <- getParam(object)
    return(pars$MAD)
}

.computeSegmentation <- function(L2R, Chr, Pos, sampleName, params, nCores){
    ksmooth <- params$ksmooth
    Kmax <- params$Kmax
    Nmin <- params$Nmin
    Mwidth <- params$Mwidth
    Alpha <- params$Alpha
    UndoSD <- params$UndoSD

    X <- cbind.data.frame(L2R=L2R, Chr=Chr, Pos=Pos)
    sX <- split(X, X$Chr)
    
    # Force nCores to be 1 on windows
    if(.Platform$OS.type == "windows" && nCores > 1){
        warning("nCores > 1 is not supported on windows.", immediate. = TRUE)
        nCores <- 1
    }

    out <- mclapply(sX,
        function(tmp, ksmooth, UndoSD, Alpha, Kmax, Nmin, Mwidth, sampleName){
            if(nrow(tmp)<2)
                stop("Too few probes to run a segmentation.")
            cna.obj <- CNA(tmp$L2R, tmp$Chr, tmp$Pos, data.type = "logratio",
                        sampleid = sampleName, presorted = TRUE)
            if(!is.na(ksmooth))
                cna.obj <- smooth.CNA(cna.obj, smooth.region = ksmooth)
            seg.cna.obj <- segment(cna.obj, undo.splits = "sdundo",
                            undo.SD = UndoSD, alpha = Alpha, kmax = Kmax,
                            nmin = Nmin, min.width = Mwidth, verbose=0)        
            return(seg.cna.obj$output)        
        }, ksmooth, UndoSD, Alpha, Kmax, Nmin, Mwidth, sampleName,
        mc.cores = nCores)

    segTable <- as.data.frame(do.call(rbind, out))
    rownames(segTable) <- seq_len(nrow(segTable))
    segTable
}

.getSegLen <- function(seg){
    abs(seg$loc.end - seg$loc.start)/1e3
}
.smoothSeg <- function(segTable, minSeg){
    minSeg <- as.numeric(minSeg)
    splitSegTables <- split(segTable, segTable$chrom)
    adjustedLocs <- lapply(splitSegTables, function(sst){
        if(nrow(sst)<2)
            return(sst)
        L <- .getSegLen(sst)
        while(any(L < minSeg)){
            i <- which(L < minSeg)[1]
            j <- .getCloser(sst, i)
            sst <- .mergeSegments(sst, i, j)
            sst <- sst[-i,]
            L <- .getSegLen(sst)
            }
        return(sst)
        })

    adjustedLocs <- as.data.frame(do.call(rbind, adjustedLocs))
    rownames(adjustedLocs) <- seq(1, nrow(adjustedLocs))

    return(adjustedLocs)
}

.getCloser <- function(sst, idx){
    if(idx==1){
        return(idx+1)
    } else if (idx==nrow(sst)){
        return(idx-1)
    } else {
        delta <- abs(sst$seg.mean[c(idx-1,idx+1)] - sst$seg.mean[idx])
        i <- ifelse(which.min(delta)==1, idx-1, idx+1)
        return(i)
    }
}

.mergeSegments <- function(sst, i, j){
    if(j<i){
        sst$loc.end[j] <- sst$loc.end[i]
    } else {
        sst$loc.start[j] <- sst$loc.start[i]
    }
    sst$num.mark[j] <- sst$num.mark[j] + sst$num.mark[i]
    sst$seg.mean[j] <- ifelse(
        abs(sst$seg.mean[i]) <= abs(sst$seg.mean[j]),
        sst$seg.mean[i], sst$seg.mean[j]
    )
    
    return(sst)
}

.computeMedSegm <- function(segTable, L2R){
    nMark <- segTable$num.mark
    e = 0
    seg.med <- Sd <- c()

    for(i in seq_len(nrow(segTable))){
        s = e + 1
        e = e + nMark[i]
        tmpL2R <- L2R[s:e]
        tmpMed <- tukey.biweight(tmpL2R[!is.na(tmpL2R)])
        tmpSd <- sd(tmpL2R[!is.na(tmpL2R)], na.rm=TRUE)
        seg.med <- c(seg.med, tmpMed)
        Sd <- c(Sd, tmpSd)
    }

    segTable <- cbind.data.frame(segTable, seg.med = seg.med, probes.Sd = Sd)
    return(segTable)
}

.mergeLevels <- function(st, thresMin=0.1, ...){
    op <- options()
    options(warn = -1)

    vObs <- st$seg.mean
    vPred <- st$seg.med
    vFit <- mergeLevels(vObs, vPred, thresMin=thresMin, verbose=0, ...)
    st$seg.med <- vFit$vecMerged
    options(op)

    return(st)
}

.probeSegValue <- function(segTable){
    segValues <- segTable$seg.med
    nMarks <- segTable$num.mark

    return(rep(segValues, times = nMarks))
}

.probeCopyValue <- function(segTable){
    copies <- segTable$estimCopy
    nMarks <- segTable$num.mark

    return(rep(copies, times = nMarks))
}

.estimateRatio <- function(p, expect){

    if(max(p, na.rm = TRUE) < 1e-3)
        return(0)

    return(2^expect[which.max(p)])
}

.estimateCopy <- function(st, ploidy, expect = log2(seq(1, 60)/2)){

    P <- lapply(1:nrow(st), function(ii){
        mi <- st$seg.med[ii]
        if(mi>5)
            return( c(rep(0, length(expect) - 1), 1) )
        si <- st$probes.Sd[ii]
        p <- sapply(expect, function(e){ dnorm(mi, e, si) })
        return(p)
        })

    ratio <- sapply(P, function(p){ .estimateRatio(p, expect) })
    copies <- ifelse(ratio == 0, 0, ifelse(ratio > 0, 2*ratio, -1/ratio))
    st$estimCopy <- copies + (ploidy - 2)

    return(st)
}

###############################################
## helpers called in byGeneTable.R
###############################################

.cmValues <- function(segTable, HG){

    cmLocs <- .locateCM(segTable, HG)
    out <- lapply(cmLocs, function(locs){
        c(segTable$seg.med[locs[1]], segTable$seg.med[locs[2]])
        })
    
    return(out)
}

.locateCM <- function(segTable, HG){

    chrs <- unique(segTable$chrom)
    cmLocs <- lapply(chrs, function(chr){

        cStart <- HG$centromerStart[HG$chrom==chr]
        cEnd <- HG$centromerEnd[HG$chrom==chr]
        tmp <- segTable[segTable$chrom==chr,]

        locStart <- which(tmp$loc.start<=cStart & cStart<=tmp$loc.end)
        if(length(locStart)==0){
            locStart <- which.min(abs(tmp$loc.end - cStart))
        }

        locEnd <- which(tmp$loc.start<=cEnd & cEnd<=tmp$loc.end)
        if(length(locEnd)==0){
            locEnd <- which.min(abs(tmp$loc.start - cEnd))
        }

        as.numeric(rownames(tmp)[c(locStart, locEnd)])
        })

    return(cmLocs)
}

.relativeLog <- function(bygene, cmValues, HG){
    relativeLog <- lapply(1:length(cmValues), function(chr){
        tmp <- bygene[bygene$chr==chr, ]
        ii <- which(tmp$chrStart < HG$centromerStart[HG$chrom==chr])
        jj <- which(tmp$chrStart > HG$centromerEnd[HG$chrom==chr])
        rl <- rep(NA, nrow(tmp))
        rl[ii] <- tmp$Log2Ratio[ii] - cmValues[[chr]][1]
        rl[jj] <- tmp$Log2Ratio[jj] - cmValues[[chr]][2]
        return(rl)
        })
    return(do.call(c, relativeLog))
}

.bygeneToSegValues <- function(bygene, segTable){
    segValues <- lapply(seq_len(nrow(bygene)), function(ii){
        gene <- bygene$symbol[ii]
        chr <- bygene$chr[ii]
        Start <- bygene$chrStart[ii]
        End <- bygene$chrEnd[ii]
        if(is.na(Start) || is.na(End))
            return(NULL)
        
        ii <- which(segTable$chrom==chr &
                        segTable$loc.start<Start &
                        Start<segTable$loc.end)
        jj <- which(segTable$chrom==chr &
                        segTable$loc.start<End &
                        End<segTable$loc.end)
        idx <- union(ii, jj)
        if(is.null(idx) || length(idx) == 0)
            return(NULL)
        
        lrr <- segTable$seg.med[idx]
        l <- abs(segTable$loc.end[idx] - segTable$loc.start[idx])/1e3
        nm <- segTable$num.mark[idx]
        copy <- segTable$estimCopy[idx]
        cbind("symbol" = gene, "Log2Ratio" = lrr, "num.mark" = nm,
            "segNum" = idx, "segLength(kb)" = l, "estimCopy" = copy)

    })
    segValues <- do.call(rbind, segValues)
    if(is.null(segValues))
        return(NULL)
    
    as.data.frame(segValues)
}

.getSegFromGene <- function(segTable, symbol, HG, geneDB, columns){

    symbol <- toupper(symbol)

    cols <- c('SYMBOL','ENTREZID', 'GENENAME', 'MAP')
    if(!all(is.na(columns)))
        if(!all(columns %in% columns(org.Hs.eg.db))){
            noMatch <- setdiff(columns, columns(org.Hs.eg.db))
            stop(sprintf("Item not allowed: %s", columns[noMatch]))
        } else{
            cols <- c(cols, columns)
        }

    suppressMessages(
        bySymbol <- try(select(org.Hs.eg.db,
                        keys = symbol,
                        keytype = 'SYMBOL',
                        columns = cols
                        ), silent = TRUE)
        )

    if(inherits(bySymbol, "try-error")){
        message(sprintf("\n'%s' not found.", symbol))
        return(NULL)
    }

    bySymbol <- bySymbol[!is.na(bySymbol$ENTREZID),]    
    entrez <- as.numeric(bySymbol$ENTREZID)
    byRange <- geneDB[geneDB$gene_id %in% entrez]
    bygene <- merge(bySymbol, as.data.frame(byRange),
                    by.x = "ENTREZID", by.y = "gene_id", all = TRUE)

    bygene <- .renameGeneList(bygene)
    segValues <- .bygeneToSegValues(bygene, segTable)
    bygene <- merge(bygene, segValues, by = "symbol", all = TRUE)
    .addGenomeLoc(bygene, HG)
}

.getGenesFromSeg <- function(chr, Start, End, geneDB, columns){
    # chr: a integer, from 1 to 24
    # Start, End: numeric. Start/End segment position (from segmentation table)

    if(chr==23) chr <- "X"
    if(chr==24) chr <- "Y"
    
    chr <- sprintf("chr%s", chr)

    ii <- which(as.vector(seqnames(geneDB)) == chr)
    jj <- intersect(ii, which(Start <= start(geneDB) & start(geneDB) <= End))
    kk <- intersect(ii, which(Start <= end(geneDB) & end(geneDB) <= End))
    idx <- unique(union(jj, kk))

    if(length(idx) == 0)
        return(NULL)

    cols <- c('SYMBOL','ENTREZID', 'GENENAME', 'MAP')
    if(!all(is.na(columns)))
        if(!all(columns %in% columns(org.Hs.eg.db))){
            noMatch <- setdiff(columns, columns(org.Hs.eg.db))
            stop(sprintf("Item not allowed: %s", columns[noMatch]))
        } else{
            cols <- c(cols, columns)
        }

    suppressMessages(
        bySymbol <- try(select(org.Hs.eg.db,
                        keys=geneDB$gene_id[idx],
                        keytype='ENTREZID',
                        columns=cols
                        ), silent = TRUE)
        )
    if(inherits(bySymbol, "try-error"))
        return(NULL)

    byRange <- as.data.frame(geneDB[idx])
        
    geneList <- merge(bySymbol, byRange,
                        by.x = "ENTREZID", by.y = "gene_id", all = TRUE)

    .renameGeneList(geneList)
}

.renameGeneList <- function(geneList){
    colnames(geneList) <- tolower(colnames(geneList))
    oldNames <- c("genename", "map", "seqnames", "start", "end")
    newNames <- c("fullName", "cytoband", "chr", "chrStart", "chrEnd")
    colnames(geneList)[colnames(geneList) %in% oldNames] <- newNames
        
    geneList <- geneList[order(geneList$symbol),]
    geneList$chr <- as.character(geneList$chr)
    .chrAsNum(geneList)
}

.chrAsNum <- function(geneList){
    geneList$chr <- gsub("chr", "", geneList$chr)
    geneList$chr <- gsub("X", "23", geneList$chr)
    geneList$chr <- gsub("Y", "24", geneList$chr)
    geneList$chr <- as.numeric(geneList$chr)
    geneList
}

.addGenomeLoc <- function(bygene, HG){

    ss <- split(bygene, bygene$chr)
    bygene <- lapply(ss, function(tmp){
        chr <- unique(tmp$chr)
        tmp$genomeStart <- tmp$chrStart + HG$cumlen[chr]
        return(tmp)
    })
    bygene <- as.data.frame(do.call(rbind, bygene))
    bygene <- bygene[order(bygene$symbol),]
    rownames(bygene) <- seq_len(nrow(bygene))
    
    # Render numeric
    renderNum <- c("chr", "chrStart", "chrEnd", "width", "Log2Ratio",
                    "num.mark", "segNum", "segLength(kb)", "relativeLog",
                    "genomeStart", "estimCopy")
    idx <- which(colnames(bygene) %in% renderNum)
    for(jj in idx){
        bygene[,jj] <- as.numeric(as.character(bygene[,jj]))
    }

    bygene
}

.ByGene <- function(segTable, symbol, genome, columns, verbose){

    hg18 <- hg18; hg19 <- hg19; hg38 <- hg38

    HG <- switch(genome,
        hg18 = hg18,
        hg19 = hg19,
        hg38 = hg38)

    geneDB <- .createGeneDB(genome)

    if(!"seg.med" %in% colnames(segTable))
        segTable$seg.med <- segTable$seg.mean

    if(!is.null(symbol))
        return(.getSegFromGene(segTable, symbol, HG, geneDB, columns))

    if(verbose) message("Creating byGene table...")
    bygene <- lapply(seq_len(nrow(segTable)), function(ii){
        chr <- segTable$chrom[ii]
        Start <- segTable$loc.start[ii]
        End <- segTable$loc.end[ii]
        l <- abs(End - Start)/1e3
        lrr <- segTable$seg.med[ii]
        nm <- segTable$num.mark[ii]
        copy <- segTable$estimCopy[ii]
        g <- .getGenesFromSeg(chr, Start, End, geneDB, columns)

        if(is.null(g))
            return(NULL)

        cbind.data.frame(g, "Log2Ratio" = lrr, "num.mark" = nm,
                        "segNum" = ii, "segLength(kb)" = round(l, 2),
                        "estimCopy" = ifelse(is.null(copy), NA, copy)
                        )
    })

    bygene <- do.call(rbind, bygene)
    cmValues <- .cmValues(segTable, HG)
    bygene$relativeLog <- .relativeLog(bygene, cmValues, HG)
    .addGenomeLoc(bygene, HG)
}

.getPatientId <- function(sampleId){
    gsub("(.*)_(.*)_(.*)+", "\\2", sampleId)
}

###############################################
## helpers called in plot functions
###############################################
.simulateLRfromST <-  function(st){
    lr <- lapply(1:nrow(st), function(ii){
        mu <- st$seg.mean[ii]
        s <- min(st$probes.Sd[ii]/5, .1)
#        s <- min(st$probes.Sd[ii]/15, .07)
        n <- round(st$num.mark[ii]/10)
        if(n>25){
            rnorm(n, mu, s)
        }
    })
    sort(do.call(c, lr))
}

.addDens <- function(x, m, s, p, best, ...){
    d <- dnorm(x, m, s)
    lines(x, d*p, lwd=3, ...)
    polygon(x, d*p, ...)
    text(m, max(d*p)+.25, labels=format(m, digits = 2),
        cex=ifelse(best, 2, 1.25))
}

.makeTitle <- function(object, gain, loss, showCopy){
    if(!showCopy){
            out <- paste(getInfo(object, 'sampleName'),
                        '-', getInfo(object, 'analysisDate'),
                        '\nGain threshold: ', round(gain, 3),
                        ' Loss threshold:', round(loss, 3)
                        )
    } else{
        out <- paste(getInfo(object, 'sampleName'),
                    '-', getInfo(object, 'analysisDate')
                    )
        }
    return(out)
}

.mainPlot <- function(segTable, cumLen, cumCentr, pCol, ylim, Title){

    N <- sum(segTable$num.mark, na.rm=TRUE)
    w <- 0.1 # N/20e3
    X <- lapply(1:nrow(segTable), function(i){
        n <- ceiling(segTable$num.mark[i]*w)
        n <- max(50, n)
        x <- seq(segTable$loc.start[i], segTable$loc.end[i], len=n)
        y <- rnorm(n, segTable$seg.med[i], segTable$probes.Sd[i]/4)
        return(cbind(chr=segTable$chrom[i], loc=x, l2r=y))
        })
    X <- as.data.frame(do.call(rbind, X))

    gPlot <- ggplot(data = X, aes_string(x="loc", y="l2r")) +
            geom_point(pch = 19, cex = 0.1, col = pCol) +
            geom_hline(yintercept = 0) +
            geom_vline(xintercept = cumLen[1:23], color = 'grey30',
                linetype = 2, size = 0.25) +
            ggtitle(Title) +
            xlab('Genomic position (bp)') + 
            ylab('Log2(Ratio)') +
            theme_bw() +
            theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                plot.margin=unit(c(0,4,4,0),"mm"),
                plot.title = element_text(lineheight=.8, size = rel(2.0),
                    face="bold"),
                axis.title.x = element_text(size = rel(1.8), angle = 00),
                axis.text.x = element_text(size = rel(1.5)),
                axis.title.y = element_text(size = rel(1.8), angle = 90),
                axis.text.y = element_text(size = rel(1.5))
                ) +
            coord_cartesian(ylim = ylim) +
            scale_y_continuous(breaks = seq(round(ylim[1]), round(ylim[2]), 
                by = 0.5)) +
            annotate(
                'text',
                x = c(-1e8, cumCentr[1:23]), y = rep(max(ylim)-0.2, 24),
                label = c("Chr", seq(1, 23)), size = 4, colour = 'grey40'
                )
    return(gPlot)

}

.addSegments <- function(gPlot, subTable, GLcolors){
    gPlot <- gPlot + 
            geom_segment(
                data = subTable,
                    aes_string(
                        x = "loc.start", xend = "loc.end",
                        y = "seg.med", yend = "seg.med"
                        ),
                    colour = GLcolors,
                    size = 2
                    )
    gPlot
}

.plotCopy <- function(segTable, cumLen, cumCentr, GLcol, Title = NULL){

    .getCol <- function(x, GLcol){
        ifelse(x>2, GLcol[1], ifelse(x<2, GLcol[2], "black"))
    }
    .addSegments <- function(gPlot, segTable, GLcol){

        idx <- which(abs(segTable$loc.end - segTable$loc.start) < 2e7)
        if(length(idx)>0){
            segTable$loc.start[idx] <- segTable$loc.start[idx] - 1e7
            segTable$loc.end[idx] <- segTable$loc.end[idx] + 1e7
        }

        gPlot <- gPlot + 
            geom_segment(
                data = segTable,
                aes_string(
                    x = "loc.start", xend = "loc.end",
                    y = "estimCopy", yend = "estimCopy"
                ),
                colour = .getCol(segTable$estimCopy, GLcol),
                size = 2
            )
        gPlot
    }

    s <- segTable$loc.start[1]
    e <- segTable$loc.end[nrow(segTable)]
    m <- min(segTable$estimCopy, na.rm = TRUE)
    M <- max(segTable$estimCopy, na.rm = TRUE)

    if(M > 10){
        M <- 10
        idx <- which(segTable$estimCopy > M)
        segTable$estimCopy[idx] <- M
        yticks <- seq(0, 10, by = 2)
        ytags <- c(sprintf("%s.0", seq(0, 8, by = 2)), "10+")
    } else {
        yticks <- seq(0, M + 2, by = 2)
        ytags <- sprintf("%s.0", yticks)
    }

    X <- data.frame(loc = c(s, e), cp = c(m, M))


    # if(M <= 10){
    #     maxy <- M + 2
    #     yticks <- seq(0, maxy, by = 2)
    # } else if(M <= 20){
    #     maxy <- M + 4
    #     yticks <- seq(0, maxy, by = 4)
    # } else{
    #     maxy <- M + 5
    #     yticks <- seq(0, maxy, by = 5)
    # }

    maxy <- M + 2
    miny <- 0
    ylim <- range(miny, maxy)

    gPlot <- ggplot(data = X, aes_string(x="loc", y="cp")) +
        geom_hline(yintercept = 2, size = 1) +
        geom_hline(yintercept = yticks,
                    color = 'grey30', linetype = 1, size = 0.25) +
        geom_vline(xintercept = cumLen[1:23], color = 'red',
                    linetype = 2, size = 0.25) +
        ggtitle(Title) +
        xlab('Genomic position (bp)') + 
        ylab('Copy number') +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.margin=unit(c(0,4,4,0),"mm"),
            plot.title = element_text(lineheight=.8, size = rel(2.0),
                                        face="bold"),
            axis.title.x = element_text(size = rel(1.8), angle = 00),
            axis.text.x = element_text(size = rel(1.5)),
            axis.title.y = element_text(size = rel(1.8), angle = 90),
            axis.text.y = element_text(size = rel(1.5))
        ) +
        coord_cartesian(ylim = ylim) +
        scale_y_continuous(breaks = yticks,
                            labels = ytags) +
        annotate(
            'text',
            x = c(-1e8, cumCentr[1:23]), y = rep(max(ylim)-0.5, 24),
            label = c("Chr", seq(1, 23)), size = 4, colour = 'grey40'
        )

    return(.addSegments(gPlot, segTable, GLcol))
}

.addTagToPlot <- function(gPlot, bg, showCopy){
    if(nrow(bg) == 0){
        message("No gene information available.")
        return(gPlot)
    }

    ymin <- min(gPlot$coordinates$limits$y)
    ymax <- max(gPlot$coordinates$limits$y)
    e <- (ymax - ymin)*.15

    if(showCopy){
        bg$Log2Ratio <- as.numeric(as.character(bg$estimCopy))
    }

    bg <- data.frame(symbol = bg$symbol,
                    genomeStart = bg$genomeStart,
                    Log2Ratio = bg$Log2Ratio,
                    Zero = rep(0, nrow(bg))
                    )

    genomeStart <- Log2Ratio <- Zero <- NULL

    gPlot2 <- gPlot

    for(ii in seq_len(nrow(bg)))
        gPlot2 <- gPlot2 +
            geom_segment(aes_string(x = bg$genomeStart[ii],
                                    xend = bg$genomeStart[ii],
                                    y = ifelse(showCopy, 2, 0),
                                    yend = bg$Log2Ratio[ii]
                                    ),
                        colour = "purple",
                        size = 0.5
                        )

    gPlot2 <- gPlot2 + 
        geom_point(data = bg, aes(x=genomeStart, y=Log2Ratio), size = 3,
                    color = "green")

#     gPlot2 <- gPlot + 
#         geom_point(data = bg, aes(x=genomeStart, y=Log2Ratio), size = 5) +
#         geom_point(data = bg, aes(x=genomeStart, y=Log2Ratio), size = 3,
#                    color = "red")

    gPlot2 <- gPlot2 + 
            annotate( 'text',
                    x = bg$genomeStart - 2.5e8,
                    y = ifelse(bg$Log2Ratio+e/2 < ymax-e, bg$Log2Ratio+e/2,
                                bg$Log2Ratio-e/2),
                    label = sprintf("%s: %s", bg$symbol,
                                    format(bg$Log2Ratio, digits = 2)
                                    ),
                    size = 6, colour = 'grey25'
                    ) +
            theme(legend.position="none")

    return(gPlot2)
}


# .addTagToPlot <- function(gPlot, bg){
#     if(nrow(bg) == 0){
#         message("No gene information available.")
#         return(gPlot)
#     }

#     genomeStart <- Log2Ratio <- NULL

#     ylim <- max(gPlot$coordinates$limits$y)
#     gPlot2 <- gPlot + 
#         geom_point(data = bg, aes(x=genomeStart, y=Log2Ratio), size = 5) +
#         geom_point(data = bg, aes(x=genomeStart, y=Log2Ratio), size = 3,
#             color = "red") +
#         annotate( 'text',
#             x = bg$genomeStart,
#             y = ifelse(bg$Log2Ratio+1<ylim, bg$Log2Ratio+1, bg$Log2Ratio-1),
#             label = sprintf("%s\nLog2R: %s", bg$symbol,
#                 format(bg$Log2Ratio, digits=2)
#                 ),
#             size = 7, colour = 'grey25'
#         ) +
#         theme(legend.position="none")
#     return(gPlot2)
# }

###############################################
## helpers called in view.R
###############################################
.convertLoc <- function(Table, HG){
    ss <- split(Table, Table$chrom)
    sconv <- lapply(ss, function(tmp){
        chr <- unique(tmp$chrom)
        tmp$loc.start <- tmp$loc.start + HG$cumlen[chr]
        tmp$loc.end <- tmp$loc.end + HG$cumlen[chr]
        return(tmp)
        })
    return(as.data.frame(do.call(rbind, sconv)))
}
