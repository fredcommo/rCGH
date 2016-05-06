#####################
# ALL ACCESSORS
#####################
setMethod(f="getInfo",
    signature="rCGH",
    definition=function(object, item = NULL){
    info <- object@info
    if (is.null(item)){
        return(as.data.frame(info))
    } else if (all(item %in% names(info))) {
        return(info[item])
    } else{
        idx <- !item %in% names(info)
        message(sprintf("'%s' not available.", item[idx]))
    }
})
setMethod(f="getCNset",
    signature="rCGH",
    definition=function(object){ return(object@cnSet) })
setMethod(f="getParam",
    signature="rCGH",
    definition=function(object){ return(object@param) })
setMethod(f="getSegTable",
    signature="rCGH",
    definition=function(object, minLen = NULL){
        segTable <- object@segTable
        if(!is.null(minLen))
            segTable <- .smoothSeg(segTable, minLen)
        return(segTable)
        })


#####################
## ALL METHODS
#####################

setMethod(f="adjustSignal",
    signature="rCGH",
    definition=function(object, Scale=TRUE, Cy=TRUE, GC=TRUE, Ref="cy3",
    suppOutliers=TRUE, nCores=NULL, verbose=TRUE){

        if(!.validrCGHObject(object))
            return(NULL)

        cnSet <- getCNset(object)

        if(!inherits(object, "rCGH-Agilent")){
            Cy <- FALSE
            GC <- FALSE
        }

        if(Cy){
            if(verbose){
                message("Recall you are using ", Ref, " as reference.")
                message("Cy effect adjustment...")
                }
            cnSet <- .CyAdjust(cnSet, Ref)
        }

        if(GC){
            if(verbose) message("GC% adjustment...")
            cnSet <- .GCadjust(cnSet)
        }

        object@param$CyAdjusted = Cy
        object@param$GCAdjusted = GC
        object@param$dLRs <- .dlrs(cnSet$Log2Ratio)
        object@param$MAD <- .MAD(cnSet$Log2Ratio)

        if(verbose){
            message("Log2Ratios QCs:")
            message('\tdLRs: ', round(object@param$dLRs, 3))
            message('\tMAD: ', round(object@param$MAD, 3))
            message()
        }

        if(Scale){
            if(verbose) message("Scaling...")
            cnSet$Log2Ratio <- scale(cnSet$Log2Ratio, center=FALSE)
            cnSet$Log2Ratio <- cnSet$Log2Ratio*1.2
        }
        
        nCores <- .setCores(nCores, verbose)

        if(suppOutliers){
            if(verbose) message('Signal filtering...')
            L2R <- cnSet$Log2Ratio
            Chr <- cnSet$ChrNum
            S <- NA
            cnSet$Log2Ratio <- .modelSignal(L2R, Chr, G=1:5, method="lr",
                alpha = 2e3, S, nCores, verbose)
        }

        if("Allele.Difference" %in% colnames(cnSet)){
            if(!all(is.na(cnSet$Allele.Difference))){
                if(verbose) message("Modeling allelic Difference...")
                signal <- cnSet$Allele.Difference
                chr <- cnSet$ChrNum
                S <- ifelse(inherits(object, "rCGH-oncoScan"), 0.05, 0.04)
#                a <- ifelse(inherits(object, "rCGH-oncoScan"), 1e3, 2.5e3)
                modelAllDif <- .modelSignal(signal, chr, G=1:7,
                    method="loh", alpha=2.5e3, S, nCores, verbose)
                cnSet$modelAllDif <- modelAllDif
            }
        }
        
        object@cnSet <- cnSet
        return(object)
    }
)

setMethod(f="segmentCGH",
    signature="rCGH",
    definition=function(object, Smooth=TRUE, UndoSD=NULL, minLen=10,
        nCores=NULL, verbose=TRUE){

        if(!.validrCGHObject(object))
            return(NULL)

        ploidy <- as.numeric(getInfo(object, "ploidy"))
        cnSet <- getCNset(object)
        cnSet <- cnSet[order(cnSet$ChrNum, cnSet$ChrStart),]
        params <- getParam(object)
        params$minSegLen <- minLen

        if(Smooth){
            mad <- .getMAD(object)
            params$ksmooth <- floor(150*mad)*2 + 1
        }

        if(is.null(UndoSD)){
            mad <- .getMAD(object)
            alpha <- 0.5
            if(inherits(object, "rCGH-Illumina")){
                alpha <- .95
            }
            if(inherits(object, "rCGH-oncoScan")){
                alpha <- .30
            }
            params$UndoSD <- alpha * mad^(1/2) #- .02
        } else {
            params$UndoSD <- UndoSD
        }

        L2R <- cnSet$Log2Ratio
        Chr <- cnSet$ChrNum
        Pos <- cnSet$ChrStart
        sampleName <- getInfo(object, "sampleName")

        if(is.na(sampleName)){
            sampleName <- "sample_x"
        }

        nCores <- .setCores(nCores)

        if(verbose){
            usd <- format(params$UndoSD, digits = 3)
            message("Computing LRR segmentation using UndoSD: ", usd)
        }
        segTable <- .computeSegmentation(L2R, Chr, Pos, sampleName,
            params, nCores)

        if(!is.null(minLen) && minLen < 0){
            message("'minLen', the minimal segment length can't be < 0")
            minLen <- NULL
        }

        if(!is.null(minLen)){
            if(verbose)
                message("Merging segments shorter than ", minLen, "Kb.")
            segTable <- .smoothSeg(segTable, minLen)
        }
        
        segTable <- .computeMedSegm(segTable, L2R)
        segTable <- .mergeLevels(segTable)
        segTable <- .estimateCopy(segTable, ploidy)
        probeValues <- .probeSegValue(segTable)
        if(verbose) message("Number of segments: ", nrow(segTable))

        params$nSegment <- nrow(segTable)
        object@param <- params
        object@segTable <- segTable
        object@cnSet <- cbind.data.frame(cnSet, Segm = probeValues)
        return(object)
    }
)

setMethod(f="EMnormalize",
        signature="rCGH",
        definition=function(object, G=2:6, peakThresh=0.5, mergeVal=0.1,
            Title=NA, verbose=TRUE){

        if(!.validrCGHObject(object)) return(NULL)

        op <- options()
        options(warn = -1)

        ploidy <- as.numeric(getInfo(object, "ploidy"))
        segTable <- getSegTable(object)
        if(nrow(segTable) == 0){
            stop("Please run the segmentation step before centralizing.")
        }
        simulLR <- .simulateLRfromST(segTable)
        EM <- Mclust(simulLR, G=G)
        nG <- EM$G
        m <- EM$parameters$mean
        p <- EM$parameters$pro
        s <- EM$parameters$variance$sigmasq
        if(length(s)==1)
            s <- rep(s, length(m))

        ord <- order(m)
        m <- m[ord]
        p <- p[ord]
        s <- s[ord]

        if(mergeVal>0){
            if(verbose)
                message("Merging peaks closer than ", mergeVal, " ...")
                mergedPars <- .mergePeaks(nG, simulLR, m, s, p,
                                            mergeVal, verbose)
                m <- mergedPars$m
                s <- mergedPars$s
                p <- mergedPars$p
                nG <- length(m)
            }

        peaks <- sapply(seq_len(nG), function(ii){
            d <- dnorm(simulLR, m[ii], sqrt(s[ii]))
            max(d*p[ii])
        })
        bestPeak <- which(peaks>=max(peaks)*peakThresh)[1]

        if(verbose){
            message("Gaussian mixture estimation:")
            message("n.peaks =  ", nG)
            message("\nGroup parameters:")
            for (grp in seq_len(nG)){
                msg <- sprintf(
                    "Grp %s:\nprop: %s,\tmean: %s,\tSd: %s,\tpeak height: %s",
                    grp, round(p[grp], 3), round(m[grp], 3),
                    round(sqrt(s[grp]), 3), round(peaks[grp], 3)
                )
                message(msg)
            }
            message()
        }

        correct <- m[bestPeak]
        if(verbose)
            message("Correction value:  ", round(correct, 3))

        if(is.na(Title)){
            Title <- sprintf("%s\nCorrection value = %s",
                        getInfo(object, "sampleName"), round(correct, 5))
        }

        segTable$seg.mean <- segTable$seg.mean - correct
        segTable$seg.med <- segTable$seg.med - correct
        segTable <- .estimateCopy(segTable, ploidy)

        cnSet <- getCNset(object)
        cnSet$Log2Ratio <- cnSet$Log2Ratio - correct
        cnSet$estimCopy <- .probeCopyValue(segTable)

        # Re-assign to object
        object@segTable <- segTable
        object@cnSet <- cnSet
        object@param$EMcentralized <- TRUE
        object@param$nPeak <- nG
        object@param$peakProp <- as.numeric(p)
        object@param$peakMeans <- as.numeric(m)
        object@param$peakSigmaSq <- as.numeric(s)
        object@param$centralPeak <- as.numeric(bestPeak)
        object@param$correctionValue <- as.numeric(correct)

        if(verbose)
            message("Use plotDensity() to visualize the LRR densities.")
        options(op)

        return(object)
    }
)

byGeneTable <- function(segTable, symbol = NULL,
    genome = c("hg19", "hg18", "hg38"), columns = NA, verbose = TRUE){

    genome <- match.arg(genome)
    .ByGene(segTable, symbol, genome, columns, verbose)

    }
