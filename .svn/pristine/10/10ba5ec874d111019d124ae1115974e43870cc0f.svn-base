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
# setMethod(f="getByGene",
#     signature="rCGH",
#     definition=function(object, gene=NULL){
#         bygene <- object@byGene
#         if(is.null(gene)){
#             return(bygene)
#         } else{
#             gene <- toupper(gene)
#             return(bygene[which(bygene$symbol == gene),])
#             }
#         })
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
            message('dLRs: ', round(object@param$dLRs, 3))
            message('MAD: ', round(object@param$MAD, 3))
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
            cnSet$Log2Ratio <- .modelSignal(L2R, Chr, G=1:5, method="lr",
                nCores, verbose)
        }

        if("Allele.Difference" %in% colnames(cnSet)){
            if(!all(is.na(cnSet$Allele.Difference))){
                if(verbose) message("Modeling allelic Difference...")
                AD <- cnSet$Allele.Difference
                Chr <- cnSet$ChrNum
                cnSet$modelAllDif <- .modelSignal(AD, Chr, G=1:7, method="loh",
                    nCores, verbose)
            }
        }
        
        object@cnSet <- cnSet
        return(object)
    }
)

setMethod(f="EMnormalize",
            signature="rCGH",
            definition=function(object, cut=c(0.01, 0.99), G=2:6, useN=25e3,
                peakThresh=0.5, ksmooth=NA, mergeVal=0.1, Title=NA,
                verbose=TRUE){
            
            if(!.validrCGHObject(object)) return(NULL)

            op <- options()
            options(warn = -1)

            cnSet <- getCNset(object)
            LR <- cnSet$Log2Ratio
            Q <- quantile(LR, probs = cut, na.rm=TRUE)
            LR <- LR[LR>Q[1] & LR<Q[2]]

            if(is.na(ksmooth)){
                mad <- .getMAD(object)
                ksmooth <- floor(150*mad)*2 + 1
            }

            if(verbose){
                message("Smoothing param: ", ksmooth)
                message("Analyzing mixture...")
            }
            runLR <- runmed(LR, k=ksmooth)
            idx <- round(seq(1, length(runLR), len=useN))
            EM <- Mclust(runLR[idx], G=G)
            nG <- EM$G
            m <- EM$parameters$mean
            p <- EM$parameters$pro
            s <- EM$parameters$variance$sigmasq
            if(length(s)==1)
                s <- rep(s, length(m))

            if(mergeVal>0){
                if(verbose)
                    message("Merging peaks closer than ", mergeVal, " ...")
                mergedPars <- .mergePeaks(nG, length(runLR), m, s, p,
                    mergeVal, verbose)
                m <- mergedPars$m
                s <- mergedPars$s
                p <- mergedPars$p
                nG <- length(m)
            }

            peaks <- sapply(seq_len(nG), function(ii){
                d <- dnorm(runLR[idx], m[ii], sqrt(s[ii]))
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

            cnSet$Log2Ratio <- cnSet$Log2Ratio - correct
            object@cnSet <- cnSet
            object@param$EMcentralized <- TRUE
            object@param$LRcut <- cut
            object@param$ksmooth <- ksmooth
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

setMethod(f="segmentCGH",
    signature="rCGH",
    definition=function(object, Smooth=TRUE, UndoSD=NULL, minLen=10,
        nCores=NULL, verbose=TRUE){

        if(!.validrCGHObject(object))
            return(NULL)

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
            message("Computing LRR segmentation using UnodSD: ", usd)
        }
        segTable <- .computeSegmentation(L2R, Chr, Pos, sampleName,
            params, nCores)

        if(verbose)
            message("Merging segments shorter than ", minLen, "Kb.")
        segTable <- .smoothSeg(segTable, minLen)

        segTable <- .computeMedSegm(segTable, L2R)
        segTable <- .mergeLevels(segTable)
        probeValues.left <- .probeSegValue(segTable, use.medians = TRUE)
        if(verbose) message("Number of segments: ", nrow(segTable))

        params$nSegment <- nrow(segTable)
        object@param <- params
        object@segTable <- segTable
        object@cnSet <- cbind.data.frame(cnSet, Segm = probeValues.left)
        return(object)
    }
)

byGeneTable <- function(segTable, symbol = NULL, verbose = TRUE){

    # if(grepl("rCGH", class(segTable)))
    #     segTable <- getSegTable(segTable)

    .ByGene(segTable, symbol, verbose)
    }

# setMethod(f="byGeneTable",
#     signature="rCGH",
#     definition=function(object, symbol = NULL, verbose = TRUE){

#         if(!.validrCGHObject(object))
#             return(NULL)
        
#         segTable <- getSegTable(object)
#         .ByGene(segTable, symbol, verbose)

#     }
# )

