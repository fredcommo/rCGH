#####################
# SET FUNCTION
#####################

setReplaceMethod(f="setInfo",
    signature="rCGH",
    definition=function(object, item = NULL, value = NULL){
    if (is.null(item)){
        stop("No item specified.")
    } else if (is.null(value)){
        stop("No value specified.")
    } else {
        object@info[item] <- value
    }
    return(object)
})

setReplaceMethod(f="recenter",
                signature="rCGH",
                definition=function(object, value=0){

                    if(!.validrCGHObject(object))
                        return(NULL)

                    pars <- getParam(object)
                    peaks <- pars$peakMeans
                    nPeaks <- length(peaks)
                    if (value<0 || value>nPeaks){
                        stop("\nNon valid value.\nThere is only ", nPeaks,
                            " possible peaks.\n")
                    } else {
                        m <- peaks[value]
                        cnSet <- getCNset(object)
                        correct <- pars$correctionValue
                        cnSet$Log2Ratio <- cnSet$Log2Ratio + correct - m
                        segTable <- getSegTable(object)
                        segTable$seg.mean <- segTable$seg.mean + correct - m
                        segTable$seg.med <- segTable$seg.med + correct - m
                        object@cnSet <- cnSet
                        object@segTable <- segTable
                        object@param$centralPeak <- value
                        object@param$correctionValue <- m
                    }
                    message("Profile recentered on: ", format(m, digits=3))
                    return(object)
                }
)
