#####################
## ALL GENERICS
#####################

setGeneric("setInfo<-",
    def=function(object, item = NULL, value = NULL)
    standardGeneric("setInfo<-"))
setGeneric("recenter<-",
    def=function(object, value = 0)
    standardGeneric("recenter<-"))
setGeneric("getInfo",
    def=function(object, item = NULL) standardGeneric("getInfo"))
# setGeneric("getByGene",
#     def=function(object, gene=NULL) standardGeneric("getByGene"))
setGeneric("getCNset",
    def=function(object) standardGeneric("getCNset"))
setGeneric("getParam",
    def=function(object) standardGeneric("getParam"))
setGeneric("getSegTable",
    def=function(object, minLen = NULL) standardGeneric("getSegTable"))
setGeneric("adjustSignal",
    def=function(object, Scale = TRUE, Cy = TRUE, GC = TRUE, Ref = "cy3",
    suppOutliers = TRUE, nCores = NULL, verbose = TRUE)
    standardGeneric("adjustSignal"))
setGeneric("EMnormalize",
    def=function(object, cut = c(.01, .99), G = 2:6, useN = 25e3,
        peakThresh = 0.5, ksmooth = NA, mergeVal = 0.1, Title = NA,
        verbose = TRUE) standardGeneric("EMnormalize"))
setGeneric("segmentCGH",
    def=function(object, Smooth = TRUE, UndoSD = NULL, minLen = 10,
        nCores = NULL, verbose = TRUE) standardGeneric("segmentCGH"))
# setGeneric("byGeneTable",
#     def=function(object, symbol = NULL, verbose = TRUE)
#     standardGeneric("byGeneTable"))
setGeneric("plotDensity",
    def=function(object, breaks = NULL, Title = NULL,...)
    standardGeneric("plotDensity"))
setGeneric("plotProfile",
    def=function(object, symbol = NULL, gain = .5, loss = (-.5), minLen = 10,
        Title = NULL, ylim = NULL) standardGeneric("plotProfile"))
setGeneric("plotLOH",
    def=function(object, Title = NULL) standardGeneric("plotLOH"))
setGeneric("multiplot",
    def=function(object, symbol = NULL, gain = .5, loss = (-.5), minLen = 10,
        L = matrix(seq(1, 12)), p = c(2/3, 1/3), Title = NULL, ylim = NULL)
    standardGeneric("multiplot")
    )
setGeneric("view",
    def=function(object, browser = TRUE, ...) standardGeneric("view"))
