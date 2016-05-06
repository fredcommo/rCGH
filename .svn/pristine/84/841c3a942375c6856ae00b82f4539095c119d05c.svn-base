test_validCytoScan <- function(){
    filePath <- system.file("extdata", "Affy_cytoScan.cyhd.CN5.CNCHP.txt.bz2",
    	package = "rCGH")
	checkTrue(rCGH:::.validCytoScan(filePath))	
}
test_validSNP6 <- function(){
    filePath <- system.file("extdata", "Affy_snp6_cnchp.txt.bz2",
        package = "rCGH")
    checkTrue(rCGH:::.validSNP6(filePath))  
}
test_validAgilent <- function(){
    filePath <- system.file("extdata", "Agilent4x180K.txt.bz2",
        package = "rCGH")
    checkTrue(rCGH:::.validAgilent(filePath))  
}
test_Cores <- function(){
    rCGH:::.setCores(4, FALSE) == 4
}
test_ToFewCores <- function(){
    rCGH:::.setCores(0, FALSE) == 1
}
test_ToManyCores <- function(){
    rCGH:::.setCores(1e12, FALSE) == parallel:::detectCores()
}
test_NotNumCores <- function(){
    checkException(rCGH:::.setCores("a", FALSE), silent = TRUE)
}
test_Segmentation <- function(){
    L2R <- rnorm(1e3)
    Chr <- rep(1:3, each = 1e3)
    Chr[1e3] <- 4
    Pos <- 1e3
    checkException(rCGH:::.computeSegmentation(L2R, Chr, Pos, NA, NA, 1),
        silent = TRUE)
}
test_setter <- function(){
    filePath <- system.file("extdata", "Affy_cytoScan.cyhd.CN5.CNCHP.txt.bz2",
                        package = "rCGH")
    cgh <- readAffyCytoScan(filePath, sampleName = "AffyScHD")
    setInfo(cgh, "item1") <- 35
    setInfo(cgh, "item2") <- TRUE

    checkTrue(
        getInfo(cgh, "item1") == 35 &&
        getInfo(cgh, "item2") == TRUE &&
        checkException(setInfo(cgh) <- "whatever", silent = TRUE) &&
        checkException(setInfo(cgh, "whatever"), silent = TRUE)
        )
}
test_getter <- function(){
    filePath <- system.file("extdata", "Affy_cytoScan.cyhd.CN5.CNCHP.txt.bz2",
                        package = "rCGH")
    cgh <- readAffyCytoScan(filePath, sampleName = "AffyScHD")

    checkTrue(
        getInfo(cgh, "fileName") == "Affy_cytoScan.cyhd.CN5.CNCHP.txt.bz2" &&
        getInfo(cgh, "sampleName") == "AffyScHD" &&
        getInfo(cgh, "analyseDate") == format(Sys.Date(), "%Y-%m-%d") &&
        getInfo(cgh, "rCGH_version") == as.character(packageVersion("rCGH"))
        )
}
test_Pipeline <- function(){
    # Reading
    filePath <- system.file("extdata", "Affy_cytoScan.cyhd.CN5.CNCHP.txt.bz2",
        package = "rCGH")
    cgh <- readAffyCytoScan(filePath, sampleName = "sample1", labName = "myLab")
    if(!checkTrue(inherits(cgh, "rCGH")))
        stop("Error when reading Affy_cytoScan")

    # Adjusting
    cgh <- adjustSignal(cgh, nCores=1)
    pars <- getParam(cgh)
    if(length(pars)<10)
        stop("Error when Adjusting: too few parameters.")
    if(is.null(pars$dLRs))
        stop("Error when Adjusting: dLRs is missing.")
    if(is.null(pars$MAD))
        stop("Error when Adjusting: MAD is missing.")

    # Centering
    cgh <- EMnormalize(cgh)
    pars <- getParam(cgh)
    if(length(pars)<18)
        stop("Error when centralizing: too few parameters.")
    if(is.null(pars$nPeak))
        stop("Error when centralizing: nPeak is missing.")
    if(is.null(pars$correctionValue))
        stop("Error when centralizing: correctionValue is missing.")

    # Segmenting
    cgh <- segmentCGH(cgh, nCores=1)
    pars <- getParam(cgh)
    st <- getSegTable(cgh)
    if(is.null(pars$nSegment))
        stop("Error when segmenting: nSegment is missing.")
    if(nrow(st)==0)
        stop("Error when segmenting: empty segmentation table.")

    # Extracting genes
#    bg <- byGeneTable(cgh)
    bg <- byGeneTable(st)
    if(nrow(bg)==0)
        stop("Error when extracting genes: empty gene table.")

    # Plotting profile
    p <- plotProfile(cgh)
    if(is.null(p))
        stop("Error when plotting profile.")

    # Plotting LOH
    p <- plotLOH(cgh)
    if(is.null(p))
        stop("Error when plotting LOH.")

    TRUE
}
