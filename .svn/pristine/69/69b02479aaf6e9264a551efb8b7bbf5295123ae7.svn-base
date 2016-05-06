################################
## Build a Agilent object
################################
readAgilent <- function(filePath, sampleName=NA, labName=NA, supFlags=TRUE,
    genome = c("hg19", "hg18", "hg38"), verbose=TRUE){

    if(!.validAgilent(filePath))
        return(NULL)

    fileName <- gsub("(.*)/", "", filePath)
    genome <- match.arg(genome)
    object <- new(
        "rCGH-Agilent",
        info = c(fileName=fileName, sampleName=sampleName, 
            labName=labName, platform='Agilent', suppressFlags=supFlags,
            genome = genome)
        )
    object@info <- c(object@info, .readAgilentInfo(filePath, verbose))
    object@cnSet <- .readAgilentMatrix(filePath, verbose)

    if(supFlags){
        object <- .suppressFlags(object, verbose)
    }

    object <- .suppressDuplic(object, verbose)
    object <- .preset(object)
    setInfo(object, "rCGH_version") <- as.character(packageVersion("rCGH"))

#    .createGeneDB(genome)

    if(verbose)
        message("Genome build: ", genome)

    return (object)
}

############################
## Build a SNP6 object
############################
readAffySNP6 <- function(filePath, sampleName=NA, labName=NA,
    useProbes=c("snp", "cn", "all"), genome = c("hg19", "hg18", "hg38"),
    verbose=TRUE){

    if(!.validSNP6(filePath))
        return(NULL)
    
    useProbes <- match.arg(useProbes)
    genome <- match.arg(genome)

    if(verbose)
        message(toupper(useProbes), " probes will be used.")

    fileName <- gsub("(.*)/", "", filePath)
    object <- new(
        "rCGH-SNP6",
        info = c(fileName=fileName, sampleName=sampleName,
            labName=labName, usedProbes=useProbes, genome = genome)
        )

    affyData <- .readSNP6(filePath, useProbes, verbose)
    object@info <- c(object@info, affyData$infos)
    object@cnSet <- affyData$cnSet

    if(verbose)
        message("Adding presettings...")

    object <- .preset(object)
    setInfo(object, "rCGH_version") <- as.character(packageVersion("rCGH"))

#    .createGeneDB(genome)

    if(verbose)
        message("Genome build: ", genome)

    return (object)
}

################################
## Build a AffyCytoScan object
################################

readAffyCytoScan <- function(filePath, sampleName=NA, labName=NA, 
    useProbes=c("snp", "cn", "all"), genome = c("hg19", "hg18", "hg38"),
    verbose=TRUE){

    if(!.validCytoScan(filePath))
        return(NULL)

    useProbes <- match.arg(useProbes)
    genome <- match.arg(genome)

    if(verbose)
        message(toupper(useProbes), " probes will be used.")

    fileName <- gsub("(.*)/", "", filePath)
    object <- new(
        "rCGH-cytoScan",
        info = c(fileName=fileName, sampleName=sampleName,
            labName=labName, usedProbes=useProbes, genome = genome)
        )

    affyData <- .readCytoScan(filePath, useProbes, verbose)
    object@info <- c(object@info, affyData$infos)
    object@cnSet <- affyData$cnSet

    if(verbose) 
        message("Adding presettings...")

    object <- .preset(object)
    setInfo(object, "rCGH_version") <- as.character(packageVersion("rCGH"))

#    .createGeneDB(genome)

    if(verbose)
        message("Genome build: ", genome)

    return (object)
}


################################
## Build a generic object
################################
readGeneric <- function(filePath, sampleName=NA, labName=NA,
    genome = c("hg19", "hg18", "hg38"), verbose=TRUE){

    fileName <- gsub("(.*)/", "", filePath)
    genome <- match.arg(genome)
    object <- new(
        "rCGH-generic",
        info = c(fileName=fileName, sampleName=sampleName,
                labName=labName, usedProbes=NA, genome = genome)
    )
    
    object@cnSet <- .readGeneric(filePath)

    if(verbose) 
        message("Adding presettings...")

    object <- .preset(object)
    setInfo(object, "rCGH_version") <- as.character(packageVersion("rCGH"))
    
#    .createGeneDB(genome)

    if(verbose)
        message("Genome build: ", genome)

    return(object)
}
