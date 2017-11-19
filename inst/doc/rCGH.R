## ----style-knitr, eval=TRUE, echo=FALSE, results="asis"---------------------------------
BiocStyle::latex()

## ----setup, include=FALSE-----------------------------------------------------
options(width = 80)
require(knitr)
opts_chunk$set(dev='png', prompt=TRUE, comment=NA, tidy=FALSE)

## ----readFiles, warning=FALSE, message=FALSE----------------------------------
library(rCGH)
filePath <- system.file("extdata", "Affy_cytoScan.cyhd.CN5.CNCHP.txt.bz2",
    package = "rCGH")
cgh <- readAffyCytoScan(filePath, sampleName = "CSc-Example",
                        labName = "myLab")

## ----cgh, warning=FALSE, message=FALSE----------------------------------------
cgh

## ----addInfo------------------------------------------------------------------
setInfo(cgh, "item1") <- 35
setInfo(cgh, "item2") <- TRUE
setInfo(cgh, "item3") <- "someComment"

## ----getInfo------------------------------------------------------------------
getInfo(cgh)
getInfo(cgh, c("item1", "item3"))

## ----adjustSignal-------------------------------------------------------------
cgh <- adjustSignal(cgh, nCores=1)

## ----SegmentCGH---------------------------------------------------------------
cgh <- segmentCGH(cgh, nCores=1)
segTable <- getSegTable(cgh)

## ----segTable-----------------------------------------------------------------
head(segTable)

## ----EMnormalize--------------------------------------------------------------
cgh <- EMnormalize(cgh)

## ----plotDensity, fig.width=7, fig.height=5, fig.show='hide'------------------
plotDensity(cgh)

## ----byGeneTable--------------------------------------------------------------
geneTable <- byGeneTable(segTable)
head(geneTable, n=3)

## ----byGeneTable2-------------------------------------------------------------
byGeneTable(segTable, "erbb2", genome = "hg19")[,1:6]
byGeneTable(segTable, "erbb2", genome = "hg18")[,1:6]

## ----getParams----------------------------------------------------------------
getParam(cgh)[1:3]

## ----getProfile, fig.width=7.7, fig.height=9.5, fig.show='hide'---------------
multiplot(cgh, symbol = c("egfr", "erbb2"))

## ----recenter, fig.width=7.5, fig.height=4, fig.show='hide'-------------------
# Recentering on peak #2
recenter(cgh) <- 2
plotProfile(cgh, symbol = c("egfr", "erbb2"))

## ----view, eval=FALSE, echo=TRUE----------------------------------------------
#  view(cgh)

## ----exampleFiles-------------------------------------------------------------
list.files(system.file("extdata", package = "rCGH"))

## ----session------------------------------------------------------------------
toLatex(sessionInfo())

