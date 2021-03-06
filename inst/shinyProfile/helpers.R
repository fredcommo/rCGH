###########################
# PLOT HELPER FUNCTIONS
###########################

# Segmentation table

.getName <- function(segTable){
    return( unique(as.character(segTable$ID)) )
}
.selectChr <- function(Choice){
    if(Choice == 'All') return(1:23)
    else return(as.numeric(Choice))
}
.mainPlot <- function(segTable, s=3){

    if(nrow(segTable)<1)
        return(NULL)

    N <- sum(segTable$num.mark, na.rm=TRUE)
    w <- N/20e3

    X <- lapply(1:nrow(segTable), function(i){
        n <- ceiling(segTable$num.mark[i]/w)
        n <- max(50, n)
        x <- seq(segTable$loc.start[i], segTable$loc.end[i], len=n)
        y <- rnorm(n, segTable$seg.med[i], segTable$probes.Sd[i]/s)
        return(cbind(loc=x, l2r=y))
    })
    X <- as.data.frame(do.call(rbind, X))

    gPlot <- ggplot(data=X, aes(x=loc, y=l2r)) +
        geom_point(pch = 19, cex = 0.15, col = 'grey50') +
        geom_hline(yintercept = 0) +
        xlab('Genomic position (bp)') +
        ylab('Log2(Ratio)') +
        theme_bw() +
        theme(
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.margin=unit(c(4,4,4,0),"mm"),
            axis.text=element_text(size=12),
            axis.title=element_text(size=18),
            axis.title.x=element_text(vjust=-.5),
            axis.title.y=element_text(hjust=.5, vjust = .75),
            plot.title = element_text(lineheight=1.1, size = 25, vjust = 2, face="bold")
            )
        return(gPlot)
}
.updateScale <- function(gPlot, chr, hg19, Ymin, Ymax){
    ymin <- min(-2.5, gPlot$data$l2r, na.rm=TRUE) - 1
    ymin <- ymin*Ymin
    ymax <- max(gPlot$data$l2r, na.rm=TRUE) + 1
    ymax <- ymax*Ymax
    if(chr != "All"){
        chr <- as.numeric(chr)
        cumLen <- cumsum(as.numeric(hg19$length))
        xmin <- ifelse(chr==1, 0, cumLen[chr-1])
        xmax <- cumLen[chr]
        gPlot <- gPlot + 
            coord_cartesian(xlim = range(xmin, xmax),
                ylim=range(ymin, ymax)) +
            scale_y_continuous(breaks = seq(round(ymin), round(ymax), by = 0.5)) +
            geom_point(pch = 19, cex = 1, col = 'grey50')
    }
    else
        gPlot <- gPlot + 
            coord_cartesian(xlim = range(-0.5e8, hg19$cumlen[24]+0.5e8),
                ylim=range(ymin, ymax)) +
            scale_y_continuous(breaks = seq(round(ymin), round(ymax), by = 0.5))

    return(gPlot)
}
.addSegments <- function(gPlot, segTable, chr, gain, loss, segLen, GLcols){
    if(chr=="All"){
        chr <- 1:23
    } else{
        chr <- as.numeric(chr)
    }

    if(segLen %in% c("All", "")){
        segLen <- Inf
    } else{
        segLen <- as.numeric(segLen)
    }

    gainCol <- rgb(0, 0.45, 1, 1)
    lossCol <- "red3"
    if(grepl("red/blue", GLcols)){
        gainCol <- "red3"
        lossCol <- rgb(0, 0.45, 1, 1)
    }

    segTable <- segTable[which(segTable$chrom %in% chr),]
    L <- abs(segTable$loc.end - segTable$loc.start)/1e6
    idx <- which((segTable$seg.med<= loss | segTable$seg.med>= gain) & 
        L <= segLen)

    if(length(idx)>0){
        subTable <- segTable[idx,]
        GLcolors <- ifelse(subTable$seg.med <= loss, lossCol,
            ifelse(subTable$seg.med >= gain, gainCol, "black")
            )
        gPlot <- gPlot+
            geom_segment(
                data=subTable,
                aes(x=loc.start, xend=loc.end, y=seg.med, yend=seg.med),
                colour=GLcolors, size=2
                )
        }

    return(gPlot)   
}
.addChr <- function(gPlot, chr, hg){
    if(chr=="All") chr <- as.numeric(1:23)
    else chr <- as.numeric(chr)

    ylim <- gPlot$coordinates$limits$y
    if(is.null(ylim)){
        ylim <- range(gPlot$data$l2r)
    }
    cumCentr <- 1/2*hg$length+hg$cumlen
    gPlot <- gPlot+
        geom_vline(xintercept = hg$cumlen[chr], color = 'grey30',
            linetype = 2, size = 0.25) +
        annotate('text', x=cumCentr[chr], y=rep(max(ylim, na.rm=TRUE)*.95,
            length(chr)), label=chr, size = 4, colour = 'grey40')
    return(gPlot)     
}
.addTitle <- function(gPlot, sampleName, gain, loss){
    Title = paste(sampleName, '\nGain threshold: ', round(gain, 3),
        ' Loss threshold:', round(loss, 3))
    gPlot <- gPlot + ggtitle(Title)
    return(gPlot)  
}
.addTag <- function(gPlot, geneAnnot, Yexpand, gain, loss){

    ylim <- gPlot$coordinates$limits$y
    ymin <- min(ylim)
    ymax <- max(ylim)
    symbol <- as.character(geneAnnot$symbol)
    x <- geneAnnot$genomeStart
    lr <- geneAnnot$Log2Ratio
    yLabel <- ifelse((lr+1.75)<ymax, lr+1.25, lr-1.25)

    if(is.na(lr)) return(gPlot)

    Col <- "grey25"
    gPlot <- gPlot +
        annotate("text",
            x=max(x, 2e8), y=yLabel,
            label=paste0(symbol, '\n(Log2R = ', round(lr, 3), ')'),
            fontface="bold", colour=Col) +
        geom_point(x=x, y=lr, size=6, pch=19, colour="black") +
        geom_point(x=x, y=lr, size=5, pch=19, colour="antiquewhite") +
        geom_point(x=x, y=lr, size=4, pch=19, colour="darkorchid4") +
        geom_point(x=x, y=lr, size=1, pch=19, colour="cyan")

    return(gPlot)
}
.resizeLOH <- function(lohPlot, chr="All"){
    if(is.null(lohPlot))
        return(NULL)

    lohPlot <- lohPlot +
        theme_bw() +
        theme(
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.margin=unit(c(0,4,4,0),"mm"),
            axis.text=element_text(size=12),
            axis.title=element_text(size=18),
            axis.title.x=element_text(vjust=-.5),
            axis.title.y=element_text(hjust=.5, vjust = .75)
            )

    if(chr == "All"){
        lohPlot <-  lohPlot + 
            coord_cartesian(xlim=range(-0.5e8, hg19$cumlen[24]+0.5e8),
                ylim=range(-2.0, 2.0)) +
            ggtitle("")
    } else{
        chr <- as.numeric(chr)
        locs <- c(hg19$cumlen[chr], hg19$cumlen[chr+1])
        m <- sum(locs)/2
        lohPlot <-  lohPlot + 
            geom_point(pch = 19, cex = 1, col = rgb(0,0,0,.75)) +
            coord_cartesian(xlim = locs,
                ylim = range(-2.0, 2.0)) +
            ggtitle("")
        }

    return(lohPlot)
}
.geneOfInt <- function(symbol, geneTable){
    symbol <- toupper(symbol)
    tmp <- geneTable[which(geneTable$symbol == symbol),]
    
    if(nrow(tmp)==0)
        return(NULL)
        
    tmp
}
.renderLink <- function(uid){
    sprintf("<a href=\"http://www.ncbi.nlm.nih.gov/gene/?term=%s[uid]\" 
    target=\"_blank\" style=\"font-size:18px; \">%s</a>", uid, uid)
}

###########################
# MERGING SEGMENTS
.getSegLen <- function(seg){
    abs(seg$loc.end - seg$loc.start)/1e3
}
.smoothSeg <- function(segTable, minSeg){
    minSeg <- as.numeric(minSeg)

    if(is.na(minSeg) || minSeg <= 10)
    	return(segTable)

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
.getCloser <- function(segTable, idx){
    if(idx==1){
        return(idx+1)
    } else if (idx==nrow(segTable)){
        return(idx-1)
    } else {
        delta <- abs(segTable$seg.med[c(idx-1,idx+1)] - segTable$seg.med[idx])
        i <- ifelse(which.min(delta)==1, idx-1, idx+1)
        return(i)
    }
}
.mergeSegments <- function(segTable, i, j){
    if(j<i){
        segTable$loc.end[j] <- segTable$loc.end[i]
    } else {
        segTable$loc.start[j] <- segTable$loc.start[i]
    }
    segTable$num.mark[j] <- segTable$num.mark[j] + segTable$num.mark[i]
    segTable$seg.mean[j] <- ifelse(
        abs(segTable$seg.mean[i]) <= abs(segTable$seg.mean[j]),
        segTable$seg.mean[i], segTable$seg.mean[j]
        )
    segTable$seg.med[j] <- ifelse(
        abs(segTable$seg.med[i]) <= abs(segTable$seg.med[j]),
        segTable$seg.med[i], segTable$seg.med[j]
        )

    return(segTable)
}

# End helper functions
###########################
# GENES ANNOTATIONS HELPER FUNCTIONS
###########################

require(TxDb.Hsapiens.UCSC.hg18.knownGene)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
require( org.Hs.eg.db)

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

#    assign("geneDB", geneDB, envir = .GlobalEnv)
}

.getGenesFromSeg <- function(chr, Start, End, DB){
    # chr: a integer, from 1 to 24
    # Start, End: numeric. Start/End segment position (from segmentation table)

#    geneDB <- geneDB

    if(chr==23) chr <- "X"
    if(chr==24) chr <- "Y"
    
    chr <- sprintf("chr%s", chr)

    ii <- which(as.vector(seqnames(DB)) == chr)
    jj <- intersect(ii, which(Start <= start(DB) & start(DB) <= End))
    kk <- intersect(ii, which(Start <= end(DB) & end(DB) <= End))
    idx <- unique(union(jj, kk))

    if(length(idx) == 0)
        return(NULL)

    suppressMessages(
        bySymbol <- try(select(org.Hs.eg.db,
                        keys=DB$gene_id[idx],
                        keytype='ENTREZID',
                        columns=c('SYMBOL', 'GENENAME', 'MAP')
                        ), silent=TRUE)
        )
    
    if(inherits(bySymbol, "try-error"))
        return(NULL)
    
    byRange <- as.data.frame(DB[idx])
        
    geneList <- merge(bySymbol, byRange,
                        by.x = "ENTREZID", by.y = "gene_id", all = TRUE)
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
.addGenomeLoc <- function(bygene, hg){
#    hg <- hg19
    ss <- split(bygene, bygene$chr)
    bygene <- lapply(ss, function(tmp){
        chr <- unique(tmp$chr)
        tmp$genomeStart <- tmp$chrStart + hg$cumlen[chr]
        return(tmp)
    })
    bygene <- as.data.frame(do.call(rbind, bygene))
    bygene <- bygene[order(bygene$symbol),]
    rownames(bygene) <- seq_len(nrow(bygene))
    bygene
}
.genometoChrLoc <- function(segTable, hg){

#    hg <- hg19

    splitTable <- split(segTable, segTable$chrom)
    newTable <- lapply(splitTable, function(tmp){
        chr <- unique(tmp$chrom)
        tmp$loc.start <- tmp$loc.start - hg$cumlen[chr]
        tmp$loc.end <- tmp$loc.end - hg$cumlen[chr]
        tmp
    })
    do.call(rbind.data.frame, newTable)
}
.filterBygene <- function(bg, chr, greater, lower, segLen){

    ii <- which(bg$chr %in% chr)
    jj <- which(bg$Log2Ratio>=greater | bg$Log2Ratio<=lower)
    kk <- which(bg$"segLength(kb)"/1e3 <= segLen)
    idx <- Reduce(intersect, list(ii, jj, kk))
    bg[idx,]

}
ByGene <- function(st, hg, geneDB){
    if(is.null(st) || st == 1)
        return(NULL)

#    hg <- hg19

    st <- .genometoChrLoc(st, hg)
    bygene <- lapply(seq_len(nrow(st)), function(ii){
        g <- .getGenesFromSeg(chr=st$chrom[ii], Start=st$loc.start[ii], End=st$loc.end[ii], geneDB)
        if(is.null(g))
            return(NULL)

        cbind.data.frame(g,
                        Log2Ratio=st$seg.med[ii],
                        num.mark=st$num.mark[ii], segNum=ii,
                        "segLength(kb)"=round(abs(st$loc.start[ii] - st$loc.end[ii])/1e3, 2)
                        )
    })
    bygene <- do.call(rbind, bygene)
    .addGenomeLoc(bygene, hg)
}
