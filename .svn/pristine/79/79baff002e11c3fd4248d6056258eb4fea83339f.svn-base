###########################
###########################
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
#    myBlue <- rgb(0, 0.45, 1, 1)

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
.addChr <- function(gPlot, chr, hg19){
    if(chr=="All") chr <- as.numeric(1:23)
    else chr <- as.numeric(chr)

    ylim <- gPlot$coordinates$limits$y
    if(is.null(ylim)){
        ylim <- range(gPlot$data$l2r)
    }
    cumCentr <- 1/2*hg19$length+hg19$cumlen
    gPlot <- gPlot+
        geom_vline(xintercept = hg19$cumlen[chr], color = 'grey30',
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
.addTag <- function(gPlot, geneAnnot, Yexpand, gain, loss, GLcols){

    ylim <- gPlot$coordinates$limits$y
    ymin <- min(ylim)
    ymax <- max(ylim)
    symbol <- as.character(geneAnnot$symbol)
    x <- geneAnnot$genomeStart
    lr <- geneAnnot$Log2Ratio
    yLabel <- ifelse((lr+1.75)<ymax, lr+1.25, lr-1.25)

    if(is.na(lr)) return(gPlot)

    gainCol <- rgb(0, 0.45, 1, 1)
    lossCol <- "red3"
    if(grepl("red/blue", GLcols)){
        gainCol <- "red3"
        lossCol <- rgb(0, 0.45, 1, 1)
    }
    Col <- ifelse(lr<= loss, lossCol, ifelse(lr>=gain, gainCol, 'black'))
    gPlot <- gPlot +
        annotate("text",
            x=max(x, 2e8), y=yLabel,
            label=paste0(symbol, '\n(Log2R = ', round(lr, 2), ')'),
            cex=5, colour="grey25", fontface = "bold") +
        geom_point(x=x, y=lr, cex=5, pch=19, colour="darkorchid4") +
        geom_point(x=x, y=lr, cex=4, pch=19, colour="antiquewhite") +
        geom_point(x=x, y=lr, cex=1, pch=19, colour=Col) #"darkorchid4"
#        geom_point(x=x, y=lr, cex=.25, pch=19, colour="cyan")

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
        delta <- abs(segTable$seg.mean[c(idx-1,idx+1)] - segTable$seg.mean[idx])
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
    return(segTable)
}

# End helper functions
###########################
###########################
