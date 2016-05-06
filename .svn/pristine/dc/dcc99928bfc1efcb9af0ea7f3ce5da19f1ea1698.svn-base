setMethod(f="plotDensity",
            signature="rCGH",
            definition=function(object, breaks=NULL, Title=NULL, ...){

                if(!.validrCGHObject(object)) return(NULL)

                pars <- getParam(object)
                cut <- pars$LRcut
                correct <- pars$correctionValue
                ksmooth <- pars$ksmooth
                nG <- pars$nPeak
                mu <- pars$peakMeans
                s2 <- pars$peakSigmaSq
                p <- pars$peakProp
                best <- pars$centralPeak

                cnSet <- getCNset(object)
                LR <- cnSet$Log2Ratio + correct
                Q <- quantile(LR, probs = cut, na.rm=TRUE)
                LR <- LR[LR>Q[1] & LR<Q[2]]
                runLR <- runmed(LR, k = ksmooth)

                if(is.null(Title)){
                    Title <- sprintf("%s\nCorrection value = %s",
                                getInfo(object, 'sampleName'),
                                round(correct, 3))
                }

                x <- seq(min(LR, na.rm=TRUE), max(LR, na.rm=TRUE), len=1000)
                if(is.null(breaks))
                    breaks <- floor(log10(length(runLR))*50)

                h <- hist(runLR, breaks=breaks, plot=FALSE)
                plot(h, freq=FALSE, border="grey50",
                    ylim=range(0, max(h$density)*1.35),
                    xlab=expression(Log[2](Ratio)), main=Title, ...)
                for(ii in seq_len(nG)){
                    col <- rgb(ii/nG, 0.2, (nG-ii)/nG,
                                ifelse(ii==best, .75, .25))
                    mu_i <- mu[ii]; s_i <- sqrt(s2[ii]); p_i <- p[ii]
                    .addDens(x, mu_i, s_i, p_i, best=ii==best, col=col)
                }
        }
)

setMethod(f="plotProfile",
    signature="rCGH",
    definition=function(object, symbol=NULL, gain=.5, loss=(-.5), minLen = 10,
        pCol = "grey50", GLcol = c("blue", "red3"), Title=NULL, ylim=NULL){

        if(!.validrCGHObject(object))
            return(NULL)

        if(!is.null(ylim) && !all(is.numeric(ylim)))
            stop("ylim must be numeric.\n")

        if(!is.null(ylim) && length(ylim)!=2)
            stop("'ylim' must be a vector of 2 values or let to NULL to use
                the full range of values.\n")

        if(is.null(GLcol)){
            GLcol <- c("blue", "red3")
            warning("'GLcol' cannot be NULL. 'blue/red' has been used.")
        }

        if(length(GLcol) < 2)
            warning("You specified less than 2 colors for 'GLcol'.
                Gains only have been colored.")

        if(length(GLcol) > 2){
            warning("You specified more than 2 colors for 'GLcol'.
                Only the 2 first have been used.")
        }

        hg18 <- hg18; hg19 <- hg19; hg38 <- hg38

        HG <- switch(getInfo(object, "genome"),
            hg18 = hg18, hg19 = hg19, hg38 = hg38)

        segTable <- getSegTable(object, minLen)
        if(nrow(segTable)==0){
            message("No data available, yet.")
            return(NULL)
        }
        s = 4
#        myBlue <- rgb(0, 0.45, 1, 1)

        segTable <- segTable[which(segTable$chrom != 24),]
        segTable <- .convertLoc(segTable, HG)

        cumLen <- HG$cumlen[1:23]

        if(is.null(ylim)){
            miny <- max(-2.5, min(segTable$seg.med, na.rm=TRUE)) - .75
            maxy <- max(2, max(segTable$seg.med, na.rm=TRUE)) + .75
            ylim <- range(miny, maxy)
        }

        idx <- which(segTable$seg.med<= loss | segTable$seg.med>= gain)
        subTable <- segTable[idx,]
        GLcolors <- ifelse(subTable$seg.med<= loss, GLcol[2],
                            ifelse(subTable$seg.med>= gain, GLcol[1], NA))

        if(is.null(Title)){
            Title = paste(getInfo(object, 'sampleName'),
                '-', getInfo(object, 'analyseDate'),
                '\nGain threshold: ', round(gain, 3),
                ' Loss threshold:', round(loss, 3))
        }

        platform <- getInfo(object, "platform")
        N <- sum(segTable$num.mark, na.rm=TRUE)
        w <- N/20e3
        X <- lapply(1:nrow(segTable), function(i){
            n <- ceiling(segTable$num.mark[i]/w)
            n <- max(50, n)
            x <- seq(segTable$loc.start[i], segTable$loc.end[i], len=n)
            y <- rnorm(n, segTable$seg.med[i], segTable$probes.Sd[i]/s)
            return(cbind(chr=segTable$chrom[i], loc=x, l2r=y))
            })
        X <- as.data.frame(do.call(rbind, X))

        cumCentr <- 1/2*HG$length[1:23] + cumLen

        # main plot
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

        if(nrow(subTable)>0)
            gPlot <- .addSegments(gPlot, subTable, GLcolors)

        if(!is.null(symbol)){
            bg <- byGeneTable(getSegTable(object, minLen),
                symbol, getInfo(object, "genome"), TRUE)
            return(.addTagToPlot(gPlot, bg))
            }

        return(gPlot)
    }
)

setMethod(f="plotLOH",
    signature="rCGH",
    definition=function(object, Title=NULL){


        if(!.validrCGHObject(object))
            return(NULL)

        hg18 <- hg18; hg19 <- hg19; hg38 <- hg38
        
        HG <- switch(getInfo(object, "genome"),
            hg18 = hg18, hg19 = hg19, hg38 = hg38)

        cnSet <- getCNset(object)

        if(!"modelAllDif" %in% colnames(cnSet)){
            message("No data available for plotting LOH")
            return(NULL)
        }

        # Use snp probes only
        snpSet <- cnSet[grep("^S|^rs", cnSet$ProbeName), ]

        if(is.null(Title)){
            Title = paste(getInfo(object, 'sampleName'),
                '-', getInfo(object, 'analyseDate'))
        }

        ss <- split(snpSet, snpSet$ChrNum)
        gLocs <- lapply(ss, function(tmp){
            n <- nrow(tmp)
            chr <- unique(tmp$ChrNum)
            return(tmp$ChrStart + HG$cumlen[chr])
        })
        gLocs <- do.call(c, gLocs)
        marks <- sapply(2:nrow(HG), function(ii)
            (HG$cumlen[ii-1] + HG$cumlen[ii])/2
        )

        if(inherits(object, "rCGH-Illumina")){
            values <- snpSet$Allele.Difference
        } else{
            values <- snpSet$modelAllDif
            }

        idx <- sample(1:length(values), min(length(values), 50e3))
        X <- data.frame(loc=gLocs[idx], AD=values[idx])
        
        cumCentr <- 1/2*HG$length + HG$cumlen
        gPlot <- ggplot(data = X, aes_string(x="loc", y="AD")) +
            geom_point(pch = 19, cex = 0.2, col = rgb(0,0,0,.75)) +
            geom_hline(yintercept = 0) +
            geom_hline(yintercept = c(-1, 1),
                color = "blue", linetype=2) +
            geom_hline(yintercept = seq(-1.5, 1.5, by = 1),
                color = "lightblue", linetype=2) +
            geom_vline(xintercept = HG$cumlen[1:23], color = 'red',
                linetype = 2, size = 0.25) +
            ggtitle(Title) +
            xlab('Genomic position (bp)') +
            ylab('Allelic Difference') +
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
            )
        if(inherits(object, "rCGH-Illumina")){
                gPlot <- gPlot +
                        coord_cartesian(ylim = range(0, 1.25)) +
                        scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
                        annotate(
                            'text',
                            x = c(-1e8, cumCentr[1:23]), y = rep(1.20, 24),
                            label = c("Chr", seq(1, 23)),
                            size = 4, colour = "grey30")
            }else{
                gPlot <- gPlot +
                        coord_cartesian(ylim = range(-1.99, 1.99)) +
                        scale_y_continuous(breaks = seq(-1.5, 1.5, by = 0.5)) +
                        annotate(
                            'text',
                            x = c(-1e8, cumCentr[1:23]), y = rep(1.75, 24),
                            label = c("Chr", seq(1, 23)),
                            size = 4, colour = "grey30")
            }
        return(gPlot)
    }
)

setMethod(f="multiplot",
    signature="rCGH",
    definition=function(object, symbol=NULL, gain=.5, loss=(-.5), minLen = 10,
        pCol = "grey50", GLcol = c("blue", "red3"),
        L=matrix(seq(1, 12)), p=c(2/3, 1/3), Title=NULL, ylim=NULL){

    if(sum(p)!=1)
        stop("Proportions in 'p' must sum to 1.\n")

    # To initialize the plot window.
#    plot.new()
#    dev.off()

    n <- nrow(L)

    plot1 <- plotProfile(object, symbol, gain, loss, minLen, pCol, GLcol,
        Title, ylim)
    if(is.null(plot1))
        stop("Nothing to plot.\n")

    if(!is.null(Title))
        plot1 <- plot1 + ggtitle(Title)
    
    plot2 <- plotLOH(object)
    
    if(is.null(plot2))
        print(plot1)
    else{
        plot2 <- plot2 + ggtitle("")

        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(L), ncol(L))))
        print(plot1, vp=viewport(layout.pos.row = 1:floor(n*p[1]),
            layout.pos.col = 1))
        print(plot2, vp=viewport(layout.pos.row = floor(n*p[1]+1):n,
            layout.pos.col = 1))
    }
}
)

setMethod(f="view",
    signature="rCGH",
    definition=function(object, browser = TRUE, ...) {

        if(!.validrCGHObject(object))
            return(NULL)

        if (!requireNamespace("shiny", quietly=TRUE)){
            stop("You may install shiny (>= 0.11.1) first")
        }

        hg18 <- hg18; hg19 <- hg19; hg38 <- hg38
        
        genome <- getInfo(object, "genome")
        HG <- switch(genome,
            hg18 = hg18, hg19 = hg19, hg38 = hg38)

        path <- system.file("shinyProfile", package="rCGH")

        lf <- list.files(file.path(path, "data"),
                        full.names = TRUE, pattern = c(".rda|.rds"))
        if(length(lf)>0)
            for(f in lf)
                system(sprintf("rm %s", f))

        if(inherits(object, "rCGH-Agilent")){
            w <- 10
        } else{
            w <- 60
        }

        if(!browser){
            browser <- getOption("shiny.launch.browser", interactive())
        }

        segTable <- .convertLoc(getSegTable(object), HG)
        if(nrow(segTable)==0){
            stop(
                "No segmentation table available to generate a genomic profile. 
                    Please run segmentCGH() first.\n"
            )
        } else{
            segTable <- segTable[which(segTable$chrom %in% 1:23),]
            segTable$num.mark <- round(segTable$num.mark/w)
            saveRDS(segTable, file=file.path(path, "data/st.rds"))
            cat(genome, file = file.path(path, "data/hg.txt"), sep = "\n")
        }

        cnSet <- getCNset(object)
        if("modelAllDif" %in% colnames(cnSet)){
            loh <- plotLOH(object)
            save(loh, file=file.path(path, "data/loh.rda"))
            }

        runApp(path, launch.browser = browser, ...)
        }
)
