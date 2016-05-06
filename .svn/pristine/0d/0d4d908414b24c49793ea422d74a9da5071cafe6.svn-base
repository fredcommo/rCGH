###########################
###########################
require(ggplot2)
require(grid)

source('helpers.R')
#source('byGene.R')
cat("Loading files...")
genome <- readLines("data/hg.txt")
hg <- switch(genome, hg18 = hg18, hg19 = hg19, hg38 = hg38)

geneDB <- .createGeneDB(genome)
segTable <- readRDS(file.path(getwd(), "data/st.rds"))
err <- try(load(file.path(getwd(), "data/loh.rda")), silent = TRUE)
if(inherits(err, "try-error"))
    loh <- NULL
cat("Done.\n")

print(geneDB)

###########################

shinyServer(function(input, output, session) {

    options(shiny.deprecation.messages=FALSE)

    Input <- reactiveValues(
    	segTable = data.frame(),
    	geneTable = data.frame()
    	)
    observe({
    	Input$segTable <- .smoothSeg(segTable, input$minSeg)
    	Input$geneTable <- ByGene(Input$segTable, hg, geneDB)
    	})

    reCenterSeg <- reactive({
        st <- Input$segTable
        if(!is.null(st)){
            st$seg.med <- st$seg.med + input$center
            return(st)
            }
        return(NULL)
        })

    reCenterGenes <- reactive({
        if(is.null(Input$geneTable))
            return(NULL)

        geneTable <- Input$geneTable
        geneTable$Log2Ratio <- geneTable$Log2Ratio + input$center
		lrr <- geneTable$Log2Ratio
		fc <- ifelse(lrr>=0, round(2^lrr, 1), round(-1/2^lrr, 1))
	    geneTable$"ApproxCN" <- 2*fc
        return(geneTable)
        })

    CGHplot <- reactive({
        seg <- reCenterSeg()
        gPlot <- .mainPlot(seg)
        gPlot <- .addSegments(gPlot, seg, input$chr, input$gain, input$loss, input$segLen, input$GLcols)
        gPlot <- .updateScale(gPlot, input$chr, hg19, input$Ymin, input$Ymax)
        gPlot <- .addChr(gPlot, input$chr, hg)
        gPlot <- .addTitle(gPlot, unique(segTable$ID), input$gain, input$loss)
        
        return(gPlot)
    })

    createCGHplot <- reactive({
        gPlot <- CGHplot()
        
        geneTable <- reCenterGenes()
        gene <- toupper(input$geneSymbol)
        if(!gene %in% c('NONE', '')){
            if(gene %in% geneTable$symbol){
            	geneAnnot <- geneTable[which(geneTable$symbol == gene),]
                gPlot <- .addTag(gPlot, geneAnnot)
            }
        }
        print(gPlot)
    })

    createLOH <- reactive({
            .resizeLOH(loh, input$chr)
        })

    createSummary <- reactive({

        gene <- toupper(input$geneSymbol)
        if(gene %in% c('NONE', ''))
            return(NULL)

		geneTable <- reCenterGenes()
		if(!gene %in% geneTable$symbol){
            msg <- sprintf("'%s' may not be a valid HUGO symbol", gene)
            return( data.frame(Error = msg) )			
		}
		else{
			geneAnnot <- geneTable[which(geneTable$symbol == gene),]
	        chr <- as.integer(geneAnnot$chr)
	        chrStart <- as.integer(geneAnnot$chrStart)
	        chrEnd <- as.integer(geneAnnot$chrEnd)
	        geneAnnot$position <- sprintf("chr%s:%s-%s", chr, chrStart, chrEnd)
	        geneAnnot$segNum <- as.integer(geneAnnot$segNum)
	        geneAnnot$"segLength(kb)" <- as.integer(geneAnnot$"segLength(kb)")
	        geneAnnot <- geneAnnot[,c("symbol", "entrezid", "fullName",
	            "position", "segNum", "segLength(kb)", "Log2Ratio", "ApproxCN")]

	        return(geneAnnot)
			}

        })

    createFullTable <- reactive({

        if(is.null(Input$geneTable))
            return(NULL)

        geneTable <- reCenterGenes()
        geneTable <- geneTable[,c("symbol", "entrezid", "fullName",
            "chr", "cytoband", "chrStart", "chrEnd",
            "segNum", "segLength(kb)", "Log2Ratio", "ApproxCN")]
        geneTable$Log2Ratio <- round(geneTable$Log2Ratio, 2)
        geneTable$entrezid <- .renderLink(geneTable$entrezid)
        geneTable
        })

    filterGeneTable <- reactive({
        geneTable <- createFullTable()

        if(input$chr=="All"){
            chr <- 1:23
        } else{
            chr <- as.numeric(input$chr)
        }

        if(input$segLen %in% c("All", "")){
            segLen <- Inf
        } else{
            segLen <- as.numeric(input$segLen)
        }

        greater <- as.numeric(input$gain)
        lower <- (-abs(as.numeric(input$loss)))

        geneTable <- .filterBygene(geneTable, chr, greater, lower, segLen)
        return(geneTable)
        })

    createTitle1 <- reactive({
        return( unique(segTable$ID) )
        })

    createTitle2 <- reactive({
        hi <- input$gain
        lo <- input$loss
        return( sprintf("Gain threshold: %s, Loss threshold: %s", hi, lo) )
        })

    output$Profile <- renderPlot({ createCGHplot() },
        res=120, width=1500, height=520)
    output$LOH <- renderPlot({ createLOH() }, res=120, width=1500, height=320)
    output$geneSummary <- renderTable({ createSummary() })
    output$tableTitle1 <- renderText({ createTitle1() })
    output$tableTitle2 <- renderText({ createTitle2() })
    output$fullTable <- renderDataTable({ filterGeneTable() },
        options=list(lengthMenu=c(25, 50, 100)), escape = FALSE)

    # Download functions
    output$downloadPlot <- downloadHandler(
        filename <- function(){
            sprintf("%s_genomic_profile.png", createTitle1())
            },
        content <- function(file) {
            png(file, width=1200, height=500)
            print( createCGHplot() )
            dev.off()
            }, contentType = "image/png"
            )

    output$downloadLOH <- downloadHandler(
        filename <- function(){
            sprintf("%s_LOH_profile.png", createTitle1())
            },
        content <- function(file) {
            png(file, width=1200, height=350)
            print( createLOH() )
            dev.off()
            }, contentType = "image/png"
            )

    output$downloadData <- downloadHandler(
        filename <- function(){
            sprintf("%s_geneTable_chr%s_gain%s_loss%s.xls",
            createTitle1(), input$chr, input$gain, input$loss)
            },
        content <- function(file){
            out <- createFullTable()
            out$entrezid <- gsub(".*term=|\\[uid\\].*", "",
                out$entrezid)
            write.table(out, file, sep="\t", row.names=FALSE)
        }
        )
    session$onSessionEnded(function() { stopApp() })
    })
