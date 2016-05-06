###########################
###########################
require(ggplot2)
require(grid)

source('helpers.R')
source('byGene.R')
cat("Loading files...")
segTable <- readRDS(file.path(getwd(), "data/st.rds"))
#geneTable <- readRDS(file.path(getwd(), "data/bg.rds"))
err <- try(load(file.path(getwd(), "data/loh.rda")), silent = TRUE)
if(inherits(err, "try-error"))
    loh <- NULL
cat("Done.\n")
###########################

shinyServer(function(input, output, session) {

    options(shiny.deprecation.messages=FALSE)
    
    gene <- reactiveValues(symbol=character())
    observe({gene$symbol <- toupper(input$geneSymb)})

    reCenterSeg <- reactive({
        seg <- segTable
        if(as.numeric(input$minSeg) > 10)
            seg <- .smoothSeg(seg, input$minSeg)
        seg$seg.med <- seg$seg.med + input$center
        return(seg)
        })

    reCenterGenes <- reactive({
        geneTable <- ByGene(reCenterSeg())
        geneTable$Log2Ratio <- geneTable$Log2Ratio + input$center
        return(geneTable)
        })

    createCGHplot <- reactive({
        seg <- reCenterSeg()
        gPlot <- .mainPlot(seg)
        gPlot <- .addSegments(gPlot, seg, input$chr, input$gain, input$loss, input$segLen, input$GLcols)
        gPlot <- .updateScale(gPlot, input$chr, hg19, input$Ymin, input$Ymax)
        gPlot <- .addChr(gPlot, input$chr, hg19)
        gPlot <- .addTitle(gPlot, unique(segTable$ID), input$gain, input$loss)

        if(!gene$symbol %in% c('NONE', '')){
            geneAnnot <- try(.geneOfInt(gene$symbol, reCenterGenes()),
                silent = TRUE)
            if(class(geneAnnot)[1] != 'try-error' & !is.null(geneAnnot)){
                gPlot <- .addTag(gPlot, geneAnnot, input$Yexpand, input$gain,
                    input$loss, input$GLcols)
            }
        }
    print(gPlot)
    })

    createLOH <- reactive({
            .resizeLOH(loh, input$chr)
        })

    createSummary <- reactive({
        if(gene$symbol %in% c('NONE', '')) return(NULL)

        geneAnnot <- try(.geneOfInt(gene$symbol, reCenterGenes()),
            silent = TRUE)
        if(class(geneAnnot)[1] != 'try-error' & !is.null(geneAnnot)){
            selected <- geneAnnot[,c("symbol", "entrezid", "fullName",
                "cytoband", "Log2Ratio", "segNum", "segLength(kb)")]
            selected$Log2Ratio <- round(selected$Log2Ratio, 3)
            for(ii in seq_len(ncol(selected)))
                selected[,ii] <- as.character(selected[,ii])
        } else{
            selected <- data.frame(message=sprintf("\"%s\" does not seem to be 
                an official symbol.", gene$symbol))
        }
        return(selected)
        })

    createFullTable <- reactive({
        gt <- reCenterGenes()
        gt$Log2Ratio <- round(gt$Log2Ratio, 3)

        if(input$chr=="All"){
            chr <- 1:23
        } else{
            chr <- as.numeric(input$chr)
        }

        gt <- gt[gt$chr %in% chr,]
        greater <- as.numeric(input$gain)
        lower <- (-abs(as.numeric(input$loss)))

        if(input$segLen %in% c("All", "")){
            segLen <- Inf
        } else{
            segLen <- as.numeric(input$segLen)
        }

        gt <- gt[(gt$Log2Ratio>=greater | gt$Log2Ratio<=lower) & 
            gt$"segLength(kb)"/1e3 <= segLen,]
        if(nrow(gt)>0){
            gt[, "entrezid"] <- .renderLink(gt[, "entrezid"])
            return(gt[,c("symbol", "fullName", "chr", "cytoband", 
                "entrezid", "Log2Ratio", "segNum", "segLength(kb)")])
        } else return(NULL)
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
    output$fullTable <- renderDataTable({ createFullTable() },
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
