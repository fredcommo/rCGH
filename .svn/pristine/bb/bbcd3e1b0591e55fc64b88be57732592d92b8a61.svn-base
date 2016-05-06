
shinyUI(
    pageWithSidebar(

        headerPanel("Interactive rCGH Viewer"),

        sidebarPanel(
        #includeCSS("shinySafir.css"),
        tags$head( tags$link(rel="stylesheet", type="text/css", href="style.css") ),

        h4('Gene symbol'),
        div(id="geneSymbol", style="margin-top: -20px;", textInput("geneSymb",
            '', 'NONE')),
        tags$hr(),
                      
        h4('Show chromosome'),
        div(id="showChrom", style="margin-top: -20px;",
            selectInput(inputId = "chr", label = "", choices = c('All', 1:23),
                selected = 'All')),
        tags$hr(),
        
        column(width=12,
               div(class = "form-group",
                   radioButtons("GLcols", "Gain/Loss colors",
                                choices = c("blue/red", "red/blue"), inline = TRUE)
               )
        ),
        
        textInput("minSeg", "Merging segments shorter than (Kb)", 25),
        sliderInput("center", "Recenter profile", min=-1.5, max=1.5, value=0,
            step = .1),
        sliderInput("Ymax", "Rescale max(y)", min=.1, max=1, value=1, step=.1),
        sliderInput("Ymin", "Rescale min(y)", min=.1, max=1, value=1, step=.1),
        sliderInput("gain", "Gain threshold (Log2ratio)", min=0, max=2,
            value=.5, step = .25),
        sliderInput("loss", "Loss threshold (Log2ratio)", min=-2, max=0,
            value=-.5, step = .25),
        withTags(div(class='row-fluid', style="margin-bottom: 75px;",
            align="center",
                div(class="span3", style="float: left; width: 40% !important;",
                    p("Segment length (< Mb)", style="font-size: 14px; 
                        font-weight: bold;")),
                div(class='span3', style="float: right; width: 40% !important; 
                    margin-top: -20px;", textInput("segLen", "", "All") )
                )
        ),

        tags$hr(),

        h4("Download"),
        withTags(
            div(
                class='row-fluid', style="margin-bottom: 65px;", align="center",
                div(class='span6 offset4 text-center', style="float: left; width: 35% !important",
                    downloadButton('downloadPlot', 'Profile') ),
                div(class='span6 offset4 text-center', style="float: left; width: 35% !important",
                    downloadButton('downloadLOH', 'LOH') ),
                div(class='span6 offset4 text-center', style="float: left; width: 25% !important",
                    downloadButton('downloadData', 'Table') )
                )
            ),

        tags$hr(),
        withTags(
            div(
            class='row-fluid', style="margin-top: 20px;", align="left",
            a("@Contact us", style="font-size: 14px;",
            href="mailto:frederic.commo@gustaveroussy.fr?Subject=rCGH%20Viewer",
            target="_top")
            )
            )
        ),

        mainPanel(
            tabsetPanel(
                tabPanel("CGH profile",
                    plotOutput("Profile", width = "100%", height = "100%"),
                    tags$hr(),
                    plotOutput("LOH", width = "100%", height = "100%"),
                    tags$hr(),
                    div(
                        class='row-fluid',
                        style="width: 100%;",
                        align="center", tableOutput("geneSummary")
                        )
                    ),
                tabPanel("Genes table",
                    h4(textOutput("tableTitle1"), align="center"),
                    h4(textOutput("tableTitle2"), align="center"),
                    tags$hr(),
                    dataTableOutput("fullTable")
                    )
                )
            )
        )
)
