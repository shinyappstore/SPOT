options(repos = BiocManager::repositories())
options(shiny.maxRequestSize = 100*1024^2)
library(rsconnect)
library(DT)
library(shiny)
library(shinyWidgets)
library(matrixStats)
library(plotly)
library(openxlsx)
library(rlist)
library(shinyEffects)
library(shinycssloaders)
library(scales)
library(tidyr)
library(plyr)
library(EnvStats)
library(stringr)
library(sortable)
library(shinybusy)
library(reshape2)
library(edgeR)
library(shinyLP)


params <- list(allsc_genes        = "Files/sc_P_berghei_averaged.csv",
               allk_genes         = "Files/bulk_H_sapiens_averaged.csv",
               allsa_genes        = "Files/bulk_Calu3_A.csv",
               dotplot_genes      = "Files/sc_P_berghei_dotplot.csv",
               Counts             = "Files/sc_P_berghei_counts.csv"
)

source("helper_module.R")

sc_genes <- read.csv2(params$allsc_genes, stringsAsFactors = FALSE)

human_genes <- read.csv2(params$allk_genes, stringsAsFactors = FALSE)

cov_genes <- read.csv2(params$allsa_genes, stringsAsFactors = FALSE)

sc_dot_plot <- read.csv2(params$dotplot_genes, stringsAsFactors = FALSE)

Sc_counts <- as.matrix(read.csv2(params$Counts, sep = ";", stringsAsFactors = FALSE))

################################################################################
##                                   UI                                       ##
################################################################################

ui <- navbarPage(
  span(HTML(paste0(tags$sup(tags$i("X"), style = "font-size : 18px; color:#383D3B ; font-weight:bold"),"SPOT")), style = " color: #c00000; font-size : 32px; font-weight: bold"),
  id = "navbar", selected = "value", windowTitle = "SPOT - Swift Profiling Of Transcriptomics",

  tags$style(HTML("  *{font-family: Verdana;}
                      .navbar { background-color: #f1f1f1}
                      .navbar-default .navbar-nav > li > a {color:#c00000;}
                      .navbar-default .navbar-nav > .active > a,
                      .navbar-default .navbar-nav > .active > a:focus,
                      .navbar-default .navbar-nav > .active > a:hover {color: white;background-color: #c00000;}
                      .navbar-default .navbar-nav > li > a:hover {color: white;background-color:#c00000;text-decoration:underline;text-color:white}
                      .button_DEA .bttn-danger{background-color: #c00000; border-color: black; margin-bottom: 10px }
                      .button_dd .bttn-danger{color: white; border-color: grey; background-color: #c00000; margin-bottom: 10px}
                      .text_about {color:black; font-size: 14px;text-align: justify; margin-top: 10px; margin-bottom: 20px; margin-left: 70px; margin-right: 70px}
                      .text_cite {color:black; font-size: 14px; font-weight: bold ;text-align: justify; margin-left: 70px; margin-right: 70px}
                      .text_ref {color:black; font-size: 14px; font-weight: bold ;text-align: justify; margin-top: 35px; margin-left: 70px; margin-right: 70px}
                      .text_dwld {color:black; font-size: 14px; font-weight: bold;text-align: justify; margin-top: 5px; margin-bottom: 5px}
                      .text_it {color:black; font-size: 14px; font-style: italic;text-align: justify; margin-top: 10px; margin-bottom: 20px; margin-left: 70px; margin-right: 70px}
                      .pic {margin-top: 35px; margin-left: 0px; margin-right: 70px}
                  ")),
  tabPanel(span("SPOT expression profiles", style = "font-size : 20px; font-weight: bold"), value = "value",
           
           fluidPage(
             setShadow(id = "Component2a"),
             setShadow(id = "Component2b"),
             setShadow(id = "Component2c"),
             
             conditionalPanel(
               condition = "input.dataSwitch2 == 'Plasmodium (single cell)'",
               withSpinner(plotlyOutput("Component2a"), type = "2", color.background = "white", color = "#c00000"),
             ),
             
             conditionalPanel(
               condition = "input.dataSwitch2 == 'Human organs (bulk)'",
               withSpinner(plotlyOutput("Component2b"), type = "2", color.background = "white", color = "#c00000"),
               
             ),
             
             conditionalPanel(
               condition = "input.dataSwitch2 == 'COVID cell line (bulk)'",
               withSpinner(plotlyOutput("Component2c"), type = "2", color.background = "white", color = "#c00000"),
               
             ),
             chooseSliderSkin("Flat"),
             
             setSliderColor(rep("#c00000", 100), c(1:100)),
             
             hr(),
             
             fluidRow(
               
               column(2, 
                      
                      radioGroupButtons(
                        inputId = "Algos",
                        label = "Choose algorithm",
                        choices = c("SPOT", "Correlation"),
                        individual = TRUE,
                        checkIcon = list(
                          yes = tags$i(class = "fa fa-circle", 
                                       style = "color: #c00000"),
                          no = tags$i(class = "fa fa-circle-o", 
                                      style = "color: #c00000"))
                      ),
                      
                      conditionalPanel(
                        condition = "input.dataSwitch2 == 'Human organs (bulk)' ",
                        tags$h5("Subset samples", style = "font-weight: bold "),
                        #br(),
                        div(class = "button_dd",
                            dropdown(
                              
                              bucket_list(
                                header = "Drag items from available datasets and drop one or several for inclusion in analysis.\n wpc -> week post conception",
                                add_rank_list(input_id = "bucket_in",
                                              text = "Available datasets",
                                              labels = c("4wpc", "10wpc", "20wpc", "infant", "toddler", "school", 
                                                         "teenager", "youngAdult", "youngMidAge", "olderMidAge", "senior")
                                              ),
                                add_rank_list( input_id = "bucket_out",
                                               text = "Selected for analysis",
                                               labels = c("newborn")
                                               ),
                                orientation = "horizontal"
                              ),
                              
                              style = "bordered", icon = icon("bars"), label = "Subset",
                              status = "danger", width = "1000px", 
                              height = "200px",
                              animate = animateOptions(
                                enter = animations$fading_entrances$fadeInLeftBig,
                                exit = animations$fading_exits$fadeOutRightBig
                              )
                            )),
                      ),
                      
                      pickerInput(
                        inputId = "dataSwitch2",
                        label = "Select Dataset", 
                        choices = c("Plasmodium (single cell)", "Human organs (bulk)", "COVID cell line (bulk)")
                      ),
                    
               ),
               
               conditionalPanel(
                 condition = "input.dataSwitch2 == 'Plasmodium (single cell)' ",
                 slidersUI("1", colnames(sc_genes[4:13])),
               ),
               
               conditionalPanel(
                 condition = "input.dataSwitch2 == 'Human organs (bulk)'",
                 uiOutput("sliders_K"),
               ),
               
               conditionalPanel(
                 condition = "input.dataSwitch2 == 'COVID cell line (bulk)'",
                 slidersUI("2", colnames(cov_genes[3:13])),
               ),
               
               column(2,
                      
                      conditionalPanel(
                        condition = "input.dataSwitch2 == 'Plasmodium (single cell)' ",
                        selectInput("spec_gene2", "Pick genes to visualize:", choices = sc_genes[,1], NULL, multiple = TRUE), #sc_genes$PB_ID
                      ),
                      
                      radioGroupButtons(
                        inputId = "Radio2",
                        label = "Alter Visualization",
                        choices = c("Table", "Bar chart", "Dot plot"),
                        individual = TRUE,
                        checkIcon = list(
                          yes = tags$i(class = "fa fa-circle", 
                                       style = "color: #c00000"),
                          no = tags$i(class = "fa fa-circle-o", 
                                      style = "color: #c00000"))
                      ),
                      div(class = "text_dwld",
                          "Download:"
                      ),
                      downloadButton('downloadData1', 'as .xlsx'),
                      downloadButton('downloadData4', 'as .csv')
                      
               )))),
  
  tabPanel(span("Differential Expression Analysis", style = "font-size : 20px; font-weight: bold"),
           fluidPage(
             
             plotlyOutput("Component3"),
             
             hr(),
             
             fluidRow(
               
               column(2, 
                      
                      radioGroupButtons(
                        inputId = "algorithm",
                        label = "Choose algorithm",
                        choices = c("Wilcox", "MAST", "DESeq2"),
                        individual = TRUE,
                        checkIcon = list(
                          yes = tags$i(class = "fa fa-circle", 
                                       style = "color: #c00000"),
                          no = tags$i(class = "fa fa-circle-o", 
                                      style = "color: #c00000"))
                      ),
                      
                      conditionalPanel(
                        condition = "input.algorithm == 'DESeq2'",
                        helpText("Attention: this calculation will take", 
                                 "several minutes. For more specific questions",
                                 "we suggest working directly with DESeq2."),
                        
                      ),
                      div(class = "text_dwld",
                          "Download:"
                      ),
                      downloadButton('downloadData2', 'as .xlsx'),
                      downloadButton('downloadData5', 'as .csv')
                      
               ),
               slidersUI("5", colnames(sc_genes[4:13])),
               column(2,
                      
                      div(class = "button_DEA",
                          actionBttn("action_DEA", label = "Start DEA!", icon = NULL, style = "material-flat", color = "danger",
                                     size = "lg", block = FALSE, no_outline = TRUE
                          )),
                      
                      radioGroupButtons(
                        inputId = "Radio3",
                        label = "Alter Visualization",
                        choices = c("Table", 
                                    "Bar chart", "Dot plot"),
                        individual = TRUE,
                        checkIcon = list(
                          yes = tags$i(class = "fa fa-circle", 
                                       style = "color: #c00000"),
                          no = tags$i(class = "fa fa-circle-o", 
                                      style = "color: #c00000"))
                      ),
                      
               ),
             ))),
  
  tabPanel(span("Compare profiles", style = " font-size : 20px; font-weight: bold "), 
           
           fluidPage(
             setShadow(id = "Component1"),
             withSpinner(plotlyOutput("Component1"), type = "2", color.background = "white", color = "#c00000"),
             fluidRow(
               column(10,
                      DTOutput("tab_gen1")
                      ),
               
               column(2,
                      div(class = "text_dwld",
                          pickerInput(
                            inputId = "dataSwitch4",
                            label = "Select Dataset", 
                            choices = c("Plasmodium (single cell)", 
                                        "Human organs (bulk)", 
                                        "COVID cell line (bulk)")
                            )
                          ),
                      ),
               )
             )
           ),
  
  tabPanel(span("Upload your own data", style = "font-size : 20px; font-weight: bold "),
           fluidPage( 
             
             conditionalPanel(
               condition = "input.Radio6 == 'Normalized counts'",
               withSpinner(plotlyOutput("Component4a"), type = 2, color = "#c00000", color.background = "white"),
             ),
             
             conditionalPanel(
               condition = "input.Radio6 == 'Raw counts' & input.method == 'SPOT'",
               plotlyOutput("Component4b"),
             ),
             
             conditionalPanel(
               condition = "input.Radio6 == 'Raw counts' & input.method == 'DEA'",
               plotlyOutput("Component4ba"),
             ),
             
             hr(),
             
             fluidRow(
               column(2,
                      radioGroupButtons(
                        inputId = "Radio6",
                        label = "Choose count format",
                        choices = c("Normalized counts", 
                                    "Raw counts"),
                        individual = TRUE,
                        checkIcon = list(
                          yes = tags$i(class = "fa fa-circle", 
                                       style = "color: #c00000"),
                          no = tags$i(class = "fa fa-circle-o", 
                                      style = "color: #c00000"))
                      ),
                      
                      fileInput('file1', 'Choose your dataset',
                                accept = c('text/csv','text/comma-separated-values',
                                           'text/tab-separated-values', 'text/plain',
                                           '.csv','.tsv', 'text/xlsx','.xlsx'),
                      ),
                      
                      conditionalPanel(
                        condition = "input.Radio6 == 'Raw counts'",
                      radioGroupButtons(
                        inputId = "method",
                        label = "Choose method",
                        choices = c("SPOT", "DEA"),
                        individual = TRUE,
                        checkIcon = list(
                          yes = tags$i(class = "fa fa-circle", 
                                       style = "color: #c00000"),
                          no = tags$i(class = "fa fa-circle-o", 
                                      style = "color: #c00000"))
                      )),
                      
                      conditionalPanel(
                        condition = "input.method == 'DEA'",
                        radioGroupButtons(
                          inputId = "algorithm_upload",
                          label = "Choose algorithm",
                          choices = c("Wilcox", "MAST", "DESeq2"),
                          individual = TRUE,
                          checkIcon = list(
                            yes = tags$i(class = "fa fa-circle", 
                                         style = "color: #c00000"),
                            no = tags$i(class = "fa fa-circle-o", 
                                        style = "color: #c00000"))
                        ))
                  
                     ),
               
               uiOutput("reactive_sliders4a"),
               
               column(2,
                      
                      conditionalPanel(
                        condition = "input.method == 'DEA'",
                        div(class = "button_DEA",
                            actionBttn("action_DEA_upload", label = "Start DEA!", icon = NULL, style = "material-flat", color = "danger",
                                       size = "lg", block = FALSE, no_outline = TRUE
                            )),
                        ),
                      
                      radioGroupButtons(
                        inputId = "Radio5",
                        label = "Alter Visualization",
                        choices = c("Table", 
                                    "Bar chart",
                                    "Dot plot"),
                        individual = TRUE,
                        checkIcon = list(
                          yes = tags$i(class = "fa fa-circle", 
                                       style = "color: #c00000"),
                          no = tags$i(class = "fa fa-circle-o", 
                                      style = "color: #c00000"))
                      ),
                      div(class = "text_dwld",
                          "Download:"
                      ),
                      downloadButton('downloadData3', 'as .xlsx'),
                      downloadButton('downloadData6', 'as .csv')
               )
             ))),
  
  tabPanel(span("About", style = "font-size : 20px; font-weight: bold "),
           setShadow(id = "UMAP"),
           fluidRow(
             column(5,
                    div(class = "text_ref",
                        "SPOT: enabling Swift Profiling Of Transcriptomes"
                    ),
                    div(class = "text_about",
                        "The increasing number of single cell and bulk RNAseq data sets describing complex 
                        gene expression profiles in different organisms, organs or cell types calls for an 
                        intuitive tool allowing rapid comparative analysis. Here we present Swift Profiling 
                        Of Transcriptomes (SPOT) as a web tool that allows not only differential expression
                        analysis but also fast ranking of genes fitting transcription profiles of interest.
                        Based on a heuristic approach the spot algorithm ranks the genes according to their
                        proximity to the user-defined gene expression profile of interest. The best hits are 
                        visualized as a table, bar chart or dot plot and can be exported as an Excel file.
                        While the tool is generally applicable, we tested it on RNAseq data from malaria parasites
                        that undergo multiple stage transformations during their complex life cycle as well as
                        on data from multiple human organs during development and cell lines infected by the
                        SARS-CoV-2 virus. SPOT should enable non-bioinformaticians to easily analyse their own 
                        and any available dataset. "
                    ),
                    div(class = "text_it",
                        "Elias Farr, Julia M Sattler, Friedrich Frischknecht, SPOT: a web-tool enabling Swift Profiling Of Transcriptomes; Biorxiv 2021",
                        tags$a(href="https://www.biorxiv.org/content/10.1101/2021.03.03.433767v1", "Link")
                    )
             ),
             column(5, 
                     div(class= "pic",
                         #img(src='Figure1.PNG', align = "bottom", height = 420)
                         HTML('<iframe width="800" height="430" src="//www.youtube.com/embed/3c38RN-ef5Y" frameborder="20" allowfullscreen></iframe>')
                         )
                    )
             ),
           fluidRow(
             
             column(10, 
                    
                    div(class = "text_ref",
                        "Further information and code: "
                    ),
                    div(class = "text_about",
                        tags$a(href="https://github.com/EliasFarr/SPOT.git", "Github"),
                        tags$a(href="https://www.biorxiv.org/content/10.1101/2021.03.03.433767v1", "Publication"),
                        tags$a(href="https://www.biorxiv.org/content/biorxiv/early/2021/03/04/2021.03.03.433767/DC1/embed/media-1.zip?download=true", "Documentation")
                      
                    ), 
                    div(class = "text_cite",
                        "Data was obtained from: "
                    ),
                    div(class = "text_about",
                        "Howick, V.M. et al., 2019. The Malaria Cell Atlas: Single parasite transcriptomes across the complete Plasmodium life cycle. Science, 365(6455), S.774.",
                        tags$a(href="https://science.sciencemag.org/content/365/6455/eaaw2619", "Paper"),
                        tags$a(href="https://github.com/vhowick/MalariaCellAtlas/tree/master/Expression_Matrices/Smartseq2", "Data")
                    ), 
                    div(class = "text_about",
                        "Cardoso-Moreira, M. et al., 2019. Gene expression across mammalian organ development. Nature, 571(7766), S.505-509.",
                        tags$a(href="https://www.nature.com/articles/s41586-019-1338-5", "Paper"),
                        tags$a(href="https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6814/", "Data")
                    ),
                    div(class = "text_about",
                         "Wyler E., et al., 2021. Transcriptomic profiling of SARS-CoV-2 infected human cell lines identifies HSP90 as target for COVID-19 therapy. iScience, 24(3), 102151",
                         tags$a(href="https://www.sciencedirect.com/science/article/pii/S258900422100119X?via%3Dihub#undfig1", "Paper"),
                         tags$a(href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148729", "Data")
                    )
                    
                    
             )
             
           )
  )
  )