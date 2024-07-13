
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

server <- function(input, output) {
################################################################################ 
##                              Component 1: SPOT                             ##
################################################################################
  # color vector for bar charts
  colors <- c("#618C84","#726E75","#948B89","#D0D1AC")
  
  # create table for section 3 "compare"
  table1 <- sc_genes[,c(1:13)]
  # table2 <- human_genes[,c(1:10)]
  # table3 <- cov_genes[,c(1:10)]
  # 
  # 
  # delete names for spot representation
  sc_genes = sc_genes[,-3]
  
  
  observe({
    # get slider values 
    Variables = callModule(sliders_mod, "1", colnames(sc_genes[3:12]))
    # define length of output table 
    output_length = 50
    
    if(input$Algos == "SPOT"){
      #save results from spot calculation
      top_df <<- spot(sc_genes, Variables, columns = c(3:12), Candidate_Number = output_length)
      colnames(top_df)[ncol(top_df)] = "Gene ID"
      
    }else if(input$Algos == "Correlation"){
      
      top_df <<- correlation(sc_genes, Variables, columns = c(3:12), Candidate_Number = output_length)
      colnames(top_df)[ncol(top_df)] = "Gene ID"
    }
    top_df_t <- as.data.frame(t(top_df))

    output$Component2a <- renderPlotly({ 
      
      if(input$Radio2 == "Table"){
        
        columnwidth = c(85, 250, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70)
        generate_table(top_df, c(ncol(top_df), 2:13), columnwidth, Table_rows = 50)
        
      } else if(input$Radio2 == "Bar chart"){
        
        if(is.null(input$spec_gene2) == FALSE){
          
          Gene_indices = get_indices(data = top_df, TextInput = input$spec_gene2)
          subplot_profiles(dat = top_df, top_df_t, Gene_indices, col_vec = c("#618C84","#726E75","#948B89","#D0D1AC"), range_ = c(0,12))
          
        }else{ 
          subplot_profiles(top_df, top_df_t, 1:4, range_ = c(0,12), profile_columns = c(4:13))
          
        }
      }else if(input$Radio2 == "Dot plot"){
        
        dotplot_high_pct1 = subset(sc_dot_plot, sc_dot_plot$Gene %in% top_df[1:10,1])
        dotplot_high_pct1[,3] = dotplot_high_pct1[,3]*20
        
        sorter_df = data.frame(0:9)
        rownames(sorter_df) = top_df[1:10,1]
        Rank = c(1:length(dotplot_high_pct1$Gene))
        j = 1
        for(i in dotplot_high_pct1$Gene){
          Rank[j] = sorter_df[i,1]
          j = j + 1
        }
        dotplot_high_pct1$Gene = paste(Rank, dotplot_high_pct1$Gene)
        sc_dotplot(dotplot_high_pct1)
      }
    })
    
  })
  
  
  
  observe({
    output$sliders_K <- renderUI({
      
      gene_subset <<- human_genes[, c(1,2,(which(unlist(str_split(colnames(human_genes)[3:(ncol(human_genes)-1)], pattern = "_"))[seq(2, 270, 2)] %in% input$bucket_out) + 2))]
      slidersUI("2", colnames(gene_subset)[3:ncol(gene_subset)])
    })
    
    
  
    output$Component2b <- renderPlotly({
      
      gene_subset <<- human_genes[, c(ncol(human_genes),2,(which(unlist(str_split(colnames(human_genes)[3:(ncol(human_genes)-1)], pattern = "_"))[seq(2, 270, 2)] %in% input$bucket_out) + 2))]
      
      if(ncol(gene_subset) > 2){
          Variables = callModule(sliders_mod, "2", colnames(gene_subset)[3:ncol(gene_subset)])
          
          if(input$Algos == "SPOT" & dim(gene_subset)[2] > 2){
            
            top_df <<- spot(gene_subset, Variables, columns = c(3:ncol(gene_subset)))
            
          }else if(input$Algos == "Correlation"){
            
            top_df <<- correlation(gene_subset, Variables, columns = c(3:ncol(gene_subset)))
          }
          
          top_df[,3:ncol(top_df)] = round(top_df[,3:ncol(top_df)], 2)#round values
          
          top_df_t = as.data.frame(t(top_df))
          
          if(input$Radio2 == "Table"){
            
            columnwidth = c( 85, 70)
            
            generate_table(top_df, columns = c(1:(ncol(top_df)-1)), columnwidth, header_size = c(12,10))
            
          } else if(input$Radio2 == "Bar chart"){
            
            if(is.null(input$spec_gene2) == FALSE){
              #change ids und create ranks to get the column index of the genes chosen
              Gene_indices = get_indices(data = top_df, TextInput = input$spec_gene2)
              subplot_profiles(dat = top_df, top_df_t, Gene_indices, col_vec = c("#618C84","#726E75","#948B89","#D0D1AC"))
              
            }else{ 
              
              subplot_profiles(top_df, top_df_t, 1:4, profile_columns = c(4:(ncol(top_df)-1)), range_ = c(0, max(top_df[4:ncol(top_df)])))
              
            }
          }else{
            bulk_dotplot(top_df[1:10,1:(ncol(top_df)-1)], ytitle = "Ensembl ID", xtitle = F)
            
          }
       }else{
         "auto"
       }
    })
  
  })
  
  observe({
    
    # get slider values 
    Variables = callModule(sliders_mod, "2", colnames(cov_genes[3:13]))
    
    # define length of output table 
    output_length = 50

    if(input$Algos == "SPOT"){
      #save results from spot calculation
      top_df <<- spot(cov_genes, Variables, columns = c(3:13), Candidate_Number = output_length)
      colnames(top_df)[ncol(top_df)] = "Gene ID"
      
    }else if(input$Algos == "Correlation"){
      
      top_df <<- correlation(cov_genes, Variables, columns = c(3:13), Candidate_Number = output_length)
      colnames(top_df)[ncol(top_df)] = "Gene ID"
    }
    top_df_t <- as.data.frame(t(top_df))

    output$Component2c <- renderPlotly({ 
      
      if(input$Radio2 == "Table"){
        
        columnwidth = c(95, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70)
        generate_table(top_df, c(ncol(top_df), 2:14), columnwidth, Table_rows = 50)
      
      } else if(input$Radio2 == "Bar chart"){
        
        subplot_profiles(top_df, top_df_t, 1:4, range_ = c(0, 1000), profile_columns = c(4:12))
          
      }else if(input$Radio2 == "Dot plot"){
        
        bulk_dotplot(top_df[1:10,1:(ncol(top_df)-1)], ytitle = "Ensembl ID", xtitle = F)
      }
    })
    
  })
################################################################################
##                              Component 2: DEA                              ##
###############################################################################
  
  observeEvent(input$action_DEA,
               
               { 
                 show_modal_progress_line() # show the modal window
                 library(Seurat)
                # import::from(Seurat, CreateSeuratObject, NormalizeData, FindMarkers)
                 Variables <- callModule(sliders_mod, "5", colnames(sc_genes[3:12]))
                 States = c("Ookinete", "Oocyst", "Sporozoite", "Liver", "Merozoite", "Ring", "Trophozoite", "Schizont",
                            "Female", "Male")
                 
                 rownames(Sc_counts) = Sc_counts[,1]
                 Sc_counts = Sc_counts[,-1]
                 Indents = Sc_counts[1,]
                 Sc_counts = Sc_counts[-1,]
                 
                 SS3_plasmo = CreateSeuratObject(counts = Sc_counts, project = "SPOT", min.cells = 3, min.features = 200)
                 Idents(SS3_plasmo) = as.character(Indents)
                 SS3_plasmo <- NormalizeData(SS3_plasmo)
                 update_modal_progress(0.2) # update progress bar value
                 if(input$algorithm == "Wilcox"){
                   library(limma)
                   #import::from(limma, rankSumTestWithCorrelation)
                   Markers <- FindMarkers(SS3_plasmo, ident.1 = States[which(Variables > 0)], ident.2 = States[which(Variables == 0)], test.use = "wilcox")
                   #detach("imports")
                   #detach(package:limma, unload = T)
                 }else if(input$algorithm == "MAST"){
                   library(MAST)
                   #import::from(MAST, FromMatrix, zlm)
                   Markers <- FindMarkers(SS3_plasmo, ident.1 = States[which(Variables > 0)], ident.2 = States[which(Variables == 0)], test.use = "MAST")
                   #detach("imports")
                   detach(package:MAST, unload = T)
                 }else if(input$algorithm == "DESeq2"){
                   library(DESeq2)
                   #import::from(DESeq2, DESeqDataSetFromMatrix, estimateSizeFactors, estimateDispersions, nbinomWaldTest, results)
                   Markers <- FindMarkers(SS3_plasmo, ident.1 = States[which(Variables > 0)], ident.2 = States[which(Variables == 0)], test.use = "DESeq2")
                   #detach("imports")
                   detach(package:DESeq2, unload = T)
                 }
                 update_modal_progress(0.9)
                 remove_modal_progress()
                 Markers_up = subset(Markers, Markers[,"avg_log2FC"] > 0)
                 
                 Markers_up$PB_ID = str_replace(rownames(Markers_up), "-", "_")
                 sc_dea = subset(sc_genes, sc_genes$PB_ID %in% Markers_up$PB_ID)
                 rownames(sc_dea) = sc_dea$PB_ID
                 merge_df = merge(Markers_up, sc_dea, sort = F)
                 
                 top_df <<- reactive({as.data.frame(merge_df %>% dplyr::arrange(desc("p_val_adj")))})
                 
                 top_df_t = reactive({as.data.frame(t(top_df()))})
                 
                 output$Component3 <- renderPlotly({ 
                  
                   if(input$Radio3 == "Table"){
                     
                     generate_table(top_df(), columns = c(ncol(top_df()), 7, 3, 6, 8:17), header_size = c(12), columnwidth = c(90, 200, 80))
                     
                   }else if(input$Radio3 == "Bar chart"){
                     
                     subplot_profiles(top_df(), top_df_t(), 1:4, range_ = c(0,20) , profile_columns = c(8:17))
                     
                   }else if(input$Radio3 == "Dot plot"){
                     
                     dotplot_high_pct1 = subset(sc_dot_plot, sc_dot_plot$Gene %in% top_df()[1:10,1])
                     dotplot_high_pct1[,3] = dotplot_high_pct1[,3]*20
                     sorter_df = data.frame(0:9)
                     rownames(sorter_df) = top_df()[1:10,1]
                     Rank = c(1:length(dotplot_high_pct1$Gene))
                     j = 1
                     for(i in dotplot_high_pct1$Gene){
                       Rank[j] = sorter_df[i,1]
                       j = j + 1
                     }
                     dotplot_high_pct1$Gene = paste(Rank, dotplot_high_pct1$Gene)
                     sc_dotplot(dotplot_high_pct1)
                   }
                 })
               }
  )
  
  
################################################################################
##                              Component 3: Compare                          ##
################################################################################
  observe({
    
    if(input$dataSwitch4 == "Plasmodium (single cell)"){
      
    }else if(input$dataSwitch4 == "Human organs (bulk)"){
      
      table1 = human_genes[,c(1:(ncol(human_genes)-1))]
      
    }else if(input$dataSwitch4 == "COVID cell line (bulk)"){
      
      table1 = cov_genes[,c(1:(ncol(cov_genes)-1))]
      
    }
    
    output$tab_gen1 <- DT::renderDT(datatable(table1, style = "bootstrap4", filter = 'top',
                                             options = list(pageLength = 5,
                                             initComplete = DT::JS(
                                               "function(settings, json) {",
                                               "$(this.api().table().header()).css({'background-color': '#c00000', 'color': '#fff'});",
                                               "}"), scrollX = T)
                                             ) %>% DT::formatStyle(1:13, color = "DimGrey"),
                                   server = T)
  
    output$Component1 <- renderPlotly({
      
      table1_t = as.data.frame(t(table1))
      
      if(is.null(input$tab_gen1_rows_selected)){
        
        if(input$dataSwitch4 == "Plasmodium (single cell)"){
          
          plot_profiles(table1, table1_t, ID = 240, c(4:ncol(table1)), color = colors[1], range_w = c(0, max(table1[4:ncol(table1)])))
          
        }else {
          
          plot_profiles(table1, table1_t, ID = 513, c(3:ncol(table1)), color = colors[1], range_w = c(0, max(table1[3:ncol(table1)])))
          
        }
        
      }else if(length(input$tab_gen1_rows_selected) > 0){
        
        if(input$dataSwitch4 == "Plasmodium (single cell)"){
          
          subplot_profiles(table1, table1_t, input$tab_gen1_rows_selected, profile_columns = c(4:ncol(table1)), range_ = c(0, max(table1[4:ncol(table1)])))
          
        }else {
          
          subplot_profiles(table1, table1_t, input$tab_gen1_rows_selected, profile_columns = c(3:ncol(table1)), range_ = c(0, max(table1[3:ncol(table1)])))
        
        }
      }
    })
  })
################################################################################
##                              Component 4: Upload                           ##
################################################################################ 
  
  options(shiny.maxRequestSize = 10*1024^2)
  
  observe(
    output$reactive_sliders4a <- renderUI({
      
        colFile <- input$file1[1,]
        
        if (is.null(colFile))
          return(NULL)
        
        if(endsWith(colFile$datapath, ".xlsx" )){
          
          colputFile <- read.xlsx(colFile$datapath)
          
        }else{
          
          colputFile <- read.csv2(colFile$datapath)
          
        }
        if(input$Radio6 == "Normalized counts"){
          
          slidersUI("6", colnames(colputFile)[2:ncol(colputFile)])
          
        }else{
          
          rownames(colputFile) = colputFile[,1]
          colputFile = colputFile[,-1]
          celltypes = unique(as.character(colputFile[1,]))
          slidersUI("8", celltypes)
          
        }
        
    })
  )
 
  observe(
    output$Component4a <- renderPlotly({
  
          inFile <- input$file1

          if (is.null(inFile))
            return(NULL)
          
          if(endsWith(inFile$datapath, ".xlsx" )){
            
            inputFile <- read.xlsx(inFile$datapath)
            
          }else{
            inputFile <- read.csv2(inFile$datapath)
          }
          
          Variables = callModule(sliders_mod, "6", colnames(inputFile))
          output_length = 100
          top_df <<- spot(inputFile, Variables, colnames(inputFile)[2:length(inputFile)], Candidate_Number = output_length, preamble = 1)
          colnames(top_df)[1] <<- "Genes"
          top_df[2:length(top_df)] = round(top_df[,2:length(top_df)], 2)
          top_df_t = as.data.frame(t(top_df))
            
          if(input$Radio5 == "Table"){
            
            generate_table(top_df[,1:(ncol(top_df)-1)], c_orientation = c("left", "center"), h_orientation = c("left", "center"))
            
          }else if(input$Radio5 == "Bar chart"){
            
            subplot_profiles(top_df, top_df_t, 1:4, profile_columns = c(3:(ncol(top_df)-1)), range_ = c(0, max(top_df[2:ncol(top_df)])), Norm = "[TPM]", title_ = "Expression Profiles")
          
          }else{
            
            bulk_dotplot(top_df[1:10,1:(ncol(top_df)-1)], label = top_df[1:10, 1], boarder1 = 1, boarder2 = 3)
          }
    })
  )
  
  observe(
    output$Component4b <- renderPlotly({
      
      inFile2 <- input$file1
      
      if (is.null(inFile2))
        return(NULL)
      
      if(endsWith(inFile2$datapath, ".xlsx" )){
        show_modal_progress_line()
        inputFile2 <- read.xlsx(inFile2$datapath)
        
      }else{
        show_modal_progress_line()
        inputFile2 <- read.csv2(inFile2$datapath)
      }
      update_modal_progress(0.2)
      rownames(inputFile2) = inputFile2[,1]
      inputFile2 = inputFile2[,-1]
      celltypes = unique(as.character(inputFile2[1,]))
      update_modal_progress(0.4)
      m <- matrix(0, ncol = length(celltypes), nrow = (nrow(inputFile2) - 1))
      pool_df = data.frame(m)
      counter = 1
      for(i in celltypes){
        pool = inputFile2[2:nrow(inputFile2),(inputFile2[1,] == i)]
        pool_numeric = apply(pool, 2, as.numeric)
        pool_sums = rowSums(pool_numeric)
        pool_df[,counter] = pool_sums
        colnames(pool_df)[counter] = i
        counter = counter + 1
      }
      update_modal_progress(0.6)
      y = calcNormFactors(pool_df)
      
      counts_norm = inputFile2
      for(j in 1:ncol(inputFile2)){
        cell_type = as.character(inputFile2[1,j])
        counts_norm[2:nrow(counts_norm),j] = as.numeric(counts_norm[2:nrow(counts_norm),j]) * as.numeric(y[cell_type])
      }
      update_modal_progress(0.8)
      m_norm <- matrix(0, ncol = length(celltypes), nrow = (nrow(counts_norm)-1))
      pool_df_norm = data.frame(m_norm)
      counter = 1
      for(i in celltypes){
        pool = counts_norm[-1,(counts_norm[1,]==i)]
        pool_numeric = apply(pool,2, as.numeric)
        pool_sums = rowSums(pool_numeric)
        pool_df_norm[,counter] = pool_sums
        colnames(pool_df_norm)[counter] = i
        counter = counter + 1
      }
      
      update_modal_progress(0.9)
      up_genes = round(pool_df_norm, 2)
      rownames(up_genes) = rownames(inputFile2)[2:nrow(inputFile2)]
      up_genes$Gene_ID = rownames(up_genes)
      up_genes_sort = up_genes[,c(6,1:5)]
      remove_modal_progress()
      
      if(input$method == "SPOT"){
        
        Variables = callModule(sliders_mod, "8", colnames(up_genes_sort))
        output_length = 100
        top_df <<- spot(up_genes_sort, Variables, colnames(up_genes_sort)[2:length(up_genes_sort)], Candidate_Number = output_length, preamble = 1)
        colnames(top_df)[1] <<- "Genes"
        top_df[2:length(top_df)] = round(top_df[,2:length(top_df)], 2)
        top_df_t = as.data.frame(t(top_df))
        
        if(input$Radio5 == "Table"){
          
          generate_table(top_df[,1:(ncol(top_df)-1)], c_orientation = c("left", "center"), h_orientation = c("left", "center"))
          
        }else if(input$Radio5 == "Bar chart"){
          
          subplot_profiles(top_df, top_df_t, 1:4, profile_columns = c(3:(ncol(top_df)-1)), range_ = c(0, max(top_df[2:ncol(top_df)])), Norm = "[TPM]", title_ = "Expression Profiles")
          
        }else{
          
          bulk_dotplot(top_df[1:10,1:(ncol(top_df)-1)], label = top_df[1:10, 1], boarder1 = 1, boarder2 = 3)
        }
      }else if(input$method == "DEA"){
        observeEvent(input$action_DEA_upload,
                     
                     { 
                       show_modal_progress_line() # show the modal window
                       library(Seurat)
                       # import::from(Seurat, CreateSeuratObject, NormalizeData, FindMarkers)
                       Variables <- callModule(sliders_mod, "8", colnames(up_genes_sort))
                       States <- celltypes
                       Indents <- as.matrix(inputFile2[1,])
                       Sc_counts <- inputFile2[-1,]
                       SS3_plasmo = CreateSeuratObject(counts = Sc_counts, project = "SPOT", min.cells = 3, min.features = 20)
                       Idents(SS3_plasmo) = as.character(Indents)
                       SS3_plasmo <- NormalizeData(SS3_plasmo)
                       SS3_plasmo
                       
                       update_modal_progress(0.2)# update progress bar value
                       if(input$algorithm == "Wilcox"){
                         library(limma)
                         #import::from(limma, rankSumTestWithCorrelation)
                         Markers <- FindMarkers(SS3_plasmo, ident.1 = States[which(Variables > 0)], ident.2 = States[which(Variables == 0)], test.use = "wilcox")
                         #detach("imports")
                         #detach(package:limma, unload = T)
                       }else if(input$algorithm == "MAST"){
                         library(MAST)
                         #import::from(MAST, FromMatrix, zlm)
                         Markers <- FindMarkers(SS3_plasmo, ident.1 = States[which(Variables > 0)], ident.2 = States[which(Variables == 0)], test.use = "MAST")
                         #detach("imports")
                         detach(package:MAST, unload = T)
                       }else if(input$algorithm == "DESeq2"){
                         library(DESeq2)
                         #import::from(DESeq2, DESeqDataSetFromMatrix, estimateSizeFactors, estimateDispersions, nbinomWaldTest, results)
                         Markers <- FindMarkers(SS3_plasmo, ident.1 = States[which(Variables > 0)], ident.2 = States[which(Variables == 0)], test.use = "DESeq2")
                         #detach("imports")
                         detach(package:DESeq2, unload = T)
                       }
                       update_modal_progress(0.9)
                       remove_modal_progress()
                       Markers_up = subset(Markers, Markers[,"avg_log2FC"] > 0)
                       Markers_up$Gene_ID = str_replace(rownames(Markers_up), "-", "_")
                       sc_dea = subset(up_genes_sort, rownames(up_genes_sort) %in% Markers_up$Gene_ID)
                       rownames(sc_dea) = sc_dea$Gene_ID
                       merge_df = merge(Markers_up, sc_dea, sort = F)
                       top_df <<- reactive({as.data.frame(merge_df %>% dplyr::arrange(desc("p_val_adj")))})

                       top_df_t = reactive({as.data.frame(t(top_df()))})

                       output$Component4ba <- renderPlotly({

                         if(input$Radio5 == "Table"){

                           generate_table(top_df()[,1:ncol(top_df())], c_orientation = c("left", "center"), h_orientation = c("left", "center"))
                           
                         }else if(input$Radio5 == "Bar chart"){

                           subplot_profiles(top_df(), top_df_t(), 1:4, profile_columns = c(7:ncol(top_df())), range_ = c(0, max(top_df()[2:ncol(top_df())])), Norm = "[TPM]", title_ = "Expression Profiles")
                           
                         }else if(input$Radio5 == "Dot plot"){

                           dotplot_high_pct1 = subset(sc_dot_plot, sc_dot_plot$Gene %in% top_df()[1:10,1])
                           dotplot_high_pct1[,3] = dotplot_high_pct1[,3]*20
                           sorter_df = data.frame(0:9)
                           rownames(sorter_df) = top_df()[1:10,1]
                           Rank = c(1:length(dotplot_high_pct1$Gene))
                           j = 1
                           for(i in dotplot_high_pct1$Gene){
                             Rank[j] = sorter_df[i,1]
                             j = j + 1
                           }
                           dotplot_high_pct1$Gene = paste(Rank, dotplot_high_pct1$Gene)
                           sc_dotplot(dotplot_high_pct1)
                         }
                       })
                     }
        )
      }
    })
    
  )

  output$downloadData1 <- downloadHandler(
    filename = function() {
      paste0("SPOT results ", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      write.xlsx(top_df[,1:(ncol(top_df)-1)], file)
    }
  )
  output$downloadData2 <- downloadHandler(
    filename = function() {
      paste("SPOT results ", Sys.Date(), ".xlsx", sep="")
    },
    content = function(file) {
      write.xlsx(top_df()[,1:(ncol(top_df())-1)], file)
    }
  )
  output$downloadData3 <- downloadHandler(
    filename = function() {
      paste("SPOT results ", Sys.Date(), ".xlsx", sep="")
    },
    content = function(file) {
      write.xlsx(top_df[,1:(ncol(top_df)-1)], file)
    }
  )
  output$downloadData4 <- downloadHandler(
    filename = function() {
      paste0("SPOT results ", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv2(top_df[,1:(ncol(top_df)-1)], file)
    }
  )
  output$downloadData5 <- downloadHandler(
    filename = function() {
      paste("SPOT results ", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv2(top_df()[,1:(ncol(top_df())-1)], file)
    }
  )
  output$downloadData6 <- downloadHandler(
    filename = function() {
      paste("SPOT results ", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv2(top_df[,1:(ncol(top_df)-1)], file)
    }
  )
  
################################################################################
##                              Component 5: About                            ##
################################################################################
  
}