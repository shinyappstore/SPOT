################################################################################
##                                   UI Functions                             ##
################################################################################
#describe positions for sliders
usher <- function(data){
  a = 1
  b = a + floor((length(data)-1)/4)
  c = b + 1
  d = c + floor((length(data)-2)/4)
  e = d + 1
  f = e + floor((length(data)-3)/4)
  g = f + 1
  h = g + floor((length(data)-4)/4)
  return(c(a, b, c, d, e, f, g, h))
}

#create sliders with columnsname as label
create_sliders <- function(data){
  
  sliders = lapply(data, FUN = function(i){
    sliderTextInput(
      inputId = i,
      label = paste("Weight:", strsplit(i, split = "-")[[1]][2]),
      choices = c(0, 0.5, 1, 1.5, 2),
      grid = TRUE,
      selected = 0
      )
    
    })
  boarder = usher(data)
  list(column(2, sliders[boarder[1]:boarder[2]]),
       column(2, sliders[boarder[3]:boarder[4]]),
       column(2, sliders[boarder[5]:boarder[6]]),
       column(2, sliders[boarder[7]:boarder[8]])
       
  )
}
################################################################################
##                                  Modules                                   ##
################################################################################

#create namespace for every set of sliders
slidersUI <- function(id, data){
  ns <- NS(id)
  tagList(
    create_sliders(ns(data))
    )
}

sliders_mod <- function(input, output, session, data){
  slider_vals <- c()
  for(i in 1:length(data)) {
    slider_vals <- c(slider_vals, input[[data[i]]])
  }
  return(slider_vals)
}
################################################################################
##                                  Server functions                          ##
################################################################################

#plot profiles as barplot
plot_profiles <- function(data, data1 = as.data.frame(t(data)), ID, profile_columns, color = "#c00000",
                          grouped = FALSE, profile_columns2 = NULL, range_w = c(0,12), Norm_str = "[TMM]", title = F) {
  
  data1$names = rownames(data1)
  data1$names <- factor(data1$names, levels = data1[["names"]])
  
  profile <- plot_ly(data1,
                     x = ~names[profile_columns],
                     y = t(data[ID, profile_columns]),
                     name = data1[1, ID],
                     type = "bar",
                     marker = list(color = color,
                                   line = list(color = "#506784",
                                              width = 1.5))
  ) %>% plotly::layout( title = title,  
                xaxis = list(title = data1[1, ID]),
                yaxis = list(title = "Expression", range = range_w)
  )  
  if(grouped == TRUE){
    profile <- profile %>% add_trace(x = data1$names[profile_columns2],
                                     y = t(data[ID, profile_columns2]),
                                     name = paste(data1[1, ID], "B"),
                                     marker = list(color = c('blue'))
    )                      
  }
  profile                     
} 

#create subplots
subplot_profiles <- function(dat, dat1, Gene_names, profile_columns = c(5:14), col_vec = c( "#618C84","#726E75","#948B89","#D0D1AC"), range_ = c(0:12), Norm = "[TMM]", title_ = "RNAseq Profile"){
  plots = list()
  j = 0
  for(i in Gene_names){
    j = j + 1
    plot = plot_profiles(data = dat, data1 = dat1, ID = as.numeric(i) , profile_columns = profile_columns, TRUE, color = col_vec[j], range_w = range_ , Norm_str = Norm, title = title_)
    plots = list.append(plots, plot)
  }
  subplot(plots, nrows = ceiling(length(Gene_names)/2), margin = 0.06, titleX = TRUE, titleY = TRUE)
}  

#generate tables 
generate_table <- function(data, columns = c(1:ncol(data)), columnwidth = 80, header_size = 11, font_size = 10, 
                           c_orientation = c("left", "left", "center"), 
                           h_orientation = c("left", "left", "center"),
                           Table_rows = 50 ){
  m <- list(l = 3, r = 30, b = 20, t = 5, pad = 4)
  
  plot_ly(type = "table",
          columnwidth = columnwidth,
          header = list(
            values = paste("<b>",t(as.matrix(colnames(data)[columns])), "</b>"),
            line = list(color = '#506784'),
            fill = list(color = '#c00000'),
            align = h_orientation,
            font = list(color = 'white', size = c(header_size))),
          cells = list(
            values = t(as.matrix(unname(data[1:Table_rows, columns]))),
            line = list(color = '#506784'),
            fill = list(color = '#f1f1f1'),
            align = c_orientation,
            font = list(color = '#506784', size = font_size))) %>% 
    plotly::layout(autosize = FALSE, width = 1500, margin = m)
}

#get indices for specific column search
get_indices <- function(data = top_df, TextInput){
  data$rank = rownames(data)
  rownames(data) = data[,1]
  Gene_indices = data[TextInput, "rank"]
  rownames(data) = data$rank 
  data = data[, - ncol(data)]
  return(Gene_indices)
}

#create single cell dotplot
sc_dotplot <- function(data){
  fig <- plot_ly(data, x = ~State, y = ~Gene, text = ~paste("Detected in ", Pct*5, "% of cells", "<br>Average expression in detected cells: ", Expression), type = 'scatter', mode = 'markers',
               marker = list(size = ~Pct, opacity = 0.8), color = ~Expression, colors = c( "blue", "#c00000"),
               hoverinfo = "text")
  fig <- fig %>% plotly::layout(title = 'Gene expression through the lifecycle',
                      xaxis = list(showgrid = T ,
                                   tickvals = list(0,1,2,3,4,5,6,7,8,9),
                                   tickmode = "array"),
                      yaxis = list(showgrid = T, type = "category",
                                   categoryorder = "category descending"))
   fig <- fig %>% plotly::layout(showlegend = F)
}

#create unicolor barplot 
bulk_dotplot <- function(data, ytitle = "Genes", xtitle = "Variables", label = data$Ensembl_ID, boarder1 = 2, boarder2 = 4){
  data$Ensembl_ID = paste(0:9, label)
  data = melt(data[,c(1:boarder1, boarder2:ncol(data))])
  
  fig <- plot_ly(data, x = ~variable, y = ~Ensembl_ID, text = paste("Expression: TPM", data$value), type = 'scatter', mode = 'markers',
                 marker = list(size = ~(value/max(value))*20, opacity = 0.8), color = "#c00000", colors = c( "#c00000"),
                 hoverinfo = "text")
  fig <- fig %>% plotly::layout(title = 'Gene expression through the lifecycle',
                        xaxis = list(showgrid = T , title = xtitle,
                                     tickvals = list(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22),
                                     tickmode = "array"),
                        yaxis = list(showgrid = T, type = "category", title = ytitle,
                                     categoryorder = "category descending"))
  
  
  fig <- fig %>% plotly::layout(showlegend = F)
}
################################################################################
##                                  Ranking                                   ##
################################################################################

#calculate spot score
spot <- function(data, Variables, columns = c(1:ncol(data)), preamble = c(1:2), Candidate_Number = 50){
    data1 = as.matrix(scale(data[,columns]))
  
    data1[which(data1 > 4)] = 4
    if(length(which(Variables > 0)) == 1){
      mean_unwanted_cols = rowMeans(data1[, which(Variables == 0)])
      mean_wanted_cols = data1[, which(Variables > 0)]
    }else if(length(which(Variables > 0)) == 0){
      mean_unwanted_cols = rowMeans(data1[, which(Variables == 0)])
      mean_wanted_cols = 0
    }else if(length(which(Variables == 0)) == 1){
      mean_unwanted_cols = data1[, which(Variables == 0)]
      mean_wanted_cols = data1[, which(Variables > 0)] %*% Variables[which(Variables > 0)]/sum(Variables[which(Variables > 0)])
    }else if(length(which(Variables == 0)) == 0){
      mean_unwanted_cols = 0
      mean_wanted_cols = data1[, which(Variables > 0)] %*% Variables[which(Variables > 0)]/sum(Variables[which(Variables > 0)])
    }else{
      mean_unwanted_cols = rowMeans(data1[, which(Variables == 0)])
      mean_wanted_cols = data1[, which(Variables > 0)] %*% Variables[which(Variables > 0)]/sum(Variables[which(Variables > 0)])
    }
    spot_Score = (mean_wanted_cols - mean_unwanted_cols)*(1 - mean_unwanted_cols)
    merged_df <- cbind(data[,preamble], spot_Score, data[,columns], data[,ncol(data)])
    merged_df = subset(merged_df, (mean_wanted_cols - mean_unwanted_cols) > 0)
    merged_df_sorted <- merged_df %>% dplyr::arrange(desc(spot_Score))
    top_df <- merged_df_sorted[1:Candidate_Number,]
}

#calculate correlation
correlation <- function(data, Variables, columns = c(1:ncol(data)), preamble = c(1:2), Candidate_Number = 50){
  
  data1 = as.matrix(scale(data[,columns]))
  data1[which(data1 > 4)] = 4
  Correlation = cor(t(data1), Variables)
  merged_df <- cbind(data[,preamble], Correlation, data[,columns], data[,ncol(data)])
  merged_df_sorted <- merged_df %>% dplyr::arrange(desc(Correlation))
  top_df <<- merged_df_sorted[1:Candidate_Number,]
}
