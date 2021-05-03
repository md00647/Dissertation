###LIBRARIES REQUIRED FOR THE APPLICATION TO RUN--------------------------------------------------------
library(DT)
library(dplyr)
library(data.table)  
library(R.utils)
library(stringr)
library(shinyjs)
library(tidyr)
library(plyr)
library(data.table)
library(corrplot)
library(RColorBrewer)
library(colorspace)
library(fpc)
library(factoextra)
library(shiny)

###LOAD IN CSV CONTAINING THE DATASETS SAVED ON COMPUTER AND THEIR PATHS
Dataset_Paths <- read.csv("Dataset_Paths.csv") 

# Sort them alphabetically for ease of selection in the application
Dataset_Paths <- Dataset_Paths[order(Dataset_Paths[,'Dataset_Names']), ]

#UI PAGE ------------------------------------------------------------------------------------------------
ui <- fluidPage(
  useShinyjs(),
  # App title ----
  titlePanel("SNP Comparison"),
  
  # Sidebar layout----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      tags$head(tags$style(".shiny-plot-output{height:100vh !important;}")), 
      # Inputs for SNPs 
      textInput("text_in", "Type in your desired SNPs, using their rsid, separated by commas with no space"),
      fileInput("file_in", "Upload SNP list as a csv file", multiple = FALSE),
      # Key Dataset selection
      selectInput("data_in", "Key Dataset for Phenotype", choices = Dataset_Paths$Dataset_Names, multiple = FALSE),
      # Datasets to compare against
      checkboxGroupInput("data_comp", "Comparison Datasets", choices = Dataset_Paths$Dataset_Names),
      # Dealing with missing data
      selectInput("miss", "How do you want to deal with Missingness?", choices = c("Mean_Imputation", "SNP_Removal"), selected = NULL),
      # Generate Outputs
      actionButton("go", "Generate Output"),
      # Download Table Button
      downloadButton("dl", "Download your Summary Table"),
      downloadButton("dl2", "Download your HeatMap"),
      # Reset Button
      actionButton("ref", "Refresh Inputs")
    ), 
    
    # Main panel for displaying SNP Outputs ----
    mainPanel(tabsetPanel(type = "tabs", 
                          tabPanel("Summary Table", dataTableOutput("SNP_out")),
                          tabPanel("Heat Map", plotOutput("map"), style='height: 100%')
    )
    )   
    
  )
) 

#SERVER PAGE ------------------------------------------------------------------------------------------------
server <- function(input, output,session) {
  
  ###LOADING INPUTS### ----------------------------------------------------------------------------------------
  
  # Use Text SNP Input to generate a table containing all of the imputed ones
  snp_text <- reactive({
    
    inText <- input$text_in
    
    if (is.null(inText))
      return(NULL)
    
    textIn <- as.data.frame(unlist(inText, ",")) 
    names(textIn)[1]<- "SNP"
    return(textIn)
  })
  
  # Use File SNP Input to generate a table containing all of the imputed ones
  snp_file <- reactive({
    
    inFile <- input$file_in
    
    if (is.null(inFile))
      return(NULL)
    
    fullcsv <- read.csv(inFile$datapath, header = FALSE, stringsAsFactors = FALSE)
    names(fullcsv)[1]<- "SNP"
    return(fullcsv)
  })
  
  # Merge SNP inputs so we have all of the desired SNPs
  snp_all <- reactive({
    rbind(snp_text(), snp_file()) %>% unique()
  })
  
  ### LOAD DATASETS FROM COMPUTER USING THE INPUTS ### ---------------------------------------------------
  
  # Filter Paths csv based on user selected inputs
  load <- reactive({
    subset(Dataset_Paths, Dataset_Names %in% input$data_in | Dataset_Names %in% input$data_comp)
  })
  
  # Import the required datasets as a list
  datasets_list <- reactive({
    lapply(load()$Path, read.csv, header = TRUE, sep=",", stringsAsFactors = FALSE)
  }) 
  
  # Convert list to a dataset so it can be filtered
  data_req  <- reactive({
    do.call(rbind.data.frame, datasets_list())
  })
  
  ###USING SNP INPUT TO FILTER THE  DATASETS###----------------------------------------------------------------
  
  # Filter the key dataset by the SNP input
  searches <- reactive({
    filter(data_req(), Dataset == input$data_in) %>%
      subset(SNP %in% snp_all()$SNP)  
  })
  
  # Filter all the datasets based on the SNPs desired which are in the key dataset
  all_snps <- reactive({
    subset(data_req(), SNP %in% searches()$SNP)
  })
  
  ### CREATE Z-SCORE COLUMN AND FIX THE COLUMNS ### ----------------------------------------------------------------------
  
  # Z-Score will be used to compare the association of each phenotype with the key phenotype
  data_z <- reactive({
    mutate(all_snps(), "Z-Score" = as.numeric(BETA/SE))
  })
  
  # Remove the rows which have multiple alleles in their EA and NEA as it would not work to compare them
  #Paste the columns together and then keep rows where the number of characters is equal to 3, one for each column and a space
  data_alleles <- reactive({
    data_z()[nchar(paste(data_z()$EA, data_z()$NEA)) == 3,]
  })
  
  # Upper case the Alleles so they are harmonised and easier to compare
  data_upper <- reactive({
    mutate(data_alleles(), EA = toupper(data_alleles()$EA), NEA = toupper(data_alleles()$NEA))
  })
  
  # Ensure the resulting table is a dataframe 
  data_dt <- reactive({
    as.data.table(data_upper())
  })
  
  ### ALLELE ALIGNMENT AND Z-SCORE GENERATION### ---------------------------------------------------------------------------------------
  
  ##Pivot table so Z-Scores are along the top and remove the SE, B and P variables
  pivot <- reactive({
    subset(data_dt(), select = -c(4, 5, 6))
  })
  
  # Filter by the key data column, transpose and merge back onto original dataset so we can compare its alleles against the others
  key_data <- reactive({
    subset(pivot(), pivot()$Dataset %in% input$data_in)
  })
  
  key_data_wide <- reactive({
    pivot_wider(key_data(), names_from = "Dataset", values_from = c("Z-Score", "EA", "NEA", "Gene")) 
  })
  
  all <- reactive({
    merge(pivot(), key_data_wide(), by = "SNP")
  })
  
  # Allele Alignment 
  # If key dataset Z-Scores are negative, flip the alleles and the sign of the Z-Score to make it positive
  data_align <- reactive({ all() %>% 
              mutate(NewZ = ifelse((.[[7]] < 0) == TRUE,  .[[7]]*-1, .[[7]]),
                     NewEA = ifelse((.[[7]] < 0) == TRUE,  .[[9]], .[[8]]),
                     NewNEA = ifelse((.[[7]] < 0) == TRUE,  .[[8]], .[[9]])
           )
        })
  
  
  # Replace the original Z-Score, EA and NEA column values for key datasets with the newly generated ones 
  align_rename <- reactive({
     data_align() %>% replace(7, data_align()$NewZ) %>% replace(8, data_align()$NewEA) %>% replace(9, data_align()$NewNEA) %>% 
      select(-c(11, 12, 13))
  })
  
  # Filter the dataset to not include the key dataset so you do not compare it against itself
  # Waste of processing power otherwise.
  no_key <- reactive({
        filter(align_rename(), Dataset != input$data_in)
        })
  
  # Format the dataset as a list to allow for comparison of alleles. 
  # Convert it to a dataframe
  all2 <- reactive({
    data.frame(lapply(align_rename(), as.character), stringsAsFactors=FALSE)
  })
  
  # Create new Z.Score columns to use for ifelse statement to modify Z-scores according to EA/NEA values of key dataset
  all_mutate <- reactive({ all2() %>%
      dplyr::mutate(
        `Z.Score` = as.numeric(all2()$`Z.Score`),
        `Z.Score.Updated` = as.numeric(0),
        `Z.Score.Updated2` = as.numeric(0)
      )
  })
  
  # An if else statement in programming is a conditional statement that runs a different set of statements depending on whether an expression is true or false
  # If they are the same, return the original Z.Score, otherwise change value to NA
  Z_Upd <- reactive({ all_mutate() %>%
      dplyr::mutate(`Z.Score.Updated` =
                      ifelse(all_mutate()$EA == all_mutate()[,8] &  all_mutate()$NEA ==  all_mutate()[,9],  all_mutate()$`Z.Score`, "NA"))
  })
  
  # If NEA=EA of comparison dataset and EA=NEA of comparison dataset, then multiple the Z-Score by -1, otherwise assign it to NA, as before
  Z_Upd2 <- reactive({ Z_Upd() %>%
      dplyr::mutate(`Z.Score.Updated2` =
                      ifelse(Z_Upd()$NEA == Z_Upd()[,8] & Z_Upd()$EA == Z_Upd()[,9], Z_Upd()$`Z.Score`*-1, Z_Upd()$`Z.Score.Updated`)
      )  
  })
  
  # Modify columns: EA and NEA become the values of the key dataset 
  all_final <- reactive({
    Z_Upd2() %>%
      select(-c(2, 3, 5, 6, 7, 11)) %>% dplyr::rename(`Z-Score` = `Z.Score.Updated2`, 
                                                  `EA` = 3,
                                                  `NEA` = 4,
                                                  `Gene` = 5)
  })
  
  # Key dataset 
  all_key <- reactive({
    align_rename() %>% 
      select(1, 7, 8, 9, 10) %>% 
      unique() %>% 
      mutate(Dataset = input$data_in) 
  }) 
 
  # Rename key dataset so that its values can be used in output table 
   key_rename <- reactive({
     all_key() %>% dplyr::rename(`Z-Score` = 2,
                               `EA` = 3,
                               `NEA` = 4,
                               `Gene` = 5)
   })
   
  # Create final dataset, merging the key dataset and comparison dataset dataframes which contain the updated Z-scores and alleles of key dataset
  total <- reactive({
    rbind(all_final(), key_rename())
  })
  
  # Reshape so all of the Z-Scores are along the top row
  output_z <- reactive({
    reshape(total(), idvar=c("SNP", "EA", "NEA", "Gene"), timevar = "Dataset", direction="wide") %>%
      relocate(paste("Z-Score.",input$data_in, sep = ""), .after = "Gene")
  }) 
  
  ### REMOVE ROWS ### ------------------------------------------------------------------------------------------------------------------------
  #Remove rows with more than 40% NA to reduce the amount of bias if user selects mean imputation
  # These would be removed in SNP removal anyway
  keep_rows <- reactive({
    as.data.frame(output_z()[which(rowMeans(!is.na(output_z())) > 0.6), ])
  })
  
  ### MISSINGNESS ### ------------------------------------------------------------------------------------------------------------------------
  
  all_complete <- reactive({
    
    data2 <- keep_rows()
    
    if (input$miss == "Mean_Imputation") {  ##finds mean of each row
      
      for(i in 4:ncol(data2)) {
        data2[, i ][is.na(data2[, i ])] <- mean(as.numeric(data2[,i ]), na.rm = TRUE) #use mean imputation
        data3 <- data2
        
      }
    } else if (input$miss == "SNP_Removal")   {
      data2 <-  data2[complete.cases(data2),]    ##select the rows which do not contain NA
      data3 <- data2
    }  
    
    return(data3)
  }) 

  ### SUMMARY TABLE OUTPUT ### -----------------------------------------------------------------------------------------------------
  output$SNP_out <- 
    
    DT::renderDataTable({
      
      input$go
      
      shiny::validate(
        need(length(input$data_in) > 0 && length(input$data_comp) > 0 && (length(input$text_in) > 0 | length(input$file_in) > 0), 
             "Please input your SNPs, Key Dataset and Comparison Datasets.")
      )
      
      all_complete()
      
    })
  
  ###HIERACHICAL CLUSTERING### ----------------------------------------------------------------------------------------------
  
  # Call in functions script to perform clustering and heatmap
  # Code from supervisors team
  source('MD_step4_functions.R')  
  
  # Create ID label of rsid and EA/NEA that will be on the y-axis of the heatmap  
  # Make the id column the rownames
  all_complete_id <- reactive({
    id_col <- all_complete() %>%
      dplyr::mutate(id = paste0(all_complete()$SNP,"_",all_complete()$Gene))
    rownames(id_col) = id_col$id
    id_col
  })
  
  # Remove the SNP, EA, NEA, id columns as these are not needed for clustering
  all_cols <- reactive({
    all_complete_id() %>% select(-c("SNP", "EA", "NEA", "Gene", "id"))
  })
  
  # Make all columns numeric 
  all_cols2 <- reactive({data.frame(lapply(all_cols(), function(x) as.numeric(as.character(x))),
                                    check.names=F, row.names = rownames(all_cols()))
  })
  
  # Create empty array based on size of previous datatable
  z_mask <- reactive({
    makeStarMask(all_cols2())
  })
  
  # Create dataframe for matrix
  z_matrix <- reactive({
    all_cols2()
  })
  
  # Truncates to 2sd
  z_matrix_truncated <- reactive({ 
    truncate_2sd(z_matrix())
  })
  
  # Run bootstrap analysis, and then use it for row dendrogram
  
  bootmeans <- list()
  bootbrds <- list()
  clusters <- list()
  
  for (i in 1:2){
    kbest.p <- i
    cboot.hclust <- reactive({ clusterboot(z_matrix_truncated(),
                                           clustermethod=hclustCBI,
                                           method="ward.D", 
                                           k=kbest.p, 
                                           B=100, 
                                           bootmethod = "boot",
                                           count = FALSE)
    })
    bootmeans[[i]] <- reactive({cboot.hclust()$bootmean})
    bootbrds[[i]] <- reactive({cboot.hclust()$bootbrd})
    clusters[[i]] <- list()
    for (j in 1:i){
      clusters[[i]][[j]] <- list()
      clusters[[i]][[j]][[1]] <- reactive({bootmeans()[[i]][j]})
      cluster <- reactive({names(cboot.hclust()$partition)[cboot.hclust()$partition == j]})
      clusters[[i]][[j]][[2]] <- reactive({cluster()})
    }
  }
  
  ###DENDOGRAM GENERATION###-------------------------------------------------------------------------------------------------
  
  clust <- reactive({
    cboot.hclust()$result$result
  })
  
  #rownames
  clust_dendro <- reactive({
    as.dendrogram(clust(), hang=-1)
  })
  
  
  col_clust <- reactive({
    clusterCols(z_matrix_truncated())
  })
  
  #columns
  col_clust_dendro <- reactive({
    as.dendrogram(col_clust(), hang=-1) 
  }) 
  
  
  ###HEATMAP GENERATION### --------------------------------------------------------------------------------------------------
  
  # Use truncated matrix to produce heatmap
  z_matrix_truncated_forplot <- reactive({
    truncate(z_matrix_truncated())
  })
  
  breaks = reactive({
    seq(min(z_matrix_truncated_forplot()), max(z_matrix_truncated_forplot()), 0.01)
  }) 
  col <- reactive({
    colorRampPalette(rev(brewer.pal(8, "RdYlBu")), interpolate = "linear")(length(breaks())-1)
  })
  
  # Heatmap Output
  
  output$map <- shiny::renderPlot({
    
    shiny::validate(
      need(length(input$data_in) > 0 && length(input$data_comp) > 0 && (length(input$text_in) > 0 || length(input$file_in) > 0), "Please input your SNPs, Key Dataset and Comparison Datasets.")
    )
    
    gplots::heatmap.2(as.matrix(z_matrix_truncated_forplot()), 
                      main = "A Heat Map to show the association of SNPs against several phenotypes", 
                      col=col(), 
                      scale="row",
                      key=TRUE, 
                      symkey=FALSE, 
                      density.info="density",
                      Rowv = clust_dendro(),
                      Colv = col_clust_dendro(), 
                      margins=c(15,15),
                      trace="none",
                      srtCol=30, 
                      cexRow = 0.3 + 1/log10(nrow(z_matrix())),
                      cexCol = 0.35 + 1/log10(ncol(z_matrix())),  
                      ylab = "SNPs with key phenotype alleles",
                      xlab = "Phenotypes",  
                      key.title = "Z-score Value") 
    
    
  })
  
  ###DOWNLOAD OUTPUTS###-----------------------------------------------------------------------------------------------
  
  output$dl <- downloadHandler(
    filename = function() {
      paste("SNP comparison across phenotypes","csv", sep = ".")
    },
    content = function(fname) {
      write.csv(all_complete(), fname, row.names = TRUE)
    }
  )
  
  
  output$dl2 <- downloadHandler(
    filename <- function() {
      paste("HeatMap","png", sep = ".") },
    
    content <- function(file) {
      png(file)
      
      shiny::renderPlot({
        gplots::heatmap.2(as.matrix(z_matrix_truncated_forplot()), 
                          main = "The association of SNPs against several phenotypes", 
                          col=col(), 
                          scale="row",
                          key=TRUE, 
                          symkey=FALSE, 
                          density.info="density",
                          Rowv = clust_dendro(),
                          Colv = col_clust_dendro(), 
                          margins=c(15,35),
                          trace="none",
                          srtCol=15, 
                          cexRow = 0.5 + 1/log10(nrow(z_matrix())),
                          cexCol = 0.5 + 1/log10(ncol(z_matrix())),  
                          ylab = "SNPs with comparison phenotype alleles",
                          xlab = "Phenotypes",  
                          key.title = "Z-score Value")
        
      })
      dev.off()
      
    }  
  )
  
  ###RESET the inputs after clicking the refresh button###----------------------------------------------------------------
  
  observeEvent(input$ref,
               {
                 shinyjs::reset("text_in")
                 shinyjs::reset("file_in")
                 shinyjs::reset("data_in")
                 shinyjs::reset("data_comp")
                 shinyjs::reset("miss")
               }) 
}

#ShinyApp --------------------------------------------------------------------------------------------------------
shinyApp(ui = ui, server = server)
