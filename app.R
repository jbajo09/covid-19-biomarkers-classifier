# Before the deploy to shinyapps.io, run:
# library(BiocManager)
# options(repos = BiocManager::repositories())
# See here:
# https://community.rstudio.com/t/failing-to-deploy-shinyapp-depending-on-bioconductor-packages/6970/5


# Packages
if (!require("shiny")) install.packages("shiny")
library(shiny)
if (!require("shinydashboard")) install.packages("shinydashboard")
library(shinydashboard)
if (!require("dplyr")) install.packages("dplyr")
library(dplyr)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("KnowSeq")) BiocManager::install("KnowSeq")
library(KnowSeq)
if (!require("reshape2")) install.packages("reshape2")
library(reshape2)
if (!require("caret")) install.packages("caret")
library(caret)
if (!require("ggplot2")) install.packages("ggplot2")
library(ggplot2)
if (!require("ggalluvial")) install.packages("ggalluvial")
library(ggalluvial)
if (!require("DT")) install.packages('DT')
library(DT)
library(waiter)
if (!require("ROCR")) install.packages('ROCR')
library(ROCR)
if (!require("limma")) install.packages('limma')
library(limma)
if (!require("M3C")) install.packages("M3C")
library(M3C)
options(repos = BiocManager::repositories())

setwd(getwd())
load("data.RData")






# File to slightly modify dataPlot function
#source("www/dataPlot.R")

# Define some spinners
spinner_abrir <- tagList(
  spin_folding_cube(),
  span(br(), h4("Loading application..."), style="color:white;")
)

spinner <- tagList(
  spin_chasing_dots(),
  span(br(), h4("Loading..."), style="color:white; display: inline-block;")
)

ui <- dashboardPage(title = "COVID-19 Gene expression analysis", # Title in web browser
                    ## Theme
                    skin = "black",
                    ## Header
                    dashboardHeader(title = span(
                      "",
                      style = "font-family: Lucida Console; font-weight: bold"
                    )),
                    ## Sidebar
                    dashboardSidebar(
                      tags$head(
                        tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
                      ),
                      
                      sidebarMenu(
                        menuItem("Introduction", tabName = "intro", icon = icon("file-alt")),
                        menuItem('Preprocessing', tabName = 'preprocessing', icon = icon('tools')),
                        menuItem("Genes selection", tabName = "genes", icon = icon("dna")),
                        menuItem("Model training", tabName = "entrenamiento", icon = icon("play")),
                        menuItem("Model validation", tabName = "validation", icon = icon("check-circle")),
                        menuItem("Code", tabName = "codigo", icon = icon("code"))
                      )
                    ),
                    ## Body
                    dashboardBody(
                      use_waiter(),
                      # Spinners to show on load, or when the application is busy
                      #waiter_show_on_load(spinner_abrir, color = "#027368"),
                      #waiter_on_busy(spinner, color = "#027368"),
                      tabItems(
                        # Tab 1
                        tabItem(tabName = "intro",
                                
                                h1("About this web application"),
                                tags$i("COVID-19 biomarkers classifier"), "is a web application that allows users to carry out COVID-19 gene expression analysis and develop Machine Learning classifier using COVID-19 biomarkers",
                                br(), br(),
                                
                                "The ", tags$i("COVID-19 biomarkers classifier"), "application is part of the", 
                                tags$a(
                                  " work entitled COVID-19 Biomarkers Recognition & Classification Using Intelligent Systems",
                                  href = "https://github.com/jbajo09/covid-19-biomarkers-classifier",
                                  target="_blank"
                                ),
                                "It's developed in R-Shiny and the code is ",
                                tags$a(
                                  "open source.",
                                  href = "https://github.com/jbajo09/covid-19-biomarkers-classifier",
                                  target="_blank"
                                ),
                                
                                h2("Abstract "),
                                
                                h3(tags$b("Background")),
                                "SARS-CoV-2 has paralyzed mankind due to its high transmissibility along with its associated mortality, causing millions of infections and deaths worldwide. The search for gene expression biomarkers from the host transcriptional responses to infection may help understand the underlying mechanisms by which the virus generates COVID-19. This research proposes a smart methodology integrating different RNA-Seq datasets from SARS-CoV-2, other respiratory diseases, and healthy patients.",
                                h3(tags$b("Methods")),
                                "The proposed pipeline exploits the functionality of the 'KnowSeq' R/Bioc package, attaining a significantly large gene expression dataset, integrating different data sources by increasing the number of samples, and endowing the results with higher statistical significance and robustness in comparison with previous studies in the literature. A detailed preprocessing step was carried out to homogenize the samples and build a clinical decision system for SARS-CoV-2 using machine learning techniques such as feature selection algorithm and supervised classification system. This Clinical Decision System uses the most differentially expressed genes among different diseases (including SARS-Cov-2) to develop a four-class classifier. ",
                                h3(tags$b("Results")),
                                "The multiclass classifier designed can discern SARS-CoV-2 samples, reaching an accuracy equal to 91.5%, a mean F1-Score equal to 88.5%, and a SARS-CoV-2 AUC equal to 94% by using only 15 genes as predictors. A biological interpretation of the gene signature extracted reveals relations with processes involved in viral responses. ",
                                
                                h3(tags$b("Conclusions")),
                                "This work proposes a COVID-19 gene signature composed of 15 genes selected after applying the  feature selection ' minimum Redundancy Maximum Relevance' algorithm. The integration between several RNA-Seq datasets was a success, allowing for a considerable number of samples and therefore providing greater statistical significance to the results than previous studies. Both classification accurate and biological interpretation were achieved by using and studying the selected genes.",
                                
                                
                                # Images
                                br(), br(), br(),
                                fluidRow(column(6, tags$img(src = "ugr.png", height = "100px")),
                                         column(6, tags$img(src = "knowseq.png", height = "120px")))
                                
                        ),
                        
                        # Tab 2
                        tabItem(tabName = "preprocessing",
                                h1("Normalize between arrays"),
                                'If your data is a collection of different series, you can carry out a normalization between arrays, scaling the log-ratios
to have the same median-absolute-deviation (MAD) across arrays',
                                br(),
                                br(),
                                actionButton(inputId = "boton_normalization",
                                             label = "Normalize between arrays",
                                             icon = icon("dna", lib = "font-awesome"),
                                             width = "50%"),
                                br(),
                                br(),
                                
                                conditionalPanel(condition = "input.boton_normalization!=0",
                                                 
                                                 'Normalization between arrays finished',
                                ),
                                
                                br(),
                                br(),
                                
                                h1("Batch effect treatment using SVA"),
                                'In order to treat Batch effect Surrogate Variable Analysis can be carried out',
                                br(),
                                br(),
                                actionButton(inputId = "boton_batch",
                                             label = "Batch effect SVA",
                                             icon = icon("dna", lib = "font-awesome"),
                                             width = "50%"),
                                br(),
                                br(),
                                conditionalPanel(condition = "input.boton_batch!=0",
                                                 
                                                 'Batch effect treatment with SVA completed',
                                ),
                                br(),
                                br(),
                                
                                h1('Outlier removal'),
                                'A majority voting system can be carried out in order to remove outliers. Three different tests are carried out. An interquantile range method, a Kolmogorov Smirnov (K-S) test and MA-plot. ',
                                br(),
                                br(),
                                actionButton(inputId = "boton_outlier",
                                             label = "Remove outlier",
                                             icon = icon("dna", lib = "font-awesome"),
                                             width = "50%"),
                                br(),
                                br(),
                                conditionalPanel(condition = "input.boton_outlier!=0",
                                                 
                                                 'Outliers removed',
                                                 
                                                 textOutput('outlier'),
                                ),
                                
                                br(),
                                
                                conditionalPanel(condition = "input.boton_outlier!=0",
                                                 
                                                 h2("Distribution of classes after preprocessing"),
                                                 
                                                 tableOutput("tabla2"),
                                                 
                                                 textOutput('dim')),
                                
                        ),
                      
                        # Tab 3
                        tabItem(tabName = "genes",
                                'At this point data is splitted in training and test samples with a 80%-20% ratio. In order to replicate the same results the same partition used in our paper is used. Training and test samples are indicated in the repository with the code under training-test file.',
                                h1("DEGs extraction 5-folds CV"),
                                textInput(inputId = "LFC", label = "Select the threshold LFC", value = 1, width = "50%"),
                                
                                sliderInput(inputId = "COV", label = "Select the threshold COV", value = 2, min = 1, max = 3, step = 1, width = "50%"),
                                
                                textInput(inputId = "pvalue", label = "Select the threshold p-value", value = 0.01, width = "50%"),
                                
                                actionButton(inputId = "boton_genes",
                                             label = "Extract DEGs",
                                             icon = icon("dna", lib = "font-awesome"),
                                             width = "50%"),
                                br(),
                                br(),
                                conditionalPanel(condition = "input.boton_genes!=0",
                                                 
                                                 h3("Number of DEGs extracted"),
                                                 br(),
                                                 textOutput("numberDEGs")
                                                 
                                                 
                                ),
                                br(),
                                br(),
                                h1('Feature selection algorithm mRMR'),
                                actionButton(inputId = "boton_ranking",
                                             label = "Apply mRMR",
                                             icon = icon("dna", lib = "font-awesome"),
                                             width = "50%"),
                                br(),
                                
                                conditionalPanel(condition = "input.boton_ranking!=0",
                                                 br(),
                                                 'Ranking obtained with mRMR',
                                                 br(),
                                                 #fluidRow(
                                                 #column(4, h4(tags$b("  MRMR")), tableOutput("genes_mrmr")),
                                                 textOutput('genes_mrmr')
                                                 
                                ),
                        ),   
                        
                        # Tab 4
                        tabItem(tabName = "entrenamiento",
                                h1("Model training"),
                                
                                # Choose number of folds
                                selectInput("number_folds",
                                            label = "Number of folds",
                                            choices = c(3, 5, 10),
                                            selected = 5,
                                            width = "50%"),
                                
                                # Train model button
                                actionButton(inputId = "boton_model_training",
                                             label 
                                             = "Train model",
                                             icon = icon("play", lib = "font-awesome"),
                                             width = "50%"),
                                
                                br(),
                                br(),
                                # OPtimal k
                                textOutput("optimal_svm"),
                                
                                br(),
                                h2("Training CV plot"),
                                br(),
                                plotOutput("train", width = "90%", height = '400px'),
                                br()
                                
                        ),
                        
                        # Tab 5
                        tabItem(tabName = "validation",
                                h1("Model validation"),
                                
                                
                                sliderInput(inputId = "numero_genes_validation", label = "Select the number of genes to use:",
                                            value = 15, min = 1, max = 49, step = 1, width = "50%"),
                                
                                actionButton(inputId = "boton_model_validation",
                                             label = "Validate model in test",
                                             icon = icon("play", lib = "font-awesome"),
                                             width = "50%"),
                                
                                br(),
                                br(),
                                
                                plotOutput("results_validation",
                                           width = "60%"),
                                br(),
                                br(),
                                
                                h2('Area Under the Curve'),
                                actionButton(inputId = "boton_calculate_AUC",
                                             label = "Calculate AUC",
                                             icon = icon("play", lib = "font-awesome"),
                                             width = "50%"),
                                textOutput("meanAUC"),
                                plotOutput("AUC",
                                           width = "60%"),
                                br(),
                                br(),
                                h1('T-Distributed Stochastic  Neighbour Embedding'),
                                br(),
                                sliderInput(inputId = "numero_genes_TSNE", label = "Select the number of genes to use:",
                                            value = 4, min = 1, max = 50, step = 1, width = "50%"),
                                br(),
                                actionButton(inputId = "boton_tsne",
                                             label = "Carry out T-SNE without any preprocessing step",
                                             icon = icon("play", lib = "font-awesome"),
                                             width = "60%"),
                                br(),
                                conditionalPanel(condition = "input.boton_tsne!=0",
                                                 plotOutput('tsne2', width = '60%')),
                            
                                
                        ),
                        
                        # Tab 6
                        tabItem(tabName = "codigo",
                                h1("Code"),
                                tags$h4(
                                  "In ", tags$a(href = "https://github.com/jbajo09/covid-19-biomarkers-classifier", "this repository"),
                                  "you can find the code of this web application.")
                        )
                      ) # Close tabs
                    ) # Close dashboard body
) # Close dashboard page

# Extend size of accepted files (40MB instead of the 5MB - default)
options(shiny.maxRequestSize = 600*1024^2)

server <- function(input, output){
  
  values <- reactiveValues(ranking = NULL, matrix = NULL , labels=NULL, optimalSVM_train = NULL, DEGs = NULL)
  
  # Server of tab : Preprocessing
  
  
  observeEvent(input$boton_normalization, {
    
    w <- Waiter$new(html = tagList(spin_folding_cube(),
                                   span(br(), br(), br(), h4("Carrying out normalization between arrays..."),
                                        style="color:white;")))  
    w$show()
    
    GSE156063_m <- t(GSE156063[,1:15959])
    GSE156063_lab <- GSE156063$GSE156063_lab
    for (i in 1:length(GSE156063_lab)){
      if (GSE156063_lab[i]=='no virus'){
        GSE156063_lab[i] <- 'NON.VIRAL.ARI'
      } else if (GSE156063_lab[i]=='other virus'){
        GSE156063_lab[i] <- 'VIRAL.ARI'
      }
    }
    
    GSE149273_m <- t(GSE149273[,1:17605])
    GSE149273_lab <- GSE149273$GSE149273_lab
    for (i in 1:length(GSE149273_lab)){
      if (GSE149273_lab[i]=='RVA' | GSE149273_lab[i]=='RVC' ){
        GSE149273_lab[i] <- 'VIRAL.ARI'
      }
    }
    
    GSE152075_m <- t(GSE152075[,1:35486])
    GSE152075_lab <- GSE152075$GSE152075_lab
    
    GSE163151_m <- t(GSE163151[,1:21909])
    GSE163151_lab <- GSE163151$GSE163151_lab
    for (i in 1:length(GSE163151_lab)){
      if (GSE163151_lab[i]=='COVID-19'){
        GSE163151_lab[i] <- 'SC2'
      } else if (GSE163151_lab[i]=='Donor control'){
        GSE163151_lab[i] <- 'Control'
      } else if (GSE163151_lab[i]=='Non-viral acute respiratory illness'){
        GSE163151_lab[i] <- 'NON.VIRAL.ARI'
      } else if (GSE163151_lab[i]=='Viral acute respiratory illness'){
        GSE163151_lab[i] <- 'VIRAL.ARI'
      }
    }
    
    GSE162835_m <- t(GSE162835[,1:20601])
    GSE162835_lab <- GSE162835$GSE162835_labels
    
    # Superserie 5 bacth
    I1 <- intersect(rownames(GSE156063_m), rownames(GSE149273_m))
    I2 <- intersect(I1, rownames(GSE152075_m))
    I3 <- intersect(I2,  rownames(GSE163151_m))
    I4 <- intersect(I3, rownames(GSE162835_m))
    GSE156063_I <- GSE156063_m[which(rownames(GSE156063_m) %in%  I4),]
    GSE149273_I <- GSE149273_m[which(rownames(GSE149273_m) %in%  I4),]
    GSE152075_I <- GSE152075_m[which(rownames(GSE152075_m) %in%  I4),]
    GSE163151_I <- GSE163151_m[which(rownames(GSE163151_m) %in%  I4),]
    GSE162835_I <- GSE162835_m[which(rownames(GSE162835_m) %in%  I4),]
    GSE156063_I <- GSE156063_I[order(rownames(GSE156063_I)),] 
    GSE149273_I <- GSE149273_I[order(rownames(GSE149273_I)),] 
    GSE152075_I <- GSE152075_I[order(rownames(GSE152075_I)),] 
    GSE163151_I <- GSE163151_I[order(rownames(GSE163151_I)),] 
    labels <- c(GSE156063_lab,GSE149273_lab,GSE152075_lab,GSE163151_lab,GSE162835_lab)
    matrix <- cbind(GSE156063_I,GSE149273_I,GSE152075_I,GSE163151_I,GSE162835_I)


    matrix_scale <- normalizeBetweenArrays(matrix, method = 'scale')
    
    w$hide()
    
    values$matrix <- matrix_scale
    values$labels <- labels
  })
  
  observeEvent(input$boton_batch, {
    
    w <- Waiter$new(html = tagList(spin_folding_cube(),
                                   span(br(), br(), br(), h4("Carrying out SVA..."),
                                        style="color:white;")))  
    w$show()
    
    
    
    matrix <- values$matrix
    labels <- values$labels

    matrix_batch <- batchEffectRemoval(matrix,labels, method = 'sva')
    
    
    w$hide()
    
    values$matrix <- matrix_batch
  })
  
  observeEvent(input$boton_outlier, {
    
    w <- Waiter$new(html = tagList(spin_folding_cube(),
                                   span(br(), br(), br(), h4("Removing outlier..."),
                                        style="color:white;")))  
    w$show()
    
    matrix <- values$matrix
    labels <- values$labels
    outliers <- RNAseqQA(matrix,toRemoval = TRUE, toPNG = FALSE, toPDF = FALSE) 
    w$hide()
    
    output$outlier <- renderText(outliers$outliers)
    
    values$labels <- labels[-which(colnames(matrix) %in% outliers$outliers)]
    values$matrix <- outliers$matrix
    
    # Table
    output$tabla2 <- renderTable({
      
      tabla_aux1 <- as.data.frame(table(values$labels))
      return(tabla_aux1)
    })
    
    output$dim <- renderText(paste0('The final dimension matrix contain information of  ', dim(values$matrix)[1], ' genes for ', dim(values$matrix)[2], ' samples'))
  })
  
  # Server of tab: Genes selection ------
  
  w <- Waiter$new(html = tagList(spin_folding_cube(),
                                 span(br(), br(), br(), h4("Running DEGs 5-folds CV..."),
                                      style="color:white;")))
  
  observeEvent(input$boton_genes, {
    
    w$show()
    
    w$hide()
    
    #values$DEGsMatrix_rna_batch <- DEGsMatrix_rna
    
    #DEGs extraction
    w <- Waiter$new(html = tagList(spin_folding_cube(),
                                   span(br(), br(), br(), h4("Running DEGs 5-folds CV..."),
                                        style="color:white;")))
    w$show
    indices = index
    matrix <- values$matrix
    labels <- values$labels
    particion.entrenamiento <- matrix[,indices]
    particion.test <- matrix[,-indices]

    # Labels
    labels_train <- labels[indices]
    labels_test  <- labels[-indices]
    
    set.seed(2332)
    DEGs_CV <- DEGsExtraction(particion.entrenamiento, as.factor(labels_train), lfc=as.numeric(input$LFC), cov = as.numeric(input$COV), pvalue = as.numeric(input$pvalue), number = Inf, CV = TRUE)
    DEGs <- DEGs_CV$Common_DEGs
    output$numberDEGs<- renderText(paste0( length(DEGs), " DEGs were extracted "))
    w$hide()
    values$DEGs <- DEGs
    
  }) # close button
    # mRMR method
  
  observeEvent(input$boton_ranking, {
    w <- Waiter$new(html = tagList(spin_folding_cube(),
                                   span(br(), br(), br(), h4("Running mRMR algorithm..."),
                                        style="color:white;")))
    w$show()
    indices = index
    matrix <- values$matrix
    labels <- values$labels
    particion.entrenamiento <- matrix[,indices]
    particion.test <- matrix[,-indices]

    # Labels
    labels_train <- labels[indices]
    labels_test  <- labels[-indices]
    
    DEGs <- values$DEGs 
    
    set.seed(2332)
    mrmrRanking <- featureSelection(t(particion.entrenamiento), as.factor(labels_train), DEGs,
                                    mode = "mrmr")
    mrmrRanking <- names(mrmrRanking)
    w$hide()
    
    values$ranking <- mrmrRanking
    
    
    # Ranking tables
    
    output$genes_mrmr <- renderText(mrmrRanking)
    
    
  }) # Close button
  
  # Server of tab: Model training ------
  
  w2 <- Waiter$new(html = tagList(spin_folding_cube(),
                                  span(br(), br(), br(), h4(""),
                                       style="color:white;")))  
  
  observeEvent(input$boton_model_training, {
    
    w2$show()
    
    indices = index
    matrix <- values$matrix
    labels <- values$labels
    particion.entrenamiento <- matrix[,indices]
    particion.test <- matrix[,-indices]
    
    # Labels
    labels_train <- labels[indices]
    labels_test  <- labels[-indices]
    
    
    w2$hide()
    
    
    ranking <- values$ranking
    
    w3 <- Waiter$new(html = tagList(spin_folding_cube(),
                                    span(br(), br(), br(), h4("Training kNN algorithm..."),
                                         style="color:white;")))  
    w3$show()
    set.seed(2332)
    results_cv <- svm_trn(t(as.matrix(particion.entrenamiento)), labels_train, ranking,
                          numFold = as.numeric(input$number_folds))
    values$optimalSVM_train <- results_cv$bestParameters
    
    w3$hide()
    
    output$optimal_svm <- renderText(paste0("\nOptimal number of neighbours k = ", results_cv$bestK))
    
    output$train <- renderPlot({
      plot(results_cv$accuracyInfo$meanAccuracy, type = 'l', col= 'black', ylab='Metric Performance', xlab='Genes', lwd=2, ylim = c(0,1), panel.first = grid(col='gray45'), cex.axis=1.3,cex.lab=1.3)
      lines(results_cv$sensitivityInfo$meanSensitivity, col='blue', lwd=2, lty=2)
      lines(results_cv$specificityInfo$meanSpecificity, col='#FF8B00', lwd=2, lty=4)
      lines(results_cv$F1Info$meanF1, col='red', lwd=2, lty=4)
      legend(x=15.9 ,y =0.8405, c('Accuracy', 'Sensitivity','Specificity','F1-Score'), lty = c(1,2,4,5), col = c('black','blue','#FF8B00','red'), cex=1.3)
      
      
    })
    
  }) 
  
  # Server of tab: Model validation  ------
  
  w3 <- Waiter$new(html = tagList(spin_folding_cube(),
                                  span(br(), br(), br(), h4("Validating model..."),
                                       style="color:white;")))  
  
  observeEvent(input$boton_model_validation, {
    
    w3$show()
    
    set.seed(2332)
    indices = index
    matrix <- values$matrix
    labels <- values$labels
    particion.entrenamiento <- matrix[,indices]
    particion.test <- matrix[,-indices]
    
    # Labels
    labels_train <- labels[indices]
    labels_test  <- labels[-indices]
    
    w3$hide()
    
    
    ranking <- values$ranking
    
    
    w3 <- Waiter$new(html = tagList(spin_folding_cube(),
                                    span(br(), br(), br(), h4("Validating SVM algorithm..."),
                                         style="color:white;")))  
    w3$show()
    set.seed(2332)
    results_validation <- svm_test(train = t(particion.entrenamiento), as.factor(labels_train),
                                   test = t(particion.test), as.factor(labels_test),
                                   ranking, bestParameters  = values$optimalSVM_train)
    w3$hide()
    
    
    output$results_validation <- renderPlot({
      tabla <- results_validation$cfMats[[input$numero_genes_validation]]$table
      plotConfMatrix(tabla)
    })
    
  })
  
  w4 <- Waiter$new(html = tagList(spin_folding_cube(),
                                  span(br(), br(), br(), h4("Calculating AUC..."),
                                       style="color:white;")))  
  
  
  #AUC
  observeEvent(input$boton_calculate_AUC, {
    
    w4$show()
    
    set.seed(2332)
    indices = index
    matrix <- values$matrix
    labels <- values$labels
    particion.entrenamiento <- matrix[,indices]
    particion.test <- matrix[,-indices]
    
    # Labels
    labels_train <- labels[indices]
    labels_test  <- labels[-indices]
    w4$hide()
    
    
    ranking <- values$ranking
    
    w4 <- Waiter$new(html = tagList(spin_folding_cube(),
                                    span(br(), br(), br(), h4("Calculating AUC..."),
                                         style="color:white;")))
    
    w4$show()
    set.seed(2332)
    response <- as.factor(labels_test)
    aucs <- rep(NA, length(levels(response))) # store AUCs
    legendLabels <- as.character()
    colours <- c('red','blue','green','black')
    a <- list()
    
    for (i in seq_along(levels(response))) {
      cur.class <- levels(response)[i]
      binaryTraining.labels <- as.factor(labels_train == cur.class)
      binaryTest.labels <- as.factor(labels_test == cur.class)
      
      tri_svm_cv_mrmr_results <- svm_trn(t(as.matrix(particion.entrenamiento)), binaryTraining.labels, ranking, numFold = 5)
      
      tri_svm_test_mrmr_results <- svm_test(t(as.matrix(particion.entrenamiento)), binaryTraining.labels,t(as.matrix(particion.test)), binaryTest.labels, ranking, bestParameters = tri_svm_cv_mrmr_results$bestParameters)
      
      score <- tri_svm_test_mrmr_results$predictions[[as.numeric(input$numero_genes_validation)]]
      score <- as.vector(score)
      score[score=='FALSE'] <- 0
      score[score=='TRUE'] <- 1
      binaryTest.labels <- as.vector(binaryTest.labels)
      binaryTest.labels[binaryTest.labels=='FALSE'] <- 0
      binaryTest.labels[binaryTest.labels=='TRUE'] <- 1
      pred <- prediction(as.numeric(score), as.numeric(binaryTest.labels))
      perf <- performance(pred, "tpr", "fpr")
      roc.x <- unlist(perf@x.values)
      roc.y <- unlist(perf@y.values)
      a[[i]] <- roc.x
      a[[i+4]] <- roc.y
      # store AUC
      auc <- performance(pred, "auc")
      auc <- unlist(slot(auc, "y.values"))
      aucs[i] <- auc
      legendLabels[i] <- paste(levels(response)[i], " AUC: ",format(round(aucs[i], 4), nsmall = 3),sep = "")
    }
    
    output$meanAUC <- renderText((paste0("Mean AUC under the precision-recall curve is: ", round(mean(aucs), 2)))) 
    
    w4$hide()
    
    output$AUC <- renderPlot({
      plot(x=NA, y=NA, xlim=c(0,1), ylim=c(0,1),ylab="Sensitivity", xlab="1 - Specificity", bty='n', cex.lab=1.3, cex.axis=1.3)
      lines(x=c(0,1), c(0,1))
      lines(a[[5]] ~ a[[1]], col = colours[1], lwd = 2)
      lines(a[[6]] ~ a[[2]], col = colours[2], lwd = 2)
      lines(a[[7]] ~ a[[3]], col = colours[3], lwd = 2)
      lines(a[[8]] ~ a[[4]], col = colours[4], lwd = 2)
      legend(x=0.38 ,y =0.305, legendLabels, lty=1, ncol= 1,inset = c(0,0),  col = colours, cex = 1.3,lwd=3)
      
      
      
    })
    
    
  })
  
  observeEvent(input$boton_tsne, {
    
    w5 <- Waiter$new(html = tagList(spin_folding_cube(),
                              span(br(), br(), br(), h4("Carrying out T-SNE ..."),
                                   style="color:white;")))
    matrix <- values$matrix
    labels <- values$labels
    ranking <- values$ranking
    
    w5$show()
    
    output$tsne2 <- renderPlot({
      tsne(matrix[which(rownames(matrix)%in% ranking[1:as.numeric(input$numero_genes_TSNE)]),],labels=as.factor(labels),controlscale=TRUE, scale=3, colvec = c('red','blue','green','black'),seed = 1,axistextsize=0)
    })
    w5$hide()
    
  })
}

shinyApp(ui, server)