# ==== Loading library ==== 

library(shiny)
library(plotly)
library(shinythemes)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(DT)

# ==== Global variables ====

load("dea.NT.rda")
load("dea.TP.rda")
load("dea.NT.TP.AA.rda")
load("dea.NT.TP.EA.rda")

genesNT <- read.table("inc/genes.DEA.NT.lst")
genesNT <- as.factor(genesNT$V1)

genesTP <- read.table("inc/genes.DEA.TP.lst")
genesTP <- as.factor(genesTP$V1)

genesAA <- read.table("inc/genes.DEA.NT.vs.TP.lst.AA")
genesAA <- as.factor(genesAA$V1)

genesEA <- read.table("inc/genes.DEA.NT.vs.TP.lst.EA")
genesEA <- as.factor(genesEA$V1)

# ==== ui.R ==== 

ui <- navbarPage(theme=shinytheme("paper"),"Gene Expression Analysis",
  
  fluidRow(
    column(12,
      h6("Breast cancer is one of the most common cancers with more than 2,300,000 cases and 684,996 deaths each year worldwide. In the United States it is the fourth leading cause of death from cancer. (INCA, 2021)")
    ),
    column(12,
      h6("The objective of this study is to investigate the relationship between the ancestry of African American (AA) and European-American (EA) patients with breast cancer and the possible related clinical and genetic changes.")
    ),
    column(12,
      h6("Clinical and genomic data of 1085 patients from the Cancer Genome Atlas project (TCGA-BRCA) were available on the Genomic Data Commons portal (portal.gdc.cancer.gov).")
    ),
    column(12,
      h6("Ancestry analysis and classification of these same patients were obtained from the TCGAA (52.25.87.215/TCGAA).")
    ),
    column(12,
      h6("Next, in each tab, you can check the differential expression analysis on a Volcano Plot and the list of genes differentially expressed in a table for the following categories:")),
    column(11, h6(strong(">> African American: Normal Samples x Tumor Samples"),  br())),
    column(11, h6(strong(">> European American: Normal Samples x Tumor Samples"),  br())),
    column(11, h6(strong(">> Tumor Samples: African American x European American"),  br()))
    ),
  
  tabPanel("AA - NT x TP",
      sidebarPanel(
        fluidRow(
          h6("Volcano Plot is a type of plot that shows statistical significance (P value) versus magnitude of change (fold change). It enables quick visual identification of genes that are statistically significant.")
          ),
        conditionalPanel(
          'input.tab === "Table"',
          checkboxGroupInput("vars_tbl_dea.NT.TP.AA", "Select columns to show:",
                             names(dea.NT.TP.AA2), selected = names(dea.NT.TP.AA2)[c(1:4)]) 
        ),
        width = 4),
    
      mainPanel(
        tabsetPanel(
          id = 'tab',
          tabPanel("Volcano plot",  plotOutput('volcano1'),
                   sidebarPanel(
                     fluidRow(
                       h6("198 samples analyzed"))
                     )),
          tabPanel("Table", DT::dataTableOutput("table1"))
        ),
      width = 8
      )
    ),
  
  tabPanel("EA - NT x TP",
           sidebarPanel(
             fluidRow(
               h6("Volcano Plot is a type of plot that shows statistical significance (P value) versus magnitude of change (fold change). It enables quick visual identification of genes that are statistically significant.")
             ),
             conditionalPanel(
               'input.tab === "Table',
               checkboxGroupInput("vars_tbl_dea.NT.TP.EA", "Select columns to show:",
                                  names(dea.NT.TP.EA2), selected = names(dea.NT.TP.EA2)[c(1:4)]) 
             ),
             width = 4),
           
           mainPanel(
             tabsetPanel(
               id = 'tab',
               tabPanel("Volcano plot",  plotOutput('volcano2'),
                        sidebarPanel(
                          fluidRow(
                            h6("913 samples analyzed"))
                        )),
               tabPanel("Table", DT::dataTableOutput("table2"))
             ),
             width = 8
           )
  ),
  
  tabPanel("TP - AA x EA",
           sidebarPanel(
             fluidRow(
               h6("Volcano Plot is a type of plot that shows statistical significance (P value) versus magnitude of change (fold change). It enables quick visual identification of genes that are statistically significant.")
             ),
             conditionalPanel(
               'input.tab === "Table"',
               checkboxGroupInput("vars_tbl_dea.TP", "Select columns to show:",
                                  names(dea.TP2), selected = names(dea.TP2)[c(1:4)]) 
             ),
             width = 4),
           
           mainPanel(
             tabsetPanel(
               id = 'tab',
               tabPanel("Volcano plot",  plotOutput('volcano3'),
                        sidebarPanel(
                          fluidRow(
                            h6("989 samples analyzed"))
                        )),
               tabPanel("Table", DT::dataTableOutput("table3"))
             ),
             width = 8
           )
  )
  
)

# ==== server.R ====

server <- function(input, output) {
  
  cutoffAA_down <- sort(dea.NT.TP.AA[dea.NT.TP.AA$logFC < -3 , "FDR"])[20]
  cutoffAA_up <- sort(dea.NT.TP.AA[dea.NT.TP.AA$logFC > 3 , "FDR"])[20]
  
  output$volcano1 <- renderPlot({
    if (is.null(data())) {return(NULL)}
    ggplot(data = dea.NT.TP.AA, aes(x = logFC, y = -log10(FDR)))+
             geom_point(aes(colour=FDR < 0.01), pch=20, size=2) +
             geom_vline(xintercept=c(-3,3), linetype="dotted") +
             geom_hline(yintercept=c(-log10(0.01)), linetype="dotted") +
             ggtitle("Differential gene expression of Normal vs Tumor samples of African American ancestry") +
             xlab("Gene expression change\n log2(FC)") + 
             ylab("Significance\n -log10(FDR)") +
             xlim(c(-10,10)) +
             theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
             geom_label_repel(data = dea.NT.TP.AA[dea.NT.TP.AA$logFC < -3 & dea.NT.TP.AA$FDR < cutoffAA_down, ], 
                              aes(label=symbol),
                              seed = 123,
                              xlim = c(NA, -5),  # down regulated genes labels before -5 logFC
                              nudge_x = -5,
                              hjust = 1,
                              direction = "y",
                              max.iter = 1e5,
                              force = 0.5,
                              size = 3,
                              max.overlaps = Inf) +          
             geom_label_repel(data = dea.NT.TP.AA[dea.NT.TP.AA$logFC > 3 & dea.NT.TP.AA$FDR < cutoffAA_up, ],
                              aes(label=symbol),
                              seed = 123,
                              xlim = c(5, NA),  # up regulated genes labels after 5 logFC
                              nudge_x = 6,
                              hjust = 0,
                              direction = "y",
                              max.iter = 1e5,
                              force = 0.5,
                              size = 3,
                              max.overlaps = Inf)
  })
  
  dea.NT.TP.AA2 = as.data.frame(dea.NT.TP.AA[dea.NT.TP.AA$symbol %in% genesAA, ])
  dea.NT.TP.AA2 <- dea.NT.TP.AA2 %>%
    select(symbol, description, logFC, FDR, pvalue, ensgene, baseMean, lfcSE, stat)
  row.names(dea.NT.TP.AA2) <- NULL
  
  dea.NT.TP.AA2$logFC = round(dea.NT.TP.AA2$logFC, 5)
  dea.NT.TP.AA2$FDR = round(dea.NT.TP.AA2$FDR, 5)
  dea.NT.TP.AA2$pvalue = round(dea.NT.TP.AA2$pvalue, 5)
  dea.NT.TP.AA2$baseMean = round(dea.NT.TP.AA2$baseMean, 5)
  dea.NT.TP.AA2$lfcSE = round(dea.NT.TP.AA2$lfcSE, 5)
  dea.NT.TP.AA2$stat = round(dea.NT.TP.AA2$stat, 5)
  
  dea.NT.TP.AA3 = dea.NT.TP.AA2[, 1:9]
  output$table1 = renderDT(
    dea.NT.TP.AA3[, input$vars_tbl_dea.NT.TP.AA, drop = FALSE] , 
    server = TRUE)
  
  cutoffEA_down <- sort(dea.NT.TP.EA[dea.NT.TP.EA$logFC < -3 , "FDR"])[20]
  cutoffEA_up <- sort(dea.NT.TP.EA[dea.NT.TP.EA$logFC > 3 , "FDR"])[20]
  
  output$volcano2 <- renderPlot({
    if (is.null(data())) {return(NULL)}
    ggplot(data = dea.NT.TP.EA, aes(x = logFC, y = -log10(FDR))) + 
      geom_point(aes(colour=FDR < 0.01), pch=20, size=2) +
      geom_vline(xintercept=c(-3,3), linetype="dotted") +
      geom_hline(yintercept=c(-log10(0.01)), linetype="dotted") +
      ggtitle("Differential gene expression of Normal vs Tumor samples of European American ancestry") +
      xlab("Gene expression change\n log2(FC)") + 
      ylab("Significance\n -log10(FDR)") +
      xlim(c(-10,10)) +
      ylim(c(-10,350)) +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
      geom_label_repel(data = dea.NT.TP.EA[dea.NT.TP.EA$logFC < -3 & dea.NT.TP.EA$FDR < cutoffEA_down, ], 
                       aes(label=symbol),
                       seed = 123,
                       xlim = c(NA, -5),  # down regulated genes labels before -5 logFC
                       nudge_x = -5,
                       hjust = 1,
                       direction = "y",
                       max.iter = 1e5,
                       force = 0.5,
                       size = 3,
                       max.overlaps = Inf) +          
      geom_label_repel(data = dea.NT.TP.EA[dea.NT.TP.EA$logFC > 3 & dea.NT.TP.EA$FDR < cutoffEA_up, ],
                       aes(label=symbol),
                       seed = 123,
                       xlim = c(5, NA),  # up regulated genes labels after 5 logFC
                       nudge_x = 6,
                       hjust = 0,
                       direction = "y",
                       max.iter = 1e5,
                       force = 0.5,
                       size = 3,
                       max.overlaps = Inf)
  })
  
  dea.NT.TP.EA2 = as.data.frame(dea.NT.TP.EA[dea.NT.TP.EA$symbol %in% genesEA, ])
  dea.NT.TP.EA2 <- dea.NT.TP.EA2 %>%
    select(symbol, description, logFC, FDR, pvalue, ensgene, baseMean, lfcSE, stat)
  row.names(dea.NT.TP.EA2) <- NULL
  
  dea.NT.TP.EA2$logFC = round(dea.NT.TP.EA2$logFC, 5)
  dea.NT.TP.EA2$FDR = round(dea.NT.TP.EA2$FDR, 5)
  dea.NT.TP.EA2$pvalue = round(dea.NT.TP.EA2$pvalue, 5)
  dea.NT.TP.EA2$baseMean = round(dea.NT.TP.EA2$baseMean, 5)
  dea.NT.TP.EA2$lfcSE = round(dea.NT.TP.EA2$lfcSE, 5)
  dea.NT.TP.EA2$stat = round(dea.NT.TP.EA2$stat, 5)
  
  dea.NT.TP.EA3 = dea.NT.TP.EA2[, 1:9]
  output$table2 = renderDT(
    dea.NT.TP.EA3[, input$vars_tbl_dea.NT.TP.EA, drop = FALSE] , 
    server = TRUE)
  
  cutoffTP_down <- sort(dea.TP[dea.TP$logFC < -3 , "FDR"])[20]
  cutoffTP_up <- sort(dea.TP[dea.TP$logFC > 3 , "FDR"])[20]
  
  output$volcano3 <- renderPlot({
    if (is.null(data())) {return(NULL)}
    ggplot(data = dea.TP, aes(x = logFC, y = -log10(FDR))) + 
      geom_point(aes(colour=FDR < 0.01), pch=20, size=2) +
      geom_vline(xintercept=c(-3,3), linetype="dotted") +
      geom_hline(yintercept=c(-log10(0.01)), linetype="dotted") +
      ggtitle("Differential gene expression of African American vs European American ancestry of Tumor samples") +
      xlab("Gene expression change\n log2(FC)") + 
      ylab("Significance\n -log10(FDR)") +
      xlim(c(-10,10)) +
      ylim(c(-10,350)) +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
      geom_label_repel(data = dea.TP[dea.TP$logFC < -3 & dea.TP$FDR < cutoffTP_down, ], 
                       aes(label=symbol),
                       seed = 123,
                       xlim = c(NA, -5),  # down regulated genes labels before -5 logFC
                       nudge_x = -5,
                       hjust = 1,
                       direction = "y",
                       max.iter = 1e5,
                       force = 0.5,
                       size = 3,
                       max.overlaps = Inf) +          
      geom_label_repel(data = dea.TP[dea.TP$logFC > 3 & dea.TP$FDR < cutoffTP_up, ],
                       aes(label=symbol),
                       seed = 123,
                       xlim = c(5, NA),  # up regulated genes labels after 5 logFC
                       nudge_x = 6,
                       hjust = 0,
                       direction = "y",
                       max.iter = 1e5,
                       force = 0.5,
                       size = 3,
                       max.overlaps = Inf)
  })
  
  dea.TP2 = as.data.frame(dea.TP[dea.TP$symbol %in% genesTP, ])
  dea.TP2 <- dea.TP2 %>%
    select(symbol, description, logFC, FDR, pvalue, ensgene, baseMean, lfcSE, stat)
  row.names(dea.TP2) <- NULL
  
  dea.TP2$logFC = round(dea.TP2$logFC, 5)
  dea.TP2$FDR = round(dea.TP2$FDR, 5)
  dea.TP2$pvalue = round(dea.TP2$pvalue, 5)
  dea.TP2$baseMean = round(dea.TP2$baseMean, 5)
  dea.TP2$lfcSE = round(dea.TP2$lfcSE, 5)
  dea.TP2$stat = round(dea.TP2$stat, 5)
  
  dea.TP3 = dea.TP2[, 1:9]
  output$table3 = renderDT(
    dea.TP3[, input$vars_tbl_dea.TP, drop = FALSE] , 
    server = TRUE)
  
}

# Run the application 
shinyApp(ui = ui, server = server)