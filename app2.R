# ==== Loading library ==== 

library(shiny)
library(plotly)
library(shinythemes)
library(ggplot2)
library(ggrepel)
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
  
  tabPanel("AA - TP x NT",
    sidebarLayout(
      sidebarPanel(
        selectInput(
          inputId ="G_groups",
          label = "A- Choose Group to plot:",
          choices = c("1- Genes down regulated and
                      up regulated " = "genesAA1",
                      "2- Genes down regulated " = "genesAA2",
                      "3- Genes up regulated  " = "genesAA3")), 
        selectInput(
          inputId = "My_dataset",
          label = "B- Choose Gene ID to show it's full name:",
          choices = levels(genesAA)),
      width = 3),
    
      mainPanel(
        tabsetPanel(
          id = 'tab',
          tabPanel("Volcano plot",  plotOutput('volcano1') ),
          tabPanel("Table", DT::dataTableOutput("table1"))
        ),
      width = 9)
    )),
  
  tabPanel("EA - TP x NT",
    sidebarLayout(
      sidebarPanel(
        selectInput(
          inputId ="G_groups",
          label = "A- Choose Group to plot:",
          choices = c("1- Genes down regulated and
                       up regulated " = "genesEA1",
                       "2- Genes down regulated " = "genesEA2",
                       "3- Genes up regulated  " = "genesEA3")), 
        selectInput(
          inputId = "My_dataset",
          label = "B- Choose Gene ID to show it's full name:",
          choices = levels(genesEA)),
        width = 3),
             
      mainPanel(
        tabsetPanel(
          id = 'tab',
          tabPanel("Volcano plot",  plotlyOutput('volcano2') ),
          tabPanel("Table", DT::dataTableOutput("table2"))
        ),
      width = 9)
    )),
  
  tabPanel("TP - AA x EA",
    sidebarLayout(
      sidebarPanel(
        selectInput(
          inputId ="G_groups",
          label = "A- Choose Group to plot:",
          choices = c("1- Genes down regulated and
                       up regulated " = "genesTP1",
                       "2- Genes down regulated " = "genesTP2",
                       "3- Genes up regulated  " = "genesTP3")), 
        selectInput(
          inputId = "My_dataset",
          label = "B- Choose Gene ID to show it's full name:",
          choices = levels(genesTP)),
        width = 3),
             
        mainPanel(
           tabsetPanel(
             id = 'tab',
             tabPanel("Volcano plot",  plotlyOutput('volcano3') ),
             tabPanel("Table", DT::dataTableOutput("table3"))
           ),
           width = 9)
        )),
  
  tabPanel("NT - AA x EA",
    sidebarLayout(
      sidebarPanel(
        selectInput(
          inputId ="G_groups",
          label = "A- Choose Group to plot:",
          choices = c("1- Genes down regulated and
                          up regulated " = "genesNT1",
                             "2- Genes down regulated " = "genesNT2",
                             "3- Genes up regulated  " = "genesNT3")), 
       selectInput(
         inputId = "My_dataset",
         label = "B- Choose Gene ID to show it's full name:",
         choices = levels(genesNT)),
      width = 3),
      
      mainPanel(
        tabsetPanel(
          id = 'tab',
          tabPanel("Volcano plot",  plotlyOutput('volcano4') ),
          tabPanel("Table", DT::dataTableOutput("table4"))
        ),
        width = 9)
      ))
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
  
  output$table1 <- renderDT(datatable(dea.NT.TP.AA[dea.NT.TP.AA$symbol %in% genesAA, ], 
                                      rownames = FALSE))
  
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
  
  output$table2 <- renderDT(datatable(dea.NT.TP.EA[dea.NT.TP.EA$symbol %in% genesEA, ], 
                                      rownames = FALSE))
  
  cutoffTP <- sort(dea.TP$pvalue)[10]
  shrink.deseq.cutTP <- dea.TP %>% 
    mutate(TopGeneLabel=ifelse(pvalue<=cutoffTP, symbol, ""))
  
  output$volcano3 <- renderPlot({
    if (is.null(data())) {return(NULL)}
    ggplot(shrink.deseq.cutTP, aes(x = logFC, y= -log10(FDR))) + 
    geom_point(aes(colour=FDR < 0.01), pch=20, size=2) +
    labs(x="log Fold Change", y="-log10(FDR)") + 
    geom_label_repel(aes(label=TopGeneLabel), 
                     seed = 123,
                     max.time = 3,
                     max.iter = Inf,
                     size = 3,
                     box.padding = 2, 
                     max.overlaps = Inf)
  })
  
  output$table3 <- renderDT(datatable(dea.TP[dea.TP$symbol %in% genesTP, ], 
                                      rownames = FALSE))
  
  cutoffNT <- sort(dea.NT$pvalue)[10]
  shrink.deseq.cutNT <- dea.NT %>% 
    mutate(TopGeneLabel=ifelse(pvalue<=cutoffNT, symbol, ""))
  
  output$volcano4 <- renderPlot({
    if (is.null(data())) {return(NULL)}
    ggplot(shrink.deseq.cutNT, aes(x = logFC, y= -log10(FDR))) + 
    geom_point(aes(colour=FDR < 0.01), pch=20, size=2) +
    labs(x="log Fold Change", y="-log10(FDR)") + 
    geom_label_repel(aes(label=TopGeneLabel), 
                     seed = 123,
                     max.time = 3,
                     max.iter = Inf,
                     size = 3,
                     box.padding = 2, 
                     max.overlaps = Inf)
  })
  
  output$table4 <- renderDT(datatable(dea.NT[dea.NT$symbol %in% genesNT, ], 
                                      rownames = FALSE))
  
}

# Run the application 
shinyApp(ui = ui, server = server)