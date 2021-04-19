# ==== Loading library ==== 

library(shiny)
library(plotly)
library(shinythemes)
library(ggplot2)
library(DT)

# ==== Global variables ====

genesAA <- read.table("C:/Users/raque/Downloads/2020.5/Shiny-and-DEA---Estagio/inc/genes.DEA.NT.vs.TP.lst.AA")
#genesEA
#genesNT
#genesTP

# ==== ui.R ==== 

ui <- navbarPage(theme=shinytheme("superhero"),"Gene Expression Analysis",
  
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
          choices = levels(genesAA$V1)),
      width = 3),
    
      mainPanel(
        tabsetPanel(
          id = 'tab',
          tabPanel("Volcano plot",  plotlyOutput('volcano1') ),
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
          choices = levels(genesAA$V1)),
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
          choices = levels(genesAA$V1)),
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
         choices = levels(genesAA$V1)),
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
  
  output$table1 <- renderDT(datatable(genesAA))
  
}
  
  # output$odataset <- renderPrint({
  #   paste(input$My_dataset," = ", gn$Gene[gn$GeneID==input$My_dataset])
  # })
  # 
  # abbreviation <- reactive((GeneCard_ID_Convert(input$My_dataset)))
  # 
  # # output for the odataset_link
  # output$odataset_link <- renderPrint({
  #   tags$a(
  #     href = paste(
  #       "https://www.genecards.org/cgi-bin/carddisp.pl?gene=",
  #       as.character(abbreviation()[1]),
  #       sep = ''
  #     ),
  #     as.character(abbreviation()[1])
  #   )
  # })
  # 
  # 
  # full_file_name <-reactive(paste("./inc/", input$G_groups, ".csv", sep = ""))
  # 
  # output$downloadData <- downloadHandler(
  #   
  #   filename = full_file_name,
  #   
  #   content = function(file){
  #     write.csv(read.csv(full_file_name()), quote = FALSE,file)
  #   } )
  # 
  # output$myplot = renderPlot({
  #   g_x <- read.csv(full_file_name())
  #   
  #   p <- ggplot(g_x, aes(x=Gene_ID, y=log(Relative_expression_levels),
  #                        fill=Resistant_or_Susceptible_strains)) +
  #     
  #     geom_bar(stat="identity", position=position_dodge()) +
  #     geom_errorbar(aes(ymin=log(Relative_expression_levels)-(SD/10),
  #                       ymax=log(Relative_expression_levels)+(SD/10)),width=.3,
  #                   position=position_dodge(.9))
  #   p + scale_fill_brewer(palette="Paired")+
  #     ggtitle(paste("Relative expression levels of candidate gene list","\n",
  #                   "expressed as mean fold difference between pre- and",
  #                   "\n", "post-infection ± standard deviation (SD) ")) +
  #     guides(fill=guide_legend(title=NULL))
  #   
  #   p$theme <- theme(axis.text.x = element_text(angle = 90, hjust = 1))
  #   p$labels$x <- "Gene ID"
  #   p$labels$y <- "Log (base 10) Relative Expression Levels"
  #   p$labels$fill <- NULL
  #   
  #   return(p)
  #   
  # })
  # 
  # 
  # # renderDT() from DT library is a replacement for Shiny renderDataTable()
  # output$datatable1 <- renderDT(datatable(g1))
  # output$datatable2 <- renderDT(datatable(g2))
  # output$datatable3 <- renderDT(datatable(g3))
  # 
  # output$text1 <- renderUI({
  #   if(input$More_info=="Introduction"){
  #     includeHTML("inc/introduction.html")
  #   } else if(input$More_info=="Information"){
  #     includeHTML("inc/information.html")
  #   } else if(input$More_info=="Help"){
  #     includeHTML("inc/help.html")
  #   } else if(input$More_info=="Table-1"){
  #     DTOutput('datatable1')
  #   } else if(input$More_info=="Table-2"){
  #     DTOutput('datatable2')
  #   } else if(input$More_info=="Table-3"){
  #     DTOutput('datatable3')
  #   } else if(input$More_info=="References"){
  #     includeHTML("inc/references.html")
  #   }
  # })

# Run the application 
shinyApp(ui = ui, server = server)