#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# ==== Loading library ===============================================================

library(shiny)
library(ggplot2)
library(DT)
library(GeneBook)
library(shinythemes)

# ==== Global variables ===============================================================

gn <- read.csv(paste(getwd(),'/inc/genes_names.csv',sep = ''))
g1 <- read.csv(paste(getwd(),"/inc/g1.csv",sep = ''))
g2 <- read.csv(paste(getwd(),"/inc/g2.csv",sep = ''))
g3 <- read.csv(paste(getwd(),"/inc/g3.csv",sep = ''))

# ==== ui.R ===============================================================

ui <- fluidPage(shinythemes::themeSelector(),
   
  # Application title
  titlePanel("Gene Expression"),
  
  # Sidebar 
  sidebarLayout(
    # Sidebar Panel
    sidebarPanel(
      selectInput(
        inputId ="G_groups",
        label = "A- Choose Group to plot:",
        choices = c("1- Genes down regulated and
                    up regulated " = "g1",
                    "2- Genes down regulated " = "g2",
                    "3- Genes up regulated  " = "g3")), 
      
      selectInput(
        inputId = "My_dataset",
        label = "B- Choose Gene ID to show it's full name:",
        choices = levels(gn$GeneID)),
      
      selectInput(
        inputId = "More_info",
        label = "C- Documentation:",
        choices = c('Introduction', 'Information', 'Help',
                    'References', 'Table-1','Table-2','Table-3'),
        selected = "Information"),
    
      submitButton(
           text = "Apply Changes",
           icon = icon("sliders"), 
           width ='200px')
    ),
    
    # Main Panel 
    mainPanel(
      downloadButton(
        outputId = "downloadData",
        label = "Download Data"),
      
      plotOutput(
        outputId = "myplot",
        width = "100%",
        height = "400px"),
      
      verbatimTextOutput(
        outputId = "odataset"),
      
      uiOutput(
        outputId = "odataset_link"),
      
      uiOutput(
        outputId = "text1")
      
    )
  )
)

# ==== server.R ===============================================================
server <- function(input, output) {
   
  output$odataset <- renderPrint({
    paste(input$My_dataset," = ", gn$Gene[gn$GeneID==input$My_dataset])
  })
  
  abbreviation <- reactive((GeneCard_ID_Convert(input$My_dataset)))
  
  # output for the odataset_link
  output$odataset_link <- renderPrint({
    tags$a(
      href = paste(
        "https://www.genecards.org/cgi-bin/carddisp.pl?gene=",
        as.character(abbreviation()[1]),
        sep = ''
      ),
      as.character(abbreviation()[1])
    )
  })
  
  
  full_file_name <-reactive(paste("./inc/", input$G_groups, ".csv", sep = ""))
  
  output$downloadData <- downloadHandler(
    
    filename = full_file_name,
    
    content = function(file){
      write.csv(read.csv(full_file_name()), quote = FALSE,file)
    } )
  
  output$myplot = renderPlot({
    g_x <- read.csv(full_file_name())
    
    p <- ggplot(g_x, aes(x=Gene_ID, y=log(Relative_expression_levels),
                         fill=Resistant_or_Susceptible_strains)) +
      
      geom_bar(stat="identity", position=position_dodge()) +
      geom_errorbar(aes(ymin=log(Relative_expression_levels)-(SD/10),
                        ymax=log(Relative_expression_levels)+(SD/10)),width=.3,
                    position=position_dodge(.9))
    p + scale_fill_brewer(palette="Paired")+
      ggtitle(paste("Relative expression levels of candidate gene list","\n",
                    "expressed as mean fold difference between pre- and",
                    "\n", "post-infection ± standard deviation (SD) ")) +
      guides(fill=guide_legend(title=NULL))
    
    p$theme <- theme(axis.text.x = element_text(angle = 90, hjust = 1))
    p$labels$x <- "Gene ID"
    p$labels$y <- "Log (base 10) Relative Expression Levels"
    p$labels$fill <- NULL
    
    return(p)
    
  })
  
  
  # renderDT() from DT library is a replacement for Shiny renderDataTable()
  output$datatable1 <- renderDT(datatable(g1))
  output$datatable2 <- renderDT(datatable(g2))
  output$datatable3 <- renderDT(datatable(g3))
  
  output$text1 <- renderUI({
    if(input$More_info=="Introduction"){
      includeHTML("inc/introduction.html")
    } else if(input$More_info=="Information"){
      includeHTML("inc/information.html")
    } else if(input$More_info=="Help"){
      includeHTML("inc/help.html")
    } else if(input$More_info=="Table-1"){
      DTOutput('datatable1')
    } else if(input$More_info=="Table-2"){
      DTOutput('datatable2')
    } else if(input$More_info=="Table-3"){
      DTOutput('datatable3')
    } else if(input$More_info=="References"){
      includeHTML("inc/references.html")
    }
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)

