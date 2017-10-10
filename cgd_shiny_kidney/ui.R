require(shiny)
library(qtlcharts)

# Define UI 
shinyUI(pageWithSidebar(
  
  # no title
  headerPanel(""),
  
  # Sidebar with a slider input for number of observations
  sidebarPanel(
    # MGI symbol
    htmlOutput("selectTextInput"),
    # mRNA or protein?
    htmlOutput("selectLevel"),
    # type of plot
    htmlOutput("selectType"),
    # interactive
    radioButtons("interactive", "Interactive plots", c("Off"=0, "On"=1), selected = 0, 
                 inline = TRUE),
    # chromosome (for LOD and mediation plots)
    conditionalPanel(condition="input.plotType < 5 & input.interactive == 0 | input.plotType == 9",
                     htmlOutput("selectChr")),
    downloadButton('downloadRDS', 'Download RDS'),
    downloadButton('downloadCSV', 'Download CSV')  
  ),
  
  # Show a plot of the generated distribution
  mainPanel(
    conditionalPanel(condition="input.interactive == 0", plotOutput("lodPlot", height=530)),
    conditionalPanel(condition="input.interactive == 1", iplot_output("interactivePlot", height=530)),
#    plotOutput("lodPlot", height=530),
#    iplot_output("interactivePlot", height=530),
    htmlOutput("description")  
  )    
  
))
