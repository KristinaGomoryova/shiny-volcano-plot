library(shiny)
library(plotly)
library(tidyverse)

ui <- fluidPage(
  titlePanel("Interactive shiny volcano plot"),
  fluidRow(
    column(
      width = 7,
      plotlyOutput("volcanoPlot", height = "500px")
    ),
    column(
      width = 5,
      dataTableOutput("selectedProteinsTable")
    )
  )
)

server <- function(input, output) {
  
  differentialExpressionResults <- read.delim("input_proteins.txt", stringsAsFactors = FALSE) 
  differentialExpressionResults["group"] <- "NS" 
  differentialExpressionResults[which(differentialExpressionResults['adj_pval']< 0.05 & abs(differentialExpressionResults["diff"]) < 1), "group"] <- "p val < 0.05" 
  differentialExpressionResults[which(differentialExpressionResults['adj_pval']> 0.05 & abs(differentialExpressionResults["diff"]) > 1), "group"] <- "|FC| > 1" 
  differentialExpressionResults[which(differentialExpressionResults['adj_pval']< 0.05 & abs(differentialExpressionResults["diff"]) > 1), "group"] <- "p val < 0.05 & |FC| > 1" 
  differentialExpressionResults["minusLog10Pvalue"] = -log10(differentialExpressionResults$adj_pval)
  differentialExpressionResults["tooltip"] = differentialExpressionResults$name
  
  output$volcanoPlot <- renderPlotly({
    
    plot <- differentialExpressionResults %>%
      ggplot(aes(x = diff,
                 y = minusLog10Pvalue,
                 colour = group,
                 text = tooltip,
                 key = row.names(differentialExpressionResults))) +
      geom_point() +
      xlab("log fold change") +
      ylab("-log10(P-value)")
    
    plot %>%
      ggplotly(tooltip = "tooltip") %>%
      layout(dragmode = "select")
  })
  
  output$selectedProteinsTable <- renderDataTable({
    
    eventData <- event_data("plotly_selected")
    
    selectedData <- differentialExpressionResults %>% slice(0)
    if (!is.null(eventData)) selectedData <- differentialExpressionResults[eventData$key,]
    
    selectedData %>%
      transmute(
        protein = name,
        `log fold change` = signif(diff, digits = 2),
        `p-value` = signif(adj_pval, digits = 2)
      )
  },
  options = list(dom = "tip", pageLength = 10, searching = FALSE)
  )
}

shinyApp(ui, server, options = list(height = 600))
