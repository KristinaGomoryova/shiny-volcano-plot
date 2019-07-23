library(shiny)
library(plotly)
library(tidyverse)

ui <- fluidPage(
  titlePanel("Interactive shiny volcano plot"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("def_logFC", 
                  label = "Log Fold Change:",
                  min = 0, max = 5, value = 1),
      sliderInput("def_adj_pval", 
                  label = "Adjusted pvalue:",
                  min = 0, max = 2, value = 0.05, round = FALSE, step = 0.05
                  ), width = 2),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Plot", 
                           plotlyOutput("VolcanoPlot"),
                           dataTableOutput("selectedProteinsTable"),
                           downloadButton("downloadData", "Download selected proteins")
                           ),
                  tabPanel("Table", 
                           DT::dataTableOutput("allProteinsTable"),
                           downloadButton("downloadDataAll", "Download table"),
                           plotOutput("static_plot"))
      )
    )
  )
)

server <- function(input, output) {
  
  differentialExpressionResults <- read.delim("input_table.txt", stringsAsFactors = FALSE) 
  
  output$VolcanoPlot <- renderPlotly({
    
  differentialExpressionResults["group"] <- "NS" 
  differentialExpressionResults[which(differentialExpressionResults['adj_pval']< input$def_adj_pval & abs(differentialExpressionResults["diff"]) < input$def_logFC), "group"] <- paste("adj.pval < ", input$def_adj_pval)
  differentialExpressionResults[which(differentialExpressionResults['adj_pval']> input$def_adj_pval & abs(differentialExpressionResults["diff"]) > input$def_logFC), "group"] <- paste("FC > ", input$def_logFC) 
  differentialExpressionResults[which(differentialExpressionResults['adj_pval']< input$def_adj_pval & abs(differentialExpressionResults["diff"]) > input$def_logFC), "group"] <- paste("adj.pval < ", input$def_adj_pval, " & FC >", input$def_logFC) 
  
  differentialExpressionResults["minusLog10Pvalue"] = -log10(differentialExpressionResults$adj_pval)
  differentialExpressionResults["tooltip"] = differentialExpressionResults$name
  
    plot <- differentialExpressionResults %>%
      ggplot(aes(x = diff,
                 y = minusLog10Pvalue,
                 colour = group,
                 text = tooltip,
                 key = row.names(differentialExpressionResults))) +
      geom_point() +
      xlab("log fold change") +
      ylab("-log10(adj p-value)")
    
    plot %>%
      ggplotly(tooltip = "tooltip") %>%
      layout(dragmode = "select")
  })
  
  selprots <- reactive({
    eventData <- event_data("plotly_selected")
    
    selectedData <- differentialExpressionResults %>% slice(0)
    if (!is.null(eventData)) selectedData <- differentialExpressionResults[eventData$key,]
    
    selectedData %>%
      transmute(
        protein = name,
        `log fold change` = signif(diff, digits = 2),
        `p-value` = signif(adj_pval, digits = 2)
      )
  })
  
  output$selectedProteinsTable <- renderDataTable({
    selprots() 
  },
  options = list(dom = "tip", pageLength = 10, searching = FALSE)
  )
  
  table_full <- reactive({
    differentialExpressionResults["group"] <- "NS" 
    differentialExpressionResults[which(differentialExpressionResults['adj_pval']< input$def_adj_pval & abs(differentialExpressionResults["diff"]) < input$def_logFC), "group"] <- paste("adj.pval < ", input$def_adj_pval) 
    differentialExpressionResults[which(differentialExpressionResults['adj_pval']> input$def_adj_pval & abs(differentialExpressionResults["diff"]) > input$def_logFC), "group"] <- paste("FC > ", input$def_logFC)
    differentialExpressionResults[which(differentialExpressionResults['adj_pval']< input$def_adj_pval & abs(differentialExpressionResults["diff"]) > input$def_logFC), "group"] <- paste("adj.pval < ", input$def_adj_pval, " & FC >", input$def_logFC) 
    differentialExpressionResults
  })
  
  output$allProteinsTable <- DT::renderDataTable({
    help_full <- table_full()
    help_full['unlist_names'] <- vapply(strsplit(help_full$name,";"), `[`, 1, FUN.VALUE=character(1))
    help_full['Publications'] <- paste("https://www.uniprot.org/uniprot/",help_full$unlist_names,"/publications", sep="")
    help_full['Publications'] <- paste0("<a href='",help_full$Publications,"'>",help_full$Publications,"</a>")
    help_full
  }, escape = FALSE)
  
  output$static_plot <- renderPlot({
    s = input$allProteinsTable_rows_selected
    help <- as.data.frame(differentialExpressionResults[,c("diff","adj_pval")])
    plot(x=help$diff, y= -log10(help$adj_pval), cex = 1, pch=16, col = "black", xlab = "logFC", ylab = "-log10(adj.pval)")
    if (length(s)) points(help[s, , drop = FALSE], pch = 19, cex = 1.4, col = "blue")
  }
  )
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("dataset", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(selprots(), file, row.names = FALSE)
    }
  )
  
  
  output$downloadDataAll <- downloadHandler(
    filename = function() {
      paste("dataset", ".csv", sep = "")
    },
    content = function(file) {
    write.csv(table_full, file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server, options = list(height = 600))
