library(shiny)
library(plotly)
library(data.table)
library(DT)

load("initial_switch_gene_list.RData")

ui <- shinyUI(fluidPage(
  # titlePanel("Heatmaps"),
  fluidRow( 
    div(style="display: inline-block;vertical-align:top;", sliderInput("treshold", "Exp treshold", min = 0, max = 100, value = 5)),
    div(style="display: inline-block;vertical-align:top;", sliderInput("sample_size", "Sample number higher than exp treshold", min = 0, max = 73, value = 2)),
    div(style="display: inline-block;vertical-align:top;", sliderInput("correlation", "Correlation", min = 0, max = 1, value = 0.7))
    # div(style="display: inline-block;vertical-align:top;", selectizeInput(inputId = "source_data_type", label = "Select initial data", selected = "Fraction_derived", choices = c("Fraction_derived", "DESeq_normalized")) )
    # br(),
    # imageOutput("paper_map")
  ),
  fluidRow(
    tabsetPanel(
      tabPanel(title = "Gene search", 
               div(style="display: inline-block;vertical-align:top;", selectizeInput(inputId = "selected_gene", label = "Select gene", selected = "NRAS", choices = c("NRAS", "BRAF", "TERT", "MET", "AXL", "KIT", "MITF", "AURKB", "CDK6", "ETV5", "FANCB", "PCNA", "EZH2", "E2F1", "BRCA1", "CDK4", "CXCL14", "IFNGR1", "CXCL12", "CCL21", "TGFBR2", "NOTCH3", "PDGFRB"))),
               div(style="display: inline-block;vertical-align:top; top: 8px; width: 150", textInput("input_gene", label = "Insert gene", placeholder = "Gene symbol")),
               div(style="display: inline-block;vertical-align:top; top: 8px;", actionButton("gene_search", label = "Serach", icon = icon("search", lib = "glyphicon"))),
               br(),
               div(style="display: inline-block;vertical-align:top;", plotlyOutput("searchable_map", height = "500px", width = "600px")),
               div(style="display: inline-block;vertical-align:top;", dataTableOutput("isoform_table")),
               div(style="display: inline-block;vertical-align:top;", htmlOutput("error_text")),
               plotlyOutput("fraction_plot")
      ),
      tabPanel(title = "Metagene browsing",
               plotlyOutput("metagene_click_map", height = "500px", width = "600px", inline = T),
               plotlyOutput("isoforms_barplot")
               # plotlyOutput("isoform_quantity", height = "500px", width = "600px", inline = T)
      ),
      # tabPanel(title = "Spot selection",
      #          plotlyOutput("isoform_map", height = "500px", width = "600px", inline = T),
      #          plotlyOutput("mappings", height = "500px", width = "600px", inline = T),
      #          plotlyOutput("gene_map", height = "500px", width = "600px", inline = T),
      #          div(style="display: inline-block;vertical-align:top;", plotlyOutput("correlation_quantities_map", height = "500px", width = "600px")),
      #          div(style="display: inline-block;vertical-align:top;", dataTableOutput("enrichment_table"))
      # ),
      tabPanel(title = "Isoform switch visualizer",
               div(style="display: inline-block;vertical-align:top;", selectizeInput(inputId = "design_type", label = "Comparison groups", selected = "Mel_vs_Nv", choices = c("M1_vs_M2", "M1_vs_N1", "M2_vs_N2", "Nv_vs_Mel", "N1_vs_N2", "Type1_vs_Type2"))),
               div(style="display: inline-block;vertical-align:top;", selectizeInput(inputId = "selected_switched_gene", label = "Select gene", selected = "MSR1", choices = initial_switch_gene_list)),
               # div(style="display: inline-block;vertical-align:top;", sliderInput("diff_cut", "Diff cut off", min = 0, max = 0.57, value = 0.3)),
               plotOutput("switch_plot", height = "850px")
      )
    )
  )
))

