library(shiny)
library(plotly)
library(data.table)
library(DT)
library(dplyr)
library(shinyjs)
library(IsoformSwitchAnalyzeR)

coordintae_transformer <- function(grid_size, metagene_num) {
  y <- metagene_num %/% grid_size
  x <- metagene_num %% grid_size
  
  zeros <- which(x == 0)
  non_zeros <- which(x != 0)
  
  x[non_zeros] <- x[non_zeros] - 1
  x[zeros] <- grid_size - 1
  
  y[which(x == grid_size - 1)] <- y[which(x == grid_size - 1)] - 1
  
  return(cbind(x = x, y = y)) 
}

spot_mapper <- function(metagenes, gene_som_metagene_data, isoform_som_metagene_data, isoform_to_gene, heatmap_colors) {
  
  gene_list <- substr(unique(isoform_to_gene[names(isoform_som_metagene_data)[which(isoform_som_metagene_data %in% metagenes)]]), 1, 15)
  
  matched_metagenes <- gene_som_metagene_data[gene_list]
  
  matched_metagenes <- matched_metagenes[which(!is.na(matched_metagenes))]
  
  n.map <- matrix(0, 40, 40)
  
  n.map[as.numeric(names(table(unname(matched_metagenes))))] <- table(unname(matched_metagenes))
  n.map[which(n.map==0)] <- NA
  n.map <- matrix(n.map, 40)
  
  plot_ly(z = t(n.map), colors = heatmap_colors, type = "heatmap")
  
}

correaltion_cluster_enricher <- function(isoform_summary_data, selected_metagenes, corr_threshold, som_gene_list, opossom_genesets) {
  
  isoform_summary_data <- isoform_summary_data[which(isoform_summary_data$max_cor >= corr_threshold),]
  
  isoform_summary_data <- isoform_summary_data[which(isoform_summary_data$target_metagene %in% selected_metagenes),]
  
  cor_isoforms_gensets <- oposSOM:::GeneSet.Fisher(substr(isoform_summary_data$gene, 1, 15), som_gene_list, opossom_genesets, sort = T)
  
  genset_data <- data.frame(geneset = names(cor_isoforms_gensets),
                            pval = unname(cor_isoforms_gensets))
  
  return(genset_data)
}

server <- shinyServer(function(input, output, session) {
  
  load("new_app_data.RData")
  ## load("/home/ubuntu/shiny_server_apps/switch_app/app_data.RData")
  ## load("/home/ubuntu/shiny_server_apps/switch_app/shiny_switch_data.RData")
  
  #### reactive values ####
  v <- reactiveValues(overview_map = overview_map, gene_search_error = NULL, max_cor_isoforms = NULL, isoform_summary_data = NULL, filtered_isoforms = NULL,
                      transcript_annotation = NULL, fraction_barplot = NULL, annotated_map = overview_map, 
                      correlation_quantities_map = correlation_quantities_map, input_gene = "NRAS", 
                      switch_data_frame = initial_switch_table, switch_list = fully_analyzid_switch_lists$Nv_vs_Mel,
                      enrichment_table = NULL,
                      gene_max_cor_isoforms = NULL, metagen_max_cor_isoforms = NULL, mapped_heatmap = isoform_overview_map)
  
  #### gene selection observers ####
  observeEvent(input$selected_gene, {
    v$input_gene <- input$selected_gene
  })
  
  observeEvent(input$gene_search, {
    v$input_gene <- toupper(input$input_gene)
  })
  
  #### gene isoform observer ####
  observeEvent(c(v$input_gene, input$treshold, input$sample_size, input$correlation, input$source_data_type),{
    
    if(TRUE) {
      v$filtered_isoforms <-  rownames(fraction_derived_isoform_normalized_counts)[which(rowSums2(fraction_derived_isoform_normalized_counts > input$treshold) >= input$sample_size)]
      v$isoform_summary_data <- isoform_summary_data_fraction_derived
      v$max_cor_isoforms <- max_cor_isoforms_fraction_derived
    } else {
      v$filtered_isoforms <-  rownames(normalized_isoform_counts)[which(rowSums2(normalized_isoform_counts > input$treshold) >= input$sample_size)]
      v$isoform_summary_data <- isoform_summary_data_normalized_isoform_counts_derived
      v$max_cor_isoforms <- max_cor_isoforms_normalized
    } 
     
    
    ### low expression isoform filtering
    v$isoform_summary_data <- v$isoform_summary_data[which(v$isoform_summary_data$isoform %in% v$filtered_isoforms),]
    
    if(v$input_gene %in% v$isoform_summary_data$symbol) {
      
      
      v$gene_search_error <- NULL
      
      selected_gene_data <- v$isoform_summary_data[v$isoform_summary_data$symbol == v$input_gene,]
      
      transcript_annotation <- cbind(selected_gene_data[,c("isoform", "transcript_type", "is_functional")], selected_gene_data$max_cor >= input$correlation, selected_gene_data$max_cor < input$correlation)
      
      colnames(transcript_annotation)[4:5] <- c("correlated", "non_correlated")
      
      transcript_annotation$correlated[which(is.na(transcript_annotation$correlated))] <- FALSE
      transcript_annotation$non_correlated[which(is.na(transcript_annotation$non_correlated))] <- TRUE
      
      v$gene_max_cor_isoforms <- v$max_cor_isoforms[as.character(transcript_annotation$isoform)]
      
      total <- c(nrow(transcript_annotation),
                 "",
                 sum(transcript_annotation$Is_functional),
                 sum(transcript_annotation$correlated),
                 sum(transcript_annotation$non_correlated))
      
      v$transcript_annotation <- rbind(total, apply(transcript_annotation, 2, as.character))
      
      
      if(nrow(selected_gene_data) > 0) {
        
        selected_gene_data$isoform <- factor(selected_gene_data$isoform, levels = selected_gene_data[["isoform"]])
        
        v$fraction_barplot <- plot_ly(x = selected_gene_data$isoform, y = selected_gene_data$mean_fraction, type = "bar", name = 'Mean fraction', source = "barplot") %>%
                # add_trace(y = selected_gene_data$fraction_sd^2, name = "Fraction var") %>%
                add_trace(y = selected_gene_data$max_cor, name = "Max cor") %>%
                layout(yaxis = list(title = 'Fraction, correlation'), barmode = 'group')
      }
      
      gene_metagene <- unique(selected_gene_data$source_metagene)
      
      selected_gene_data <- selected_gene_data[which(selected_gene_data$max_cor >= input$correlation),]
      
      if(nrow(selected_gene_data) > 0) {
        aggregated <- aggregate(selected_gene_data$symbol_vis, by = list(selected_gene_data$target_metagene), FUN = function(x) {paste0(x, collapse = "\n")})
        
        annotation_text <- c(v$input_gene, aggregated$x)
        
        annotation_coords <- coordintae_transformer(grid_size = 40, metagene_num = c(gene_metagene, aggregated$Group.1))
        
        v$annotated_map <- plot_ly(z = t(matrix(overview_metadata, 40, 40)), colors = heatmap_colors, type = "heatmap") %>% 
          add_markers(x = annotation_coords[-1,1], y = annotation_coords[-1,2], text = annotation_text[-1],
                      inherit = F, name = unique(selected_gene_data$symbol), 
                      marker = list(color = ~"white", line=list(width=1, color='black'))) %>%
          add_annotations(text = annotation_text[1], x = annotation_coords[1,1], y = annotation_coords[1,2], 
                          bgcolor = "#b3d9ff", color = "black", arrowsize = 1, opacity = 0.7)
      } else {
        v$annotated_map <- plot_ly(z = t(matrix(overview_metadata, 40, 40)), colors = heatmap_colors, type = "heatmap") %>% 
          add_annotations(text = v$input_gene, x = coordintae_transformer(grid_size = 40, metagene_num = gene_metagene)[,1], y = coordintae_transformer(grid_size = 40, metagene_num = gene_metagene)[,2], 
                          bgcolor = "#b3d9ff", color = "black", arrowsize = 1, opacity = 1)
      }
      
         
      
    } else {
      v$transcript_annotation <- NULL
      
      v$gene_search_error <- "Please insert correct gene name"
      
    }
    
  })
  
  #### metagene unique correlation heatmam observer ####
  observeEvent(c(input$treshold, input$sample_size, input$correlation), {
    metagene_unique_corr <- sapply(1:1600, function(x) {sum(v$isoform_summary_data[which(v$isoform_summary_data$max_cor >= input$correlation), "target_metagene"] == x)})
    
    metagene_unique_corr <- log(metagene_unique_corr)
    
    metagene_unique_corr[is.infinite(metagene_unique_corr)] <- 0
    
    v$enrichment_table <- NULL
    
    # v$correlation_quantities_map <- plot_ly(z = t(matrix(metagene_unique_corr, 40, 40)), colors = heatmap_colors, type = "heatmap", source = "erichment_metagenes") %>%
    #   add_markers(x = coordintae_transformer(grid_size = 40, metagene_num = 1:1600)[,1], y = coordintae_transformer(grid_size = 40, metagene_num = 1:1600)[,2], text = "",
    #               inherit = F, name = "", 
    #               marker = list(color = ~"#ffffff00", line=list(width=1, color='black')))
  })
  
  #### gene searching error ####
  output$error_text <- renderText({
    paste("<font color=\"#ff0000\"><b>", v$gene_search_error, "</b></font>")
  })
  
  #### isofrom info data table ####
  output$isoform_table <- DT::renderDataTable({datatable(
    v$transcript_annotation,
    class   = 'cell-border compact hover',
    fillContainer = T,
    escape = F,
    selection = list(mode = 'single'),
    options = list(dom = 'rtp',
                   pageLength = 50,
                   searchHighlight = TRUE,
                   initComplete = JS("function(settings, json) {",
                                     "$(this.api().table().header()).css({'background-color': '#3474B7', 'color': '#fff'});",
                                     "}"),
                   scrollY = 500,scrollX = TRUE
    ),
    style = 'bootstrap', editable = FALSE)
  })
  
  #### gene mapping fractional barplot ####
  output$fraction_plot <- renderPlotly({
    v$fraction_barplot
  })
  
  #### gene isoforms mapping map ####
  output$searchable_map <- renderPlotly({
    hover_isoform <- event_data("plotly_hover", source = "barplot")[,"x"]
    
    # input$isoform_table_rows_selected
    
    coords <- coordintae_transformer(grid_size = 40, metagene_num = v$gene_max_cor_isoforms[hover_isoform])
    if(!is.null(hover_isoform)) {
      add_annotations(v$annotated_map, text = hover_isoform, x = coords[,1], y = coords[,2], bgcolor = "#99cc00")
    } else {
      v$annotated_map
    }
    
  })
  
  #### metagene click observer ####
  observeEvent(c(event_data("plotly_click", source = "click_map"), input$correlation, input$source_data_type, input$treshold, input$sample_size), {
    
    # ifelse(input$source_data_type == "Fraction_derived", v$isoform_summary_data <- isoform_summary_data_fraction_derived, v$isoform_summary_data <- isoform_summary_data_normalized_isoform_counts_derived)
    
    click_metagene <- as.integer(event_data("plotly_click", source = "click_map")[,c("x", "y")] + 1)
    
    if(length(click_metagene) > 0) {
      
      metagene_num <- (click_metagene[1]) + (click_metagene[2] - 1)*40
      
      clicked_metagene_data <- v$isoform_summary_data[which(v$isoform_summary_data$source_metagene == metagene_num),]
      
      clicked_metagene_data <- clicked_metagene_data[which(clicked_metagene_data$max_cor > input$correlation),]
      
      v$metagen_max_cor_isoforms <- v$max_cor_isoforms[as.character(clicked_metagene_data$isoform)]
      
      if(nrow(clicked_metagene_data) > 0) {
        
        # View(clicked_metagene_data)
        
        aggregated <- aggregate(clicked_metagene_data$symbol_vis, by = list(clicked_metagene_data$target_metagene), FUN = function(x) {paste0(x, collapse = "\n")})
        
        
        coords <- coordintae_transformer(grid_size = 40, metagene_num = c(aggregated$Group.1, metagene_num))
        
        col <- c(rep("#FFFFFF", length(aggregated$x)), "black")
        
        v$click_heatmap <- plot_ly(z = t(matrix(overview_metadata, 40, 40)), colors = heatmap_colors, type = "heatmap",  source = "click_map") %>% 
          add_markers(x = coords[,1], y = coords[,2], text = c(aggregated$x, paste0("Selected ", metagene_num)),
                      inherit = F, name = metagene_num, 
                      marker = list(color = ~col, line=list(width=1, color='black'))) %>% event_register("plotly_click")
        
        
        clicked_metagene_data$symbol_vis <- factor(clicked_metagene_data$symbol_vis, levels = clicked_metagene_data[["symbol_vis"]])
        
        v$isoforms_barplot <- plot_ly(x = clicked_metagene_data$symbol_vis, y = clicked_metagene_data$mean_fraction, type = "bar", name = 'Mean fraction', source = "metagene_barplot") %>%
          # add_trace(y = clicked_metagene_data$fraction_sd^2, name = "Fraction var") %>%
          add_trace(y = clicked_metagene_data$max_cor, name = "Max cor") %>%
          layout(yaxis = list(title = 'Fraction, correlation'), barmode = 'group')
        
        
        # add_annotations(text = click_metagene_corr_isoforms, x = coords[,1], y = coords[,2], bgcolor = "#99cc00")
      } else {
        v$isoforms_barplot <- NULL
        v$click_heatmap <- plot_ly(z = t(matrix(overview_metadata, 40, 40)), colors = heatmap_colors, type = "heatmap",  source = "click_map") %>% 
          add_markers(x = coordintae_transformer(grid_size = 40, metagene_num = metagene_num)[,1], y = coordintae_transformer(grid_size = 40, metagene_num = metagene_num)[,2], text = paste0("Selected ", metagene_num), 
                      inherit = F, name = metagene_num, 
                      marker = list(color = ~"black")) %>% event_register("plotly_click")
      }
      
      } else {
        v$isoforms_barplot <- NULL
        v$click_heatmap <- plot_ly(z = t(matrix(overview_metadata, 40, 40)), colors = heatmap_colors, type = "heatmap",  source = "click_map")
      }
      
  })
  
  #### switched genes selection observers ####
  observeEvent(c(input$design_type, input$diff_cut), {
    
    if(input$design_type %in% c("M1_vs_M2", "M1_vs_N1", "M2_vs_N2",  "N1_vs_N2") ) {
      v$switch_list <- fully_analyzid_switch_lists$joined_subgroup
      
      v$switch_data_frame <- extractTopSwitches(v$switch_list, n=NA, dIFcutoff = 0.1)
      
      v$switch_data_frame <- v$switch_data_frame[which(v$switch_data_frame$condition_1 %in% condition_list[[input$design_type]][1] & v$switch_data_frame$condition_2 %in% condition_list[[input$design_type]][2]),]
      
      updateSelectizeInput(session, 'selected_switched_gene', choices = v$switch_data_frame$gene_name, selected = v$switch_data_frame$gene_name[1], server = TRUE)
    } else {
      v$switch_list <- fully_analyzid_switch_lists[[input$design_type]]
      
      v$switch_data_frame <- extractTopSwitches(v$switch_list, n=NA, dIFcutoff = 0.1)
      
      updateSelectizeInput(session, 'selected_switched_gene', choices = v$switch_data_frame$gene_name, selected = v$switch_data_frame$gene_name[1], server = TRUE)
    }
    
    # updateSliderInput(session, "diff_cut", max = round(max(abs(v$switch_list$isoformSwitchAnalysis$dIF)), digits = 2) - 0.01, value = 0.3)
  })
  
  #### isoform switch plotter ####
  output$switch_plot <- renderPlot({
    if(input$selected_switched_gene %in% v$switch_data_frame$gene_name) {
      switchPlot(
        v$switch_list,
        gene = input$selected_switched_gene,
        localTheme = theme_bw(base_size = 15),
        IFcutoff = 0.05, 
        condition1 = condition_list[[input$design_type]][1], 
        condition2 = condition_list[[input$design_type]][2]
      )
    }
  })  
  
  #### metagene click map ####
  output$metagene_click_map <- renderPlotly({
    
    hover_isoform <- event_data("plotly_hover", source = "metagene_barplot")[,"x"]
    
    if(!is.null(hover_isoform)) {
      coords <- coordintae_transformer(grid_size = 40, metagene_num = v$metagen_max_cor_isoforms[substr(hover_isoform, 1, unlist(gregexpr(pattern = "(", text = hover_isoform, fixed = T)) - 1)])
      add_annotations(v$click_heatmap, text = hover_isoform, x = coords[,1], y = coords[,2], bgcolor = "#99cc00")
    } else {
      v$click_heatmap
    }
    
  })
  
  output$isoform_quantity <- renderPlotly({
    isoform_quantity_map
  })
  
  output$isoforms_barplot <- renderPlotly({
    v$isoforms_barplot
  })
  
  # observeEvent(event_data("plotly_selected", source = "brush_map"), {
  #   selected_metagenes <- event_data("plotly_selected", source = "brush_map")$pointNumber + 1
  #   
  #   v$mapped_heatmap <- spot_mapper(metagenes = selected_metagenes, gene_som_metagene_data = gene_som_metagene_data, 
  #                                   isoform_som_metagene_data = isoform_som_metagene_data, 
  #                                   isoform_to_gene = isoform_to_gene, heatmap_colors = heatmap_colors)
  #   
  # })
  
  # observeEvent(event_data("plotly_selected", source = "erichment_metagenes"), {
  #   v$enrichment_table <- NULL
  #   selected_metagenes <- event_data("plotly_selected", source = "erichment_metagenes")$pointNumber + 1
  #   
  #   v$enrichment_table <- correaltion_cluster_enricher(isoform_summary_data = v$isoform_summary_data,
  #                                                      selected_metagenes = selected_metagenes, corr_threshold = input$correlation, 
  #                                                      som_gene_list = som_gene_list, opossom_genesets = opossom_genesets)
  #   
  # })
  
  #### enrichment table ####
  # output$enrichment_table <- DT::renderDataTable({datatable(
  #   v$enrichment_table,
  #   class   = 'cell-border compact hover',
  #   fillContainer = T,
  #   escape = F,
  #   selection = list(mode = 'single'),
  #   options = list(dom = 'rtp',
  #                  pageLength = 50,
  #                  searchHighlight = TRUE,
  #                  initComplete = JS("function(settings, json) {",
  #                                    "$(this.api().table().header()).css({'background-color': '#3474B7', 'color': '#fff'});",
  #                                    "}"),
  #                  scrollY = 400,scrollX = TRUE
  #   ),
  #   style = 'bootstrap', editable = FALSE)
  # })
  
  #### isoform level som map ####
  # output$isoform_map <- renderPlotly({
  #   isoform_overview_map %>%
  #     add_markers(x = coordintae_transformer(grid_size = 40, metagene_num = 1:1600)[,1], y = coordintae_transformer(grid_size = 40, metagene_num = 1:1600)[,2], text = "",
  #                 inherit = F, name = "", 
  #                 marker = list(color = ~"#ffffff00", line=list(width=1, color='black')))
  # })

  #### mapping cluster between isoform and gene som ####
  # output$mappings <- renderPlotly({
  #   v$mapped_heatmap
  # })

  #### gene overview map ####
  # output$gene_map <- renderPlotly({
  #   overview_map
  # })
  
  #### correlation quantity map ####
  # output$correlation_quantities_map <- renderPlotly({
  #   v$correlation_quantities_map
  # })
  
  output$paper_map <- renderImage({
    return(list(
      src = "~/shiny_server_apps/heatmap/paper_map.png",
      contentType = "image/png",
      alt = "map", width = 450,
      height = 300
    ))
  }, deleteFile = F)
})
