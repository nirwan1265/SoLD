library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(dplyr)
library(vroom)
library(plotly)
library(Rtsne)
library(umap)
library(heatmaply)
library(RColorBrewer)
library(tibble)

## Read the files
lowinput <- vroom("data/lowinput_exp.csv")
control <- vroom("data/control_exp.csv")

dataset_list <- list(
  lowinput = vroom::vroom("data/lowinput_exp.csv", show_col_types = FALSE),
  control = vroom::vroom("data/control_exp.csv", show_col_types = FALSE)
)

# Function to convert dataset to z-scores
# zscore_transform <- function(df) {
#   df %>%
#     column_to_rownames("Genotype") %>%
#     as.matrix() %>%
#     scale(center = TRUE, scale = TRUE) %>%
#     as.data.frame() %>%
#     rownames_to_column("Genotype")
# }
# 
# # Apply to both datasets
# lowinput <- as.tibble(zscore_transform(lowinput))
# control  <- as.tibble(zscore_transform(control))
# 
# # Create a list of z-scored datasets for use in the app
# dataset_list <- list(
#   lowinput = lowinput,
#   control  = control
# )

# Class
lipid_class_info <- vroom::vroom("data/lowinput_class.csv", show_col_types = FALSE) %>%
  filter(!is.na(Class))

# Load gene annotation
gene_annotation <- vroom::vroom("data/gene_annotation.txt", delim = "\t", show_col_types = FALSE)

# Load RDS files
annotation_data <- list(
  lowinput = readRDS("data/all_annotations_lowinput.rds"),
  control = readRDS("data/all_annotations_control.rds")
)

# Path to Manhattan figures
#manhattan_dir <- "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/SoLD/figures/manhattans/lowinput"
manhattan_dir <- "manhattans/lowinput"
## If you're using prcomp, you can skip factoextra. If you want fancy visuals, use library(factoextra).
## If you're doing UMAP, library(umap) is also helpful.

## -- Load your data. Adjust paths if needed --
# e.g., lowinput <- vroom::vroom("../data/lowinput_exp.csv")
#       control  <- vroom::vroom("../data/control_exp.csv")

# For demonstration, we assume `control` and `lowinput` are in memory.

ui <- dashboardPage(
  dashboardHeader(title = "SoLD"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Introduction", tabName = "intro_tab", icon = icon("info-circle")),
      menuItem("Data Preview",   tabName = "data_tab",  icon = icon("table")),
      menuItem("Boxplots",       tabName = "plot_tab",  icon = icon("chart-bar")),
      menuItem("Correlation",    tabName = "correlation_tab", icon = icon("project-diagram")),
      menuItem("PCA",            tabName = "pca_tab",   icon = icon("project-diagram")),
      menuItem("t-SNE / UMAP",   tabName = "tsne_tab",  icon = icon("braille")),
      #menuItem("Volcano Plot",   tabName = "volcano_tab", icon = icon("fire")),
      #menuItem("Heatmap",        tabName = "heatmap_tab", icon = icon("th")),
      menuItem("GWAS",           tabName = "gwas_tab",    icon = icon("dna")),
      menuItem("Gene Hits",      tabName = "gene_hits_tab", icon = icon("chart-bar")),
      menuItem("Genotype", tabName = "genotype_tab", icon = icon("filter"))
    )
  ),
  dashboardBody(
    tabItems(
      
      
      ##------------------------------------------------------------------
      ##  TAB 1: INTRODUCTION
      ##------------------------------------------------------------------
      tabItem(
        tabName = "intro_tab",
        fluidRow(
          box(
            width = 12,
            title = "Welcome to the SOrghum Lipidomics Database (SoLD)",
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            p("This Shiny app presents lipidomics data generated from an experiment using the Sorghum Association Panel (SAP), a genetically diverse set of sorghum lines."),
            p("Our goal is to provide a reference lipid profile under two distinct field conditions, enabling researchers to explore, compare, and utilize lipid data from sorghum grown under contrasting environments."),
            br(),
            h4("Experimental Setup"),
            p("We planted SAP lines in two field conditions over two growing seasons (2019 and 2022):"),
            tags$ul(
              tags$li(strong("Control Condition:"), " Standard agricultural input with adequate nitrogen (N), phosphorus (P), and normal planting time."),
              tags$li(strong("Low Input Condition:"), " Characterized by low nitrogen, low phosphorus, and early planting to mimic cooler, more stressful conditions.")
            ),
            br(),
            p("After the plants matured, we performed lipidomics analysis, profiling more than 800 lipid compounds from each condition. These include:"),
            tags$ul(
              tags$li("Condition-specific lipids unique to either control or low-input fields"),
              tags$li("Shared lipids found in both conditions but varying in abundance")
            ),
            br(),
            h4("What You Can Do with This App"),
            tags$ul(
              tags$li("Explore lipid profiles across different genotypes and conditions"),
              tags$li("Visualize condition-specific lipid signatures"),
              tags$li("Download curated datasets"),
              tags$li("Use the data as a reference for your own sorghum lipidomics experiments")
            ),
            br(),
            p("Our long-term goal is to establish a publicly accessible sorghum lipid database to support the research community in advancing sorghum studies, and metabolic studies.")
          )
        )
      ),
      
      
      ##------------------------------------------------------------------
      ##  TAB 1: Data Preview (same as before, with DT to preview data)
      ##------------------------------------------------------------------
      tabItem(
        tabName = "data_tab",
        fluidRow(
          box(
            width = 4,
            title = "Choose Dataset",
            selectInput(
              inputId = "data_choice",
              label = "Dataset:",
              choices = c("control", "lowinput")
            ),
            uiOutput("class_filter_ui")
          )
        ),
        fluidRow(
          box(
            width = 12,
            title = "Table Preview",
            DTOutput("data_table")
          )
        )
      ),
      
      ##------------------------------------------------------------------
      ##  TAB 2: Boxplots (same as before)
      ##------------------------------------------------------------------
      
      ##----------------------------------------------------------------------
      ##  TAB 2: Boxplots with 4 Sub-tabs
      ##----------------------------------------------------------------------
      
      tabItem(
        tabName = "plot_tab",
        tabsetPanel(
          tabPanel("Absolute",
                   fluidRow(
                     box(
                       width = 4,
                       title = "Plot Controls",
                       selectInput(
                         inputId = "plot_dataset",
                         label = "Choose Dataset for Boxplot:",
                         choices = c("control", "lowinput")
                       ),
                       uiOutput("class_filter_plot_ui"),
                       uiOutput("compound_ui")
                     ),
                     box(
                       width = 8,
                       title = "Boxplot(s)",
                       plotOutput("boxplot_single"),
                       plotOutput("boxplot_both")
                     )
                   )
          ),
          tabPanel("Relative",
                   fluidRow(
                     box(
                       width = 4,
                       title = "Ratio",
                       selectInput(
                         inputId = "plot_dataset",  # reused from earlier
                         label = "Choose Dataset:",
                         choices = c("control", "lowinput")
                       ),
                       uiOutput("ratio_ui_1"),
                       uiOutput("ratio_ui_2")
                     ),
                     box(
                       width = 8,
                       title = "Ratio Boxplot",
                       plotOutput("boxplot_ratio"),
                       plotOutput("boxplot_ratio_combined")
                       )
                     )
                   )
          )
      ),
      
      
      
      ##------------------------------------------------------------------
      ##  TAB 3: PCA
      ##------------------------------------------- -----------------------
      tabItem(
        tabName = "pca_tab",
        fluidRow(
          box(
            width = 4,
            title = "PCA Controls",
            selectInput(
              inputId = "pca_dataset",
              label = "Select Dataset:",
              choices = c("control", "lowinput")
            ),
            checkboxInput("pca_center", "Center", TRUE),
            checkboxInput("pca_scale", "Scale", TRUE),
            selectInput(
              inputId = "pca_mode",
              label = "PCA Mode:",
              choices = c("Samples", "Lipids"),
              selected = "Samples"
            ),
            helpText("This tab runs a simple PCA on samples or lipids."),
            verbatimTextOutput("pca_selected")  # ✅ this goes here
          ),
          box(
            width = 8,
            title = "PCA Plot",
            plotlyOutput("pca_plot")  # ✅ make sure it's plotlyOutput, not plotOutput
          )
        )
      ),
      
      
      
      ##------------------------------------------------------------------
      ##  TAB 4: t-SNE / UMAP
      ##------------------------------------------------------------------
      
      tabItem(
        tabName = "tsne_tab",
        fluidRow(
          box(
            width = 4,
            title = "t-SNE / UMAP Controls",
            selectInput(
              inputId = "dr_dataset",
              label = "Select Dataset:",
              choices = c("control", "lowinput")
            ),
            selectInput(
              inputId = "dr_method",
              label = "Method:",
              choices = c("t-SNE", "UMAP")
            ),
            selectInput(
              inputId = "dr_mode",
              label = "Dimension Mode:",
              choices = c("Samples", "Lipids"),
              selected = "Samples"
            ),
            numericInput("perplexity", "t-SNE perplexity:", 10, min = 1, max = 50),
            helpText("For t-SNE, perplexity is used. For UMAP, it's ignored."),
            verbatimTextOutput("dr_selected")  # Show selected labels
          ),
          box(
            width = 8,
            title = "t-SNE / UMAP Plot",
            plotlyOutput("dr_plot")  # Use plotly for interactivity
          )
        )
      ),
      
      
      ##------------------------------------------------------------------
      ##  TAB 5: Volcano Plot
      ##------------------------------------------------------------------
      tabItem(
        tabName = "volcano_tab",
        fluidRow(
          box(
            width = 12,
            title = "Volcano Plot",
            helpText("We run a simple t-test for each lipid in the intersection of columns between 'control' and 'lowinput'. 
                      X-axis = log2 Fold Change, Y-axis = -log10 p-value. 
                      *Use caution: real analyses may require more rigorous stats.*"),
            plotOutput("volcano_plot")
          )
        )
      ),
      
      ##------------------------------------------------------------------
      ##  TAB 6: Heatmap
      ##------------------------------------------------------------------
      tabItem(
        tabName = "heatmap_tab",
        fluidRow(
          box(
            width = 4,
            title = "Heatmap Controls",
            selectInput(
              inputId = "heatmap_dataset",
              label = "Select Dataset:",
              choices = c("control", "lowinput")
            ),
            helpText("This tab shows a simple heatmap of the chosen dataset, 
                     using all numeric columns (2+). 
                     You might want to subset or transform further in real usage.")
          ),
          box(
            width = 8,
            title = "Heatmap",
            plotlyOutput("heatmap_plot", height = "700px")
          )
        )
      ),
      
      ##------------------------------------------------------------------
      ##  TAB 7: GWAS
      ##------------------------------------------------------------------
        # GWAS tab
        tabItem(tabName = "gwas_tab",
                fluidRow(
                  box(title = "GWAS Controls", width = 4, status = "primary", solidHeader = TRUE,
                      selectInput("gwas_dataset", "Select Dataset:", choices = names(annotation_data), selected = "lowinput"),
                      uiOutput("gwas_lipid_selector")
                  ),
                  box(title = "GWAS Manhattan Plot", width = 8, status = "info", solidHeader = TRUE,
                      uiOutput("manhattan_image")
                  )
                ),
                fluidRow(
                  box(title = "GWAS Annotation Table", width = 12, status = "success", solidHeader = TRUE,
                      DTOutput("gwas_table")
                  )
                )
        ),
      
      ##------------------------------------------------------------------
      ##  TAB 8: Gene Hits
      ##------------------------------------------------------------------
      
      tabItem(tabName = "gene_hits_tab",
              fluidRow(
                box(title = "Gene Hit Filters", width = 4, status = "primary", solidHeader = TRUE,
                    selectInput("hit_dataset", "Select Dataset:", choices = names(annotation_data), selected = "lowinput"),
                    selectInput("hit_threshold", "Select -log10(p) threshold:", choices = c(7, 5), selected = 7),
                    uiOutput("hit_class_filter")  # 👈 Add this for lipid class filter
                ),
                box(title = "Top Genes by Hit Count", width = 8, status = "success", solidHeader = TRUE,
                    DTOutput("gene_hit_table")
                )
              )
      ),
      
      
      ##------------------------------------------------------------------
      ##  TAB 9: Correlations
      ##------------------------------------------------------------------
      
      tabItem(tabName = "correlation_tab",
              fluidRow(
                box(title = "Correlation Controls", width = 4, status = "primary", solidHeader = TRUE,
                    selectInput("cor_dataset", "Select Dataset:", choices = names(dataset_list), selected = "lowinput"),
                    uiOutput("lipid1_selector"),
                    uiOutput("lipid2_selector")
                ),
                box(title = "Lipid Correlation Plot", width = 8, status = "info", solidHeader = TRUE,
                    plotOutput("correlation_plot")
                )
              )
      ),
      
      ##------------------------------------------------------------------
      ## TAB 9: Select Genotype
      ##------------------------------------------------------------------
      
      tabItem(
        tabName = "genotype_tab",
        tabsetPanel(
          tabPanel("Genotype Selection",
                   fluidRow(
                     box(
                       width = 4,
                       title = "Genotype Selection Controls",
                       uiOutput("geno_dataset_ui"),
                       uiOutput("geno_lipid_ui"),
                       selectInput("geno_filter_method", "Filter Type:",
                                   choices = c("Custom Range", "Lowest Quartile", "Highest Quartile")),
                       conditionalPanel(
                         condition = "input.geno_filter_method == 'Custom Range'",
                         uiOutput("geno_slider_ui")
                       ),
                       actionButton("filter_genotypes", "Submit", icon = icon("check"))
                     ),
                     box(
                       width = 8,
                       title = "Matching Genotypes",
                       verbatimTextOutput("genotype_output")
                     )
                   )
          ),
          
          tabPanel("Genotype Viewer",
                   fluidRow(
                     box(
                       width = 4,
                       title = "View Genotype Lipids",
                       uiOutput("geno_dataset_ui"),     # Reuse dataset selector
                       uiOutput("geno_viewer_ui"),      # Genotype selector
                       uiOutput("geno_view_lipid_ui")   # Lipid selector
                     ),
                     box(
                       width = 8,
                       title = "Lipid Level Across Conditions",
                       plotOutput("geno_view_plot")
                     )
                   )
          )
        )
      )
      
      
      )
    )
)


server <- function(input, output, session) {
  
  ##----------------------------------------------------------------------
  ##  0. Put your data in a named list for easy reference
  ##----------------------------------------------------------------------
  dataset_list <- list(
    control  = control,
    lowinput = lowinput
  )
  
  ##======================================================================
  ##  TAB 1: Data Preview
  ##======================================================================
  
  output$class_filter_ui <- renderUI({
    classes <- sort(unique(lipid_class_info$Class))
    selectInput(
      inputId = "class_choice",
      label = "Filter by Lipid Class:",
      choices = c("All", classes),
      selected = "All"
    )
  })
  
  output$data_table <- renderDT({
    df <- dataset_list[[input$data_choice]]
    
    # Transpose for display
    df_t <- as.data.frame(t(df[, -1]))
    colnames(df_t) <- df[[1]]
    df_t$LipidName <- rownames(df_t)
    rownames(df_t) <- NULL
    df_t <- df_t[, c(ncol(df_t), 1:(ncol(df_t)-1))]
    
    # Apply class filter if selected
    if (!is.null(input$class_choice) && input$class_choice != "All") {
      lipid_subset <- lipid_class_info %>%
        filter(Class == input$class_choice) %>%
        pull(Lipids)
      
      df_t <- df_t %>% filter(LipidName %in% lipid_subset)
    }
    
    datatable(
      df_t,
      options = list(scrollX = TRUE, pageLength = 10, lengthMenu = c(5, 10, 25, 50)),
      rownames = FALSE
    ) %>% formatRound(columns = 2:ncol(df_t), digits = 2)
  })
  
  
  
  
  ##======================================================================
  ##  TAB 2: Boxplots
  ##======================================================================
  
  # # UI for lipid class dropdown (with "All")
  # output$class_filter_plot_ui <- renderUI({
  #   classes <- sort(unique(na.omit(lipid_class_info$Class)))
  #   selectInput(
  #     inputId = "class_filter_plot",
  #     label = "Filter by Lipid Class:",
  #     choices = c("All", classes),
  #     selected = "All"
  #   )
  # })
  # 
  # # UI for check box to sum lipids (only if not "All")
  # output$sum_lipids_ui <- renderUI({
  #   req(input$class_filter_plot)
  #   if (input$class_filter_plot != "All") {
  #     checkboxInput("sum_lipids", "Sum all lipids in this class?", value = FALSE)
  #   }
  # })
  # 
  # # UI for lipid compound dropdown, includes "All (Sum)"
  # output$compound_ui <- renderUI({
  #   req(input$plot_dataset)
  #   df <- dataset_list[[input$plot_dataset]]
  #   all_lipids <- colnames(df)[-1]
  #   
  #   if (!is.null(input$class_filter_plot) && input$class_filter_plot != "All") {
  #     selected_lipids <- lipid_class_info$Lipids[lipid_class_info$Class == input$class_filter_plot]
  #     available_lipids <- intersect(selected_lipids, all_lipids)
  #     choices <- c("All (Sum)" = "ALL_LIPIDS_SUM", available_lipids)
  #   } else {
  #     available_lipids <- all_lipids
  #     choices <- available_lipids
  #   }
  #   
  #   selectInput(
  #     inputId = "compound_choice",
  #     label = "Choose a Lipid Compound:",
  #     choices = choices,
  #     selected = available_lipids[1]
  #   )
  # })
  # 
  # # Single dataset plot
  # output$boxplot_single <- renderPlot({
  #   req(input$plot_dataset, input$compound_choice)
  #   df <- dataset_list[[input$plot_dataset]]
  #   
  #   if (input$compound_choice == "ALL_LIPIDS_SUM") {
  #     lipids_in_class <- lipid_class_info$Lipids[lipid_class_info$Class == input$class_filter_plot]
  #     valid_lipids <- intersect(lipids_in_class, colnames(df))
  #     
  #     df_plot <- df %>%
  #       select(Sample = 1, all_of(valid_lipids)) %>%
  #       mutate(Value = rowSums(across(all_of(valid_lipids)), na.rm = TRUE)) %>%
  #       select(Sample, Value) %>%
  #       mutate(Condition = input$plot_dataset)
  #     
  #     plot_title <- paste("Summed abundance of", input$class_filter_plot, "class")
  #     
  #   } else {
  #     df_plot <- df %>%
  #       select(Sample = 1, Value = all_of(input$compound_choice)) %>%
  #       mutate(Condition = input$plot_dataset)
  #     
  #     plot_title <- paste("Distribution of", input$compound_choice, "in", input$plot_dataset)
  #   }
  #   
  #   ggplot(df_plot, aes(x = Condition, y = Value, fill = Condition)) +
  #     geom_violin(alpha = 0.5, width = 0.8) +
  #     geom_boxplot(width = 0.1, outlier.shape = NA) +
  #     labs(title = plot_title, x = "Condition", y = "Value") +
  #     theme_minimal(base_size = 16) +
  #     theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  # })
  # 
  # # Combined dataset plot
  # output$boxplot_both <- renderPlot({
  #   req(input$compound_choice, input$plot_dataset, input$class_filter_plot)
  #   
  #   other_dataset <- setdiff(names(dataset_list), input$plot_dataset)
  #   df_main <- dataset_list[[input$plot_dataset]]
  #   df_other <- dataset_list[[other_dataset]]
  #   
  #   if (input$compound_choice == "ALL_LIPIDS_SUM") {
  #     lipids_in_class <- lipid_class_info$Lipids[lipid_class_info$Class == input$class_filter_plot]
  #     valid_main <- intersect(lipids_in_class, colnames(df_main))
  #     valid_other <- intersect(lipids_in_class, colnames(df_other))
  #     
  #     if (length(valid_main) == 0 || length(valid_other) == 0) return(NULL)
  #     
  #     df_main_sum <- df_main %>%
  #       select(Sample = 1, all_of(valid_main)) %>%
  #       mutate(Value = rowSums(across(all_of(valid_main)), na.rm = TRUE),
  #              Condition = input$plot_dataset) %>%
  #       select(Sample, Value, Condition)
  #     
  #     df_other_sum <- df_other %>%
  #       select(Sample = 1, all_of(valid_other)) %>%
  #       mutate(Value = rowSums(across(all_of(valid_other)), na.rm = TRUE),
  #              Condition = other_dataset) %>%
  #       select(Sample, Value, Condition)
  #     
  #     df_combined <- bind_rows(df_main_sum, df_other_sum)
  #     
  #     plot_title <- paste("Summed abundance of", input$class_filter_plot, "across conditions")
  #     
  #   } else {
  #     compound_col <- input$compound_choice
  #     
  #     if (!(compound_col %in% colnames(df_main)) || !(compound_col %in% colnames(df_other))) {
  #       return(NULL)
  #     }
  #     
  #     df_main_plot <- df_main %>%
  #       select(Sample = 1, Value = all_of(compound_col)) %>%
  #       mutate(Condition = input$plot_dataset)
  #     
  #     df_other_plot <- df_other %>%
  #       select(Sample = 1, Value = all_of(compound_col)) %>%
  #       mutate(Condition = other_dataset)
  #     
  #     df_combined <- bind_rows(df_main_plot, df_other_plot)
  #     
  #     plot_title <- paste("Comparison of", compound_col, "across conditions")
  #   }
  #   
  #   ggplot(df_combined, aes(x = Condition, y = Value, fill = Condition)) +
  #     geom_violin(alpha = 0.5, width = 0.8) +
  #     geom_boxplot(width = 0.1, outlier.shape = NA) +
  #     labs(
  #       title = plot_title,
  #       x = "Condition",
  #       y = "Value"
  #     ) +
  #     theme_minimal(base_size = 16) +
  #     theme(
  #       axis.text.x = element_text(angle = 45, hjust = 1),
  #       legend.position = "none"
  #     )
  # })
  # 
  
  
  ##======================================================================
  ##  TAB 2: Boxplots with Sub-tabs (Individual, Class, Ratio Individual, Ratio Class)
  ##======================================================================
  
  # UI for lipid class dropdown (with "All")
  output$class_filter_plot_ui <- renderUI({
    classes <- sort(unique(na.omit(lipid_class_info$Class)))
    selectInput(
      inputId = "class_filter_plot",
      label = "Filter by Lipid Class:",
      choices = c("All", classes),
      selected = "All"
    )
  })
  
  # UI for lipid compound dropdown, includes "All (Sum)"
  output$compound_ui <- renderUI({
    req(input$plot_dataset)
    df <- dataset_list[[input$plot_dataset]]
    all_lipids <- colnames(df)[-1]
    
    if (!is.null(input$class_filter_plot) && input$class_filter_plot != "All") {
      selected_lipids <- lipid_class_info$Lipids[lipid_class_info$Class == input$class_filter_plot]
      available_lipids <- intersect(selected_lipids, all_lipids)
      choices <- c("All (Sum)" = "ALL_LIPIDS_SUM", available_lipids)
    } else {
      available_lipids <- all_lipids
      choices <- available_lipids
    }
    
    selectInput(
      inputId = "compound_choice",
      label = "Choose a Lipid Compound:",
      choices = choices,
      selected = available_lipids[1]
    )
  })
  
  
  # # Single dataset plot
  output$boxplot_single <- renderPlot({
    req(input$plot_dataset, input$compound_choice)
    df <- dataset_list[[input$plot_dataset]]

    if (input$compound_choice == "ALL_LIPIDS_SUM") {
      lipids_in_class <- lipid_class_info$Lipids[lipid_class_info$Class == input$class_filter_plot]
      valid_lipids <- intersect(lipids_in_class, colnames(df))

      df_plot <- df %>%
        select(Sample = 1, all_of(valid_lipids)) %>%
        mutate(Value = rowSums(across(all_of(valid_lipids)), na.rm = TRUE)) %>%
        select(Sample, Value) %>%
        mutate(Condition = input$plot_dataset)

      plot_title <- paste("Summed abundance of", input$class_filter_plot, "class")

    } else {
      df_plot <- df %>%
        select(Sample = 1, Value = all_of(input$compound_choice)) %>%
        mutate(Condition = input$plot_dataset)

      plot_title <- paste("Distribution of", input$compound_choice, "in", input$plot_dataset)
    }

    ggplot(df_plot, aes(x = Condition, y = Value, fill = Condition)) +
      geom_violin(alpha = 0.5, width = 0.8) +
      geom_boxplot(width = 0.1, outlier.shape = NA) +
      labs(title = plot_title, x = "Condition", y = "Value") +
      theme_minimal(base_size = 16) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  })

  # Combined dataset plot
  output$boxplot_both <- renderPlot({
    req(input$compound_choice, input$plot_dataset, input$class_filter_plot)

    other_dataset <- setdiff(names(dataset_list), input$plot_dataset)
    df_main <- dataset_list[[input$plot_dataset]]
    df_other <- dataset_list[[other_dataset]]

    if (input$compound_choice == "ALL_LIPIDS_SUM") {
      lipids_in_class <- lipid_class_info$Lipids[lipid_class_info$Class == input$class_filter_plot]
      valid_main <- intersect(lipids_in_class, colnames(df_main))
      valid_other <- intersect(lipids_in_class, colnames(df_other))

      if (length(valid_main) == 0 || length(valid_other) == 0) return(NULL)

      df_main_sum <- df_main %>%
        dplyr::select(Sample = 1, all_of(valid_main)) %>%
        dplyr::mutate(Value = rowSums(across(all_of(valid_main)), na.rm = TRUE),
               Condition = input$plot_dataset) %>%
        dplyr::select(Sample, Value, Condition)

      df_other_sum <- df_other %>%
        dplyr::select(Sample = 1, all_of(valid_other)) %>%
        dplyr::mutate(Value = rowSums(across(all_of(valid_other)), na.rm = TRUE),
               Condition = other_dataset) %>%
        dplyr::select(Sample, Value, Condition)

      df_combined <- bind_rows(df_main_sum, df_other_sum)

      plot_title <- paste("Summed abundance of", input$class_filter_plot, "across conditions")

    } else {
      compound_col <- input$compound_choice

      if (!(compound_col %in% colnames(df_main)) || !(compound_col %in% colnames(df_other))) {
        return(NULL)
      }

      df_main_plot <- df_main %>%
        select(Sample = 1, Value = all_of(compound_col)) %>%
        mutate(Condition = input$plot_dataset)

      df_other_plot <- df_other %>%
        select(Sample = 1, Value = all_of(compound_col)) %>%
        mutate(Condition = other_dataset)

      df_combined <- bind_rows(df_main_plot, df_other_plot)

      plot_title <- paste("Comparison of", compound_col, "across conditions")
    }

    ggplot(df_combined, aes(x = Condition, y = Value, fill = Condition)) +
      geom_violin(alpha = 0.5, width = 0.8) +
      geom_boxplot(width = 0.1, outlier.shape = NA) +
      labs(
        title = plot_title,
        x = "Condition",
        y = "Value"
      ) +
      theme_minimal(base_size = 16) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
  })
  
  #=======================
  # Ratio Tabs (UI)
  #=======================
  
  output$ratio_ui_1 <- renderUI({
    req(input$plot_dataset)
    df <- dataset_list[[input$plot_dataset]]
    all_lipids <- colnames(df)[-1]
    sum_lipids <- paste0("SUM_", sort(unique(lipid_class_info$Class)))
    selectInput("ratio_lipid_1", "Numerator Lipid:", choices = c(sum_lipids, all_lipids))
  })
  
  output$ratio_ui_2 <- renderUI({
    req(input$plot_dataset)
    df <- dataset_list[[input$plot_dataset]]
    all_lipids <- colnames(df)[-1]
    sum_lipids <- paste0("SUM_", sort(unique(lipid_class_info$Class)))
    selectInput("ratio_lipid_2", "Denominator Lipid:", choices = c(sum_lipids, all_lipids))
  })
  
  #=======================
  # Shared Function for Boxplot Rendering
  #=======================
  
  plot_lipid_boxplot <- function(df, values, title, condition) {
    df_plot <- data.frame(Sample = df[[1]], Value = values, Condition = condition)
    
    ggplot(df_plot, aes(x = Condition, y = Value, fill = Condition)) +
      geom_violin(alpha = 0.5, width = 0.8) +
      geom_boxplot(width = 0.1, outlier.shape = NA) +
      labs(title = title, x = "Condition", y = "log10(Ratio)") +
      theme_minimal(base_size = 16) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  }
  
  get_ratio_values <- function(df, lipid_input) {
    if (startsWith(lipid_input, "SUM_")) {
      class_name <- sub("^SUM_", "", lipid_input)
      lipids <- lipid_class_info$Lipids[lipid_class_info$Class == class_name]
      valid <- intersect(lipids, colnames(df))
      rowSums(df[, valid, drop = FALSE], na.rm = TRUE)
    } else {
      df[[lipid_input]]
    }
  }
  
  #=======================
  # Ratio Individual Dataset Plot
  #=======================
  output$boxplot_ratio <- renderPlot({
    req(input$plot_dataset, input$ratio_lipid_1, input$ratio_lipid_2)
    df <- dataset_list[[input$plot_dataset]]
    
    numerator <- get_ratio_values(df, input$ratio_lipid_1)
    denominator <- get_ratio_values(df, input$ratio_lipid_2)
    ratio_vals <- log10(numerator) - log10(denominator)
    
    plot_lipid_boxplot(df, ratio_vals,
                       paste("log10(Ratio):", input$ratio_lipid_1, "/", input$ratio_lipid_2),
                       input$plot_dataset)
  })
  
  #=======================
  # Ratio Combined Plot (Control vs LowInput)
  #=======================
  # output$boxplot_ratio_combined <- renderPlot({
  #   req(input$ratio_lipid_1, input$ratio_lipid_2)
  #   datasets <- names(dataset_list)
  #   
  #   df_combined <- do.call(rbind, lapply(datasets, function(ds) {
  #     df <- dataset_list[[ds]]
  #     numerator <- get_ratio_values(df, input$ratio_lipid_1)
  #     denominator <- get_ratio_values(df, input$ratio_lipid_2)
  #     ratio_vals <- log10(numerator) - log10(denominator) ### ASK FAUSTO
  #     data.frame(Sample = df[[1]], Value = ratio_vals, Condition = ds)
  #   }))
  #   
  #   ggplot(df_combined, aes(x = Condition, y = Value, fill = Condition)) +
  #     geom_violin(alpha = 0.5, width = 0.8) +
  #     geom_boxplot(width = 0.1, outlier.shape = NA) +
  #     labs(
  #       title = paste("log10(Ratio of", input$ratio_lipid_1, "/", input$ratio_lipid_2, ") across conditions"),
  #       y = "log10(Ratio)",
  #       x = "Condition"
  #     ) +
  #     theme_minimal(base_size = 16) +
  #     theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  # })
  # 
  
  output$boxplot_ratio_combined <- renderPlot({
    req(input$ratio_lipid_1, input$ratio_lipid_2)
    datasets <- names(dataset_list)
    
    df_combined <- do.call(rbind, lapply(datasets, function(ds) {
      df <- dataset_list[[ds]]
      numerator <- get_ratio_values(df, input$ratio_lipid_1)
      denominator <- get_ratio_values(df, input$ratio_lipid_2)
      ratio_vals <- log10(numerator / denominator)
      data.frame(Sample = df[[1]], Value = ratio_vals, Condition = ds)
    }))
    
    # Perform t-test
    ttest <- t.test(Value ~ Condition, data = df_combined)
    pval <- ttest$p.value
    
    # Format significance
    sig_label <- if (pval < 0.001) {
      "***"
    } else if (pval < 0.01) {
      "**"
    } else if (pval < 0.05) {
      "*"
    } else {
      "ns"
    }
    
    ggplot(df_combined, aes(x = Condition, y = Value, fill = Condition)) +
      geom_violin(alpha = 0.5, width = 0.8) +
      geom_boxplot(width = 0.1, outlier.shape = NA) +
      geom_text(
        aes(x = 1.5, y = max(df_combined$Value, na.rm = TRUE) * 1.05, label = sig_label),
        inherit.aes = FALSE,
        size = 6
      ) +
      labs(
        title = paste("log10(Ratio of", input$ratio_lipid_1, "/", input$ratio_lipid_2, ") across conditions"),
        y = "log10(Ratio)",
        x = "Condition"
      ) +
      theme_minimal(base_size = 16) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  })
  
  ##======================================================================
  ##  TAB 3: PCA
  ##======================================================================
 
  output$pca_selected <- renderPrint({
    eventdata <- event_data("plotly_selected")
    if (is.null(eventdata)) {
      cat("Select points on the plot to see their labels here.")
    } else {
      unique(eventdata$text)
    }
  })
  
  output$pca_plot <- renderPlotly({  # plotly for interaction
    df <- dataset_list[[input$pca_dataset]]
    
    if (input$pca_mode == "Samples") {
      numeric_df <- df[, -1]
      numeric_df <- as.data.frame(sapply(numeric_df, as.numeric))
      pca_res <- prcomp(numeric_df, center = input$pca_center, scale. = input$pca_scale)
      pca_data <- as.data.frame(pca_res$x)
      pca_data$ID <- df[[1]]  # sample ID
    } else {
      numeric_df <- df[, -1]
      numeric_df <- t(numeric_df)
      colnames(numeric_df) <- df[[1]]
      numeric_df <- as.data.frame(numeric_df)
      pca_res <- prcomp(numeric_df, center = input$pca_center, scale. = input$pca_scale)
      pca_data <- as.data.frame(pca_res$x)
      pca_data$ID <- rownames(pca_data)  # lipid name
    }
    
    # Interactive plot with plotly
    plot_ly(
      data = pca_data,
      x = ~PC1,
      y = ~PC2,
      text = ~ID,
      type = "scatter",
      mode = "markers",
      marker = list(size = 7)
    ) %>%
      layout(
        title = paste("PCA -", input$pca_mode, "(", input$pca_dataset, ")"),
        xaxis = list(title = "PC1"),
        yaxis = list(title = "PC2"),
        dragmode = "select"
      )
  })
  
  
  ##======================================================================
  ##  TAB 4: t-SNE / UMAP
  ##======================================================================
  
  output$dr_selected <- renderPrint({
    eventdata <- event_data("plotly_selected")
    if (is.null(eventdata)) {
      cat("Select points on the plot to see labels here.")
    } else {
      unique(eventdata$text)
    }
  })
  
  output$dr_plot <- renderPlotly({
    df <- dataset_list[[input$dr_dataset]]
    method <- input$dr_method
    
    # Handle sample vs lipid mode
    if (input$dr_mode == "Samples") {
      numeric_df <- df[, -1]
      numeric_df <- as.data.frame(sapply(numeric_df, as.numeric))
      row_labels <- df[[1]]
    } else {
      numeric_df <- df[, -1]
      numeric_df <- t(numeric_df)
      colnames(numeric_df) <- df[[1]]
      numeric_df <- as.data.frame(numeric_df)
      row_labels <- rownames(numeric_df)
    }
    
    # Dimensionality reduction
    if (method == "t-SNE") {
      library(Rtsne)
      dr_out <- Rtsne(numeric_df, perplexity = input$perplexity, check_duplicates = FALSE)
      dim_df <- data.frame(Dim1 = dr_out$Y[,1], Dim2 = dr_out$Y[,2], Label = row_labels)
    } else {
      library(umap)
      dr_out <- umap(numeric_df)
      dim_df <- data.frame(Dim1 = dr_out$layout[,1], Dim2 = dr_out$layout[,2], Label = row_labels)
    }
    
    plot_ly(
      data = dim_df,
      x = ~Dim1,
      y = ~Dim2,
      type = "scatter",
      mode = "markers",
      text = ~Label,
      marker = list(size = 6)
    ) %>%
      layout(
        title = paste(method, "-", input$dr_mode, "(", input$dr_dataset, ")"),
        xaxis = list(title = "Dim 1"),
        yaxis = list(title = "Dim 2"),
        dragmode = "select"
      )
  })
  
  
  
  ##======================================================================
  ##  TAB 5: Volcano Plot
  ##======================================================================
  output$volcano_plot <- renderPlot({
    # We'll do a simple approach:
    # 1) Find columns in both 'control' and 'lowinput'
    # 2) For each of those columns, run a t-test
    # 3) log2FC = mean(lowinput) - mean(control) in log2
    #    pval from t-test
    
    common_cols <- intersect(colnames(control), colnames(lowinput))
    # remove the first column (Sample ID) from analysis
    common_cols <- setdiff(common_cols, colnames(control)[1])
    
    # if no common lipid columns, just bail
    validate(
      need(length(common_cols) > 0, "No overlapping lipid columns found!")
    )
    
    # Build a results data frame
    results <- data.frame(Lipid = common_cols, log2FC = NA, pval = NA)
    
    for (i in seq_along(common_cols)) {
      col_i <- common_cols[i]
      # Extract numeric vectors
      v_control  <- as.numeric(control[[col_i]])
      v_lowinput <- as.numeric(lowinput[[col_i]])
      
      # T-test
      # For fold change: typically (LowInput / Control). We'll do log2 of means
      t_out <- t.test(v_lowinput, v_control)
      mean_control  <- mean(v_control,  na.rm = TRUE)
      mean_lowinput <- mean(v_lowinput, na.rm = TRUE)
      # Avoid zero or negative means in log (in real data handle carefully)
      # We'll do log2FC = log2( mean_lowinput / mean_control )
      fc       = mean_lowinput / mean_control
      log2fc   = log2(fc)
      
      results$log2FC[i] <- log2fc
      results$pval[i]   <- t_out$p.value
    }
    
    # -log10 p-value
    results$negLog10p <- -log10(results$pval)
    
    # Basic volcano
    ggplot(results, aes(x = log2FC, y = negLog10p)) +
      geom_point() +
      labs(
        title = "Volcano Plot (LowInput vs Control)",
        x = "log2 Fold Change",
        y = "-log10 p-value"
      ) +
      theme_minimal()
  })
  
  ##======================================================================
  ##  TAB 6: Heatmap
  ##======================================================================
  # output$heatmap_plot <- renderPlot({
  #   # We'll do a very minimal heatmap using pheatmap
  #   # For a more advanced approach, see 'ComplexHeatmap'
  #   library(pheatmap)
  #   
  #   df <- dataset_list[[input$heatmap_dataset]]
  #   
  #   # Use numeric columns only (excluding sample ID)
  #   mat <- df[, -1]
  #   mat <- as.data.frame(sapply(mat, as.numeric))
  #   
  #   # pheatmap can show row names, so let's set them to sample IDs
  #   rownames(mat) <- df[[1]]
  #   
  #   # Simple pheatmap
  #   pheatmap(mat,
  #            scale = "row",           # optional, normalizes each row
  #            show_rownames = TRUE,
  #            show_colnames = TRUE)
  # })
  
  output$heatmap_plot <- renderPlotly({
    df <- dataset_list[[input$heatmap_dataset]]
    
    # Strip first column and ensure numeric
    mat <- df[, -1]
    mat <- as.data.frame(sapply(mat, as.numeric))
    rownames(mat) <- df[[1]]
    
    # Remove any rows/columns that are all NA (just in case)
    mat <- mat[rowSums(is.na(mat)) != ncol(mat), ]
    mat <- mat[, colSums(is.na(mat)) != nrow(mat)]
    
    # Validate size
    validate(
      need(nrow(mat) > 1 && ncol(mat) > 1, "Not enough data for heatmap.")
    )
    
    # Plot heatmaply
    heatmaply(
      mat,
      scale = "row",
      dendrogram = "both",
      showticklabels = c(TRUE, TRUE),
      colors = colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlBu"))(255),
      xlab = "Lipids",
      ylab = "Samples",
      labRow = rownames(mat),
      labCol = colnames(mat),
      margins = c(40, 120, 40, 120)
    )
  })
  
  
  ##======================================================================
  ##  TAB 7: GWAS
  ##======================================================================
  
  # Lipid dropdown based on dataset
  output$gwas_lipid_selector <- renderUI({
    req(input$gwas_dataset)
    lipid_list <- names(annotation_data[[input$gwas_dataset]])
    selectInput("selected_lipid", "Select Lipid:", choices = lipid_list)
  })
  
  # Manhattan plot (dynamic path based on dataset)
  output$manhattan_image <- renderUI({
    req(input$gwas_dataset, input$selected_lipid)
    lipid <- input$selected_lipid
    
    # Build the appropriate directory path
    dir_name <- paste0("manhattans/", input$gwas_dataset)
    img_path <- file.path(dir_name, paste0(lipid, ".jpg"))
    
    # Check if file exists inside www/
    if (!file.exists(file.path("www", img_path))) {
      return(tags$p("No Manhattan plot available for this lipid."))
    }
    
    tags$img(
      src = img_path,
      width = "100%",
      style = "border:1px solid #ddd; padding:10px;"
    )
  })
  
  
  # GWAS Table
  output$gwas_table <- renderDT({
    req(input$gwas_dataset, input$selected_lipid)
    df <- annotation_data[[input$gwas_dataset]][[input$selected_lipid]]
    validate(need(!is.null(df), "No annotation table found for this lipid."))
    
    df_annotated <- df %>%
      left_join(gene_annotation, by = "GeneID")
    
    datatable(df_annotated, options = list(scrollX = TRUE, pageLength = 10))
  })
  
  ##======================================================================
  ##  TAB 8: Gene Hits
  ##======================================================================
  
  # Dropdown for class filter
  output$hit_class_filter <- renderUI({
    classes <- sort(unique(lipid_class_info$Class))
    selectInput("hit_class_choice", "Filter by Lipid Class:", choices = c("All", classes), selected = "All")
  })
  # 
  # # Gene Hits Table
  # output$gene_hit_table <- renderDT({
  #   req(input$hit_dataset, input$hit_threshold)
  #   data_list <- annotation_data[[input$hit_dataset]]
  #   threshold <- as.numeric(input$hit_threshold)
  #   
  #   # Combine all entries from all lipids into one table
  #   combined_df <- do.call(rbind, lapply(data_list, function(df) {
  #     if (!is.null(df) && all(c("GeneID", "log(p)") %in% colnames(df))) {
  #       df[df$`log(p)` >= threshold, c("GeneID", "log(p)")]
  #     }
  #   }))
  #   
  #   # Count gene occurrences
  #   gene_summary <- as.data.frame(table(combined_df$GeneID))
  #   colnames(gene_summary) <- c("GeneID", "Occurrences")
  #   
  #   # Sort by highest count and join with gene annotation
  #   gene_summary <- gene_summary %>%
  #     arrange(desc(Occurrences)) %>%
  #     left_join(gene_annotation, by = "GeneID")
  #   
  #   datatable(gene_summary, options = list(scrollX = TRUE, pageLength = 10))
  # })
  # 
  output$gene_hit_table <- renderDT({
    req(input$hit_dataset, input$hit_threshold)
    data_list <- annotation_data[[input$hit_dataset]]
    threshold <- as.numeric(input$hit_threshold)
    
    # Filter lipids by selected class
    lipid_subset <- names(data_list)
    if (!is.null(input$hit_class_choice) && input$hit_class_choice != "All") {
      class_lipids <- lipid_class_info %>%
        filter(Class == input$hit_class_choice) %>%
        pull(Lipids)
      lipid_subset <- intersect(lipid_subset, class_lipids)
    }
    
    # Combine hits only for selected lipids
    combined_df <- do.call(rbind, lapply(data_list[lipid_subset], function(df) {
      if (!is.null(df) && all(c("GeneID", "log(p)") %in% colnames(df))) {
        df[df$`log(p)` >= threshold, c("GeneID", "log(p)")]
      }
    }))
    
    # Count occurrences
    gene_summary <- as.data.frame(table(combined_df$GeneID))
    colnames(gene_summary) <- c("GeneID", "Occurrences")
    
    gene_summary <- gene_summary %>%
      arrange(desc(Occurrences)) %>%
      left_join(gene_annotation, by = "GeneID")
    
    datatable(gene_summary, options = list(scrollX = TRUE, pageLength = 10))
  })
  
  
  ##======================================================================
  ##  TAB 9: Correlations
  ##======================================================================
  
  output$lipid1_selector <- renderUI({
    df <- dataset_list[[input$cor_dataset]]
    lipid_names <- colnames(df)[-1]
    selectInput("lipid1", "Select Lipid 1:", choices = lipid_names)
  })
  
  output$lipid2_selector <- renderUI({
    df <- dataset_list[[input$cor_dataset]]
    lipid_names <- colnames(df)[-1]
    selectInput("lipid2", "Select Lipid 2:", choices = lipid_names)
  })
  
  output$correlation_plot <- renderPlot({
    req(input$cor_dataset, input$lipid1, input$lipid2)
    
    df <- dataset_list[[input$cor_dataset]]
    
    validate(
      need(all(c(input$lipid1, input$lipid2) %in% colnames(df)), "Selected lipids not found in dataset")
    )
    
    # Clean: filter out rows where either value is NA
    df_clean <- df %>%
      dplyr::filter(!is.na(.data[[input$lipid1]]) & !is.na(.data[[input$lipid2]]))
    
    # Dynamically assign variables for plotting
    ggplot(df_clean, aes(x = .data[[input$lipid1]], y = .data[[input$lipid2]])) +
      geom_point(alpha = 0.7) +
      geom_smooth(method = "lm", color = "blue", se = FALSE) +
      labs(
        title = paste("Correlation between", input$lipid1, "and", input$lipid2),
        x = input$lipid1,
        y = input$lipid2
      ) +
      theme_minimal() +
      theme_minimal(base_size = 20) +
      theme(
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_line(),          # Add base axes
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
  })
  
  # UI for lipid choices in the selected dataset
  output$geno_lipid_ui <- renderUI({
    req(input$geno_dataset)
    df <- dataset_list[[input$geno_dataset]]
    lipid_choices <- colnames(df)[-1]
    selectInput("geno_lipid", "Choose Lipid:", choices = lipid_choices)
  })
  
  # UI for slider based on selected lipid range
  output$geno_slider_ui <- renderUI({
    req(input$geno_dataset, input$geno_lipid)
    df <- dataset_list[[input$geno_dataset]]
    vals <- df[[input$geno_lipid]]
    sliderInput("geno_range", "Select Lipid Value Range:",
                min = floor(min(vals, na.rm = TRUE)),
                max = ceiling(max(vals, na.rm = TRUE)),
                value = c(floor(min(vals, na.rm = TRUE)), ceiling(max(vals, na.rm = TRUE))),
                step = 0.01)
  })
  
  ##======================================================================
  ##  TAB 9: Select Genotype
  ##======================================================================
  
  # UI: Dataset and Lipid Selectors
  output$geno_dataset_ui <- renderUI({
    selectInput("geno_dataset", "Choose Dataset:", choices = names(dataset_list))
  })
  
  output$geno_lipid_ui <- renderUI({
    req(input$geno_dataset)
    df <- dataset_list[[input$geno_dataset]]
    all_lipids <- colnames(df)[-1]
    class_options <- paste0("SUM_", sort(unique(lipid_class_info$Class)))
    selectInput("geno_lipid", "Choose Lipid:", choices = c(class_options, all_lipids))
  })
  
  # UI: Slider for Custom Range
  output$geno_slider_ui <- renderUI({
    req(input$geno_dataset, input$geno_lipid)
    df <- dataset_list[[input$geno_dataset]]
    
    values <- if (startsWith(input$geno_lipid, "SUM_")) {
      class_name <- sub("SUM_", "", input$geno_lipid)
      lipids <- lipid_class_info$Lipids[lipid_class_info$Class == class_name]
      valid <- intersect(lipids, colnames(df))
      rowSums(df[, valid, drop = FALSE], na.rm = TRUE)
    } else {
      df[[input$geno_lipid]]
    }
    
    sliderInput("geno_range", "Select Lipid Value Range:",
                min = floor(min(values, na.rm = TRUE)),
                max = ceiling(max(values, na.rm = TRUE)),
                value = quantile(values, c(0.25, 0.75), na.rm = TRUE))
  })
  
  # Logic: Handle Filter Button
  observeEvent(input$filter_genotypes, {
    req(input$geno_dataset, input$geno_lipid)
    df <- dataset_list[[input$geno_dataset]]
    
    values <- if (startsWith(input$geno_lipid, "SUM_")) {
      class_name <- sub("SUM_", "", input$geno_lipid)
      lipids <- lipid_class_info$Lipids[lipid_class_info$Class == class_name]
      valid <- intersect(lipids, colnames(df))
      rowSums(df[, valid, drop = FALSE], na.rm = TRUE)
    } else {
      df[[input$geno_lipid]]
    }
    
    selected_genos <- NULL
    
    if (input$geno_filter_method == "Lowest Quartile") {
      q1 <- quantile(values, 0.25, na.rm = TRUE)
      selected_genos <- df[[1]][values <= q1]
      
    } else if (input$geno_filter_method == "Highest Quartile") {
      q3 <- quantile(values, 0.75, na.rm = TRUE)
      selected_genos <- df[[1]][values >= q3]
      
    } else {
      req(input$geno_range)
      selected_genos <- df[[1]][values >= input$geno_range[1] & values <= input$geno_range[2]]
    }
    
    output$genotype_output <- renderPrint({
      if (length(selected_genos) == 0) {
        "No genotypes found for the selected filter."
      } else {
        cat("Genotypes in selected filter:\n\n")
        print(selected_genos)
      }
    })
  })
  
  ##======================================================================
  ##  TAB 10: Genotype Viewer
  ##======================================================================
  
  # UI: Genotype selector and lipid selector
  output$geno_viewer_ui <- renderUI({
    req(input$geno_dataset)
    df <- dataset_list[[input$geno_dataset]]
    selectInput("geno_view_id", "Select Genotype:", choices = df[[1]])
  })
  
  output$geno_view_lipid_ui <- renderUI({
    req(input$geno_dataset)
    df <- dataset_list[[input$geno_dataset]]
    all_lipids <- colnames(df)[-1]
    class_options <- paste0("SUM_", sort(unique(lipid_class_info$Class)))
    selectInput("geno_view_lipid", "Select Lipid:", choices = c(class_options, all_lipids))
  })
  
  # Plot viewer (barplot for both conditions)
  output$geno_view_plot <- renderPlot({
    req(input$geno_view_id, input$geno_view_lipid)
    
    get_lipid_value <- function(df, genotype, lipid) {
      if (startsWith(lipid, "SUM_")) {
        class_name <- sub("SUM_", "", lipid)
        lipids <- lipid_class_info$Lipids[lipid_class_info$Class == class_name]
        valid <- intersect(lipids, colnames(df))
        rowSums(df[df[[1]] == genotype, valid, drop = FALSE], na.rm = TRUE)
      } else {
        df[df[[1]] == genotype, lipid, drop = TRUE]
      }
    }
    
    values <- sapply(names(dataset_list), function(dataset) {
      get_lipid_value(dataset_list[[dataset]], input$geno_view_id, input$geno_view_lipid)
    })
    
    df_bar <- data.frame(Condition = names(values), Value = as.numeric(values))
    
    ggplot(df_bar, aes(x = Condition, y = Value, fill = Condition)) +
      geom_col(width = 0.6) +
      labs(title = paste("Lipid Level for", input$geno_view_id), y = "Value", x = "Condition") +
      theme_minimal(base_size = 16) +
      theme(legend.position = "none")
  })
  


  
}

shinyApp(ui, server)