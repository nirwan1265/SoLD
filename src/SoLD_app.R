library(vroom)
## Read the files
lowinput <- vroom("../data/lowinput_exp.csv")
control <- vroom("../data/control_exp.csv")
head(lowinput)
head(control)

library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(dplyr)


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
      menuItem("Data Preview",   tabName = "data_tab",  icon = icon("table")),
      menuItem("Boxplots",       tabName = "plot_tab",  icon = icon("chart-bar")),
      menuItem("PCA",            tabName = "pca_tab",   icon = icon("project-diagram")),
      menuItem("t-SNE / UMAP",   tabName = "tsne_tab",  icon = icon("braille")),
      menuItem("Volcano Plot",   tabName = "volcano_tab", icon = icon("fire")),
      menuItem("Heatmap",        tabName = "heatmap_tab", icon = icon("th"))
    )
  ),
  dashboardBody(
    tabItems(
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
            )
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
      tabItem(
        tabName = "plot_tab",
        fluidRow(
          box(
            width = 4,
            title = "Plot Controls",
            selectInput(
              inputId = "plot_dataset",
              label = "Choose Dataset for Boxplot:",
              choices = c("control", "lowinput")
            ),
            uiOutput("compound_ui")  # We'll build the compound dropdown dynamically
          ),
          box(
            width = 8,
            title = "Boxplot(s)",
            plotOutput("boxplot_single"), 
            plotOutput("boxplot_both")
          )
        )
      ),
      
      ##------------------------------------------------------------------
      ##  TAB 3: PCA
      ##------------------------------------------------------------------
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
            checkboxInput("pca_scale",  "Scale", TRUE),
            helpText("This tab runs a simple PCA (using columns 2+). 
                      Samples are rows. The first column is sample ID.")
          ),
          box(
            width = 8,
            title = "PCA Plot",
            plotOutput("pca_plot")
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
            numericInput("perplexity", "t-SNE perplexity:", 10, min = 1, max = 50),
            helpText("For t-SNE, perplexity is relevant. For UMAP, we simply ignore it.")
          ),
          box(
            width = 8,
            title = "t-SNE / UMAP Plot",
            plotOutput("dr_plot")
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
            plotOutput("heatmap_plot")
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
  output$data_table <- renderDT({
    df <- dataset_list[[input$data_choice]]
    
    # 1) Remove the first column (Sample ID) and transpose the rest
    df_t <- as.data.frame(t(df[ , -1]))
    
    # 2) Use the original sample IDs as column names
    colnames(df_t) <- df[[1]]
    
    # 3) Create a "LipidName" column from row names (then remove row names)
    df_t$LipidName <- rownames(df_t)
    rownames(df_t) <- NULL  # Remove row names completely
    
    # 4) Move LipidName to the first column
    df_t <- df_t[, c(ncol(df_t), 1:(ncol(df_t)-1))]
    
    # Round numeric columns to two decimals
    num_cols <- 2:ncol(df_t)
    
    datatable(
      df_t,
      options = list(
        scrollX = TRUE,      # Horizontal scrolling for many samples
        pageLength = 10,     # Default rows per page
        lengthMenu = c(5,10,25,50)  # Row pagination options
      ),
      rownames = FALSE  # This removes duplicate row names!
    ) %>%
      formatRound(columns = num_cols, digits = 2)
  })
  
  
  
  
  ##======================================================================
  ##  TAB 2: Boxplots
  ##======================================================================
  # 1) Dropdown to pick the compound
  output$compound_ui <- renderUI({
    df <- dataset_list[[input$plot_dataset]]
    cols <- colnames(df)[-1]  # skip sample column
    selectInput(
      inputId = "compound_choice",
      label = "Choose a Lipid Compound:",
      choices = cols,
      selected = cols[1]
    )
  })
  
  # 2) Single-boxplot: Violin + Boxplot (One dataset)
  output$boxplot_single <- renderPlot({
    req(input$compound_choice)
    df <- dataset_list[[input$plot_dataset]]
    compound_col <- input$compound_choice
    
    validate(
      need(compound_col %in% colnames(df),
           "Selected compound not found in this dataset.")
    )
    
    # Build data frame for plotting
    df_plot <- df %>%
      select(Sample = 1, Value = all_of(compound_col))
    
    # Add Condition column
    df_plot$Condition <- input$plot_dataset
    
    ggplot(df_plot, aes(x = Condition, y = Value, fill = Condition)) +
      geom_violin(alpha = 0.5, width = 0.8) +   # Violin for density
      geom_boxplot(width = 0.1, outlier.shape = NA) +  # Boxplot inside
      labs(
        title = paste("Distribution of", compound_col, "in", input$plot_dataset),
        x = "Condition",
        y = "Value"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
  })
  
  # 3) Combined Boxplot with Violin (Control vs LowInput)
  output$boxplot_both <- renderPlot({
    req(input$compound_choice)
    
    # Identify the other dataset
    other_dataset <- setdiff(names(dataset_list), input$plot_dataset)
    compound_col <- input$compound_choice
    
    # If the compound doesn't exist in the other dataset, skip
    if (!compound_col %in% colnames(dataset_list[[other_dataset]])) {
      return(NULL)
    }
    
    # Main dataset
    df_main <- dataset_list[[input$plot_dataset]] %>%
      select(Sample = 1, Value = all_of(compound_col)) %>%
      mutate(Condition = input$plot_dataset)
    
    # Other dataset
    df_other <- dataset_list[[other_dataset]] %>%
      select(Sample = 1, Value = all_of(compound_col)) %>%
      mutate(Condition = other_dataset)
    
    # Combine both datasets
    df_combined <- bind_rows(df_main, df_other)
    
    ggplot(df_combined, aes(x = Condition, y = Value, fill = Condition)) +
      geom_violin(alpha = 0.5, width = 0.8) +   # Violin for density
      geom_boxplot(width = 0.1, outlier.shape = NA) +  # Boxplot inside
      labs(
        title = paste("Comparison of", compound_col, "across conditions"),
        x = "Condition",
        y = "Value"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right"
      )
  })
  
  
  
  
  ##======================================================================
  ##  TAB 3: PCA
  ##======================================================================
  output$pca_plot <- renderPlot({
    df <- dataset_list[[input$pca_dataset]]
    
    # Real code might need filtering out any non-numeric columns, NA handling, etc.
    numeric_df <- df[, -1]  # remove sample column
    numeric_df <- as.data.frame(sapply(numeric_df, as.numeric))  # ensure numeric
    
    # run PCA via base R's prcomp
    pca_res <- prcomp(numeric_df, center = input$pca_center, scale. = input$pca_scale)
    
    # Make a quick data frame for plotting first two PCs
    pca_data <- as.data.frame(pca_res$x)
    pca_data$Sample <- df[[1]]  # sample names
    
    ggplot(pca_data, aes(x = PC1, y = PC2, label = Sample)) +
      geom_point() +
      labs(
        title = paste("PCA -", input$pca_dataset),
        x = "PC1",
        y = "PC2"
      ) +
      theme_minimal()
  })
  
  ##======================================================================
  ##  TAB 4: t-SNE / UMAP
  ##======================================================================
  output$dr_plot <- renderPlot({
    df <- dataset_list[[input$dr_dataset]]
    numeric_df <- df[, -1]
    numeric_df <- as.data.frame(sapply(numeric_df, as.numeric))
    
    method <- input$dr_method
    
    if (method == "t-SNE") {
      # Run t-SNE
      library(Rtsne)
      tsne_out <- Rtsne(numeric_df, perplexity = input$perplexity, check_duplicates = FALSE)
      
      # Make a data frame
      dr_df <- data.frame(
        Dim1   = tsne_out$Y[,1],
        Dim2   = tsne_out$Y[,2],
        Sample = df[[1]]
      )
      
      ggplot(dr_df, aes(x = Dim1, y = Dim2, label = Sample)) +
        geom_point() +
        labs(
          title = paste("t-SNE -", input$dr_dataset),
          x = "Dim1",
          y = "Dim2"
        ) +
        theme_minimal()
      
    } else {
      # Run UMAP
      library(umap)
      umap_out <- umap(numeric_df)
      layout   <- umap_out$layout
      
      dr_df <- data.frame(
        Dim1   = layout[,1],
        Dim2   = layout[,2],
        Sample = df[[1]]
      )
      
      ggplot(dr_df, aes(x = Dim1, y = Dim2, label = Sample)) +
        geom_point() +
        labs(
          title = paste("UMAP -", input$dr_dataset),
          x = "UMAP1",
          y = "UMAP2"
        ) +
        theme_minimal()
    }
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
  output$heatmap_plot <- renderPlot({
    # We'll do a very minimal heatmap using pheatmap
    # For a more advanced approach, see 'ComplexHeatmap'
    library(pheatmap)
    
    df <- dataset_list[[input$heatmap_dataset]]
    
    # Use numeric columns only (excluding sample ID)
    mat <- df[, -1]
    mat <- as.data.frame(sapply(mat, as.numeric))
    
    # pheatmap can show row names, so let's set them to sample IDs
    rownames(mat) <- df[[1]]
    
    # Simple pheatmap
    pheatmap(mat,
             scale = "row",           # optional, normalizes each row
             show_rownames = TRUE,
             show_colnames = TRUE)
  })
}

shinyApp(ui, server)
