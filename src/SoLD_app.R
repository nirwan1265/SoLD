library(vroom)
## Read the files
lowinput <- vroom("../data/lowinput_exp.csv")
control <- vroom("../data/control_exp.csv")
head(lowinput)
head(control)

# install.packages("shiny")
# install.packages("shinydashboard")
# install.packages("DT")
# install.packages("ggplot2")
library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(dplyr)

## -- Example data: already loaded as 'control' and 'lowinput'
## -- Make sure you have them in your environment or source them beforehand.

ui <- dashboardPage(
  dashboardHeader(title = "SoLD"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data Preview", tabName = "data_tab", icon = icon("table")),
      menuItem("Boxplots", tabName = "plot_tab", icon = icon("chart-bar"))
    )
  ),
  dashboardBody(
    tabItems(
      
      ##--- TAB 1: Data Preview
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
      
      ##--- TAB 2: Boxplots
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
            # We will render one or two plots depending on compound presence
            plotOutput("boxplot_single"), 
            plotOutput("boxplot_both")
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  # Put the datasets in a named list to make them easy to switch
  dataset_list <- list(
    control  = control,
    lowinput = lowinput
  )
  
  #--------------------#
  #  TAB 1: Data Table #
  #--------------------#
  
  output$data_table <- renderDT({
    # Pick the chosen dataset
    df <- dataset_list[[input$data_choice]]
    # Show up to the first 15 columns
    df_subset <- df[, seq_len(min(ncol(df), 15))]
    datatable(
      df_subset,
      options = list(
        scrollX = TRUE,      # horizontal scroll
        pageLength = 10      # show 10 rows by default
      )
    )
  })
  
  #-------------------#
  #  TAB 2: Boxplots  #
  #-------------------#
  
  # Dynamically populate the "compound" dropdown with all column names
  # except the first (which is presumably samples)
  output$compound_ui <- renderUI({
    df <- dataset_list[[input$plot_dataset]]
    cols <- colnames(df)[-1]  # exclude the first column if it's "Samples"
    selectInput(
      inputId = "compound_choice",
      label = "Choose a Lipid Compound:",
      choices = cols,
      selected = cols[1]
    )
  })
  
  # Single boxplot for the chosen dataset
  output$boxplot_single <- renderPlot({
    req(input$compound_choice)  # make sure the user has chosen something
    
    # Grab the relevant dataset & columns
    df <- dataset_list[[input$plot_dataset]]
    compound_col <- input$compound_choice
    
    # Make sure the chosen compound exists in the dataset
    validate(
      need(compound_col %in% colnames(df),
           "Selected compound not found in this dataset.")
    )
    
    # Build a small data frame for plotting
    df_plot <- df %>%
      select(Sample = 1, Value = all_of(compound_col))
    
    # Single boxplot
    ggplot(df_plot, aes(x = "", y = Value)) +
      geom_boxplot(fill = "#69b3a2", alpha = 0.7) +
      labs(
        title = paste("Boxplot of", compound_col, "in", input$plot_dataset),
        x = "",
        y = "Value"
      ) +
      theme_minimal()
  })
  
  # Combined boxplot if the same compound also exists in the *other* dataset
  output$boxplot_both <- renderPlot({
    req(input$compound_choice)
    
    # Determine the "other" dataset
    other_dataset <- setdiff(names(dataset_list), input$plot_dataset)
    # e.g., if user picks "control", other_dataset is "lowinput"
    
    compound_col <- input$compound_choice
    # Check if 'compound_col' is in the *other* dataset
    if (!compound_col %in% colnames(dataset_list[[other_dataset]])) {
      # If not found, return NULL so no second boxplot is shown
      return(NULL)
    }
    
    # If found, combine data into one data frame
    df_main  <- dataset_list[[input$plot_dataset]] %>%
      select(Sample = 1, Value = all_of(compound_col)) %>%
      mutate(Condition = input$plot_dataset)
    
    df_other <- dataset_list[[other_dataset]] %>%
      select(Sample = 1, Value = all_of(compound_col)) %>%
      mutate(Condition = other_dataset)
    
    df_combined <- bind_rows(df_main, df_other)
    
    # Make a boxplot grouped by Condition
    ggplot(df_combined, aes(x = Condition, y = Value, fill = Condition)) +
      geom_boxplot(alpha = 0.7) +
      scale_fill_manual(values = c("#F08080", "#87CEFA")) +
      labs(
        title = paste("Comparison of", compound_col, "across both conditions"),
        x = "Condition",
        y = "Value"
      ) +
      theme_minimal()
  })
}

shinyApp(ui, server)
