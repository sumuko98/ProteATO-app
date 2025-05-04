source("global.R")
ui <- fluidPage(
  # CSS for button styling
  tags$style(HTML("
    .btn {
      border-radius: 5px;
      margin-bottom: 10px;
      font-weight: bold;
      max-width: 100%;      /* Ensure buttons fit within the column */
      white-space: normal;  /* Allow text to wrap if it's too long */
    }
    
    .sidebar {
      padding: 15px;
      border: 1px solid #ddd;
      background-color: #f9f9f9;
      border-radius: 5px;
    }

    .btn-info {
      background-color: #5bc0de;  /* Custom color for info buttons */
    }

    .btn-info:hover {
      background-color: #31b0d5;  /* Darker hover effect */
    }
  ")),
  
  titlePanel("Protein Analysis Tool"),
  
  navbarPage("Protein Analysis",
             id = "navbarPage",
             
             #  Tab for Data Input and Processing
             tabPanel("Data Input and Processing",
                      sidebarLayout(
                        sidebarPanel(
                          helpText("=> Please upload a .txt file containing proteomics data."),
                          helpText("=> Ensure that your data includes 'Gene names', 'Protein IDs', and 'LFQ intensity' columns."),
                          helpText("=> To do statistical comparisons, use the same condition name for all replicates for the experimental design file as in the example ."),
                          fileInput("txt_file", "Upload .txt Dataset", accept = ".txt"),
                          fileInput("design_file", "Upload Experimental Design (.csv)", accept = ".csv"),
                          actionButton("loadData", "Load Data"),  # Button to load data
                          tags$hr(),
                          
                          helpText("Number of blank & duplicated gene names found in your data:"),  # Simultaneous display
                          tags$div(
                            textOutput("duplicatesCount"),  # Display the count of duplicated gene names
                            style = "border: 1px solid #ddd; padding: 9px; margin-top: 5px; background-color: #f9f9f9; border-radius: 5px;"
                          ),
                          
                          # Sliding panel for filtering options using bsCollapse
                          bsCollapse(
                            id = "filterPanel",
                            
                            # Existing Filtering Panel
                            bsCollapsePanel("Filtering",
                                            radioButtons("filterChoice", "Choose Missing Value Filtering Method:",
                                                         choices = list("No Filter" = "no_filter",
                                                                        "Stringent: Proteins identified in all replicates of at least one condition" = "stringent",
                                                                        "Less Stringent: Proteins identified in all but one replicate for at least one condition" = "less_stringent"),
                                                         selected = "less_stringent"),
                                            actionButton("filterMissing", "Apply Missing Value Filtering", class = "btn-warning"),
                                            style = "info"
                            ),
                            
                            
                            # Normalization Panel with Yes/No option
                            bsCollapsePanel("Variance Stabilizing Normalization (VSN)",
                                            radioButtons("normalizationChoice", "Would you like to apply normalization?",
                                                         choices = list("Yes" = "yes", "No" = "no"),
                                                         selected = "no"),
                                            actionButton("applyNormalization", "Apply Normalization", class = "btn-info"),
                                            style = "warning"
                            ),
                            
                            # Imputation Panel
                            bsCollapsePanel("Imputation",
                                            # Replace radio buttons with searchable selectize input
                                            selectizeInput("imputationMethod", "Choose Imputation Method:",
                                                           choices = c("bpca", "knn", "QRILC", "MLE", "MinDet", "MinProb", 
                                                                       "man", "min", "zero", "mixed", "nbavg"),
                                                           selected = "MinProb",
                                                           options = list(
                                                             placeholder = 'Type to search for a method',
                                                             allowEmptyOption = FALSE
                                                           )),
                                            actionButton("applyImputation", "Apply Imputation", class = "btn-danger"),
                                            style = "danger"
                            )
                          ),
                          
                          # Next (Dynamic) Button:
                          uiOutput("nextButtonUI"),
                          tags$hr()
                        ),
                        
                        mainPanel(
                          textOutput("fileInfo"),  # Display file info or messages
                          verbatimTextOutput("lfqColumns"),  # Show LFQ column names
                          
                          # Conditional panel for displaying the image only if data is not yet loaded
                          conditionalPanel(
                            condition = "output.dataLoaded == false",  # Condition to show image when data is not loaded
                            img(src = "image_SE.png", height = "auto", width = "80%", alt = "Reference Table")  # Reference image
                          ),
                          
                          tableOutput("duplicatesTable"),  # Show duplicated gene names
                          textOutput("duplicatesInfo"),  # Display duplicate info
                          uiOutput("uniqueNamesInfo"),  # Show uniqueness conversion info
                          verbatimTextOutput("seSummary")  # Show summarized experiment
                        )
                      )
             ),
             
             tabPanel(
               "Venn Diagram",
               sidebarLayout(
                 sidebarPanel(
                   # Dropdown for selecting a condition for replicate comparison
                   selectInput("conditionVenn", "Select Condition for Replicate Comparison:", choices = NULL),
                   
                   # Dropdown for selecting multiple conditions for condition comparison
                   selectInput("conditionCompare", "Select Conditions for Condition Comparison:", 
                               choices = NULL, multiple = TRUE)  # Enable multiple condition selection
                 ),
                 mainPanel(
                   fluidRow(
                     # First column for Replicate Comparison Venn Diagram
                     column(
                       width = 6,
                       h4("Venn Diagram for Replicate Comparison"),
                       plotOutput("vennPlot", width = "100%", height = "600px"),
                       downloadButton("downloadVennReplicate", "Download Replicate Venn Diagram")
                     ),
                     # Second column for Condition Comparison Venn Diagram
                     column(
                       width = 6,
                       h4("Venn Diagram for Condition Comparison"),
                       plotOutput("vennPlotCondition", width = "100%", height = "600px"),
                       downloadButton("downloadVennCondition", "Download Condition Venn Diagram")
                     )
                   )
                 )
               )
             )
             
             ,
             
             
             # Tab for Overview
             tabPanel("Overview",
                      sidebarLayout(
                        sidebarPanel(
                          helpText("Check the Protein ID overlap between your samples."),
                          actionButton("IDCheck", "Check for Protein ID Overlap"),  # Check Protein ID overlap
                          tags$hr(),
                          helpText("Plot the number of identified proteins per sample."),
                          actionButton("plotProteinCounts", "Plot Protein Counts"),  # Button to plot protein counts
                          tags$hr(),
                          helpText("Normalize the data to reduce technical variability."),
                          actionButton("plotNormalizeData", "Show Normalization Plot", class = "btn-info"),
                          tags$hr(),
                          helpText("Visualize the missing values across samples."),
                          actionButton("plotMissingValues", "Show Missing Values Heatmap", class = "btn-info"),
                          tags$hr(),
                          helpText("Intensity distributions of the missing & non-missing values "),
                          actionButton("plotMissingValueIntensityDistribution", "Missing Value Intensity Distribution", class = "btn-info"),
                          tags$hr(),
                          helpText("Intensity distributions before and after imputation"),
                          actionButton("plotImputIntensityDistribution", "Imputation Intensity Distribution", class = "btn-info"),
                          tags$hr()
                        ),
                        mainPanel(
                          tabBox(
                            id = "plotsTab", width = 12,
                            
                            tabPanel("Overlap Plot", 
                                     plotOutput("OverlapPlot", width = "400px"),  # Display Protein ID overlap plot
                                     downloadButton("downloadOverlapPlot", "Download Overlap Plot as PDF")
                            ),
                            
                            tabPanel("Protein Counts", 
                                     plotOutput("ProteinCountsPlot", width = "600px"),  # Display number of proteins per sample
                                     downloadButton("downloadProteinCountsPlot", "Download Protein Counts Plot as PDF")
                            ),
                            
                            tabPanel("Normalized Data", 
                                     plotOutput("NormalizedDataPlot", width = "600px"),  # Pre and post normalization visualization
                                     downloadButton("downloadNormalizedDataPlot", "Download Normalized Data Plot as PDF")
                            ),
                            
                            tabPanel("Missing Values", 
                                     plotOutput("MissingValuesHeatmap", width = "600px"),  # Display missing values heatmap
                                     downloadButton("downloadMissingValuesHeatmap", "Download Missing Values Heatmap as PDF")
                            ),
                            
                            tabPanel("Intensity Distribution", 
                                     fluidRow(
                                       column(width = 8,  
                                              plotOutput("MissingValueIntensityDistributionPlot", height = "600px"),
                                              downloadButton("downloadMissingValueIntensityPlot", "Download Intensity Distribution as PDF")
                                       )
                                     )
                            ),
                            
                            
                            tabPanel("Imputation Effect", 
                                     fluidRow(
                                       column(width = 8,
                                              plotOutput("ImputationEffectPlot", height = "600px"),  # Display imputation effect plot
                                              downloadButton("downloadImputationEffectPlot", "Download Imputation Effect Plot as PDF")
                                       )
                                     )
                            )
                          )
                        )
                      )
             ),
             
             
             
             tabPanel("Differential Enrichment Analysis",
                      sidebarLayout(
                        sidebarPanel(
                          # Contrast Type Selection
                          radioButtons("contrastType", "Select Contrast Type:",
                                       choices = list("Control" = "control", 
                                                      "All Comparisons" = "all", 
                                                      "Manual" = "manual"),
                                       selected = "control"),
                          
                          # Control Sample Dropdown (conditionally shown)
                          conditionalPanel(
                            condition = "input.contrastType == 'control'",
                            selectInput("controlSample", "Select Control Condition:",
                                        choices = NULL)  # Choices will be dynamically populated
                          ),
                          
                          # Manual Contrast Input (conditionally shown)
                          conditionalPanel(
                            condition = "input.contrastType == 'manual'",
                            textInput("manualContrast", "Specify Contrasts to Test (comma-separated):",
                                      placeholder = "e.g., Ctrl_vs_Ubi4, Ubi4_vs_Ubi6")
                          ),
                          
                          # Significance Thresholds
                          numericInput("pValueCutoff", "Adjusted p-value cutoff (Alpha Level):", 
                                       value = 0.05, min = 0, max = 1, step = 0.01),
                          
                          numericInput("log2FCCutoff", "Log2 Fold Change Cutoff:", 
                                       value = 1, min = 0),
                          
                          # Submit Button for Differential Analysis
                          actionButton("submitDifferentialAnalysis", "Submit", class = "btn-primary"),
                          
                          tags$hr(),  # Horizontal line to separate sections
                          
                          # Section for user-defined cutoffs for the protein comparison
                          h4("Proteins to Compare"),
                          
                          
                          # Section for protein input
                          h4("Select Proteins for Bar Plot"),
                          
                          # Input field for entering protein names separated by commas
                          textInput("proteinNames", "Enter at least two different Protein Names:", placeholder = "e.g., USP15, IKBKG"),
                          
                          # Submit button for generating bar plots for selected proteins
                          actionButton("submitProteins", "Generate Bar Plot", class = "btn-success"),
                          
                          
                          h4("Select P-Value Column for Histogram"),
                          uiOutput("pvalSelector"),
                        ),
                        mainPanel(
                          # Placeholder for the results
                          h4("Results will be displayed here after analysis"),
                          
                          # UI output for the contrast selector 
                          uiOutput("contrastSelector"),
                          
                          # Volcano plot output
                          withSpinner(plotOutput("volcanoPlot")),
                          
                          downloadButton("downloadVolcanoData", "Download Volcano Plot Data"),
                          
                          tags$hr(),  # Horizontal line to separate sections
                          
                          # Plot output for the bar plot of selected proteins
                          # Placeholder for the results
                          h4("Bar Plot for Selected Proteins"),
                          
                          # Display the loading spinner only when rendering the bar plot after the button click
                          uiOutput("proteinBarplotUI"),  # Use uiOutput for conditional display
                          
                          
                          #### 
                          tags$hr(),
                          
                          # **Added Histogram for P-Value Distribution**
                          h4("P-Value Distribution Histogram"),  # Highlighted change
                          withSpinner(plotOutput("pvalHistogram")),  # Display the histogram dynamically (new plot output)
                          
                          tags$hr(),
                          
                          tags$hr(),
                          # Placeholder for results count
                          h4(""),
                          textOutput("numSignificantProteins"),
                          
                          tags$hr(),  # Optional separator line
                          
                          # Placeholder for results table
                          h4("Results Table"),
                          DT::dataTableOutput("significantProteinsTable")  # Display results as a DataTable
                          
                        )
                      )
             ),
             tabPanel("PCA Plot",
                      sidebarLayout(
                        sidebarPanel(
                          h4("PCA Analysis"),
                          helpText("The PCA plot provides a high-level overview of the variance in the dataset."),
                          numericInput("n_pca", "Top N most variable proteins:", value = 500, min = 1, max = 1000),
                          numericInput("x_pca", "Select PC for X-axis:", value = 1, min = 1),
                          numericInput("y_pca", "Select PC for Y-axis:", value = 2, min = 1),
                          actionButton("submitPCA", "Generate PCA Plot", class = "btn-primary"),
                          
                          # Enhanced informative message below the Generate PCA Plot button
                          htmlOutput("numProteinsInfo"),  # Use htmlOutput to allow HTML styling
                          
                          tags$hr(),  # Separator line
                          downloadButton("downloadDep", "Download DEP file as Excel")  # Download button for dep object
                        ),
                        mainPanel(
                          uiOutput("pcaPlotUI"),  # Placeholder for dynamic plot with spinner
                          downloadButton("downloadPCAPlot", "Download PCA Plot as PDF"),
                          
                          #Frequency Plot of Significant Proteins
                          tags$hr(),
                          h4("Frequency Plot of Significant Proteins by Condition"),
                          plotOutput("freqPlot"),  # Output for the frequency plot
                          downloadButton("downloadFreqPlot", "Download Frequency Plot as PDF")  # Download button for frequency plot
                        )
                      )
             ),
             
             tabPanel("Gene Ontology (GO) Analysis",
                      sidebarLayout(
                        sidebarPanel(
                          # Dropdown for selecting the ontology with GO IDs
                          selectInput("ontologyType", "Select Ontology Type:",
                                      choices = c("Biological Process" = "GO:0008150", 
                                                  "Molecular Function" = "GO:0003674", 
                                                  "Cellular Component" = "GO:0005575"),
                                      selected = "GO:0008150"),
                          
                          # Dropdown for selecting the organism using long_name and taxon_id
                          selectInput(
                            inputId = "organism",
                            label = "Select Organism:",
                            choices = NULL,  # Will be dynamically populated with organism long names and taxon IDs
                            selected = "3702"  # Default taxon ID (e.g., Arabidopsis) can be adjusted as needed
                          ),
                          
                          # Dropdown for selecting the comparison (condition/contrast)
                          selectInput("conditionContrast", "Select Condition/Contrast:", 
                                      choices = NULL),  # Dynamically populated based on design file
                          
                          # Action button to submit GO analysis
                          actionButton("submitGO", "Run GO Analysis"),
                          
                          # Button to download the GO results
                          downloadButton("downloadGOResults", "Download GO Results"),
                          
                          hr(),
                          
                          # Output for Mapped and Unmapped IDs
                          h4("Mapped Gene IDs"),
                          verbatimTextOutput("mappedIDs"),
                          h4("Unmapped Gene IDs"),
                          verbatimTextOutput("unmappedIDs")
                        ),
                        
                        # Main panel to display results
                        mainPanel(
                          h3("GO Enrichment Results"),
                          DT::dataTableOutput("goResultsTable"),  # Table output for GO results
                          hr(),  # Add a horizontal line for separation
                          h3("GO Enrichment Plot"),
                          plotOutput("dotPlot")  # Placeholder for the dot plot
                        )
                      )
             ),
             
             tabPanel("KEGG",
                      sidebarLayout(
                        sidebarPanel(
                          h4("KEGG Analysis Options"),
                          
                          # Organism code input with instructions for lookup
                          textInput("keggOrganism", "Type in Organism Code:", 
                                    placeholder = "Enter KEGG code (e.g., 'ath' for Arabidopsis Thaliana)"),
                          tags$p("To find your organism code, copy and paste this link in your browser:",
                                 style = "font-size: 14px; color: #333; margin-top: 10px;"),
                          tags$p("https://www.genome.jp/kegg/catalog/org_list.html", 
                                 style = "font-weight: bold; color: #0073e6; font-size: 14px;"),
                          
                          # Contrast selection and analysis button
                          selectInput("conditionContrastKEGG", "Select Condition/Contrast:", choices = NULL),
                          actionButton("runKEGGAnalysis", "Run KEGG Analysis", class = "btn-primary")
                        ),
                        
                        mainPanel(
                          h4("KEGG IDs for Significant Proteins"),
                          verbatimTextOutput("keggIDs"),
                          tags$hr(),
                          
                          h4("Regulation"),
                          verbatimTextOutput("RegulationIDs"),
                          
                          h4("KEGG Pathway Enrichment Results"),
                          DT::dataTableOutput("keggEnrichmentTable"),
                          tags$hr(),
                          
                          h4("KEGG Enrichment Network Plot"),
                          uiOutput("keggCnetPlotUI")
                        )
                      )
             ),
             
             tabPanel("ANOVA",
                      sidebarLayout(
                        sidebarPanel(
                          #h4(""),
                          
                          # Checkbox Group for selecting conditions
                          checkboxGroupInput(
                            inputId = "anova_conditions",
                            label = "Select Conditions for ANOVA:",
                            choices = NULL,  # Dynamically populated in the server
                            selected = NULL  # Optionally preselect some conditions
                          ),
                          
                          # Action Button to trigger ANOVA
                          actionButton("run_anova", "Run ANOVA"),
                          
                          br(),
                          hr(),
                          
                          # Optional Message or Help Text
                          helpText(
                            HTML("Select the conditions you want to include in the ANOVA analysis.<br>
                            NOTE: You have to select at leat 3 groups to be able to run the analysis.<br>
                            Click 'Run ANOVA' to compute the results."))
                        ),
                        
                        mainPanel(
                          tabBox(
                            id = "anova_tab", width = 12,
                            # Table panel
                            tabPanel("ANOVA DEPs Table",
                                     DT::dataTableOutput("anova_deps_table"),
                                     downloadButton("download_anova_deps_table", "Download ANOVA results")
                            ),
                            # Heatmap panel
                            tabPanel("ANOVA Heatmap", 
                                     plotOutput("anovaHeatmap"),
                                     downloadButton("download_anova_heatmap", "Download Heatmap as PDF")
                            )
                          )
                        ) # Close ANOVA main panel
                      )
             ) # Close the whole ANOVA tabpanel
  )
)



server <- function(input, output, session) {
  options(shiny.maxRequestSize = 30 * 1024^2)  # Increase file size limit to 30 MB
  
  
  # Reactive value to track whether data has been loaded and whether normalization is applied
  dataLoaded <- reactiveVal(FALSE)  # Initially set to FALSE to show the image at the start
  normalizedData <- reactiveVal(FALSE)  # Tracks if normalization has been applied
  imputedDataAvailable <- reactiveVal(FALSE)  # Tracks if imputation is complete 
  
  # Helper function to choose between normalized or filtered data
  get_data_for_imputation <- function() {
    if (normalizedData()) {
      return(normalized_se())  # Use normalized data if available
    } else {
      return(filtered_se())  # Otherwise, use the filtered data
    }
  }
  
  # Event reactive for reading the uploaded dataset and filtering contaminants
  data <- eventReactive(input$loadData, {
    req(input$txt_file)
    # Read the dataset and filter contaminants
    read.delim(input$txt_file$datapath, header = TRUE) %>%
      filter(Reverse != "+", Potential.contaminant != "+", Only.identified.by.site != "+")
  })
  
  # Event reactive for reading the experimental design file
  design <- eventReactive(input$loadData, {
    req(input$design_file)
    read.csv(input$design_file$datapath, header = TRUE)
  })
  
  # Populate the control sample dropdown based on conditions from the design file
  observe({
    req(design())  # Ensure design data is available
    
    # Extract unique condition names from the design file
    condition_names <- unique(design()$condition)
    
    # Ensure condition_names is not NULL and has at least one condition
    if (!is.null(condition_names) && length(condition_names) > 0) {
      
      # Update the control sample dropdown
      updateSelectInput(session, "controlSample", choices = condition_names)
      
      # Update the replicate comparison dropdown (conditionVenn)
      updateSelectInput(session, "conditionVenn", choices = condition_names)
      
      # Update the condition comparison dropdown (conditionCompare)
      updateSelectInput(session, "conditionCompare", choices = condition_names)
      
      # Update the selectInput for 'go_condition' with the available conditions
      updateSelectInput(session, "go_condition", 
                        choices = condition_names, 
                        selected = condition_names[1])  # Default to the first condition
      
      # Update the checkboxGroupInput for ANOVA conditions
      updateCheckboxGroupInput(session, "anova_conditions", 
                               choices = condition_names, 
                               selected = condition_names)  # Default to all conditions selected
    }
  })
  
  
  # Summarized Experiment object
  data_se <- reactiveVal(NULL)
  
  # Filtered Summarized Experiment object
  filtered_se <- reactiveVal(NULL)
  
  # Normalized Summarized Experiment object
  normalized_se <- reactiveVal(NULL)
  
  # Imputed Summarized Experiment object
  imputed_se <- reactiveVal(NULL)
  
  # Differential Results object
  diff_results <- reactiveVal(NULL)
  
  # Reactive value to hold the differential analysis results
  dep <- reactiveVal(NULL)  # Holds the significant proteins after marking
  
  # Reactive values to store the cutoffs (p-value and log2 fold change) provided by the user
  userCutoffs <- reactiveValues(
    alpha = NULL,
    lfc = NULL)
  
  # Reactive value for storing ANOVA results
  anova_results <- reactiveVal(NULL)
  
  
  # Observe when Load Data button is clicked
  observeEvent(input$loadData, {
    req(input$txt_file, input$design_file)
    
    # Verify experimental design columns
    required_cols <- c("label", "condition", "replicate")
    if (!all(required_cols %in% colnames(design()))) {
      shinyalert::shinyalert(
        title = "Error",
        text = "Experimental design must contain 'label', 'condition', and 'replicate'.",
        type = "error"
      )
      return()
    }
    
    # Extract LFQ intensity column names
    lfq_columns <- colnames(data())[grepl("^LFQ.intensity", colnames(data()), ignore.case = TRUE)]
    if (length(lfq_columns) == 0) {
      shinyalert::shinyalert(
        title = "Error",
        text = "No LFQ intensity columns found in the dataset.",
        type = "error"
      )
      return()
    }
    
    
    
    ##### Venn #####
    
    ##### Venn for Condition Comparison #####
    getDetectedProteins <- function(data, lfq_columns, selected_labels) {
      # Construct the expected LFQ column names using the labels
      expected_lfq_cols <- paste0("LFQ.intensity.", selected_labels)
      
      # Find LFQ columns that match the constructed expected names
      selected_lfq_cols <- lfq_columns[lfq_columns %in% expected_lfq_cols]
      
      # Handle case where no LFQ columns match
      if (length(selected_lfq_cols) == 0) {
        return(character(0))
      }
      
      # Extract the LFQ data for the selected columns
      lfq_data <- data[, selected_lfq_cols, drop = FALSE]
      
      # Identify proteins with at least one non-zero value across selected LFQ columns
      detected_proteins <- rownames(lfq_data)[rowSums(lfq_data > 0, na.rm = TRUE) > 0]
      
      return(detected_proteins)
    }
    
    ##### Generalized Function to Render the Venn Diagram #####
    renderVennDiagram <- function(protein_lists, category_names) {
      # Remove empty lists (if any groups have no detected proteins)
      valid_indices <- sapply(protein_lists, function(x) length(x) > 0)
      protein_lists <- protein_lists[valid_indices]
      category_names <- category_names[valid_indices]
      
      # Ensure there are at least 2 groups with detected proteins
      if (length(protein_lists) < 2) {
        showNotification("Not enough groups with detected proteins to generate a Venn diagram.", type = "warning")
        return(NULL)
      }
      
      # Dynamically switch visualization based on the number of groups
      num_groups <- length(protein_lists)
      
      if (num_groups <= 5) {
        # Use VennDiagram for up to 5 groups
        venn_diagram <- VennDiagram::venn.diagram(
          x = protein_lists,
          category.names = category_names,
          filename = NULL,  # Render directly without saving to a file
          output = TRUE,
          disable.logging = TRUE,
          log = FALSE,
          fill = RColorBrewer::brewer.pal(num_groups, "Set3")[1:num_groups],
          alpha = rep(0.5, num_groups),
          lty = "solid",
          lwd = 2,
          col = rep("black", num_groups),
          cat.cex = 1.2,
          cat.fontface = "bold",
          cat.col = RColorBrewer::brewer.pal(num_groups, "Dark2")[1:num_groups],
          #cat.pos = seq(-45, 45, length.out = num_groups),
          cat.dist = rep(0.05, num_groups),
          cex = 1.2,
          fontface = "bold",
          fontfamily = "sans",
          background = "white"
        )
        grid::grid.draw(venn_diagram)
      } else {
        # Use UpSetR for more than 5 groups
        # Create a binary matrix for UpSetR visualization
        upset_data <- do.call(cbind, lapply(protein_lists, function(set) {
          as.integer(unique(unlist(protein_lists)) %in% set)
        }))
        colnames(upset_data) <- category_names
        rownames(upset_data) <- unique(unlist(protein_lists))
        
        # Convert to data frame and plot using UpSetR
        upset(as.data.frame(upset_data), sets = category_names, order.by = "freq", main.bar.color = "steelblue", sets.bar.color = "tomato")
      }
    }
    
    
    ##### Render Venn Diagram for Replicate Comparison #####
    output$vennPlot <- renderPlot({
      req(input$conditionVenn, data(), lfq_columns, design())  # Ensure inputs are available
      
      # Get replicate numbers for the selected condition
      condition <- input$conditionVenn
      replicate_numbers <- unique(design()$replicate[design()$condition == condition])
      
      # Check if replicates exist for the selected condition
      if (length(replicate_numbers) < 2) {
        showNotification("Not enough replicates for this condition.", type = "warning")
        return(NULL)
      }
      
      # Extract detected proteins for each replicate
      protein_lists <- lapply(replicate_numbers, function(rep) {
        selected_labels <- design()$label[design()$condition == condition & design()$replicate == rep]
        getDetectedProteins(data(), lfq_columns, selected_labels)
      })
      
      # Remove empty lists (if any replicates have no detected proteins)
      protein_lists <- Filter(function(x) length(x) > 0, protein_lists)
      
      # Name each list element for labeling in the Venn diagram
      category_names <- paste("Replicate", replicate_numbers)
      
      # Render the Venn diagram using the generalized function
      renderVennDiagram(protein_lists, category_names)
    }, width = 600, height = 600, res = 96)
    
    ##### Download Venn Diagram for Replicate Comparison #####
    output$downloadVennReplicate <- downloadHandler(
      filename = function() {
        paste("venn_replicate_comparison_", Sys.Date(), ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width = 8, height = 8)  # Open a PDF device
        condition <- input$conditionVenn
        replicate_numbers <- unique(design()$replicate[design()$condition == condition])
        
        # Ensure valid replicates
        if (length(replicate_numbers) < 2) {
          showNotification("Not enough replicates for this condition.", type = "warning")
          return(NULL)
        }
        
        # Extract detected proteins for each replicate
        protein_lists <- lapply(replicate_numbers, function(rep) {
          selected_labels <- design()$label[design()$condition == condition & design()$replicate == rep]
          getDetectedProteins(data(), lfq_columns, selected_labels)
        })
        
        # Remove empty lists
        protein_lists <- Filter(function(x) length(x) > 0, protein_lists)
        category_names <- paste("Replicate", replicate_numbers)
        
        # Draw the Venn diagram into the PDF
        renderVennDiagram(protein_lists, category_names)
        dev.off()  # Close the PDF device
      }
    )
    
    
    
    ##### Render Venn Diagram for Condition Comparison #####
    output$vennPlotCondition <- renderPlot({
      req(input$conditionCompare, data(), lfq_columns, design())  # Ensure inputs are available
      
      # Get the list of selected conditions
      selected_conditions <- input$conditionCompare
      
      # Extract detected proteins for each condition
      protein_lists <- lapply(selected_conditions, function(condition) {
        selected_labels <- design()$label[design()$condition == condition]
        getDetectedProteins(data(), lfq_columns, selected_labels)
      })
      
      # Remove empty lists (if any conditions have no detected proteins)
      valid_indices <- sapply(protein_lists, function(x) length(x) > 0)
      protein_lists <- protein_lists[valid_indices]
      valid_conditions <- selected_conditions[valid_indices]
      
      # Check if we have at least two valid conditions to compare
      if (length(valid_conditions) < 2) {
        showNotification("Not enough valid conditions with detected proteins to generate a Venn diagram.", type = "warning")
        return(NULL)
      }
      
      # Render the Venn diagram using the generalized function
      renderVennDiagram(protein_lists, valid_conditions)
      
    }, width = 600, height = 600, res = 96)
    ##### Download Venn Diagram for Condition Comparison #####
    output$downloadVennCondition <- downloadHandler(
      filename = function() {
        paste("venn_condition_comparison_", Sys.Date(), ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width = 8, height = 8)  # Open a PDF device
        selected_conditions <- input$conditionCompare
        
        # Extract detected proteins for each condition
        protein_lists <- lapply(selected_conditions, function(condition) {
          selected_labels <- design()$label[design()$condition == condition]
          getDetectedProteins(data(), lfq_columns, selected_labels)
        })
        
        # Remove empty lists (if any conditions have no detected proteins)
        valid_indices <- sapply(protein_lists, function(x) length(x) > 0)
        protein_lists <- protein_lists[valid_indices]
        valid_conditions <- selected_conditions[valid_indices]
        
        # Ensure valid conditions
        if (length(valid_conditions) < 2) {
          showNotification("Not enough valid conditions with detected proteins to generate a Venn diagram.", type = "warning")
          return(NULL)
        }
        
        # Draw the Venn diagram into the PDF
        renderVennDiagram(protein_lists, valid_conditions)
        dev.off()  # Close the PDF device
      }
    )
    ##### Venn #####
    
    
    ######
    ## Display the LFQ columns
    ##lfq_after <- sub("^LFQ.intensity.", "", lfq_columns)
    ##output$lfqColumns <- renderText({
    ##paste("LFQ Intensity Labels:\n", paste(lfq_after, collapse = "\n"), sep = "")
    ##})
    ######
    
    
    # Count and display duplicated gene names
    duplicate_count <- sum(duplicated(data()$Gene.names))
    output$duplicatesCount <- renderText({
      paste(duplicate_count)
    })
    
    # Step 1: Check for duplicates and display information
    has_duplicates <- data()$Gene.names %>% duplicated() %>% any()
    if (!has_duplicates) {
      # Shinyalert to show no duplicates found
      shinyalert::shinyalert(
        title = "No Duplicates",
        text = "No duplicated Gene names found.",
        type = "success"
      )
      output$duplicatesTable <- renderTable(NULL)
    } else {
      # Shinyalert to show duplicates found
      shinyalert::shinyalert(
        title = "Duplicates Found",
        text = "Duplicated Gene names found. Automatically converting to unique names.",
        type = "warning"
      )
      
      dup_table <- data() %>%
        group_by(Gene.names) %>%
        summarize(frequency = n()) %>%
        arrange(desc(frequency)) %>%
        filter(frequency > 1)
      
      output$duplicatesTable <- renderTable(dup_table)
    }
    
    # Signal that data has been loaded
    dataLoaded(TRUE)  # Set to TRUE to hide the image once data is loaded
    
    # Hide the image and update file information
    output$fileInfo <- renderText({
      paste("Duplicated gene names before the automatic conversion:")
    })
    output$lfqColumns <- renderText(NULL)  # Hide the LFQ columns text
    
    # Step 2: Convert Gene names to unique using Protein.IDs if duplicates exist
    data_unique <- make_unique(data(), "Gene.names", "Protein.IDs", delim = ";")
    output$uniqueNamesInfo <- renderUI({
      HTML('<p style="font-size: 18px; font-weight: bold; color: black;">
        Please proceed with the <span style="color: blue;">Filtering</span> step.
       </p>')
    })
    
    # Step 3: Create Summarized Experiment
    LFQ_columns <- grep("^LFQ.", colnames(data_unique))
    data_se(make_se(data_unique, LFQ_columns, design()))
    
    # Step 4: Show confirmation modal
    showModal(modalDialog(
      title = "Submission Successful",
      "Your data has been successfully processed and submitted!",
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  # Make dataLoaded available in the UI for conditional display of the image
  output$dataLoaded <- reactive({
    dataLoaded()
  })
  outputOptions(output, "dataLoaded", suspendWhenHidden = FALSE)
  
  # Observe when Apply Missing Value Filtering button is clicked
  observeEvent(input$filterMissing, {
    req(data_se())  # Ensure that a SummarizedExperiment object is available
    
    # Determine which filtering method was selected
    if (input$filterChoice == "stringent") {
      # Stringent filtering: Proteins identified in all replicates of at least one condition
      filtered_data <- filter_missval(data_se(), thr = 0)
    } else if (input$filterChoice == "less_stringent") {
      # Less stringent filtering: Proteins identified in 2 out of 3 replicates of at least one condition
      filtered_data <- filter_missval(data_se(), thr = 1)
    } else if (input$filterChoice == "no_filter") {
      # No filtering: use the raw data directly
      filtered_data <- data_se()
    }
    
    # Update the filtered SummarizedExperiment object
    filtered_se(filtered_data)
    
    # Show confirmation that the filtering has been applied
    showModal(modalDialog(
      title = "Missing Value Filtering Applied",
      paste("Filtering based on", input$filterChoice, "method is complete."),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  # Observe when Apply Normalization button is clicked
  observeEvent(input$applyNormalization, {
    req(filtered_se())  # Ensure filtered data is available
    
    if (input$normalizationChoice == "yes") {
      # Normalize the data using VSN method
      normalized_se(normalize_vsn(filtered_se()))
      normalizedData(TRUE)  # Set to TRUE when normalization is applied
      
      showModal(modalDialog(
        title = "Normalization Applied",
        "The data has been successfully normalized using VSN.",
        easyClose = TRUE,
        footer = NULL
      ))
    } else {
      showModal(modalDialog(
        title = "Normalization Skipped",
        "Normalization was not applied as per your selection.",
        easyClose = TRUE,
        footer = NULL
      ))
    }
  })
  # Helper function to choose between normalized or filtered data
  get_data_for_imputation <- function() {
    if (!is.null(normalized_se())) {
      return(normalized_se())  # Use normalized data if available
    } else if (!is.null(filtered_se())) {
      return(filtered_se())  # Otherwise, use the filtered data
    } else {
      stop("No data available for imputation")
    }
  }
  # Observe when Apply Imputation button is clicked
  observeEvent(input$applyImputation, {
    req(filtered_se())  # Ensure filtered data is available
    
    # Get the correct data for imputation (normalized or filtered)
    data_for_imputation <- get_data_for_imputation()
    
    # Apply imputation based on the user's choice of method
    imputed_data <- impute(data_for_imputation, input$imputationMethod)
    
    # Store the imputed data in the reactive value
    imputed_se(imputed_data)
    
    # Print column names of the imputed data to the console for debugging
    observeEvent(imputed_se(), {
      req(imputed_se())  # Ensure the imputed_se object is available
      cat("Column names of imputed_se:\n")
      cat(colnames(imputed_se()), sep = ", ")
      cat("\n")  # Add a new line for better readability
    })
    
    # Set imputation flag to TRUE
    imputedDataAvailable(TRUE)
    
    # Show confirmation that the imputation was applied
    showModal(modalDialog(
      title = "Imputation Applied",
      paste("Imputation using", input$imputationMethod, "method is complete."),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  # Render the "Next" button UI
  output$nextButtonUI <- renderUI({
    if (imputedDataAvailable()) {
      actionButton("Next", "Next", class = "btn-primary")  # Enabled "Next" button
    } else {
      actionButton("Next", "Next", class = "btn-primary", disabled = TRUE)  # Disabled "Next" button
    }
  })
  observeEvent(input$Next, {
    shinyalert::shinyalert(
      title = "Proceeding",
      text = "Now you can proceed to the Differential Enrichment Analysis section!",
      type = "info"
    )
  })
  # Helper function for downloading plots
  downloadPlot <- function(outputId, plotExpr, filePrefix) {
    output[[outputId]] <- downloadHandler(
      filename = function() {
        paste(filePrefix, Sys.Date(), ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width = 8, height = 6)  # Open a PDF device with specified dimensions
        print(plotExpr())  # Use the reactive expression passed in
        dev.off()  # Close the PDF device
      }
    )
  }
  
  
  #### Normalization plot
  
  # Reactive expression for the pre- and post-normalization plot
  normalizedDataPlot <- reactive({
    req(filtered_se(), normalized_se())  # Ensure both filtered and normalized data are available
    
    # Check if normalization has been applied
    if (normalizedData()) {
      plot_normalization(filtered_se(), normalized_se())  # Plot the normalization
    } else {
      showModal(modalDialog(
        title = "Normalization Not Applied",
        "Please apply normalization first.",
        easyClose = TRUE,
        footer = NULL
      ))
      return(NULL)  # Avoid rendering the plot if normalization hasn't been applied
    }
  })
  
  # Observe when Normalize Data button is clicked to show pre- and post-normalization plot
  observeEvent(input$plotNormalizeData, {
    req(filtered_se())  # Ensure filtered data is available
    
    # Render the normalized data plot
    output$NormalizedDataPlot <- renderPlot({
      normalizedDataPlot()  # Use the reactive expression to render the plot
    })
    
    # Switch to the "Normalized Data" tab after the plot is rendered
    updateTabsetPanel(session, "plotsTab", selected = "Normalized Data")
  })
  
  # Download handler for the normalized data plot
  downloadPlot("downloadNormalizedDataPlot", normalizedDataPlot, "normalized_data_plot")
  
  #### Normalization plot
  
  
  
  ##Overview##
  
  # Event Reactive Expression for Overlap Plot
  overlapPlot <- eventReactive(input$IDCheck, {
    req(data_se())  # Ensure data_se() is available and non-NULL
    plot_frequency(data_se())  # Generate the overlap plot only when the button is clicked
  })
  
  # Render the Overlap Plot in the UI
  output$OverlapPlot <- renderPlot({
    overlapPlot()  # Use the eventReactive expression for the plot
  })
  
  # Switch to the Overlap Plot tab
  observeEvent(input$IDCheck, {
    updateTabsetPanel(session, "plotsTab", selected = "Overlap Plot")
  })
  
  
  # Download handler for Overlap Plot
  downloadPlot("downloadOverlapPlot", overlapPlot, "overlap_plot")
  
  
  # Reactive expression for Protein Counts Plot
  proteinCountsPlot <- reactive({
    req(filtered_se())  # Ensure filtered_se() is available
    plot_numbers(filtered_se())  # Generate the plot
  })
  
  # Observe when Plot Protein Counts button is clicked
  observeEvent(input$plotProteinCounts, {
    req(filtered_se())  # Ensure data is available
    
    # Render the Protein Counts Plot in the UI
    output$ProteinCountsPlot <- renderPlot({
      proteinCountsPlot()  # Use the reactive expression for the plot
    })
    
    # Switch to the Protein Counts tab
    updateTabsetPanel(session, "plotsTab", selected = "Protein Counts")
  })
  
  # Download handler for Protein Counts Plot
  downloadPlot("downloadProteinCountsPlot", proteinCountsPlot, "protein_counts_plot")
  
  
  # Reactive expression for Missing Values Heatmap
  missingValuesHeatmapPlot <- reactive({
    req(filtered_se())  # Ensure filtered_se() is available
    plot_missval(filtered_se())  # Generate the missing values heatmap
  })
  
  # Observe when Missing Values Heatmap button is clicked
  observeEvent(input$plotMissingValues, {
    req(filtered_se())  # Ensure data is available
    
    # Render the Missing Values Heatmap in the UI
    output$MissingValuesHeatmap <- renderPlot({
      missingValuesHeatmapPlot()  # Use the reactive expression for the plot
    })
    
    # Switch to the Missing Values tab
    updateTabsetPanel(session, "plotsTab", selected = "Missing Values")
  })
  
  # Download handler for Missing Values Heatmap
  downloadPlot("downloadMissingValuesHeatmap", missingValuesHeatmapPlot, "missing_values_heatmap")
  
  
  # Reactive expression for Intensity Distribution Plot
  intensityDistributionPlot <- reactive({
    req(filtered_se())  # Ensure filtered_se() is available
    plot_detect(filtered_se())  # Generate the intensity distribution plot
  })
  
  # Observe when Intensity Distribution button is clicked
  observeEvent(input$plotMissingValueIntensityDistribution, {
    req(filtered_se())  # Ensure data is available
    
    # Render the Intensity Distribution Plot in the UI
    output$MissingValueIntensityDistributionPlot <- renderPlot({
      intensityDistributionPlot()  # Use the reactive expression for the plot
    })
    
    # Switch to the Intensity Distribution tab
    updateTabsetPanel(session, "plotsTab", selected = "Intensity Distribution")
  })
  
  # Download handler for Intensity Distribution Plot
  downloadPlot("downloadMissingValueIntensityPlot", intensityDistributionPlot, "intensity_distribution_plot")
  
  
  # Reactive expression for Imputation Effect Plot
  imputationEffectPlot <- reactive({
    req(normalized_se(), imputed_se())  # Ensure both normalized and imputed data are available
    plot_imputation(normalized_se(), imputed_se())  # Generate the imputation effect plot
  })
  
  # Observe when Imputation Effect button is clicked
  observeEvent(input$plotImputIntensityDistribution, {
    req(normalized_se(), imputed_se())  # Ensure data is available
    
    # Render the Imputation Effect Plot in the UI
    output$ImputationEffectPlot <- renderPlot({
      imputationEffectPlot()  # Use the reactive expression for the plot
    })
    
    # Switch to the Imputation Effect tab
    updateTabsetPanel(session, "plotsTab", selected = "Imputation Effect")
  })
  
  # Download handler for Imputation Effect Plot
  downloadPlot("downloadImputationEffectPlot", imputationEffectPlot, "imputation_effect_plot")
  
  ##Overview##
  
  
  
  ##dep object creation while selecting control##
  
  
  # Reactive values to hold the dep object and the number of proteins
  dep <- reactiveVal(NULL)
  numProteins <- reactiveVal(NULL)  # Holds the number of proteins detected
  
  # Function to perform differential enrichment analysis and store the dep object
  observeEvent(input$submitDifferentialAnalysis, {
    req(input$contrastType, imputed_se())  # Ensure contrast type and imputed data are available
    
    # Clear diff_results if switching from manual to another type
    if (input$contrastType != "manual") {
      diff_results(NULL)
    }
    
    # Perform the analysis based on the selected contrast type
    if (input$contrastType == "control") {
      req(input$controlSample)  # Ensure a control sample is selected
      
      # Perform differential analysis based on control vs treatment
      data_diff <- edited_test_diff(imputed_se(), type = "control", control = input$controlSample)
      
      # Create dep object based on user-defined cutoffs 
      dep_object <- add_rejections(data_diff, alpha = input$pValueCutoff, lfc = input$log2FCCutoff)
      
      # Store the dep object in reactive value (only for control vs treatment)
      dep(dep_object)
      
      # Store the number of proteins listed in the dep object
      protein_count <- nrow(assay(dep_object))  # Assuming the protein data is in the assay
      numProteins(protein_count)  # Store the count in the reactive value
      
      # Dynamically update the "Top N most variable proteins" input based on the number of proteins
      updateNumericInput(session, "n_pca", max = numProteins())
      
    } else if (input$contrastType == "all") {
      # Perform pairwise comparisons across all conditions
      data_diff <- edited_test_diff(imputed_se(), type = "all")
      
      # Apply rejections to mark significant proteins in all pairwise comparisons
      dep_object <- add_rejections(data_diff, alpha = input$pValueCutoff, lfc = input$log2FCCutoff)
      dep(dep_object)  # Store dep object in reactive value
      
    } else if (input$contrastType == "manual") {
      req(input$manualContrast)  # Ensure manual contrasts are provided
      contrasts <- unlist(strsplit(input$manualContrast, ",\\s*"))  # Split by commas and remove spaces
      data_diff <- edited_test_diff(imputed_se(), type = "manual", test = contrasts)
    }
    
    # Update the differential analysis results
    diff_results(data_diff)
    
    # Show confirmation modal that the analysis is complete
    showModal(modalDialog(
      title = "Differential Enrichment Analysis Complete",
      "Analysis complete based on the selected contrast type.",
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  
  # Generate the bar plot using plot_single() for selected proteins when control group is selected
  observeEvent(input$submitProteins, {
    # Ensure that the contrast type is "control" before generating the bar plot
    if (input$contrastType == "control") {
      req(dep())  # Ensure the dep object is available
      req(input$proteinNames)  # Ensure user has provided protein names
      
      # Split the input by commas for multiple protein names
      protein_list <- unlist(strsplit(input$proteinNames, ",\\s*")) 
      
      # Check if selected proteins are valid and exist in the dep object
      proteins_in_dep <- rownames(dep())
      
      validate(
        need(all(protein_list %in% proteins_in_dep), "One or more proteins are not available in the dataset.")
      )
      
      # Show spinner during plot rendering
      output$proteinBarplotUI <- renderUI({
        withSpinner(plotOutput("proteinBarplot"))
      })
      
      # Generate the bar plot using the plot_single function from DEP
      output$proteinBarplot <- renderPlot({
        plot_single(dep(), proteins = protein_list)
      })
      
    } else {
      # If not in control mode, inform the user that the bar plot can only be generated for control group
      showModal(modalDialog(
        title = "Bar Plot Not Available",
        "Bar plots can only be generated when the 'Control' contrast type is selected.",
        easyClose = TRUE,
        footer = NULL
      ))
      
      # Clear the previous plot if "All" or "Manual" is selected
      output$proteinBarplotUI <- renderUI({ NULL })
    }
  })
  
  
  #IMPORTANT
  # Function to get available contrasts from diff_results
  get_contrasts <- function(diff_results) {
    # Get all column names from rowData
    rd_names <- names(rowData(diff_results))
    
    # Extract contrasts by finding names ending with '_p.adj'
    contrasts <- sub("_p\\.adj$", "", rd_names[grep("_p\\.adj$", rd_names)])
    
    return(contrasts)
  }
  
  # Render the contrast selector UI element based on the contrast type
  output$contrastSelector <- renderUI({
    req(diff_results())
    
    # Retrieve available contrasts from diff_results
    contrasts <- get_contrasts(diff_results())
    
    # Display the select input based on contrast type
    if (input$contrastType == "manual" || input$contrastType == "control" || input$contrastType == "all") {
      if (length(contrasts) > 0) {
        selectInput("selectedContrast", "Select Contrast to Plot:",
                    choices = contrasts)
      } else {
        h5("No contrasts available to select.")
      }
    }
  })
  
  ### P - Value histogram ##
  
  # Function to get cleaned labels for p-value columns
  get_cleaned_pval_columns <- function(diff_results) {
    # Get all column names from rowData
    rd_names <- names(rowData(diff_results))
    
    # Extract p-value columns by finding names ending with '_p.val'
    pval_columns <- rd_names[grep("_p\\.val$", rd_names)]
    
    # Remove "_p.val" suffix for cleaner display in dropdown
    cleaned_labels <- sub("_p\\.val$", "", pval_columns)
    
    # Return as a named list: cleaned labels as display names, original columns as values
    names(pval_columns) <- cleaned_labels
    
    return(pval_columns)
  }
  
  # Render the dropdown menu for p-value columns
  output$pvalSelector <- renderUI({
    req(diff_results())  # Ensure differential results are available
    
    # Retrieve cleaned p-value column names
    cleaned_pval_columns <- get_cleaned_pval_columns(diff_results())
    
    # Display the select input for p-value columns with cleaned labels
    if (!is.null(cleaned_pval_columns) && length(cleaned_pval_columns) > 0) {
      selectInput("selectedPvalColumn", "Select P-Value Column for Histogram:",
                  choices = cleaned_pval_columns)  # Use cleaned labels
    } else {
      h5("No p-value columns available to select.")
    }
  })
  
  # Render the histogram plot for the selected p-value column
  output$pvalHistogram <- renderPlot({
    req(diff_results())  # Ensure differential results are available
    req(input$selectedPvalColumn)  # Ensure a p-value column is selected
    
    # Extract the selected p-value column
    pval_data <- rowData(diff_results())[[input$selectedPvalColumn]]
    
    # Generate the histogram using ggplot2
    ggplot(data.frame(pval_data = pval_data), aes(x = pval_data)) +
      geom_histogram(binwidth = 0.05, fill = "blue", color = "black", alpha = 0.7) +
      labs(title = paste("P-Value Distribution:", sub("_p\\.val$", "", input$selectedPvalColumn)),
           x = "P-Value", y = "Frequency") +
      theme_minimal() +
      geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +  # Optional significance threshold
      theme(
        plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)
      )
  })
  
  ### P - Value histogram ##
  
  # Render the custom ggplot2 volcano plot when a contrast is selected
  output$volcanoPlot <- renderPlot({
    req(diff_results())                  # Ensure diff_results is available before plotting
    req(input$selectedContrast)          # Ensure a contrast is selected
    
    # Apply the cutoffs and add rejections
    dep <- add_rejections(
      diff_results(), 
      alpha = input$pValueCutoff,  # Adjusted p-value cutoff
      lfc = input$log2FCCutoff     # Log2 fold change cutoff
    )
    
    # Access the contrast selected by the user
    selected_contrast <- input$selectedContrast
    
    # Generate the custom ggplot2 volcano plot
    ggplot(as.data.frame(rowData(dep)), 
           aes_string(x = paste0(selected_contrast, "_diff"), 
                      y = paste0("-log10(", selected_contrast, "_p.adj)"))) +
      geom_point(aes_string(color = paste0(selected_contrast, "_significant")), size = 2) +  # Points for proteins
      geom_vline(xintercept = c(-input$log2FCCutoff, input$log2FCCutoff), linetype = "dashed", color = "black") +  # Log2 fold change threshold
      geom_hline(yintercept = -log10(input$pValueCutoff), linetype = "dashed", color = "black") +  # p-value threshold
      geom_text_repel(aes_string(label = paste0("ifelse(", selected_contrast, "_p.adj < ", input$pValueCutoff, 
                                                " & abs(", selected_contrast, "_diff) > ", input$log2FCCutoff, 
                                                ", Gene.names, '')")), 
                      size = 3, max.overlaps = 100) +  # Increase max.overlaps to allow more labels
      scale_color_manual(name = "Significance",  # Updated legend title and labels
                         values = c("FALSE" = "grey", "TRUE" = "red"), 
                         labels = c("FALSE" = "Not Significant", "TRUE" = "Significant")) +  
      labs(title = paste("Volcano Plot:", selected_contrast),
           x = paste0("Log2 Fold Change (", selected_contrast, ")"), 
           y = "-log10 Adjusted P-value") +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 16, face = "bold"),  # Bold plot title
        legend.title = element_text(size = 12, face = "bold"),  # Bold legend title
        legend.text = element_text(size = 10),  # Adjust legend text size
        axis.title.x = element_text(size = 12, face = "bold"),  # Bold x-axis title
        axis.title.y = element_text(size = 12, face = "bold")   # Bold y-axis title
      )
  })
  
  # Download handler for volcano plot data
  output$downloadVolcanoData <- downloadHandler(
    filename = function() {
      paste("volcano_plot_data_", input$selectedContrast, Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(diff_results(), input$selectedContrast)  # Ensure data and contrast are available
      
      # Extract the data for the selected contrast
      volcano_data <- as.data.frame(rowData(diff_results()))  # Ensure this is scoped within `content`
      selected_contrast <- input$selectedContrast
      
      # Filter and structure the data for the selected contrast using "ID" as the main identifier
      volcano_data <- volcano_data %>%
        dplyr::select(
          ID,  # Use "ID" as the identifier column
          paste0(selected_contrast, "_diff"),
          paste0(selected_contrast, "_p.adj")
        ) %>%
        dplyr::rename(
          log2FoldChange = paste0(selected_contrast, "_diff"),
          adjPValue = paste0(selected_contrast, "_p.adj")
        )
      
      # Write to CSV
      write.csv(volcano_data, file, row.names = FALSE)
    }
  )
  
  data_results <- reactive({
    req(dep())  # Ensure the dep object is available
    get_results(dep())
  })
  
  # Render the table of significant proteins in the UI
  output$significantProteinsTable <- DT::renderDataTable({
    req(data_results())  # Ensure data_results is available
    data_results()
  })
  
  # Display the number of significant proteins
  output$numSignificantProteins <- renderText({
    req(data_results())  # Ensure data_results is available
    num_significant <- data_results() %>% filter(significant) %>% nrow()
    paste("Number of significant proteins:", num_significant)
  })
  
  # # Create the download handler for dep object
  # output$downloadDep <- downloadHandler(
  #   filename = function() {
  #     paste("dep_data_", Sys.Date(), ".xlsx", sep = "")
  #   },
  #   content = function(file) {
  #     req(dep())  # Ensure dep object is available
  #     
  #     tryCatch({
  #       # Extract the rowData from the dep object and convert it to a data frame
  #       dep_data <- as.data.frame(rowData(dep()))
  #       
  #       # Truncate long strings to a maximum of 32,767 characters
  #       dep_data <- dep_data %>%
  #         mutate(across(where(is.character), ~ substr(., 1, 32767)))
  #       
  #       # Write the dep data to an Excel file using writexl package
  #       writexl::write_xlsx(dep_data, path = file)
  #       
  #     }, error = function(e) {
  #       # Catch and display the error in a modal
  #       showModal(modalDialog(
  #         title = "Error in Download",
  #         paste("An error occurred while generating the download file:", conditionMessage(e)),
  #         easyClose = TRUE,
  #         footer = NULL
  #       ))
  #     })
  #   }
  # )
  
  
  # Use eventReactive to update the PCA plot only when the button is pressed
  pcaPlot <- eventReactive(input$submitPCA, {
    req(dep(), numProteins())  # Ensure both dep object and numProteins are available
    
    # Fetch inputs reactively when the button is pressed
    n <- input$n_pca
    x <- input$x_pca
    y <- input$y_pca
    
    # Validate the input: Check if the user selects more proteins than available
    validate(
      need(n <= numProteins(), 
           paste("Error: You can select up to", numProteins(), "proteins for the PCA plot. Please adjust the 'Top N most variable proteins' field."))
    )
    
    # Generate and return the PCA plot
    plot_pca(dep(), x = x, y = y, n = n, point_size = 4)
  })
  
  # Render the PCA plot in the UI
  output$pcaPlotUI <- renderUI({
    withSpinner(plotOutput("pcaPlot"))
  })
  
  output$pcaPlot <- renderPlot({
    req(pcaPlot())  # Only render if the eventReactive has generated the plot
    pcaPlot()  # Use the plot generated by eventReactive
  })
  
  # Display the number of proteins detected for the PCA plot
  output$numProteinsInfo <- renderText({
    req(numProteins())  # Ensure numProteins has a value before displaying it
    
    # Create a styled message to inform the user about the protein count
    HTML(paste(
      "<strong style='color: #0073e6;'>",  # Bold and blue color
      "Based on the analysis, you can select up to", numProteins(), 
      "proteins for the PCA plot.",
      "</strong>"
    ))
  })
  
  
  # Download handler for PCA plot
  downloadPlot("downloadPCAPlot", pcaPlot, "PCA_plot")
  
  
  #####Results#####
  observeEvent(input$submitPCA,{
    # Generate the results table
    data_results <- get_results(dep())
    # Print results table and number of significant proteins to console
    print("Results Table:")
    print(data_results)
    print("Number of significant proteins:")
    print(data_results %>% filter(significant) %>% nrow())
  })
  #######################################################################
  
  # Render the frequency plot of significant proteins by condition
  output$freqPlot <- renderPlot({
    req(dep())  # Ensure dep data is available
    
    # Silently skip rendering if there are two or fewer unique conditions
    if (length(unique(design()$condition)) <= 2) {
      return(NULL)  # Skip rendering without any message or error
    }
    
    plot_cond(dep())  # Call the plot_cond function to generate the plot
  })
  
  # Download handler for frequency plot
  downloadPlot("downloadFreqPlot", function() plot_cond(dep()), "Frequency_plot")
  
  ##GO##
  # Populate the conditionContrast dropdown based on available contrasts in rowData(dep)
  observeEvent(dep(), {
    # Detect contrast columns that end with "_significant"
    contrast_columns <- grep("_significant$", colnames(rowData(dep())), value = TRUE)
    contrast_choices <- gsub("_significant$", "", contrast_columns)  # Remove suffix for cleaner names
    
    # Update the dropdown with available contrasts
    updateSelectInput(session, "conditionContrast", choices = contrast_choices)
  })
  
  # Detect when the user switches to the GO Analysis tab and populate the organism dropdown
  observeEvent(input$navbarPage, {
    if (input$navbarPage == "Gene Ontology (GO) Analysis") {
      # Fetch organism data at app startup
      organism_data <- rba_panther_info(what = "organisms")
      
      # Convert to named list with long_name as label and taxon_id as value
      organism_choices <- setNames(as.list(organism_data$taxon_id), organism_data$long_name)
      
      # Update the selectInput for organisms
      updateSelectInput(session, "organism", choices = organism_choices, selected = "3702")
    }
  })
  
  # Reactive expression to get IDs based on selected condition vs control contrast
  significant_proteins_names <- reactive({
    req(input$conditionContrast)  # Ensure a contrast is selected
    
    # Construct the column name for the selected contrast's significance column
    significant_column <- paste0(input$conditionContrast, "_significant")
    
    # Extract IDs where the selected contrast column is TRUE
    rowData(dep())$name[rowData(dep())[[significant_column]] == TRUE]
  })
  
  # Reactive expression to get significant protein IDs only
  significant_proteins_IDs <- reactive({
    req(input$conditionContrast)  # Ensure a contrast is selected
    
    # Construct the column name for the selected contrast's significance column
    significant_column <- paste0(input$conditionContrast, "_significant")
    
    # Extract IDs where the selected contrast column is TRUE
    rowData(dep())$ID[rowData(dep())[[significant_column]] == TRUE]
  })
  
  # Run GO analysis when user clicks "Run GO Analysis"
  observeEvent(input$submitGO, {
    # Retrieve the selected IDs
    ids <- significant_proteins_names()
    
    # Ensure that all required inputs are selected before running the analysis
    req(ids, input$organism, input$ontologyType)
    
    # Retrieve the selected organism ID and ontology type (annotation dataset)
    selected_organism <- as.numeric(input$organism)
    selected_ontology <- input$ontologyType
    
    # Run the GO enrichment analysis with the selected parameters
    go <- rba_panther_enrich(ids, 
                             organism = selected_organism, 
                             annot_dataset = selected_ontology, 
                             cutoff = 0.05)
    
    # Rename columns for easier interpretation
    enriched_data <- go$result
    colnames(enriched_data) <- c("Proteins in List", "Fold Enrichment", "FDR", 
                                 "Expected Count", "Proteins in Reference", "p-Value", 
                                 "Plus/Minus", "GO Term ID", "GO Term Label")
    
    # Round Fold Enrichment and Expected Count to 3 decimal places
    enriched_data$`Fold Enrichment` <- round(enriched_data$`Fold Enrichment`, 3)
    enriched_data$`Expected Count` <- round(enriched_data$`Expected Count`, 3)
    
    # Render the GO results as a dynamic DataTable in the UI, excluding Plus/Minus column
    output$goResultsTable <- DT::renderDataTable({
      req(enriched_data)  # Ensure enriched_data has results
      DT::datatable(
        enriched_data[, c("GO Term Label", "Fold Enrichment", "FDR", "p-Value", 
                          "Proteins in List", "Expected Count", "Proteins in Reference")],
        options = list(pageLength = 10, autoWidth = TRUE)
      )
    })
    
    # Display Mapped and Unmapped IDs
    output$mappedIDs <- renderText({
      paste(go$input_list$mapped_id, collapse = ", ")
    })
    
    output$unmappedIDs <- renderText({
      paste(go$input_list$unmapped_id, collapse = ", ")
    })
    
    # Render Dot Plot for GO Term Enrichment
    output$dotPlot <- renderPlot({
      # Select top 10 terms based on FDR
      top_terms <- enriched_data[order(enriched_data$FDR), ][1:10, ]
      
      # Filter out rows with NA in the `GO Term Label`
      filtered_terms <- top_terms %>% 
        filter(!is.na(`GO Term Label`))
      
      # Create the lollipop plot
      ggplot(filtered_terms, aes(x = `Fold Enrichment`, y = reorder(`GO Term Label`, `Fold Enrichment`))) +
        geom_segment(aes(x = 0, xend = `Fold Enrichment`, yend = reorder(`GO Term Label`, `Fold Enrichment`)), 
                     color = "grey", size = 0.8) +  # Line segment (stick)
        geom_point(aes(size = as.integer(`Proteins in List`), color = -log10(`p-Value`))) +  # Point with integer size
        scale_color_gradient(low = "lightblue", high = "red") +
        scale_size_continuous(name = "Proteins in List", breaks = seq(min(filtered_terms$`Proteins in List`), 
                                                                      max(filtered_terms$`Proteins in List`), by = 4)) +
        labs(title = "GO Term Enrichment Lollipop Plot", x = "Fold Enrichment", y = "GO Term") +
        theme_minimal() +
        theme(
          panel.grid.major.y = element_blank(),  # Remove y-axis grid lines for cleaner look
          panel.grid.minor = element_blank()
        )
    })
    
  })
  
  #SKEGGs
  # Populate the conditionContrast dropdown for KEGG analysis, including "All Significant Proteins" as an option
  observeEvent(dep(), {
    # Detect contrast columns that end with "_significant"
    contrast_columns <- grep("_significant$", colnames(rowData(dep())), value = TRUE)
    contrast_choices <- c("All Significant Proteins", gsub("_significant$", "", contrast_columns))  # Add "All Significant Proteins"
    
    # Update the dropdown with available contrasts for KEGG, defaulting to "All Significant Proteins"
    updateSelectInput(session, "conditionContrastKEGG", choices = contrast_choices, selected = "All Significant Proteins")
  })
  
  #### REGULATION reactive ####
  significant_proteins_with_regulation <- eventReactive(input$runKEGGAnalysis, {
    req(data_results(), input$conditionContrastKEGG, input$keggOrganism)  # Ensure data_results, user selection, and organism input are available
    selected_contrast <- input$conditionContrastKEGG
    
    # Define a list to collect regulation data for each condition if "All Significant Proteins" is selected
    regulation_data_list <- list()
    
    # Choose all contrasts if "All Significant Proteins" is selected
    significant_columns <- if (selected_contrast == "All Significant Proteins") {
      grep("_significant$", colnames(data_results()), value = TRUE)
    } else {
      paste0(selected_contrast, "_significant")
    }
    
    for (significant_column in significant_columns) {
      # Extract the condition and its ratio column
      condition <- sub("_significant$", "", significant_column)
      ratio_column <- paste0(condition, "_ratio")
      
      # Filter for significant proteins and remove NAs in the ratio column
      significant_proteins <- data_results()[data_results()[[significant_column]] == TRUE & 
                                               !is.na(data_results()[[ratio_column]]), ]
      
      # Extract UniProt IDs and regulation (Up/Down) for each protein
      significant_ids <- significant_proteins$ID
      regulation <- ifelse(significant_proteins[[ratio_column]] > 0, "Up", "Down")
      
      # Convert UniProt IDs to KEGG IDs for the selected organism
      kegg_ids <- bitr_kegg(
        significant_ids,
        fromType = "uniprot",
        toType = "kegg",
        organism = input$keggOrganism
      )
      
      # Check if KEGG conversion returned results
      if (!is.null(kegg_ids) && nrow(kegg_ids) > 0) {
        # Match the converted KEGG IDs with the regulation values
        converted_data <- data.frame(
          Condition = condition,
          ID = kegg_ids$kegg,  # Use converted KEGG IDs
          Regulation = regulation[match(kegg_ids$uniprot, significant_ids)]  # Align regulation by UniProt ID
        )
        
        # Append the data frame for this condition to the list
        regulation_data_list[[condition]] <- converted_data
      }
    }
    
    # Combine all conditions into a single data frame
    regulation_data <- do.call(rbind, regulation_data_list)
    
    # Return the compiled data frame with KEGG IDs, regulation status, and conditions
    regulation_data
  })
  
  # Output the regulation information for the selected comparison
  output$RegulationIDs <- renderPrint({
    req(significant_proteins_with_regulation())  # Ensure the data is available
    
    # Fetch the updated regulation data with KEGG IDs
    regulation_data <- significant_proteins_with_regulation()
    
    # Print the data frame in a compact format
    if (nrow(regulation_data) == 0) {
      print("No significant proteins found for the selected contrast.")
    } else {
      print(regulation_data, row.names = FALSE)
    }
  })
  
  
  
  #############################
  
  # Reactive expression to convert UniProt IDs to KEGG IDs for the user-selected organism
  significant_proteins_IDs_for_kegg <- eventReactive(input$runKEGGAnalysis, {
    req(dep(), input$keggOrganism)  # Ensure dep object and organism input are available
    selected_contrast <- input$conditionContrastKEGG
    
    # Fetch significant protein IDs based on the selected option
    significant_ids <- if (selected_contrast == "All Significant Proteins") {
      rowData(dep())$ID[rowData(dep())$significant == TRUE]
    } else {
      significant_column <- paste0(selected_contrast, "_significant")
      rowData(dep())$ID[rowData(dep())[[significant_column]] == TRUE]
    }
    
    # Convert UniProt IDs to KEGG IDs for the selected organism
    kegg_ids <- bitr_kegg(
      significant_ids,
      fromType = "uniprot",
      toType = "kegg",
      organism = input$keggOrganism
    )
    
    # Return KEGG IDs if found, else show a warning
    if (nrow(kegg_ids) == 0) {
      showNotification("No KEGG IDs found for the selected contrast or significant proteins.", type = "warning")
      return(NULL)
    }
    
    kegg_ids$kegg
  })
  
  
  
  
  
  # Perform KEGG pathway enrichment analysis using the user-selected organism
  kegg_enrichment_results <- reactive({
    req(significant_proteins_IDs_for_kegg(), input$keggOrganism)  # Ensure KEGG IDs and organism input are available
    
    enrichKEGG(
      gene = significant_proteins_IDs_for_kegg(),
      organism = input$keggOrganism,             # User-provided KEGG organism code
      pvalueCutoff = 0.05                        # Set p-value cutoff for significance
    )
  })
  
  # Render the KEGG enrichment results as a table in the UI
  output$keggEnrichmentTable <- DT::renderDataTable({
    req(kegg_enrichment_results())  # Ensure results are available
    
    # Convert results to a data frame for display
    as.data.frame(kegg_enrichment_results())
  })
  
  # Trigger the KEGG analysis and render cnetplot only after button press
  observeEvent(input$runKEGGAnalysis, {
    output$keggCnetPlotUI <- renderUI({
      # Add the plot output with a spinner
      withSpinner(plotOutput("keggCnetPlot"))
    })
    
    # Generate the cnetplot after button click
    output$keggCnetPlot <- renderPlot({
      req(kegg_enrichment_results())  # Ensure KEGG enrichment results are available
      
      # Generate the cnetplot with the top 10 enriched categories
      cnetplot(kegg_enrichment_results(), showCategory = 10, circular = TRUE, color.params = list(edge = TRUE))
    })
  })
  
  
  #ANOVA#
  ##### Validate replicates for selected conditions #####
  observeEvent(input$run_anova, {
    selected_conditions <- input$anova_conditions
    
    # Check if there are at least two replicates per condition
    replicate_counts <- sapply(selected_conditions, function(condition) {
      sum(design()$condition == condition)
    })
    
    ##### Subset data for ANOVA #####
    anova_res <- reactive({
      req(imputed_se(), input$anova_conditions)  # Ensure data and user selection are available
      
      # Validate selection of at least 3 conditions
      validate(
        need(length(input$anova_conditions) >= 3,
             "You should choose at least 3 conditions to perform ANOVA analysis."))
      
      # Validate replicates and notify if requirements aren't met
      validate(
        need(all(replicate_counts >= 2), 
             "Each selected condition must have at least 2 replicates for ANOVA."))
      
      # Filter for selected conditions
      selected_conditions <- unique(input$anova_conditions)
      condition <- selected_conditions
      subset_indices <- colData(imputed_se())$condition %in% condition
      subset_se <- imputed_se()[, subset_indices]
      
      # Create desugn matrix & Fit the model
      design_formula <- formula(~ 0 + condition)
      col_data <- colData(subset_se)
      raw <- assay(subset_se)
      
      variables <- terms.formula(design_formula) %>%
        attr(., "variables") %>%
        as.character() %>%
        .[-1]
      
      # Obtain variable factors
      for(var in variables) {
        temp <- factor(col_data[[var]])
        assign(var, temp)
      }
      
      design <- model.matrix(design_formula, data = environment())
      colnames(design) <- gsub("condition", "", colnames(design))
      
      # Fit the model and extract significant proteins (threshold 0.01)
      fit <- lmFit(assay(subset_se), design = design)
      fit.smooth <- eBayes(fit)
      sig <- topTable(fit.smooth, number=Inf,
                      sort.by='none', coef=condition)
      sig$significant <- sig$adj.P.Val <= 0.01
      # se <- (sqrt(fit.smooth$s2.post) * fit.smooth$stdev.unscaled)[,condition,drop=FALSE]
      # colnames(se) <- paste("SE.", colnames(se), sep = "")
      # result <- cbind(sig, se)
      
      # Round numbers for being user-friendly
      rounded <- sig %>%
        mutate(across(where(is.numeric), round, 3))
      # Return results table
      return(sig)
    })
    
    # Render ANOVA table output
    output$anova_deps_table <- DT::renderDataTable(
      anova_res()
    )
    
    #download_handler
    output$download_anova_deps_table <- downloadHandler(
      filename = "DEPs_ANOVA.xlsx",
      content = function(file) {
        write_xlsx(anova_res(), file)
      }
    )
    
    # HEATMAP
    for_heatmap <- anova_res()[anova_res()$significant == TRUE,] %>%
      select(all_of(selected_conditions))
    
    output$anovaHeatmap <- renderPlot({
      pheatmap(for_heatmap, scale = "row")
    })
    
    output$download_anova_heatmap <- downloadHandler(
      filename = "heatmap_ANOVA.pdf",
      content = function(filename) {
        pdf(filename, height = 10)
        pheatmap(for_heatmap, scale = "row")
        dev.off()
      }
    )
  })
  
}
# Run the application
shinyApp(ui = ui, server = server)
