library(visNetwork)
library(shiny)
library(tidyverse)
library(stringdist)
library(igraph)


# Server---------------


server <- function(input, output) {
  # input the csv file with all the epitope data as in (Task_2_dataset.csv)
  inputdata <- reactive({
    req(input$file_input)
    
    
    df <-  read.csv(input$file_input$datapath)
    df
  })
  
  # input the query sequence as CSRWDEDGFYVFDYW
  input_seq_df <- eventReactive(input$calc_ld, {
    # convert the query sequence to empty data frame to help in joining to the input dataframe
    data.frame(
      id = "Query_seq"  ,
      cdr3_aa = input$input_seq ,
      clone_count = NA ,
      v_call = NA ,
      v_call.light = NA  ,
      KD.nM. = NA ,
      kon = NA  ,
      koff = NA ,
      size_kDa = NA ,
      Epitope = NA
    )
  })
  
  # join the input dataframe to input sequence to help build network
  new_seq_rep <-
    reactive({
      joined_df <- rbind(inputdata(), input_seq_df())
      joined_df$label <-
        paste0(joined_df$id, "_", joined_df$cdr3_aa) # join id and aa sequence columns to make unique 
      joined_df$group <-  # join columns vcall, vcall_light, epitope to help in visualization of network later
        paste(joined_df$v_call,
              joined_df$v_call.light,
              joined_df$Epitope, sep = ",")
      joined_df
    })
  
  
  
  # calculate the ld of the query seuqence
  ld_input <- reactive({
    ld_input <-   stringdistmatrix(a = new_seq_rep()$cdr3_aa,
                                   b = input$input_seq,
                                   method = "lv")
    rownames(ld_input) <- new_seq_rep()$label
    ld_input
  })
  
  
  
  # calculate  LV distance which will be converted to adjacency matrix for the network 
  adj_mat <- reactive({
    adj_mat <-  stringdistmatrix(a = new_seq_rep()$cdr3_aa,
                                 b = new_seq_rep()$cdr3_aa,
                                 method = "lv")
    colnames(adj_mat) <- new_seq_rep()$label
    rownames(adj_mat) <- new_seq_rep()$label
    adj_mat
  })
  
  
  
  # put various thresholds for parameters of the table, like LD score, Clone count, Kon, Koff, etc
  output$LD_score <- renderUI({
    sliderInput(
      "ld_score",
      "LD Score:",
      ceiling(mean(ld_input())),
      min = min(ld_input()[!ld_input() == 0]), #all ld scores except zero(self similarity) 
      max = min(c(12, max(ld_input(
      )))),
      round = TRUE,
      step = 1
    )
  })
  
  ld_slider_input <- eventReactive(input$update, {
    input$ld_score
  })
  
  
  output$Clone_Count <- renderUI({
    sliderInput(
      "clone_count",
      "Clone count",
      max(new_seq_rep()$clone_count, na.rm = TRUE),
      min = min(new_seq_rep()$clone_count, na.rm = TRUE),
      max = max(new_seq_rep()$clone_count, na.rm = TRUE),
      round = TRUE
    )
  })
  
  clone_count_slider_input <- eventReactive(input$update, {
    input$clone_count
  })
  
  
  output$Kon <- renderUI({
    sliderInput(
      "kon",
      "Kon",
      max(new_seq_rep()$kon, na.rm = TRUE),
      min = min(new_seq_rep()$kon, na.rm = TRUE),
      max = max(new_seq_rep()$kon, na.rm = TRUE),
      round = TRUE
    )
  })
  
  kon_slider_input <- eventReactive(input$update, {
    input$kon
  })
  
  
  
  output$Koff <- renderUI({
    sliderInput(
      "koff",
      "Koff",
      max(new_seq_rep()$koff, na.rm = TRUE),
      min = min(new_seq_rep()$koff, na.rm = TRUE),
      max = max(new_seq_rep()$koff, na.rm = TRUE),
      round = TRUE
    )
  })
  
  koff_slider_input <- eventReactive(input$update, {
    input$koff
  })
  
  
  output$Size_kDa <- renderUI({
    sliderInput(
      "size_kDa",
      "Size kDa",
      max(new_seq_rep()$size_kDa, na.rm = TRUE),
      min = min(new_seq_rep()$size_kDa, na.rm = TRUE),
      max = max(new_seq_rep()$size_kDa, na.rm = TRUE),
      round = TRUE
    )
  })
  
  size_kDa_slider_input <- eventReactive(input$update, {
    input$size_kDa
  })
  
  output$KDNM <- renderUI({
    sliderInput(
      "KDnM",
      "KD(nM)",
      max(new_seq_rep()$KD.nM., na.rm = TRUE),
      min = min(new_seq_rep()$KD.nM., na.rm = TRUE),
      max = max(new_seq_rep()$KD.nM., na.rm = TRUE),
      round = TRUE
    )
  })
  
  KDnM_slider_input <- eventReactive(input$update, {
    input$KDnM
  })
  
  
  
  adj_mat_threshold <- reactive({
    adj_mat_threshold <- adj_mat() # Use the LD matrix to convert to adjacency matrix 
    adj_mat_threshold[adj_mat_threshold <= ld_slider_input()] <- 1 # binarize the network based on the LD score because this is the main criteria for the whole analysis
    adj_mat_threshold[adj_mat_threshold > ld_slider_input()] <- 0 # If not binarized it creates many false and duplicate connections
    adj_mat_threshold
  })
  
  
  output$network_graph <- renderVisNetwork({ # make the network using the adjacency matrix calculated above
    network_graph <-
      graph_from_adjacency_matrix(adjmatrix = adj_mat_threshold(), mode = "undirected")
    network_graph <- simplify(network_graph)
    V(network_graph)$LD_input <- ld_input() # add ld score the vertices
    
    
    network_graph_df <-
      igraph::as_data_frame(network_graph, what = 'both') # convert igraph to dataframe to add attributes to the vertex
    
    network_graph_df$vertices <-
      merge(network_graph_df$vertices, # add all attributes related to each sequence to the vertices
            new_seq_rep(),
            by.x = "name",
            by.y = "label")
    
    network_graph_2 <- # convert dataframe to network again 
      graph_from_data_frame(network_graph_df$edges,
                            directed = F,
                            vertices = network_graph_df$vertices)
    
    
    theshold_network <- # subset the graph according to various thresholds of parameters
      induced_subgraph(network_graph_2, vids = c(names(V(
        network_graph_2
      ))[names(V(network_graph_2)) == paste0("Query_seq", "_", input$input_seq)], names(V(
        network_graph_2
      )[(V(network_graph_2)$LD_input <= ld_slider_input())
        &
          (V(network_graph_2)$name %in% V(network_graph_2)$name[V(network_graph_2)$clone_count <=
                                                                  clone_count_slider_input()])
        #
        &
          (V(network_graph_2)$name %in% V(network_graph_2)$name[V(network_graph_2)$kon <=
                                                                  kon_slider_input()])
        &
          (V(network_graph_2)$name %in% V(network_graph_2)$name[V(network_graph_2)$koff <=
                                                                  koff_slider_input()])
        &
          (V(network_graph_2)$name %in% V(network_graph_2)$name[V(network_graph_2)$size_kDa <=
                                                                  size_kDa_slider_input()])
        &
          (V(network_graph_2)$name %in% V(network_graph_2)$name[V(network_graph_2)$KD.nM. <=
                                                                  KDnM_slider_input()])])))
    
    
    
    
    
    igraph_to_visnetwork <- toVisNetworkData(theshold_network) #convert igraph to visnetwork
    
    visNetwork(nodes = igraph_to_visnetwork$nodes,
               edges = igraph_to_visnetwork$edges,
               width = "100%",) %>%
      visIgraphLayout(layout = "layout_with_fr", randomSeed = 1234) %>% visOptions( #visualize visnetwork
        highlightNearest = list( 
          enabled = TRUE,
          degree = 1,
          hover = T
        ),
        selectedBy = "group",
        nodesIdSelection = TRUE
      )
    
    
  })
  
  
  
  output$input_seq_df <- renderTable(input_seq_df())
  
  
}

# UI--------------
ui <-
  
  
  fluidPage(title = "Bioinformatic tool concept",
            wellPanel (
              h1("Bioinformatic tool concept"),
              p(
                "This is a concept tool to visualize diversity landscape of an antibody repertoire. useful to visually determine the correlation of sequence-similarity and other properties such as binding affinity,
heavy-light chain pairings etc."
              ),
            h2("Potential improvements"),
tags$ol(
  tags$li("Perform clustering on the network"),
  tags$li("Calculate network properties -- pageRank, centrality, betweenness, authority, cliques etc"),
  tags$li("Download associated DNA Sequence to help and facilitate wet lab scientist"),
  tags$li("For thresholding using thresolding on both sides instead of one sided thresholding"),
  tags$li("Network is being calculated twice, can be improved to help computation efficiency"),
  tags$li("Improve visalizations for various chain-- currently grouped visualization available"),
  tags$li("Add check boxes to change labels"),
  tags$li("May be descriptive graphs to help understand the basic composition of the sequences, like compositon of v chains, how many epitopes may the antibodies target etc"),
  
),          
  h2("Known bugs"),
tags$ol(
  tags$li("Generate graph button should come after ld score calculation"),
),
  h2("What if dealing with very big data?"),
  p("Given a dataset of 10 million sequences and 100 selected lead candidates sequences, what
would be a viable approach to identify additional 'interesting' sequences from the dataset
using similarity networks?"),
  p("Some of the ways can be:"),
tags$ul(
  tags$li("Collapse similar nodes together. similar LD scores can be clustered together? This might reduce the number of nodes"),
  tags$li("Subset the network based on the LD score -- perform th enetwork based analysis later"),
  tags$li("Personalized pagerank based on the 100 sequences and their neighbours will help utilize the full power of network but also reduce the load of using full network"),
),


            sidebarLayout(
              # sidebar panel UI-------------
              sidebarPanel(
                # input csv files
                
                fileInput(
                  inputId = "file_input",
                  "Upload the repoitoire sequences (Task_2_dataset.csv)",
                  multiple = FALSE,
                  accept = ".csv"
                ),
                
                # Input query sequence
                textInput(
                  inputId = "input_seq",
                  label = "Input the CD3 AA sequence to query",
                  value = "",
                  width = NULL,
                  placeholder = NULL
                ),
                
                # Calculate LD Score
                
                actionButton("calc_ld", "Calculate LD Scores"),
                
                br(),
                
                # sliders for thersholding of the parameters
                conditionalPanel(
                  "input.input_seq!=''",
                  uiOutput("LD_score"),
                  uiOutput("Clone_Count"),
                  uiOutput("Kon"),
                  uiOutput("Koff"),
                  uiOutput("Size_kDa"),
                  uiOutput("KDNM"),
                  actionButton("update", "Generate/Update graph")
                  
                ),
                
              ),
              # Main Panel ------------
              mainPanel(
                visNetworkOutput("network_graph",height = "800px"),)
            )))

shinyApp(ui = ui, server = server)
