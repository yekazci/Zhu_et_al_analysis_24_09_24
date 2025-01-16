
generate_gene_set_visualizations_w_extra_plot_w_cluster_option_network_focused <- function(pando_obj, all_gene_sets, 
                                                                           seurat_identity, 
                                                                           term_to_search, 
                                                                           output_dir, 
                                                                           num_of_gene_sets = NULL, 
                                                                           aesthetical = NULL,
                                                                           high_exp_plot = NULL,
                                                                           cell_type_of_interest = NULL,
                                                                           group_by = NULL,
                                                                           network_focused_on = NULL,
                                                                           plot_separate = NULL
                                                                                          ) {
    
  library(dplyr)
  library(stringr)
  library(here)
  library(ggplot2)
  library(patchwork)
  
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  sub_dir <- file.path(output_dir, paste0(term_to_search, '_GeneSets'))  
    
  dir.create(sub_dir)  
    
  # Filter gene sets based on the biological term
  filtered_gene_sets <- all_gene_sets %>%
    filter(str_detect(str_to_lower(gs_name), str_to_lower(term_to_search)))
  
  # Extract unique gene set names
  all_gs_names <- filtered_gene_sets %>% pull(gs_name) %>% unique()
  
  if (!is.null(num_of_gene_sets)) {
      
      all_gs_names <- all_gs_names[1:num_of_gene_sets]
      
      }
    
  # Load custom plotting scripts once
  source(here::here('r_scripts', 'custom_gene_plot.r'))
  source(here::here("r_scripts", "network_visualization_functions_w_degree_option.r"))
  
  # Initialize a list to store output paths
  output_files <- list()
  
  for (gene_set_name in all_gs_names) {  
    # Extract gene symbols for the specific gene set
    target_gene_set <- filtered_gene_sets %>%
      filter(gs_name == gene_set_name) %>%
      pull(gene_symbol)
    
    # Check if genes are present in the network modules
    network_genes <- pando_obj %>%
      NetworkModules %>%
      .@meta %>%
      select(tf, target) %>%
      unlist() %>%
      unique()
    
    selected_genes <- target_gene_set[target_gene_set %in% network_genes]
    
    if (!is.null(high_exp_plot)) {
        
        selected_genes <- identify_high_expression_genes_w_SCT_assay(seurat_object = pando_obj@data,  
                                                         selected_genes = selected_genes, 
                                                         cell_type_of_interest = cell_type_of_interest,
                                                         group_by = group_by) 
        
        }
      
    var_height_1 <- (length(selected_genes) / 100.1) %>% floor    
      

    if (length(selected_genes) == 0) {
      warning(paste("No matching genes found in the Pando network modules for gene set:", gene_set_name))
      next
    }
    
    # Set Seurat object identity
    seurat_object <- pando_obj@data
    Idents(seurat_object) <- seurat_identity
    
    # Generate and save expression plots for the genes overlapping with the network:
    expression_plot <- generate_chunked_plots_expression_for_SCT_assay(
      seurat_object = seurat_object,
      genes = selected_genes,
      aesthetical = aesthetical
    )
    
    expression_plot_file <- file.path(sub_dir, paste0(gene_set_name, "_overlap_expression_plot.png"))
    ggsave(expression_plot_file, plot = expression_plot, width = 30, height = 24 + 24*var_height_1, dpi = 300)
    
    # For the network_focused_argument:  
      
      if (!is.null(network_focused_on)) {
          
            if (network_focused_on == 'tf') {
                focus_args <- list(TF_names = selected_genes)
            } else if (network_focused_on == 'target') {
                focus_args <- list(target_genes = selected_genes)
            } else if (network_focused_on == 'both') {
                focus_args <- list(TF_names = selected_genes, target_genes = selected_genes)
            }
          
   }  else {
    
          focus_args <- list(node_names = selected_genes)
   }

    # Single call with dynamic arguments
    network_plot <- do.call(vis_subnetwork_w_GGNET2, 
                            c(list(meta_DF = pando_obj %>% NetworkModules %>% .@meta), focus_args, list(arrow = TRUE)) )

    
    network_plot_file <- file.path(sub_dir, paste0(gene_set_name, "_network_plot.png"))
    ggsave(network_plot_file, plot = network_plot, width = 30, height = 24, dpi = 300)
    
    # Generate and save the expression plots for the overlapping genes with all the network connections.
    # I will obtain all possible connections with one of nodes is from selected genes. Then,
    # I will filter the TFs and targets, separately, and plot thier expression values in respective plots:
    
      
if (!is.null(network_focused_on)) { 
    # Based on network_focused_on arguments, find the genes to be plotted
    if (network_focused_on == 'tf') {
        plotted_genes <- pando_obj %>% 
                         NetworkModules %>% 
                         .@meta %>% 
                         filter(tf %in% selected_genes) %>% 
                         select(tf, target) %>% 
                         unlist() %>% 
                         unique()
    } else if (network_focused_on == 'target') {
        plotted_genes <- pando_obj %>% 
                         NetworkModules %>% 
                         .@meta %>% 
                         filter(target %in% selected_genes) %>% 
                         select(tf, target) %>% 
                         unlist() %>% 
                         unique()
    } else if (network_focused_on == 'both') {
        plotted_genes <- pando_obj %>% 
                         NetworkModules %>% 
                         .@meta %>% 
                         filter(tf %in% selected_genes & target %in% selected_genes) %>% 
                         select(tf, target) %>% 
                         unlist() %>% 
                         unique()
    }
    
} else {
    # If network_focused_on is NULL, consider all genes connected to the selected_genes
    plotted_genes <- pando_obj %>% 
                     NetworkModules %>% 
                     .@meta %>% 
                     filter(tf %in% selected_genes | target %in% selected_genes) %>% 
                     select(tf, target) %>% 
                     unlist() %>% 
                     unique()
}
 
      
    print(plotted_genes)  
      
    if (!is.null(plot_separate)) {  
            
     all_TFs_in_the_network <- pando_obj %>% 
                          NetworkModules %>% 
                                 .@meta  %>%
                                pull(tf) %>% 
                                unique()    
        
     all_Targets_in_the_network <- pando_obj %>% 
                         NetworkModules %>% 
                         .@meta  %>%
                         pull(target) %>% 
                         unique()    
                            
        # use plotted_genes to create separate plots.
        
     var_height_2 <- (length(intersect(plotted_genes, all_TFs_in_the_network)) / 100.1) %>% floor  
      
     expression_plot_all_TF <- generate_chunked_plots_expression_for_SCT_assay(
      seurat_object = seurat_object,
      genes = intersect(plotted_genes, all_TFs_in_the_network),
      aesthetical = aesthetical )  
        
    expression_plot_all_TF_file <- file.path(sub_dir, paste0(gene_set_name, "_all_connections_TF_expression_plot.png"))
    ggsave(expression_plot_all_TF_file, plot = expression_plot_all_TF, width = 30, height = 24 + 24*var_height_2, dpi = 300)   
     
        
     var_height_3 <- (length(intersect(plotted_genes, all_Targets_in_the_network)) / 100.1) %>% floor  
   
        
     expression_plot_all_Target <- generate_chunked_plots_expression_for_SCT_assay(
      seurat_object = seurat_object,
      genes = intersect(plotted_genes, all_Targets_in_the_network),
      aesthetical = aesthetical   )
         
         
    expression_plot_all_Target_file <- file.path(sub_dir, paste0(gene_set_name, "_all_connections_Target_expression_plot.png"))
    ggsave(expression_plot_all_Target_file, plot = expression_plot_all_Target, width = 30, height = 24 + 24*var_height_3, dpi = 300)        
       
        
       # Save file paths to output list
    output_files[[gene_set_name]] <- list(
      expression_plot = expression_plot_file,
      network_plot = network_plot_file,
      expression_plot_all_Target = expression_plot_all_Target_file,
      expression_plot_all_TF = expression_plot_all_TF_file
    )    
        
    } else {    
        
        
    # use used genes to create a single plot.   
    
        
     var_height_4 <- (length(plotted_genes) / 100.1) %>% floor  
      
     expression_plot_all <- generate_chunked_plots_expression_for_SCT_assay(
      seurat_object = seurat_object,
      genes = plotted_genes,
      aesthetical = aesthetical
    )
    
    expression_plot_all_file <- file.path(sub_dir, paste0(gene_set_name, "_all_connections_expression_plot.png"))
    ggsave(expression_plot_all_file, plot = expression_plot_all, width = 30, height = 24 + 24*var_height_4, dpi = 300)    
      
        
       # Save file paths to output list
    output_files[[gene_set_name]] <- list(
      expression_plot = expression_plot_file,
      network_plot = network_plot_file,
      expression_plot_all = expression_plot_all_file
    )    
        
        
    }
    
      
      
 
  }
  
  return(output_files)
}


######################################################

generate_gene_set_visualizations <- function(pando_obj, all_gene_sets, seurat_identity, term_to_search, output_dir) {
  library(dplyr)
  library(stringr)
  library(here)
  library(ggplot2)
  library(patchwork)
  
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  sub_dir <- file.path(output_dir, paste0(term_to_search, '_GeneSets'))  
    
  dir.create(sub_dir)  
    
  # Filter gene sets based on the biological term
  filtered_gene_sets <- all_gene_sets %>%
    filter(str_detect(str_to_lower(gs_name), str_to_lower(term_to_search)))
  
  # Extract unique gene set names
  all_gs_names <- filtered_gene_sets %>% pull(gs_name) %>% unique()
  
  # Load custom plotting scripts once
  source(here::here('r_scripts', 'custom_gene_plot.r'))
  source(here::here("r_scripts", "network_visualization_functions_w_degree_option.r"))
  
  # Initialize a list to store output paths
  output_files <- list()
  
  for (gene_set_name in all_gs_names) {  
    # Extract gene symbols for the specific gene set
    target_gene_set <- filtered_gene_sets %>%
      filter(gs_name == gene_set_name) %>%
      pull(gene_symbol)
    
    # Check if genes are present in the network modules
    network_genes <- pando_obj %>%
      NetworkModules %>%
      .@meta %>%
      select(tf, target) %>%
      unlist() %>%
      unique()
    
    selected_genes <- target_gene_set[target_gene_set %in% network_genes]
    
    if (length(selected_genes) == 0) {
      warning(paste("No matching genes found in the Pando network modules for gene set:", gene_set_name))
      next
    }
    
    # Set Seurat object identity
    seurat_object <- pando_obj@data
    Idents(seurat_object) <- seurat_identity
    
    # Generate and save expression plots
    expression_plot <- generate_chunked_plots_expression_for_SCT_assay(
      seurat_object = seurat_object,
      genes = selected_genes
    )
    
    expression_plot_file <- file.path(sub_dir, paste0(gene_set_name, "_expression_plot.png"))
    ggsave(expression_plot_file, plot = expression_plot, width = 30, height = 24, dpi = 300)
    
    # Generate and save network visualization
    network_plot <- vis_subnetwork_w_GGNET2(
      meta_DF = pando_obj %>% NetworkModules %>% .@meta,
      node_names = selected_genes,
      arrow = TRUE
    )
    
    network_plot_file <- file.path(sub_dir, paste0(gene_set_name, "_network_plot.png"))
    ggsave(network_plot_file, plot = network_plot, width = 30, height = 24, dpi = 300)
    
    # Save file paths to output list
    output_files[[gene_set_name]] <- list(
      expression_plot = expression_plot_file,
      network_plot = network_plot_file
    )
  }
  
  return(output_files)
}

#######################

generate_gene_set_visualizations_w_extra_plot <- function(pando_obj, all_gene_sets, seurat_identity, term_to_search, output_dir, num_of_gene_sets = NULL, aesthetical = NULL) {
  library(dplyr)
  library(stringr)
  library(here)
  library(ggplot2)
  library(patchwork)
  
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  sub_dir <- file.path(output_dir, paste0(term_to_search, '_GeneSets'))  
    
  dir.create(sub_dir)  
    
  # Filter gene sets based on the biological term
  filtered_gene_sets <- all_gene_sets %>%
    filter(str_detect(str_to_lower(gs_name), str_to_lower(term_to_search)))
  
  # Extract unique gene set names
  all_gs_names <- filtered_gene_sets %>% pull(gs_name) %>% unique()
  
  if (!is.null(num_of_gene_sets)) {
      
      all_gs_names <- all_gs_names[1:num_of_gene_sets]
      
      }
    
  # Load custom plotting scripts once
  source(here::here('r_scripts', 'custom_gene_plot.r'))
  source(here::here("r_scripts", "network_visualization_functions_w_degree_option.r"))
  
  # Initialize a list to store output paths
  output_files <- list()
  
  for (gene_set_name in all_gs_names) {  
    # Extract gene symbols for the specific gene set
    target_gene_set <- filtered_gene_sets %>%
      filter(gs_name == gene_set_name) %>%
      pull(gene_symbol)
    
    # Check if genes are present in the network modules
    network_genes <- pando_obj %>%
      NetworkModules %>%
      .@meta %>%
      select(tf, target) %>%
      unlist() %>%
      unique()
    
    selected_genes <- target_gene_set[target_gene_set %in% network_genes]
    
    var_height_1 <- (length(selected_genes) / 100.1) %>% floor    
      

    if (length(selected_genes) == 0) {
      warning(paste("No matching genes found in the Pando network modules for gene set:", gene_set_name))
      next
    }
    
    # Set Seurat object identity
    seurat_object <- pando_obj@data
    Idents(seurat_object) <- seurat_identity
    
    # Generate and save expression plots for the genes overlapping with the network:
    expression_plot <- generate_chunked_plots_expression_for_SCT_assay(
      seurat_object = seurat_object,
      genes = selected_genes,
      aesthetical = aesthetical
    )
    
    expression_plot_file <- file.path(sub_dir, paste0(gene_set_name, "_overlap_expression_plot.png"))
    ggsave(expression_plot_file, plot = expression_plot, width = 30, height = 24 + 24*var_height_1, dpi = 300)
    
    # Generate and save network visualization
    network_plot <- vis_subnetwork_w_GGNET2(
      meta_DF = pando_obj %>% NetworkModules %>% .@meta,
      node_names = selected_genes,
      arrow = TRUE
    )
    
    network_plot_file <- file.path(sub_dir, paste0(gene_set_name, "_network_plot.png"))
    ggsave(network_plot_file, plot = network_plot, width = 30, height = 24, dpi = 300)
    
    # Generate and save the expression plots for the overlapping genes with all the network connections:
    
    all_genes_in_the_network <- pando_obj %>% 
                         NetworkModules %>% 
                         .@meta  %>% 
                         filter(tf %in% selected_genes | target %in% selected_genes) %>% 
                         select(tf, target) %>% 
                         unlist() %>% 
                         unique()
      
    var_height_2 <- (length(all_genes_in_the_network) / 100.1) %>% floor  
      
     expression_plot_all <- generate_chunked_plots_expression_for_SCT_assay(
      seurat_object = seurat_object,
      genes = all_genes_in_the_network,
      aesthetical = aesthetical
    )
    
    expression_plot_all_file <- file.path(sub_dir, paste0(gene_set_name, "_all_connections_expression_plot.png"))
    ggsave(expression_plot_all_file, plot = expression_plot_all, width = 30, height = 24 + 24*var_height_2, dpi = 300)
      
      
    # Save file paths to output list
    output_files[[gene_set_name]] <- list(
      expression_plot = expression_plot_file,
      network_plot = network_plot_file,
      expression_plot_all = expression_plot_all_file
    )
  }
  
  return(output_files)
}

#####################################################

# Write a function that computes from the input genes the genes that are of average expression in the desired identity:

identify_high_expression_genes_w_SCT_assay <- function(seurat_object, selected_genes, cell_type_of_interest, group_by) {
  
  # Check if the group_by argument exists in Idents(seurat_object)
  if (!group_by %in% colnames( slot(object = seurat_object, name = 'meta.data') ) ) {
    stop(glue::glue("The specified group_by '{group_by}' is not found in the meta.data of the Seurat object."))
  }
    
    
    # Check if the cell_type_of_interest exists in the group_by column of the meta.data: 
  if (!cell_type_of_interest %in% dplyr::pull(var = group_by, .data = seurat_object %>% slot(name = 'meta.data')) ) {
    stop(glue::glue("The specified cell_type_of_interest '{cell_type_of_interest}' is not found in the '{group_by}' column of meta.data of the Seurat object."))
  }  
  
  # Ensure selected_genes are available in the Seurat object
  available_genes <- intersect(selected_genes, rownames(seurat_object@assays$SCT$data))
  missing_genes <- setdiff(selected_genes, available_genes)
  if (length(missing_genes) > 0) {
    warning(glue::glue("The following genes are not found in the Seurat object: {paste(missing_genes, collapse = ', ')}"))
  }
  
  if (length(available_genes) == 0) {
    stop("No valid genes found in the Seurat object.")
  }
  
  # Calculate average expression for all cell types
  avg_expression <- Seurat::AverageExpression(
    object = seurat_object, 
    assays = "SCT", 
    layer = "data",
    group.by = group_by,
    features = available_genes,
    return.seurat = FALSE  # Return a data frame for easier handling
  )$SCT  # Extract the RNA assay average
  
  # Identify the genes with higher average expression in the cell_type_of_interest
  high_expression_genes <- rownames(avg_expression)[
    avg_expression[, cell_type_of_interest] == apply(X = avg_expression, MARGIN = 1, FUN = max)
  ]
  
  # alternative to above command:
    
  #   high_expression_genes <- rownames(avg_expression)[
  # apply(avg_expression, 1, which.max) == which(colnames(avg_expression) == cell_type_of_interest)
  #   ]

  
  return(sort(high_expression_genes))
}

##############################################################################

generate_gene_set_visualizations_w_extra_plot_w_cluster_option <- function(pando_obj, all_gene_sets, 
                                                                           seurat_identity, 
                                                                           term_to_search, 
                                                                           output_dir, 
                                                                           num_of_gene_sets = NULL, 
                                                                           aesthetical = NULL,
                                                                           high_exp_plot = NULL,
                                                                           cell_type_of_interest = NULL,
                                                                           group_by = NULL) {
  library(dplyr)
  library(stringr)
  library(here)
  library(ggplot2)
  library(patchwork)
  
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  sub_dir <- file.path(output_dir, paste0(term_to_search, '_GeneSets'))  
    
  dir.create(sub_dir)  
    
  # Filter gene sets based on the biological term
  filtered_gene_sets <- all_gene_sets %>%
    filter(str_detect(str_to_lower(gs_name), str_to_lower(term_to_search)))
  
  # Extract unique gene set names
  all_gs_names <- filtered_gene_sets %>% pull(gs_name) %>% unique()
  
  if (!is.null(num_of_gene_sets)) {
      
      all_gs_names <- all_gs_names[1:num_of_gene_sets]
      
      }
    
  # Load custom plotting scripts once
  source(here::here('r_scripts', 'custom_gene_plot.r'))
  source(here::here("r_scripts", "network_visualization_functions_w_degree_option.r"))
  
  # Initialize a list to store output paths
  output_files <- list()
  
  for (gene_set_name in all_gs_names) {  
    # Extract gene symbols for the specific gene set
    target_gene_set <- filtered_gene_sets %>%
      filter(gs_name == gene_set_name) %>%
      pull(gene_symbol)
    
    # Check if genes are present in the network modules
    network_genes <- pando_obj %>%
      NetworkModules %>%
      .@meta %>%
      select(tf, target) %>%
      unlist() %>%
      unique()
    
    selected_genes <- target_gene_set[target_gene_set %in% network_genes]
    
    if (!is.null(high_exp_plot)) {
        
        selected_genes <- identify_high_expression_genes_w_SCT_assay(seurat_object = pando_obj@data,  
                                                         selected_genes = selected_genes, 
                                                         cell_type_of_interest = cell_type_of_interest,
                                                         group_by = group_by) 
        
        }
      
    var_height_1 <- (length(selected_genes) / 100.1) %>% floor    
      

    if (length(selected_genes) == 0) {
      warning(paste("No matching genes found in the Pando network modules for gene set:", gene_set_name))
      next
    }
    
    # Set Seurat object identity
    seurat_object <- pando_obj@data
    Idents(seurat_object) <- seurat_identity
    
    # Generate and save expression plots for the genes overlapping with the network:
    expression_plot <- generate_chunked_plots_expression_for_SCT_assay(
      seurat_object = seurat_object,
      genes = selected_genes,
      aesthetical = aesthetical
    )
    
    expression_plot_file <- file.path(sub_dir, paste0(gene_set_name, "_overlap_expression_plot.png"))
    ggsave(expression_plot_file, plot = expression_plot, width = 30, height = 24 + 24*var_height_1, dpi = 300)
    
    # Generate and save network visualization
    network_plot <- vis_subnetwork_w_GGNET2(
      meta_DF = pando_obj %>% NetworkModules %>% .@meta,
      node_names = selected_genes,
      arrow = TRUE
    )
    
    network_plot_file <- file.path(sub_dir, paste0(gene_set_name, "_network_plot.png"))
    ggsave(network_plot_file, plot = network_plot, width = 30, height = 24, dpi = 300)
    
    # Generate and save the expression plots for the overlapping genes with all the network connections:
    
    all_genes_in_the_network <- pando_obj %>% 
                         NetworkModules %>% 
                         .@meta  %>% 
                         filter(tf %in% selected_genes | target %in% selected_genes) %>% 
                         select(tf, target) %>% 
                         unlist() %>% 
                         unique()
      
    var_height_2 <- (length(all_genes_in_the_network) / 100.1) %>% floor  
      
     expression_plot_all <- generate_chunked_plots_expression_for_SCT_assay(
      seurat_object = seurat_object,
      genes = all_genes_in_the_network,
      aesthetical = aesthetical
    )
    
    expression_plot_all_file <- file.path(sub_dir, paste0(gene_set_name, "_all_connections_expression_plot.png"))
    ggsave(expression_plot_all_file, plot = expression_plot_all, width = 30, height = 24 + 24*var_height_2, dpi = 300)
      
      
    # Save file paths to output list
    output_files[[gene_set_name]] <- list(
      expression_plot = expression_plot_file,
      network_plot = network_plot_file,
      expression_plot_all = expression_plot_all_file
    )
  }
  
  return(output_files)
}

#############################

