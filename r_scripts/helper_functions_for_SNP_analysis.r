
# Write a function to create cell type specific eGRNs using differentially accessible
# regions and global eGRN object as inputs, and find TF modules:

make_cell_type_eGRN <- function(global_eGRN_object, DA_regions) {
    
    celltype_eGRN <- global_eGRN_object
    
    celltype_eGRN@grn@networks$glm_network@coefs <- 
                celltype_eGRN@grn@networks$glm_network@coefs %>% 
                            filter(region %in% rownames(DA_regions))
    
    # Find TF modules:

    celltype_eGRN <- find_modules(
                                    celltype_eGRN, 
                                    p_thresh = 0.1,
                                    nvar_thresh = 2, 
                                    min_genes_per_module = 1, 
                                    rsq_thresh = 0.05 )
    
    return(celltype_eGRN)
}

# Another function to generate Granges objects from the cell type specific eGRN objects
# to be readily used in SNP enrichment analysis:

make_GRanges_for_SNP_enrich <- function(Pando_eGRN_obj) {


eGRN_regions <- 
    Pando_eGRN_obj %>% 
    NetworkModules() %>% 
    .@meta %>% 
    dplyr::select(regions)


eGRN_regions_GRanges <- 
    eGRN_regions %>% 
    pull(regions) %>% 
    str_split(pattern = ';') %>% 
    unlist %>% 
    unique %>% 
    Signac::StringToGRanges(sep = c('-', '-'))

 return(eGRN_regions_GRanges)   
    
}

#####################################################

make_cell_type_eGRN_using_min_cells <- function (global_eGRN_object, cell_type, main_identity, min_cells) 
{
    celltype_eGRN <- global_eGRN_object
	
    Idents(celltype_eGRN@data) <- main_identity       


    selected_peaks <-  
        Signac::AccessiblePeaks(object = celltype_eGRN@data, 
                                assay = 'peaks', idents = cell_type, 
                                min.cells = min_cells)
    
    celltype_eGRN@grn@networks$glm_network@coefs <- celltype_eGRN@grn@networks$glm_network@coefs %>% 
        filter(region %in% selected_peaks)
    celltype_eGRN <- find_modules(celltype_eGRN, p_thresh = 0.1, 
        nvar_thresh = 2, min_genes_per_module = 1, rsq_thresh = 0.05)
    return(celltype_eGRN)
}

######################################################

# Use the combined function to create cell type analysis with varying numbers for min_cells parameter:

prepare_cell_type_specific_eGRNs_and_do_SNP_analysis_w_min_cells <- function(global_eGRN_object, min_cells, main_identity, cell_type_names, taSNPLD) {
  

  # Initialize lists to store intermediate and final results
  cell_type_eGRNs <- list()
  cell_type_GRanges <- list()
  res_list_SNPs <- list()
      
    
  # Directory setup
    
  main_output_dir <- here::here('R_Outputs', paste0('SNPs_w_min_', min_cells, '_cells'))
  output_dir <- here::here(main_output_dir, paste0('cell_type_eGRNs_using_min', min_cells, '_cells'))
  snp_dir <- here::here(main_output_dir, paste0('cell_type_eGRNs_min', min_cells, '_cells_SNPs'))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(snp_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Step 1: Prepare cell type eGRNs
  for (cell_type in cell_type_names) {
    clean_name <- janitor::make_clean_names(string = cell_type)
    
    # Generate cell type-specific eGRN
    cell_type_eGRNs[[clean_name]] <- make_cell_type_eGRN_using_min_cells(
      global_eGRN_object = global_eGRN_object,
      min_cells = min_cells,
      main_identity = main_identity,
      cell_type = cell_type
    )
    
    # Generate and save network plot
    network_plot <- vis_network_w_GGNET2(meta_DF = cell_type_eGRNs[[clean_name]] %>% NetworkModules %>% .@meta)
    network_plot_file <- file.path(output_dir, paste0(clean_name, "_min", min_cells, "_cells_network_plot.png"))
    ggsave(network_plot_file, plot = network_plot, width = 30, height = 24, dpi = 300)
  }
  
  # Step 2: Create GRanges for SNP enrichment
  for (clean_name in names(cell_type_eGRNs)) {
    cell_type_GRanges[[clean_name]] <- make_GRanges_for_SNP_enrich(
      Pando_eGRN_obj = cell_type_eGRNs[[clean_name]]
    )
  }
  
  # Step 3: Perform SNP analysis
  for (clean_name in names(cell_type_GRanges)) {
      
    snp_sub_dir  <- here::here(snp_dir, clean_name)
      
    dir.create(path = snp_sub_dir, showWarnings = FALSE, recursive = TRUE)
      
    for (expansion_bases in seq(from = 0, to = 500, by = 50)) {
      result <- try({
        enrich_SNP_generic(
          granges_obj = cell_type_GRanges[[clean_name]], 
          expand_by = expansion_bases,
          SNPs = taSNPLD,
          suppress_messages = TRUE
        )
      }, silent = TRUE)
      
      # Skip if there's an error
      if (inherits(result, "try-error")) {
        cat(sprintf("No SNP overlap found for %d bp expansion. Skipping...\n", expansion_bases))
        next
      }
      
      # Store successful results
      res_list_SNPs[[clean_name]][[paste0('expanded_by_', expansion_bases, 'bp')]] <- result
      
      # Generate and save SNP plot
      snp_plot <- plot_for_diff_expansions_for_saving(
        enriched_SNP_obj = result,
        expansion_value_string = paste0(expansion_bases, 'bp')
      )
      plot_snp_file <- file.path(snp_sub_dir, paste0(clean_name, "_expanded_by_", expansion_bases, "bp_min", min_cells, "_cells_SNP_plot.png"))
      ggsave(plot_snp_file, plot = snp_plot, width = 30, height = 24, dpi = 300)
    }
  }
  
  # Return final results
  return(list(
    eGRNs = cell_type_eGRNs,
    GRanges = cell_type_GRanges,
    SNP_results = res_list_SNPs
  ))
}



##############################################################

# Use the combined function to create cell type analysis with varying numbers for min_cells parameter:

prepare_cell_type_specific_eGRNs_and_do_SNP_analysis_w_DA_regions <- function(global_eGRN_object, DA_regions, cell_type_names, taSNPLD) {
  

  # Initialize lists to store intermediate and final results
  cell_type_eGRNs <- list()
  cell_type_GRanges <- list()
  res_list_SNPs <- list()
      
    
  # Directory setup
    
  main_output_dir <- here::here('R_Outputs', paste0('SNPs_w_da_regions'))
  output_dir <- here::here(main_output_dir, paste0('cell_type_eGRNs_using_da_regions'))
  snp_dir <- here::here(main_output_dir, paste0('cell_type_eGRNs_using_da_regions_SNPs'))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(snp_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Step 1: Prepare cell type eGRNs
  for (cell_type in cell_type_names) {
    clean_name <- janitor::make_clean_names(string = cell_type)
    
    # Generate cell type-specific eGRN
    cell_type_eGRNs[[clean_name]] <- make_cell_type_eGRN(
      global_eGRN_object = global_eGRN_object,
      DA_regions = DA_regions  %>%  filter(cluster %in% cell_type)
    )
    
    # Generate and save network plot
    network_plot <- vis_network_w_GGNET2(meta_DF = cell_type_eGRNs[[clean_name]] %>% NetworkModules %>% .@meta)
    network_plot_file <- file.path(output_dir, paste0(clean_name, "_eGRN_network_plot.png"))
    ggsave(network_plot_file, plot = network_plot, width = 30, height = 24, dpi = 300)
  }
  
  # Step 2: Create GRanges for SNP enrichment
  for (clean_name in names(cell_type_eGRNs)) {
    cell_type_GRanges[[clean_name]] <- make_GRanges_for_SNP_enrich(
      Pando_eGRN_obj = cell_type_eGRNs[[clean_name]]
    )
  }
  
  # Step 3: Perform SNP analysis
  for (clean_name in names(cell_type_GRanges)) {
      
    snp_sub_dir  <- here::here(snp_dir, clean_name)
      
    dir.create(path = snp_sub_dir, showWarnings = FALSE, recursive = TRUE)
      
    for (expansion_bases in seq(from = 0, to = 500, by = 50)) {
      result <- try({
        enrich_SNP_generic(
          granges_obj = cell_type_GRanges[[clean_name]], 
          expand_by = expansion_bases,
          SNPs = taSNPLD,
          suppress_messages = TRUE
        )
      }, silent = TRUE)
      
      # Skip if there's an error
      if (inherits(result, "try-error")) {
        cat(sprintf("No SNP overlap found for %d bp expansion. Skipping...\n", expansion_bases))
        next
      }
      
      # Store successful results
      res_list_SNPs[[clean_name]][[paste0('expanded_by_', expansion_bases, 'bp')]] <- result
      
      # Generate and save SNP plot
      snp_plot <- plot_for_diff_expansions_for_saving(
        enriched_SNP_obj = result,
        expansion_value_string = paste0(expansion_bases, 'bp')
      )
      plot_snp_file <- file.path(snp_sub_dir, paste0(clean_name, "_expanded_by_", expansion_bases, "bp_SNP_plot.png"))
      ggsave(plot_snp_file, plot = snp_plot, width = 30, height = 24, dpi = 300)
    }
  }
  
  # Return final results
  return(list(
    eGRNs = cell_type_eGRNs,
    GRanges = cell_type_GRanges,
    SNP_results = res_list_SNPs
  ))
}

