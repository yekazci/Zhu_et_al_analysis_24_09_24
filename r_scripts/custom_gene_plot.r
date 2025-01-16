generate_chunked_plots_expression_for_SCT_assay <- function(seurat_object, genes, chunk_size = 20, aesthetical = NULL) {
  
  # Check if the genes are present in the Seurat object
  available_genes <- intersect(genes, rownames(seurat_object@assays$SCT$data))  # Note that the assay objects is v5.
  missing_genes <- setdiff(genes, available_genes)
  
  if (length(missing_genes) > 0) {
    warning(glue::glue("The following genes are not found in the Seurat object: {paste(missing_genes, collapse = ', ')}"))
  }
  
  if (length(available_genes) == 0) {
    stop("No genes found in the Seurat object.")
  }
  
  # Sort the available genes alphabetically
  available_genes <- sort(available_genes)
  
    
  ########################## ADDING the padding genes to the last panel or not ################
    
  if (!is.null(aesthetical)) {
    
    # Check if the total number of genes is a multiple of chunk_size
  remainder <- length(available_genes) %% chunk_size
  if (remainder > 0) {
    # Pad with the last gene repeated to make it a multiple of chunk_size
    padding_genes <- rep(tail(available_genes, 1), chunk_size - remainder)
    available_genes <- c(available_genes, padding_genes)
  }
  
  # Calculate the number of chunks needed (now a perfect multiple of chunk_size)
  num_chunks <- length(available_genes) / chunk_size 
    
   } else {
      
            # Calculate the number of chunks needed
  num_chunks <- ceiling(length(available_genes) / chunk_size)
      
  }
      
   ######################################################
 
  
  # Initialize an empty list to store plots for each chunk
  plot_list <- list()
  
  # Generate violin plots for each chunk of genes
  for (i in 1:num_chunks) {
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, length(available_genes))
    
    # Extract genes for the current chunk
    genes_chunk <- available_genes[start_idx:end_idx]
    
    # Generate stacked violin plot for the current chunk
    plot_list[[i]] <- scCustomize::Stacked_VlnPlot(
      seurat_object = seurat_object,
      features = genes_chunk,
      x_lab_rotate = TRUE,
      plot_spacing = 0.1,
      plot_legend = i == num_chunks  # Show legend only for the last chunk
    )
  }
    

  
  # Arrange the plots into rows using patchwork
  num_rows <- ceiling(num_chunks / 5)  # Arrange up to 5 plots per row
  row_plots <- vector("list", length = num_rows)
  
  for (i in 1:num_rows) {
    start_plot <- (i - 1) * 5 + 1
    end_plot <- min(i * 5, num_chunks)
    row_plots[[i]] <- patchwork::wrap_plots(plot_list[start_plot:end_plot], ncol = 5)
  }
  
  # Combine rows into a final plot
  final_plot <- patchwork::wrap_plots(row_plots, nrow = num_rows)
  
  # Add a title to the final plot
  final_plot <- final_plot + 
    patchwork::plot_annotation(
      title = "Expression Abundance Across Clusters",
      theme = theme(plot.title = element_text(hjust = 0.5, size = 16))
    )
  
  return(final_plot)
}

###########################

