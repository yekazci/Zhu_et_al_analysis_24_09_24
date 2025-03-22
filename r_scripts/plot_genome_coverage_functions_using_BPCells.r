# add highlight buffer as parameter:

library(ggplot2)

plot_detailed_genome_coverage_w_BULKED_EXPRESSION <- function(fragments, 
                                                              my_region, 
                                                              tf_name, 
                                                              target_name, 
                                                              bulked_rna_data, 
                                                              cell_types, 
                                                              total_reads_per_cell, 
                                                              bins_value, 
                                                              extend.bp, 
                                                              genes_info, 
                                                              transcripts_info, 
                                                              highlight_buffer = 50) {
    
    # Parse the user-provided region ('chr:start-end')
    original_region <- list(
        chr = str_split_i(my_region, pattern = ':', i = 1), 
        start = as.numeric(str_split_i(str_split_i(my_region, pattern = ':', i = 2), pattern = '-', i = 1)),
        end = as.numeric(str_split_i(str_split_i(my_region, pattern = ':', i = 2), pattern = '-', i = 2))
    )   
    
    # Expand the region
    region <- list(
        chr = original_region$chr,
        start = max(0, original_region$start - extend.bp),  # Prevent negative start values
        end = original_region$end + extend.bp
    )

    # Print the expanded and original regions for debugging
    print(paste("Expanded region:", region$chr, region$start, region$end))
    print(paste("Original region:", original_region$chr, original_region$start, original_region$end))

    # Generate coverage plot
    coverage_plot <- trackplot_coverage(
        fragments, 
        region = region, 
        groups = cell_types,
        total_reads_per_cell,
        bins = bins_value
    )

    
 # All these commands below are different working versions of the annotaion command.
    # I will keep them, however, I un-commented only the one that I currently use.

    highlight_buffer <- highlight_buffer  # Expand each side by.
coverage_plot <- coverage_plot +
    annotate(
        "rect",
        xmin = original_region$start - highlight_buffer, 
        xmax = original_region$end + highlight_buffer, 
        ymin = -Inf, 
        ymax = Inf,
        fill = "red", 
        alpha = 0.2
    )

  # Other alternative annotation commands:
    
    # Annotate the original region as a red shaded box using `annotate()`
    # coverage_plot <- coverage_plot +
    #     annotate(
    #         "rect",
    #         xmin = original_region$start, 
    #         xmax = original_region$end, 
    #         ymin = -Inf, 
    #         ymax = Inf,
    #         fill = "red", alpha = 0.2  # Semi-transparent red
    #     ) +
    #     annotate(
    #         "text",
    #         x = (original_region$start + original_region$end) / 2, 
    #         y = max(total_reads_per_cell, na.rm = TRUE) * 0.9, 
    #         label = "Input Region",
    #         color = "red", size = 5, fontface = "bold"
    #     )

########################################################
        # coverage_plot <- coverage_plot +
        # annotate(
        #     "rect",
        #     xmin = original_region$start, 
        #     xmax = original_region$end, 
        #     ymin = -Inf, 
        #     ymax = Inf,
        #     fill = "red", alpha = 0.2
        # ) +
        # annotate(
        #     "text",
        #     x = (original_region$start + original_region$end) / 2, 
        #     y = max(total_reads_per_cell, na.rm = TRUE) * 0.5,  # Reduce from 0.9 to 0.5
        #     label = "Input Region",
        #     color = "red", size = 5, fontface = "bold"
        # )

    #####################################################################
    
    # coverage_plot <- coverage_plot +
    # # Highlight the original region with a semi-transparent red rectangle
    # annotate(
    #     "rect",
    #     xmin = original_region$start, 
    #     xmax = original_region$end, 
    #     ymin = -Inf, 
    #     ymax = Inf,
    #     fill = "red", alpha = 0.2  # Light red shading for visibility
    # ) +
    # # Add a text label at a fixed y-position
    # annotate(
    #     "text",
    #     x = (original_region$start + original_region$end) / 2,  # Center of original region
    #     y = 50,  # Fixed y-position (adjust as needed)
    #     label = "Input Region",
    #     color = "red", size = 5, fontface = "bold"
    # )
    
####################################################################
    # coverage_plot <- coverage_plot +
    # geom_vline(xintercept = c(original_region$start, original_region$end), 
    #            color = "red", linetype = "dashed", size = 1) + 
    # annotate("text", 
    #          x = (original_region$start + original_region$end) / 2,  
    #          y = 50,
    #          label = "",
    #          color = "red", 
    #          size = 5, 
    #          fontface = "bold")

############################################################
    
    # Apply ggplot2 theme adjustments
    p1 <- coverage_plot + my_margin + ggplot2::theme(
        strip.text = ggplot2::element_text(size = 11, face = "bold"),
        axis.title.x = ggplot2::element_text(size = 14),  
        axis.title.y = ggplot2::element_text(size = 14),  
        axis.text.x  = ggplot2::element_text(size = 14, angle = 0, hjust = 1),  
        axis.text.y  = ggplot2::element_text(size = 14),  
        legend.text  = ggplot2::element_text(size = 15)
    ) +  ggtitle(label = paste0("TF: ", tf_name, ", target: ", target_name), subtitle = paste0("linked_region: ", my_region)) +  # Add the title
    theme(
        plot.title = element_text(
            hjust = 0.5,  # Center the title horizontally
            vjust = 1,    # Adjust the vertical position if needed
            size = 16,    # Set the font size
            face = "bold" # Make the title bold (optional)
        ),
    )

    # Generate gene and scalebar plots
                 
    gene_plot <- trackplot_gene(transcripts_info, region)
    scalebar_plot <- trackplot_scalebar(region)

    # retrieve gene expression value of a specific feature:

    gene_df <- mat_rna_norm_pseudobulk[target_name,] %>% data.frame(gene = .)

    # I needed to add the following to ensure the order of cell types match between tracks and bar plots.

    custom_order <- cell_types %>% levels() # cell_types needs to be FACTOR !
    
    gene_df$x_ordered <- factor(rownames(gene_df), levels = custom_order)

    # Now, it will plot in the order of levels:
    
    expression_plot <- ggplot(gene_df, aes(x = x_ordered, y = gene, fill = x_ordered)) +
    geom_bar(stat = "identity") +
    ggplot2::guides(y="none", fill="none") + 
    ggplot2::labs(x=NULL, y="log1p(PseudoBulked_RNA") +
    ggplot2::scale_fill_manual(values=discrete_palette("stallion"), drop=FALSE) +
    BPCells:::trackplot_theme()
    
    # Combine plots
    p2 <- trackplot_combine(
        list(
            scalebar_plot, 
            p1, 
            gene_plot + ggplot2::guides(color = "none")
        ),
        side_plot = expression_plot
    )

    return(p2)
}

###################################################

plot_gene_coverage_w_BULK_expression <- function(fragments, 
                                                 gene = NULL, 
                                                 cell_types, 
                                                 bulked_rna_data, 
                                                 total_reads_per_cell, 
                                                 bins_value, 
                                                 extend.bp, 
                                                 genes_info, 
                                                 transcripts_info) {

        
   region <- gene_region(genes_info, gene, extend_bp =  extend.bp)     

        coverage_plot <- trackplot_coverage(
  fragments, # frags
  region = region, 
  groups=cell_types,
  total_reads_per_cell,
  bins=bins_value
)

        # expand_by: 
    
    #         granges_obj <- 
    # GenomicRanges::resize(granges_obj, 
    #                       width = GenomicRanges::width(granges_obj) + 2 * expand_by, 
    #                       fix = "center")
    

    p1 <- coverage_plot + my_margin + ggplot2::theme(
    strip.text = ggplot2::element_text(size = 11, face = "bold"),  # Increase facet label (cell type) font size
    axis.title.x = ggplot2::element_text(size = 14),  
    axis.title.y = ggplot2::element_text(size = 14),  
    axis.text.x  = ggplot2::element_text(size = 14, angle = 0, hjust = 1),  
    axis.text.y  = ggplot2::element_text(size = 14),  
    legend.text  = ggplot2::element_text(size = 15)   
) +  ggtitle(gene) +  # Add the title
    theme(
        plot.title = element_text(
            hjust = 0.5,  # Center the title horizontally
            vjust = 1,    # Adjust the vertical position if needed
            size = 16,    # Set the font size
            face = "bold" # Make the title bold (optional)
        )
    )

    gene_plot <- trackplot_gene(transcripts_info, region)
    scalebar_plot <- trackplot_scalebar(region)


    # retrieve gene expression value of a specific feature:

    gene_df <- mat_rna_norm_pseudobulk[gene,] %>% data.frame(gene = .)

    # I needed to add the following to ensure the order of cell types match between tracks and bar plots.

    custom_order <- cell_types %>% levels() # cell_types needs to be FACTOR !
    
    gene_df$x_ordered <- factor(rownames(gene_df), levels = custom_order)

    # Now, it will plot in the order of levels:
    
    expression_plot <- ggplot(gene_df, aes(x = x_ordered, y = gene, fill = x_ordered)) +
    geom_bar(stat = "identity") +
    ggplot2::guides(y="none", fill="none") + 
    ggplot2::labs(x=NULL, y="log1p(PseudoBulked_RNA") +
    ggplot2::scale_fill_manual(values=discrete_palette("stallion"), drop=FALSE) +
    BPCells:::trackplot_theme()
    
    p2 <- 
trackplot_combine(
  list(
    scalebar_plot, 
    p1, 
    gene_plot + ggplot2::guides(color="none")
  ), side_plot = expression_plot
    
)

    return(p2)

}   

#########################################################

# add highlight buffer as parameter:

library(ggplot2)

plot_detailed_genome_coverage_w_EXPRESSION <- function(fragments, my_region, tf_name, target_name, rna_data, cell_types, total_reads_per_cell, bins_value, extend.bp, genes_info, transcripts_info, highlight_buffer = 50) {
    
    # Parse the user-provided region ('chr:start-end')
    original_region <- list(
        chr = str_split_i(my_region, pattern = ':', i = 1), 
        start = as.numeric(str_split_i(str_split_i(my_region, pattern = ':', i = 2), pattern = '-', i = 1)),
        end = as.numeric(str_split_i(str_split_i(my_region, pattern = ':', i = 2), pattern = '-', i = 2))
    )   
    
    # Expand the region
    region <- list(
        chr = original_region$chr,
        start = max(0, original_region$start - extend.bp),  # Prevent negative start values
        end = original_region$end + extend.bp
    )

    # Print the expanded and original regions for debugging
    print(paste("Expanded region:", region$chr, region$start, region$end))
    print(paste("Original region:", original_region$chr, original_region$start, original_region$end))

    # Generate coverage plot
    coverage_plot <- trackplot_coverage(
        fragments, 
        region = region, 
        groups = cell_types,
        total_reads_per_cell,
        bins = bins_value
    )

    
 # All these commands below are different working versions of the annotaion command.
    # I will keep them, however, I un-commented only the one that I currently use.

    highlight_buffer <- highlight_buffer  # Expand each side by.
coverage_plot <- coverage_plot +
    annotate(
        "rect",
        xmin = original_region$start - highlight_buffer, 
        xmax = original_region$end + highlight_buffer, 
        ymin = -Inf, 
        ymax = Inf,
        fill = "red", 
        alpha = 0.2
    )

  # Other alternative annotation commands:
    
    # Annotate the original region as a red shaded box using `annotate()`
    # coverage_plot <- coverage_plot +
    #     annotate(
    #         "rect",
    #         xmin = original_region$start, 
    #         xmax = original_region$end, 
    #         ymin = -Inf, 
    #         ymax = Inf,
    #         fill = "red", alpha = 0.2  # Semi-transparent red
    #     ) +
    #     annotate(
    #         "text",
    #         x = (original_region$start + original_region$end) / 2, 
    #         y = max(total_reads_per_cell, na.rm = TRUE) * 0.9, 
    #         label = "Input Region",
    #         color = "red", size = 5, fontface = "bold"
    #     )

########################################################
        # coverage_plot <- coverage_plot +
        # annotate(
        #     "rect",
        #     xmin = original_region$start, 
        #     xmax = original_region$end, 
        #     ymin = -Inf, 
        #     ymax = Inf,
        #     fill = "red", alpha = 0.2
        # ) +
        # annotate(
        #     "text",
        #     x = (original_region$start + original_region$end) / 2, 
        #     y = max(total_reads_per_cell, na.rm = TRUE) * 0.5,  # Reduce from 0.9 to 0.5
        #     label = "Input Region",
        #     color = "red", size = 5, fontface = "bold"
        # )

    #####################################################################
    
    # coverage_plot <- coverage_plot +
    # # Highlight the original region with a semi-transparent red rectangle
    # annotate(
    #     "rect",
    #     xmin = original_region$start, 
    #     xmax = original_region$end, 
    #     ymin = -Inf, 
    #     ymax = Inf,
    #     fill = "red", alpha = 0.2  # Light red shading for visibility
    # ) +
    # # Add a text label at a fixed y-position
    # annotate(
    #     "text",
    #     x = (original_region$start + original_region$end) / 2,  # Center of original region
    #     y = 50,  # Fixed y-position (adjust as needed)
    #     label = "Input Region",
    #     color = "red", size = 5, fontface = "bold"
    # )
    
####################################################################
    # coverage_plot <- coverage_plot +
    # geom_vline(xintercept = c(original_region$start, original_region$end), 
    #            color = "red", linetype = "dashed", size = 1) + 
    # annotate("text", 
    #          x = (original_region$start + original_region$end) / 2,  
    #          y = 50,
    #          label = "",
    #          color = "red", 
    #          size = 5, 
    #          fontface = "bold")

############################################################
    
    # Apply ggplot2 theme adjustments
    p1 <- coverage_plot + my_margin + ggplot2::theme(
        strip.text = ggplot2::element_text(size = 11, face = "bold"),
        axis.title.x = ggplot2::element_text(size = 14),  
        axis.title.y = ggplot2::element_text(size = 14),  
        axis.text.x  = ggplot2::element_text(size = 14, angle = 0, hjust = 1),  
        axis.text.y  = ggplot2::element_text(size = 14),  
        legend.text  = ggplot2::element_text(size = 15)
    ) +  ggtitle(label = paste0("TF: ", tf_name, ", target: ", target_name), subtitle = paste0("linked_region: ", my_region)) +  # Add the title
    theme(
        plot.title = element_text(
            hjust = 0.5,  # Center the title horizontally
            vjust = 1,    # Adjust the vertical position if needed
            size = 16,    # Set the font size
            face = "bold" # Make the title bold (optional)
        ),
    )

    # Generate gene and scalebar plots
                 
    gene_plot <- trackplot_gene(transcripts_info, region)
    scalebar_plot <- trackplot_scalebar(region)

    # retrieve genbe expression value of a specific feature:

    expression <- collect_features(rna_data, target_name)
    names(expression) <- "gene"

    expression_plot <- ggplot2::ggplot(expression, ggplot2::aes(metadata$main_cell_types, gene, fill=metadata$main_cell_types)) +
    ggplot2::geom_boxplot() + 
    ggplot2::guides(y="none", fill="none") + 
    ggplot2::labs(x=NULL, y=paste0(target_name, " RNA expression")) +
    ggplot2::scale_fill_manual(values=discrete_palette("stallion"), drop=FALSE) +
    BPCells:::trackplot_theme()
    
    # Combine plots
    p2 <- trackplot_combine(
        list(
            scalebar_plot, 
            p1, 
            gene_plot + ggplot2::guides(color = "none")
        ),
        side_plot = expression_plot
    )

    return(p2)
}

##############################################################

plot_gene_coverage_w_expression <- function(fragments, gene = NULL, cell_types, rna_data, total_reads_per_cell, bins_value, extend.bp, genes_info, transcripts_info) {

        
   region <- gene_region(genes_info, gene, extend_bp =  extend.bp)     

        coverage_plot <- trackplot_coverage(
  fragments, # frags
  region = region, 
  groups=cell_types,
  total_reads_per_cell,
  bins=bins_value
)

        # expand_by: 
    
    #         granges_obj <- 
    # GenomicRanges::resize(granges_obj, 
    #                       width = GenomicRanges::width(granges_obj) + 2 * expand_by, 
    #                       fix = "center")
    

    p1 <- coverage_plot + my_margin + ggplot2::theme(
    strip.text = ggplot2::element_text(size = 11, face = "bold"),  # Increase facet label (cell type) font size
    axis.title.x = ggplot2::element_text(size = 14),  
    axis.title.y = ggplot2::element_text(size = 14),  
    axis.text.x  = ggplot2::element_text(size = 14, angle = 0, hjust = 1),  
    axis.text.y  = ggplot2::element_text(size = 14),  
    legend.text  = ggplot2::element_text(size = 15)   
) +  ggtitle(gene) +  # Add the title
    theme(
        plot.title = element_text(
            hjust = 0.5,  # Center the title horizontally
            vjust = 1,    # Adjust the vertical position if needed
            size = 16,    # Set the font size
            face = "bold" # Make the title bold (optional)
        )
    )

    gene_plot <- trackplot_gene(transcripts_info, region)
    scalebar_plot <- trackplot_scalebar(region)


    # retrieve genbe expression value of a specific feature:

    	expression <- collect_features(rna_data, gene)
        names(expression) <- "gene"

    expression_plot <- ggplot2::ggplot(expression, ggplot2::aes(metadata$main_cell_types, gene, fill=metadata$main_cell_types)) +
    ggplot2::geom_boxplot() + 
    ggplot2::guides(y="none", fill="none") + 
    ggplot2::labs(x=NULL, y=paste0(gene, " RNA expression")) +
    ggplot2::scale_fill_manual(values=discrete_palette("stallion"), drop=FALSE) +
    BPCells:::trackplot_theme()
    
    p2 <- 
trackplot_combine(
  list(
    scalebar_plot, 
    p1, 
    gene_plot + ggplot2::guides(color="none")
  ), side_plot = expression_plot
    
)

    return(p2)

}   