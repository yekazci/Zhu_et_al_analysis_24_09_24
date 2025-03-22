# handle non-found TFs gracefully:

analyze_multiple_tf_footprints <- function (tf_names, peaks_gr, fragments, metadata, pwms_mat, genes, cell_type_col_name, line_thickness = 1.5) {
    suppressPackageStartupMessages({
        library(stringr)
        library(motifmatchr)
        library(chromVARmotifs)
        library(ggplot2)
        library(patchwork)
    })
    if (!is.character(tf_names)) {
        stop("tf_names must be a character vector of transcription factor names.")
    }

    motif_matches <- unlist(lapply(tf_names, function(tf) {
        str_subset(names(pwms_mat), pattern = tf)
    }))

    if (length(motif_matches) == 0) {
        warning(paste("No motifs found for TF:", paste(tf_names, collapse=", "), "Skipping..."))
        return(NULL)  # Return NULL instead of stopping execution
    }

    selected_motifs <- setNames(motif_matches, tf_names)
    message("Motifs selected: ", paste(names(selected_motifs), collapse = ", "))

    suppressWarnings({
        motif_positions <- motifmatchr::matchMotifs(pwms_mat[selected_motifs], peaks_gr, genome = "hg38", out = "positions")
    })

    names(motif_positions) <- names(selected_motifs)
    
    footprinting_plots <- lapply(names(selected_motifs), function(motif) {
        plot_tf_footprint(fragments, motif_positions[[motif]], 
            cell_groups = metadata[, cell_type_col_name], zero_based_coords = !is(genes, "GRanges"), 
            flank = 250, smooth = 2) + 
            ggplot2::labs(title = motif) + 
            ggplot2::theme_minimal() + 
            ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold")) + 
            ggplot2::guides(color = ggplot2::guide_legend(title = "Cell Type")) + 
            ggplot2::geom_line(linewidth = line_thickness)
    })

    if (length(footprinting_plots) == 1) {
        return(footprinting_plots[[1]])
    } else {
        return(patchwork::wrap_plots(footprinting_plots, guides = "collect", ncol = 1))
    }
}

###################################

plot_multiple_TFs_w_highlights <- function (tf_genes, 
                                            peaks_gr, 
                                            fragments, 
                                            metadata,
                                            cell_type_col_name,
                                            pwms_mat, 
                                            genes, 
                                            highlight_cell_type, 
                                            line_thickness = 1, 
                                            set_theme, 
                                            set_margin) {
    
    plot_list <- list()
    
    for (each_TF in tf_genes) {
        p1 <- analyze_multiple_tf_footprints(    tf_names = each_TF, 
                                                 peaks_gr = peaks_gr, 
                                                 fragments = fragments, 
                                                 metadata = metadata, 
                                                 pwms_mat = pwms_mat, 
                                                 cell_type_col_name = cell_type_col_name,
                                                 genes = genes, 
                                                 line_thickness = line_thickness  )

        if (is.null(p1)) {
            next  # Skip to the next TF if no match was found
        }

        highlight_colours <- c(rep(x = "#D3D3D380", times = length(unique(metadata[, cell_type_col_name]))))
        highlight_colours[sort(levels(metadata[, cell_type_col_name])) %in% highlight_cell_type] <- "#F6222EFF"
        p2 <- p1 & ggplot2::scale_color_manual(values = highlight_colours)
        
        plot_list[[each_TF]] <- p1 + p2 & set_theme & set_margin
    }
    
    return(plot_list)
}