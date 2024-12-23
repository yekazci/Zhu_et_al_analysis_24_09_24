
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
    select(regions)


eGRN_regions_GRanges <- 
    eGRN_regions %>% 
    pull(regions) %>% 
    str_split(pattern = ';') %>% 
    unlist %>% 
    unique %>% 
    Signac::StringToGRanges(sep = c('-', '-'))

 return(eGRN_regions_GRanges)   
    
}



