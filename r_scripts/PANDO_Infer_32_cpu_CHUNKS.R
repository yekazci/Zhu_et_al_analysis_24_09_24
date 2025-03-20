# GRN inference:

# Pando provides multiple modelling methods, and I used the default method: 
# "The default option when running infer_grn() is a generalized linear model 
# (GLM) with gaussian noise. Using the family argument, one can choose 
# other noise models, e.g to fit directly on counts instead of on log-normalized 
# data."


# JOB script:

main_dir <- 
"/fast/AG_Bunina/Yusuf/Project_Endothelial_and_Stroke/Datasets/Chromatin_and_Gene_Exp/2023_Zhu_et_al/Zhu_et_al_2023/24_09_24"

here::here(main_dir)

# load the R environment with the necessary packages such as Epiregulon:

my_epiregulon_lib <- "/fast/AG_Bunina/Yusuf/Project_Endothelial_and_Stroke/Datasets/Chromatin_and_Gene_Exp/2024_C_A_Mannens_C_et_al/04_02_25/renv/library/linux-rhel-9.4/R-4.4/x86_64-unknown-linux-gnu"

.libPaths(new = my_epiregulon_lib, include.site = FALSE)

.libPaths()

library(tidyverse)

Sys.sleep(150)

message("tidyverse loaded!")

library(Seurat)

Sys.sleep(150)

message("Seurat loaded!")

library(Signac)

Sys.sleep(150)

message("Signac loaded!")

library(Pando)

Sys.sleep(150)

message("Pando loaded!")

library(GenomicRanges)

Sys.sleep(150)

message("GenomicRanges loaded!")

library(doParallel)

Sys.sleep(150)

message("doParallel loaded!")

registerDoParallel(32)

file_name <- 'Zhu_et_al_Pando_w_motifs.RDS'

input_object <- 
    readRDS(here::here('r_objects', file_name))

Sys.sleep(15)

message("Input seurat object loaded!")

gene_names <-input_object@data %>% rownames()

# Function to retrieve elements in chunks with a 5-gene overlap
retrieve_in_chunks <- function(gene_names, chunk_size, overlap = 5) {
  
  # Get the total number of gene names provided in the input
  n <- length(gene_names)  

  # Define the step size (how much to move forward after creating each chunk).
  # The step size is calculated by subtracting the overlap from the chunk size.
  # This ensures that each successive chunk shares 'overlap' number of genes with the previous chunk.
  step_size <- chunk_size - overlap  

  # Initialize an empty list to store the resulting chunks of gene names.
  chunk_list <- list()  

  # Use a for-loop to iterate through the gene names vector in steps of 'step_size'.
  # The 'seq()' function generates a sequence of starting points (i) for each chunk.
  # The loop increments by 'step_size' to determine the next chunk’s starting point.
  for (i in seq(1, n, by = step_size)) {
    
    # Extract the current chunk of gene names from the vector.
    # The 'min()' function ensures that the chunk doesn't go beyond the last gene in 'gene_names'.
    # The chunk starts at 'i' and ends at 'i + chunk_size - 1' (or the last available gene).
    chunk <- gene_names[i:min(i + chunk_size - 1, n)]
    
    # Generate a name for each chunk (optional, but helpful for referencing).
    # The 'ceiling()' function calculates which chunk we are in based on 'i' and 'step_size'.
    chunk_name <- paste0("chunk_", ceiling(i / step_size))
    
    # Add the current chunk to the 'chunk_list' under its generated name.
    # This allows for easy access to each chunk by its name.
    chunk_list[[chunk_name]] <- chunk
  }
  
  # Return the full list of chunks. Each element of the list corresponds to a chunk of gene names.
  return(chunk_list)  
}


gene_chunks <- retrieve_in_chunks(gene_names = gene_names, 
                                  chunk_size = 500,
                                  overlap = 5)


 output_file_name <- file_name %>% str_split_i(pattern = 'et_al', i = 1)

for (i in 1:length(gene_chunks)) { 

    gene_names <- gene_chunks[[i]]

    input_obj_w_eGRN <- infer_grn(
                                      input_object,
                                      genes = gene_names,
                                      peak_to_gene_method = 'GREAT',
                                      parallel = T)

    # make sure the 'pando_eGRN_chunks_final' folder exists !!!.

    input_obj_w_eGRN %>% 
                saveRDS(here::here('r_objects', 'pando_eGRN_chunks_final',
                             paste0(output_file_name, 'et_al_w_eGRNs_', 
                                                    'SLURM', '_c', i , '.RDS')))

    
    message("Saved R object for chunk ", i ,": FINISHED!")
    
    rm(input_obj_w_eGRN)
    
}


cat("\n STATUS: COMPLETED !!! \n")
