message("Loading traseR package...")

library(traseR)

message("Loading taSNP database...")

data(taSNP)

message("Loading taSNPLD database which additionally contains LD SNPs...")

data(taSNPLD)

message("Loading CEU database which contains all known SNPs...")

data(CEU)

# generic function to calculate enriched SNPs with background of whole genome, and all SNPs.

enrich_SNP_generic <- function(granges_obj, SNPs, expand_by = NULL, suppress_messages = NULL) {

    if (!is.null(expand_by)) {
    
        granges_obj <- 
    GenomicRanges::resize(granges_obj, 
                          width = GenomicRanges::width(granges_obj) + 2 * expand_by, 
                          fix = "center")
    
    }
    
    
    if (!is.null(suppress_messages)) {
        
         message('Using whole genome as background in the statistical analysis....  \n')
        
         res_w_genome_as_bg  <- suppressMessages(traseR(snpdb = SNPs, 
                                                       region = granges_obj)) 
        
         cat('\n')
    
         cat('\n')
        
         message('Using all known SNPs as background in the statistical analysis....  \n')
        
         res_w_all_SNPs_as_bg  <- suppressMessages(traseR(snpdb = SNPs,  
                                                             snpdb.bg = CEU,  
                                                             region = granges_obj))  
        
    } 
        
     
    if (is.null(suppress_messages)) {
    
    message('Using whole genome as background in the statistical analysis....  \n')
        
         res_w_genome_as_bg  <- traseR(snpdb = SNPs, 
                                        region = granges_obj)
            
     cat('\n')
    
     cat('\n')
    
     message('Using all known SNPs as background in the statistical analysis....  \n')
    
         res_w_all_SNPs_as_bg  <- traseR(snpdb = SNPs,  
                                         snpdb.bg = CEU,  
                                         region = granges_obj)    
        
    
    }
    
    
    names(res_w_genome_as_bg) <- c("tb.all", "traits", "trait_classes", "ntraits", "ntraitclass")
    
    names(res_w_all_SNPs_as_bg) <- c("tb.all", "traits", "trait_classes", "ntraits", "ntraitclass")
    
    return(list(all_genome_as_bg = res_w_genome_as_bg, all_SNPs_as_bg = res_w_all_SNPs_as_bg))
  
}

#########################

# Add extra parameter to optionally exclude SNPs with hits less than a cut off:

plot_traits <- function(enriched_SNP_obj, qvalue, SNP_hit_cut_off = NULL) {
    
  if (!is.null(SNP_hit_cut_off))
      
      {
      
      df <- enriched_SNP_obj$traits %>% 
                    filter(taSNP.hits > 0)
      
      } else {
      
      df  <- enriched_SNP_obj$traits
      
      }
    
    
  df <- df %>%
  mutate(log2_odds_ratio = log2(odds.ratio),
         neg_log10_q = -log10(q.value),
         significant = ifelse(q.value < qvalue, "Significant", "Not Significant"))  

  ggplot(df, aes(x = log2_odds_ratio, y = neg_log10_q)) +
  geom_point(aes(size = taSNP.hits, color = significant)) +  # Points colored by significance
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "blue")) +
  ggrepel::geom_text_repel(aes(label = Trait), max.overlaps = 10, size = 4, box.padding = 0.5) +  # Repel to avoid overlap
  theme_minimal() +  # Clean theme
  labs(title = "SNP Overrepresentation Analysis", 
       x = "log2(Odds Ratio)", 
       y = "-log10(q-value)", 
       size = "Number of SNP Hits", 
       color = "Significance") +
  theme(plot.title = element_text(hjust = 0.5))
                       
}

plot_trait_classes <- function(enriched_SNP_obj, qvalue, SNP_hit_cut_off = NULL) {
    
  if (!is.null(SNP_hit_cut_off))
      
      {
      
      df <- enriched_SNP_obj$trait_classes %>% 
                    filter(taSNP.hits > 0)
      
      }  else {
      
      df  <- enriched_SNP_obj$trait_classes
      
      }
    
    
      df <- df %>%
  mutate(log2_odds_ratio = log2(odds.ratio),
         neg_log10_q = -log10(q.value),
         significant = ifelse(q.value < qvalue, "Significant", "Not Significant"))  
    
  ggplot(df, aes(x = log2_odds_ratio, y = neg_log10_q)) +
  geom_point(aes(size = taSNP.hits, color = significant)) +  # Points colored by significance
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "blue")) +
  ggrepel::geom_text_repel(aes(label = Trait_Class), max.overlaps = 10, size = 4, box.padding = 0.5) +  # Repel to avoid overlap
  theme_minimal() +  # Clean theme
  labs(title = "SNP Overrepresentation Analysis", 
       x = "log2(Odds Ratio)", 
       y = "-log10(q-value)", 
       size = "Number of SNP Hits", 
       color = "Significance") +
  theme(plot.title = element_text(hjust = 0.5))
                       
}

############

plot_for_diff_expansions <- function(enriched_SNP_obj, expansion_value_string = expansion_value) { 

p1 <-
    
plot_traits(enriched_SNP_obj = enriched_SNP_obj$all_genome_as_bg, qvalue = 0.05) +
    ggtitle(label = glue::glue('Genomic regions expanded by {expansion_value_string} from each end, with whole genome as background')) +
    theme(plot.margin = unit(c(2, 2, 2, 2), "lines")) # add some margins to top, right, bottom, left.
    
p2 <- 
    
plot_traits(enriched_SNP_obj = enriched_SNP_obj$all_SNPs_as_bg, qvalue = 0.05, SNP_hit_cut_off = 0) +
    ggtitle(label = glue::glue('Genomic regions expanded by {expansion_value_string} from each end, with all SNPs as background')) +
    theme(plot.margin = unit(c(2, 2, 2, 2), "lines"))

p3 <- 
    
plot_trait_classes(enriched_SNP_obj = enriched_SNP_obj$all_genome_as_bg, qvalue = 0.05) +
    ggtitle(label = glue::glue('Genomic regions expanded by {expansion_value_string} from each end. with whole genome as background')) +
    theme(plot.margin = unit(c(2, 2, 2, 2), "lines"))

p4 <- 
    
plot_trait_classes(enriched_SNP_obj = enriched_SNP_obj$all_SNPs_as_bg, qvalue = 0.05, SNP_hit_cut_off = 0) +
    ggtitle(label = glue::glue('Genomic regions expanded by {expansion_value_string} from each end, with all SNPs as background')) +
    theme(plot.margin = unit(c(2, 2, 2, 2), "lines"))

# Suppress warnings while printing each plot, to get rid off extra ggrepel warnings:
    
    
 suppressWarnings((gridExtra::grid.arrange(grobs = list(p1,p2,p3,p4), ncol = 2)))
    
 cat('\n')
    
 cat('\n')
    

#     suppressWarnings(print(p1))
#     cat('\n')
    
#     Sys.sleep(1)
    
#     suppressWarnings(print(p2))
#     cat('\n')
    
#     Sys.sleep(1)
    
#     suppressWarnings(print(p3))
#     cat('\n')
    
#     Sys.sleep(1)
    
#     suppressWarnings(print(p4))
#     cat('\n')

}

######################

plot_for_diff_expansions_for_saving <- plot_for_diff_expansions <- function(enriched_SNP_obj, expansion_value_string = expansion_value) { 

p1 <-
    
plot_traits(enriched_SNP_obj = enriched_SNP_obj$all_genome_as_bg, qvalue = 0.05) +
    ggtitle(label = glue::glue('Genomic regions expanded by {expansion_value_string} from each end, with whole genome as background')) +
    theme(plot.margin = unit(c(2, 2, 2, 2), "lines")) # add some margins to top, right, bottom, left.
    
p2 <- 
    
plot_traits(enriched_SNP_obj = enriched_SNP_obj$all_SNPs_as_bg, qvalue = 0.05, SNP_hit_cut_off = 0) +
    ggtitle(label = glue::glue('Genomic regions expanded by {expansion_value_string} from each end, with all SNPs as background')) +
    theme(plot.margin = unit(c(2, 2, 2, 2), "lines"))

p3 <- 
    
plot_trait_classes(enriched_SNP_obj = enriched_SNP_obj$all_genome_as_bg, qvalue = 0.05) +
    ggtitle(label = glue::glue('Genomic regions expanded by {expansion_value_string} from each end. with whole genome as background')) +
    theme(plot.margin = unit(c(2, 2, 2, 2), "lines"))

p4 <- 
    
plot_trait_classes(enriched_SNP_obj = enriched_SNP_obj$all_SNPs_as_bg, qvalue = 0.05, SNP_hit_cut_off = 0) +
    ggtitle(label = glue::glue('Genomic regions expanded by {expansion_value_string} from each end, with all SNPs as background')) +
    theme(plot.margin = unit(c(2, 2, 2, 2), "lines"))

#     suppressWarnings(print(p1))
#     cat('\n')
    
#     Sys.sleep(1)
    
#     suppressWarnings(print(p2))
#     cat('\n')
    
#     Sys.sleep(1)
    
#     suppressWarnings(print(p3))
#     cat('\n')
    
#     Sys.sleep(1)
    
#     suppressWarnings(print(p4))
#     cat('\n')
    
# Combine plots into a grid without rendering it

    combined_plot <- gridExtra::arrangeGrob(grobs = list(p1, p2, p3, p4), ncol = 2)
    
   return(combined_plot)

}