
message("make sure you installed the GGally, network, igraph and ggraph packages :)")

message("First option: One can directly use only the node_names parameter to specify 
          the gene names to visualize in the graph")

message("Second option: Use TF-names and/or target genes parameters 
         to visualize their connections")


library(GGally)
library(network)

vis_network_w_GGNET2 <- function(meta_DF, score_cut = NULL, 
                          target_genes = NULL,
                          TF_names = NULL, 
                          node_names = NULL,
                          seed_value = 123,
                          label_size = 3,
                          arrow = NULL,
                          override_node_size_value = NULL,
                          subnetwork_default_node_and_label_size = NULL,
                          degree = NULL) {
    

  # To filter the vertices by user-defined degree threshold:  
    
  if (!is.null(degree)) {
      
# Below, I convert the DF to a vector by unlist,
# then, tabulate the frequencies and select the genes that exist more or equal to the input number of degree.  
      
      keep_nodes <- meta_DF %>% unlist %>% table %>% { .[. >= degree] } %>% names
      
      meta_DF <- meta_DF %>% filter(tf %in% keep_nodes | target %in% keep_nodes)
      
    }  
    
  # Identify genes that appear both as transcription factors (TFs) and as targets.
  
  target_TFs <- intersect(meta_DF$tf, meta_DF$target)  
    
  # Filter the network data frame to only include rows where 'tf' is either in node_names, or 
  # 'target' is in node_names, or both. 
    
  if(!is.null(node_names)) {
    
    meta_DF <- meta_DF %>% filter(tf %in% node_names | target %in% node_names)
    
  }  
    
  # Filter the network data frame to only include rows where 'target' is in target_genes, if provided.
  if (!is.null(target_genes)) {
    meta_DF <- meta_DF %>% filter(target %in% target_genes)
  }
  
  # Further filter the network data frame to only include rows where 'tf' is in TF_names, if provided.
  if (!is.null(TF_names)) {
    meta_DF <- meta_DF %>% filter(tf %in% TF_names)
  }
  
  # Filter based on the absolute value of the 'estimate' column, if score_cut is provided.
  if (!is.null(score_cut)) {
    meta_DF <- meta_DF %>% filter(abs(estimate) >= score_cut)
    score_cut <- paste0(' >= ', score_cut)  # Record the applied score cut for annotation.
  } else {
    score_cut <- 'unset'  # If no score_cut is provided, note that it's unset.
  }

    
   # Add check for meta_df having at least one row after the optional filtering steps above:
  
  if (nrow(meta_DF) <= 0) {
    stop("Error: meta_DF should have at least one row.")
  }
      
    
  # Modify the 'tf' and 'target' columns to append ".target-TF" to genes that are both TFs and targets.
  meta_DF$tf <- ifelse(meta_DF$tf %in% target_TFs, paste0(meta_DF$tf, ".target-TF"), meta_DF$tf)
  meta_DF$target <- ifelse(meta_DF$target %in% target_TFs, paste0(meta_DF$target, ".target-TF"), meta_DF$target)
  
  # Rescale the edge weights (based on the 'estimate' column) to a range suitable for visualization (0.5 to 3).
  meta_DF$weight <- scales::rescale(abs(meta_DF$estimate), to = c(0.5, 3))

  # Create a directed network object from the 'tf' and 'target' columns, interpreted as an edge list.
  net <- network(meta_DF[, c("tf", "target")], directed = TRUE, matrix.type = "edgelist")
    
  net %>% print()   
 
  # Assign the computed weights to the edges in the network.
  set.edge.value(x = net, attrname = "weight", value = meta_DF$weight)
  
  # Note: Below, I wrote the network package explicitly, because igraph package can mask these functions,
  # when used in the same R session.
    
  # Extract the names of the vertices (nodes) in the network.
  vertex_names <- network::get.vertex.attribute(x = net, attrname = 'vertex.names')

  # Determine the identity of each node as 'target', 'TF', or 'target-TF'.
  identity <- ifelse(vertex_names %in% paste0(target_TFs, ".target-TF"), "target-TF", 
                     ifelse(vertex_names %in% meta_DF$target, "target", "TF"))
  
  # Set the class attribute of each vertex to the determined identity ('target', 'TF', or 'target-TF').
  network::set.vertex.attribute(x = net, attrname = 'class', value = identity)
    
  
  # ALTERNATIVE:

# Network object contains a list called 'val'. This list consists of
# list elements. Each of this list elements stores the vertex names 
# as their secon variable.

# For instance the fist vertex name can be accessed via:

# net$val[[1]][[2]]  

# below code uses sapply to access the second element of each nested lists and
# generates a vector consisting of the second elements from each list.

# net %v% "class" <-  
#   ifelse(sapply(net$val,"[[",2) %in% meta_DF$target, "target", "tf")

#

# After we set the classes of vertices (nodes) as target and tf,
# we can now delete the '.' that we previously added to distinguish
# target genes from TFs. Because, there can be target genes that
# are also Transcription factors, we applied this above.
    
    
  # Clean up the vertex names by removing the ".target-TF" suffix.
  network::set.vertex.attribute(net, "vertex.names", 
                                gsub(pattern = ".target-TF", 
                                     replacement = "", 
                                     x = network::get.vertex.attribute(net, "vertex.names"), fixed = TRUE))

  # Assign edge colors based on whether the correlation ('estimate') is positive (red) or negative (blue).
  network::set.edge.attribute(net, "color", ifelse(meta_DF$estimate > 0, "firebrick", "steelblue"))

  # Define colors for the nodes based on their class.
  col <- c("TF" = "gray60", "target" = "sandybrown", "target-TF" = "#AC94F4")

    
# set.seed(seed_value)

# p1 <- 
#   ggnet2(net = net,
#        color.legend = "Class",
#        size.legend = "Degree",
#        label = TRUE,
#        size="class",
#        edge.size = "weight",
#        edge.color = "color",
#        edge.alpha = 0.75,
#        alpha = 0.7,
#        color="class",
#        size.palette = c("target" = 3, "tf" = 1),
#        label.size = 3,
#        palette=col,
#        legend.position = 'right') + 
#   coord_equal() +
#   theme(legend.text = element_text(size = 20))

# p1
    
# I could not add legend for the edge weights therefore I added the following attribute
# to the vertex to use only in the legend. It does not change anything in the structure
# of the network.

for_legend <- 
  ifelse(vertex_names %in% meta_DF$target, "activator", "repressor")

network::set.vertex.attribute(x = net, attrname = 'for_legend', value = for_legend)
  

    if(!is.null(subnetwork_default_node_and_label_size)) {
        
    override_node_size_value = 28
    label_size = 6  
  }    

    
# Set a random seed for reproducibility in the layout of the network visualization.  
    
set.seed(seed_value)

if (!is.null(arrow)) {
    
    p1 <- 
  ggnet2(net = net,
       color.legend = "Class",
       size.legend = "State",
       label = TRUE,
       size="for_legend",
       edge.size = "weight",
       edge.color = "color",
       edge.alpha = 0.75,
       alpha = 0.7,
       color="class",
       size.palette = c("activator" = 1, "repressor" = 1),
       label.size = label_size,
       palette=col,
       legend.position = 'right',
       arrow.size = 12, # newly added  for edge arrows.
       arrow.gap = 0.025) + # newly added  for edge arrows.
  coord_equal() +
  theme(legend.text = element_text(size = 20)) +
  ggtitle(label = 'TF modules-meta df based network graph')
    
    
  } else {

     p1 <- 
  ggnet2(net = net,
       color.legend = "Class",
       size.legend = "State",
       label = TRUE,
       size="for_legend", # this is just a workaround to add the edge colours to the legend. 
                          # It will be overriden with the two edge  colours after generation of the plot.
       edge.size = "weight",
       edge.color = "color",
       edge.alpha = 0.75,
       alpha = 0.7,
       color="class",
       size.palette = c("activator" = 1, "repressor" = 1),
       label.size = label_size,
       palette=col,
       legend.position = 'right') + 
  coord_equal() +
  theme(legend.text = element_text(size = 20)) +
  ggtitle(label = 'TF modules-meta df based network graph')

}
    

p1 <- p1 + guides(size = guide_legend(override.aes = list(color = c("firebrick", "steelblue")))) # Override the relevant legend for the edges.
    
    # To make the node sizes bigger when visualizing the subnetworks:

    if(!is.null(override_node_size_value)) {
        
     p1 <- p1 + scale_size_manual(values = c(override_node_size_value, override_node_size_value))  # It will make the edge link legends bigger !!!.
  
  }
    
return(p1)   
    
}  

#########################################

vis_subnetwork_w_GGNET2 <- function(meta_DF, score_cut = NULL, 
                          target_genes = NULL,
                          TF_names = NULL, 
                          node_names = NULL,
                          seed_value = 123,
                          label_size = 6,
                          arrow = NULL,
                          override_node_size_value = 28,
                          degree = NULL) {
    

  # To filter the vertices by user-defined degree threshold:  
    
  if (!is.null(degree)) {
      
# Below, I convert the DF to a vector by unlist,
# then, tabulate the frequencies and select the genes that exist more or equal to the input number of degree.  
      
      keep_nodes <- meta_DF %>% unlist %>% table %>% { .[. >= degree] } %>% names
      
      meta_DF <- meta_DF %>% filter(tf %in% keep_nodes | target %in% keep_nodes)
      
    }  
  
  # Identify genes that appear both as transcription factors (TFs) and as targets.
  
  target_TFs <- intersect(meta_DF$tf, meta_DF$target)  
    
  # Filter the network data frame to only include rows where 'tf' is either in node_names, or 
  # 'target' is in node_names, or both. 
    
  if(!is.null(node_names)) {
    
    meta_DF <- meta_DF %>% filter(tf %in% node_names | target %in% node_names)
    
  }  
    
  # Filter the network data frame to only include rows where 'target' is in target_genes, if provided.
  if (!is.null(target_genes)) {
    meta_DF <- meta_DF %>% filter(target %in% target_genes)
  }
  
  # Further filter the network data frame to only include rows where 'tf' is in TF_names, if provided.
  if (!is.null(TF_names)) {
    meta_DF <- meta_DF %>% filter(tf %in% TF_names)
  }
  
  # Filter based on the absolute value of the 'estimate' column, if score_cut is provided.
  if (!is.null(score_cut)) {
    meta_DF <- meta_DF %>% filter(abs(estimate) >= score_cut)
    score_cut <- paste0(' >= ', score_cut)  # Record the applied score cut for annotation.
  } else {
    score_cut <- 'unset'  # If no score_cut is provided, note that it's unset.
  }
   
  # Add check for meta_df having at least one row after the optional filtering steps above:
  
  if (nrow(meta_DF) <= 0) {
    stop("Error: Sorry but meta_df should have at least one row.")
  }
    
  # Modify the 'tf' and 'target' columns to append ".target-TF" to genes that are both TFs and targets.
  meta_DF$tf <- ifelse(meta_DF$tf %in% target_TFs, paste0(meta_DF$tf, ".target-TF"), meta_DF$tf)
  meta_DF$target <- ifelse(meta_DF$target %in% target_TFs, paste0(meta_DF$target, ".target-TF"), meta_DF$target)
  
  # Rescale the edge weights (based on the 'estimate' column) to a range suitable for visualization (0.5 to 3).
  meta_DF$weight <- scales::rescale(abs(meta_DF$estimate), to = c(0.5, 3))

  # Create a directed network object from the 'tf' and 'target' columns, interpreted as an edge list.
  net <- network(meta_DF[, c("tf", "target")], directed = TRUE, matrix.type = "edgelist")
    
  net %>% print()   
 
  # Assign the computed weights to the edges in the network.
  set.edge.value(x = net, attrname = "weight", value = meta_DF$weight)
  
  # Note: Below, I wrote the network package explicitly, because igraph package can mask these functions,
  # when used in the same R session.
    
  # Extract the names of the vertices (nodes) in the network.
  vertex_names <- network::get.vertex.attribute(x = net, attrname = 'vertex.names')

  # Determine the identity of each node as 'target', 'TF', or 'target-TF'.
  identity <- ifelse(vertex_names %in% paste0(target_TFs, ".target-TF"), "target-TF", 
                     ifelse(vertex_names %in% meta_DF$target, "target", "TF"))
  
  # Set the class attribute of each vertex to the determined identity ('target', 'TF', or 'target-TF').
  network::set.vertex.attribute(x = net, attrname = 'class', value = identity)
    
  
  # ALTERNATIVE:

# Network object contains a list called 'val'. This list consists of
# list elements. Each of this list elements stores the vertex names 
# as their secon variable.

# For instance the fist vertex name can be accessed via:

# net$val[[1]][[2]]  

# below code uses sapply to access the second element of each nested lists and
# generates a vector consisting of the second elements from each list.

# net %v% "class" <-  
#   ifelse(sapply(net$val,"[[",2) %in% meta_DF$target, "target", "tf")

#

# After we set the classes of vertices (nodes) as target and tf,
# we can now delete the '.' that we previously added to distinguish
# target genes from TFs. Because, there can be target genes that
# are also Transcription factors, we applied this above.
    
    
  # Clean up the vertex names by removing the ".target-TF" suffix.
  network::set.vertex.attribute(net, "vertex.names", 
                                gsub(pattern = ".target-TF", 
                                     replacement = "", 
                                     x = network::get.vertex.attribute(net, "vertex.names"), fixed = TRUE))

  # Assign edge colors based on whether the correlation ('estimate') is positive (red) or negative (blue).
  network::set.edge.attribute(net, "color", ifelse(meta_DF$estimate > 0, "firebrick", "steelblue"))

    
  # Define colors for the nodes based on their class.
  col <- c("TF" = "gray60", "target" = "sandybrown", "target-TF" = "#AC94F4")

    
# set.seed(seed_value)

# p1 <- 
#   ggnet2(net = net,
#        color.legend = "Class",
#        size.legend = "Degree",
#        label = TRUE,
#        size="class",
#        edge.size = "weight",
#        edge.color = "color",
#        edge.alpha = 0.75,
#        alpha = 0.7,
#        color="class",
#        size.palette = c("target" = 3, "tf" = 1),
#        label.size = 3,
#        palette=col,
#        legend.position = 'right') + 
#   coord_equal() +
#   theme(legend.text = element_text(size = 20))

# p1
    
# I could not add legend for the edge weights therefore I added the following attribute
# to the vertex to use only in the legend. It does not change anything in the structure
# of the network.

for_legend <- 
  ifelse(vertex_names %in% meta_DF$target, "activator", "repressor")

network::set.vertex.attribute(x = net, attrname = 'for_legend', value = for_legend)

# Set a random seed for reproducibility in the layout of the network visualization.    
    
set.seed(seed_value)

if (!is.null(arrow)) {
    
    p1 <- 
  ggnet2(net = net,
       color.legend = "Class",
       size.legend = "State",
       label = TRUE,
       size="for_legend",
       edge.size = "weight",
       edge.color = "color",
       edge.alpha = 0.75,
       alpha = 0.7,
       color="class",
       size.palette = c("activator" = 1, "repressor" = 1),
       label.size = label_size,
       palette=col,
       legend.position = 'right',
       arrow.size = 12, # newly added  for edge arrows.
       arrow.gap = 0.025) + # newly added  for edge arrows.
  coord_equal() +
  theme(legend.text = element_text(size = 20)) +
  ggtitle(label = 'TF modules-meta df based network graph')
    
    
  } else {

     p1 <- 
  ggnet2(net = net,
       color.legend = "Class",
       size.legend = "State",
       label = TRUE,
       size="for_legend", # this is just a workaround to add the edge colours to the legend. 
                          # It will be overriden with the two edge  colours after generation of the plot.
       edge.size = "weight",
       edge.color = "color",
       edge.alpha = 0.75,
       alpha = 0.7,
       color="class",
       size.palette = c("activator" = 1, "repressor" = 1),
       label.size = label_size,
       palette=col,
       legend.position = 'right') + 
  coord_equal() +
  theme(legend.text = element_text(size = 20)) +
  ggtitle(label = 'TF modules-meta df based network graph')

}
    

p1 <- p1 + guides(size = guide_legend(override.aes = list(color = c("firebrick", "steelblue")))) # Override the relevant legend for the edges.
    
    # To make the node sizes bigger when visualizing the subnetworks:

    if(!is.null(override_node_size_value)) {
        
     p1 <- p1 + scale_size_manual(values = c(override_node_size_value, override_node_size_value))  # It will make the edge link legends bigger !!!.
  
  }
    
return(p1)   
    
}  

#########################################

library(igraph)

vis_network_IGRAPH <- function(meta_DF, score_cut = NULL, 
                               target_genes = NULL,
                               TF_names = NULL, 
                               node_names = NULL,
                               seed_value = 123,
                               label_size = 15,
                               arrow = NULL,
                               edge_arrow_size = 0.3,
                               edge_curved = 0.2,
                               node_size = 10,
                               add_legend = TRUE,  # New parameter to control node size
                               degree = NULL) {
    

  # To filter the vertices by user-defined degree threshold:  
    
  if (!is.null(degree)) {
      
# Below, I convert the DF to a vector by unlist,
# then, tabulate the frequencies and select the genes that exist more or equal to the input number of degree.  
      
      keep_nodes <- meta_DF %>% unlist %>% table %>% { .[. >= degree] } %>% names
      
      meta_DF <- meta_DF %>% filter(tf %in% keep_nodes | target %in% keep_nodes)
      
    }  
    
  # Identify genes that appear both as transcription factors (TFs) and as targets.
  
  target_TFs <- intersect(meta_DF$tf, meta_DF$target)  
  
  # Filter the network data frame to include rows where 'tf' or 'target' matches node_names, if provided.
  if(!is.null(node_names)) {
    meta_DF <- meta_DF %>% filter(tf %in% node_names | target %in% node_names)
  }  
    
  # Further filter the network data frame to include rows where 'target' matches target_genes, if provided.
  if (!is.null(target_genes)) {
    meta_DF <- meta_DF %>% filter(target %in% target_genes)
  }
  
  # Further filter the network data frame to include rows where 'tf' matches TF_names, if provided.
  if (!is.null(TF_names)) {
    meta_DF <- meta_DF %>% filter(tf %in% TF_names)
  }
  
  # Filter based on the absolute value of the 'estimate' column, if score_cut is provided.
  if (!is.null(score_cut)) {
    meta_DF <- meta_DF %>% filter(abs(estimate) >= score_cut)
    score_cut <- paste0(' >= ', score_cut)  # Record the applied score cut for annotation.
  } else {
    score_cut <- 'unset'  # If no score_cut is provided, note that it's unset.
  }

    
    # Add check for meta_df having at least one row after the optional filtering steps above:
  
  if (nrow(meta_DF) <= 0) {
    stop("Error: meta_DF should have at least one row.")
  }
    
  # Modify the 'tf' and 'target' columns to append ".target-TF" to genes that are both TFs and targets.
  meta_DF$tf <- ifelse(meta_DF$tf %in% target_TFs, paste0(meta_DF$tf, ".target-TF"), meta_DF$tf)
  meta_DF$target <- ifelse(meta_DF$target %in% target_TFs, paste0(meta_DF$target, ".target-TF"), meta_DF$target)
  
  # Rescale the edge weights (based on the 'estimate' column) to a range suitable for visualization (0.5 to 3).
  meta_DF$weight <- scales::rescale(abs(meta_DF$estimate), to = c(0.5, 3))

  # Create a directed igraph object from the 'tf' and 'target' columns, interpreted as an edge list.
  net <- igraph::graph_from_data_frame(d = meta_DF, directed = TRUE)
  
  # Set the 'class' attribute for vertices based on their identity ('target', 'TF', or 'target-TF').
  V(net)$class <- ifelse(V(net)$name %in% paste0(target_TFs, ".target-TF"), "target-TF", 
                         ifelse(V(net)$name %in% meta_DF$target, "target", "TF"))
  
  # Clean up the vertex names by removing the ".target-TF" suffix.
  V(net)$name <- gsub(pattern = ".target-TF", replacement = "", x = V(net)$name, fixed = TRUE)
  
  # Assign edge colors based on whether the correlation ('estimate') is positive (red) or negative (blue).
  E(net)$color <- ifelse(meta_DF$estimate > 0, "firebrick", "steelblue")
  
  # Assign node colors based on their class.
  col <- c("TF" = "gray60", "target" = "sandybrown", "target-TF" = "#AC94F4")  # the last one is a version of light purple. Purple is too dark in the igraph implementation.
  V(net)$color <- col[V(net)$class]
  
  # Print the graph summary.
  print(net)
    
  # Set the random seed for reproducibility of the layout.
  set.seed(seed_value)
  
  # Plot the network using igraph's plotting capabilities.
  plot(net, 
       vertex.label = V(net)$name,               # Label vertices with their names.
       vertex.label.cex = label_size / 10,       # Control label size.
       vertex.size = node_size,                  # Use the node_size parameter for vertex size.
       vertex.color = V(net)$color,              # Assign colors based on class.
       edge.width = E(net)$weight + 2,           # Set edge width based on rescaled weight. I ADDED PLUS TWO, BECAUSE IGRAPH MAKES LINES THINNER IN APPEARANCE THAN THE OTHER METHODS!!" 
       edge.color = E(net)$color,                # Set edge color.
       edge.curved = edge_curved,                        # Slightly curve edges.
       edge.arrow.size = ifelse(!is.null(arrow), edge_arrow_size, 0),  # Smaller arrow size; no arrows if `arrow = NULL`.
       main = 'TF modules-meta df based network graph')  # Title for the plot.
  
#   # Add a legend for node classes.
#   legend("topright", legend = names(col), col = col, pch = 21, pt.bg = col, pt.cex = 3, 
#          title = "Node Classes", cex = 1.75)
  
#   # Add a legend for edge types (activator/repressor).
#   legend("bottomright", legend = c("Activator", "Repressor"), col = c("firebrick", "steelblue"), 
#          lwd = 2, title = "Edge Types", cex = 1.75)

#   # Add a legend for node classes without rectangles, with increased text size and reduced spacing.
#   legend("topright", legend = names(col), col = col, pch = 21, pt.bg = col, pt.cex = 1.5, 
#          title = "Node Classes", cex = 1.75, 
#          bty = "n", x.intersp = 0.5, y.intersp = 0.5)
  
#   # Add a legend for edge types (activator/repressor) without rectangles, with increased text size and reduced spacing.
#   legend("bottomright", legend = c("Activator", "Repressor"), col = c("firebrick", "steelblue"), 
#          lwd = 2, title = "Edge Types", pch = 21, cex = 1.75, 
#          bty = "n", x.intersp = 0.5, y.intersp = 0.5)
    
  if(add_legend) {
      
          # Add a legend for node classes without rectangles, with increased text size and reduced spacing.
      legend("topright", legend = names(col), col = col, pch = 21, pt.bg = col, pt.cex = 3, 
             cex = 1.75, 
             bty = "n", x.intersp = 0.5, y.intersp = 0.5, 
             xpd = TRUE, inset = c(-0.3, 0))  # Move the legend further to the right

      # Add a legend for edge types (activator/repressor) without rectangles, with increased text size and reduced spacing.
      legend("bottomright", legend = c("Activator", "Repressor"), col = c("firebrick", "steelblue"), 
             lwd = 3, pch = NA, cex = 1.75, 
             bty = "n", x.intersp = 0.5, y.intersp = 0.5, 
             xpd = TRUE, inset = c(-0.3, 0), seg.len = 0.5)  # Move the legend further to the right
    }
 
}  


###################################

library(ggraph)

vis_subnetwork_w_GGRAPH <- function(meta_DF, score_cut = NULL, 
                             target_genes = NULL,
                             TF_names = NULL, 
                             node_names = NULL,
                             up_or_down = NULL,
                             seed_value = 123,
                             label_size = 6,
                             node_point_size = 25,
                             arc = NULL,
                             arrow = NULL,
                             layout = 'fr',
                             degree = NULL) {
    

  # To filter the vertices by user-defined degree threshold:  
    
  if (!is.null(degree)) {
      
# Below, I convert the DF to a vector by unlist,
# then, tabulate the frequencies and select the genes that exist more or equal to the input number of degree.  
      
      keep_nodes <- meta_DF %>% unlist %>% table %>% { .[. >= degree] } %>% names
      
      meta_DF <- meta_DF %>% filter(tf %in% keep_nodes | target %in% keep_nodes)
      
    }  
    
  # Identify genes that appear both as transcription factors (TFs) and as targets.
  
  target_TFs <- intersect(meta_DF$tf, meta_DF$target)    
    
  manual_color <- c("firebrick", "steelblue")  # it is used in the ggraph function to adjust the colour of the edges.
    
  # Filter the network data frame to only include rows where 'estimate' matches up_or_down, if provided.
  if (!is.null(up_or_down)) {
    if (up_or_down == 'up') {
      meta_DF <- meta_DF %>% filter(estimate > 0)
      
      manual_color <- c("firebrick")  
      
    } else if (up_or_down == 'down') {
      meta_DF <- meta_DF %>% filter(estimate < 0)
      
      manual_color <- c("steelblue")  
        
    }
  }
  
  # Filter the network data frame to only include rows where 'target' is in target_genes, if provided.
  if (!is.null(target_genes)) {
    meta_DF <- meta_DF %>% filter(target %in% target_genes)
  }
  
  # Further filter the network data frame to only include rows where 'tf' is in TF_names, if provided.
  if (!is.null(TF_names)) {
    meta_DF <- meta_DF %>% filter(tf %in% TF_names)
  }
  
  if (!is.null(node_names)) {
    meta_DF <- meta_DF %>% filter(tf %in% node_names | target %in% node_names)
  }
    
  # Filter based on the absolute value of the 'estimate' column, if score_cut is provided.
  if (!is.null(score_cut)) {
    meta_DF <- meta_DF %>% filter(abs(estimate) >= score_cut)
    score_cut <- paste0(' >= ', score_cut)  # Record the applied score cut for annotation.
  } else {
    score_cut <- 'unset'  # If no score_cut is provided, note that it's unset.
  }

    
  # Add check for meta_df having at least one row after the optional filtering steps above:
  
  if (nrow(meta_DF) <= 0) {
    stop("Error: meta_DF should have at least one row.")
  }


  # Modify the 'tf' and 'target' columns to append ".target-TF" to genes that are both TFs and targets.
  meta_DF$tf <- ifelse(meta_DF$tf %in% target_TFs, paste0(meta_DF$tf, ".target-TF"), meta_DF$tf)
  meta_DF$target <- ifelse(meta_DF$target %in% target_TFs, paste0(meta_DF$target, ".target-TF"), meta_DF$target)
  
  # Rescale the edge weights (based on the 'estimate' column) to a range suitable for visualization (0.5 to 3).
  meta_DF$weight <- scales::rescale(abs(meta_DF$estimate), to = c(0.5, 3))
  
  # Create an igraph object from the 'tf' and 'target' columns, interpreted as an edge list.
  graph <- graph_from_data_frame(d = meta_DF, directed = TRUE)
  
  # Clean up the vertex names by removing the ".target-TF" suffix.
  V(graph)$name <- gsub(pattern = ".target-TF", replacement = "", x = V(graph)$name, fixed = TRUE)
  
  # Determine the identity of each node as 'target', 'TF', or 'target-TF'.
  V(graph)$class <- ifelse(V(graph)$name %in% target_TFs, "target-TF", 
                           ifelse(V(graph)$name %in% meta_DF$target, "target", "TF"))
  
  # Assign edge colors based directly on the filtered meta_DF.
  E(graph)$color <- ifelse(E(graph)$estimate > 0, "firebrick", "steelblue")
    
  # Define colors for the nodes based on their class.
  col <- c("TF" = "gray60", "target" = "sandybrown", "target-TF" = "#AC94F4")
  
  # Set a random seed for reproducibility in the layout of the network visualization.
  set.seed(seed_value)
  
  # print grapph object:
    
    print(graph)
    
 
 if(layout %in% c('stress', 'fr', 'graphopt', 'kk')) {   
    
  if(!is.null(arc)) {
      
        # Create the ggraph plot
  p1 <- ggraph(graph, layout = layout) + 
    geom_edge_arc(aes(width = weight, color = color),
                  arrow = if (!is.null(arrow)) arrow(length = unit(5, 'mm')) else NULL,
                  start_cap = circle(12, 'mm'),
                  end_cap = circle(12, 'mm')) +
    geom_node_point(aes(color = class), size = node_point_size) +
    geom_node_text(aes(label = name), size = label_size, repel = FALSE) +
    scale_edge_width(range = c(0.5, 3)) +
    scale_edge_color_manual(values = manual_color) +
    scale_color_manual(values = col) +
    theme_void() +
    theme(legend.text = element_text(size = 20)) +
    ggtitle(label = 'TF modules-meta df based network graph')
      
  } else {    
    
  # Create the ggraph plot
  p1 <- ggraph(graph, layout = layout) + 
    geom_edge_fan(aes(width = weight, color = color),
                  arrow = if (!is.null(arrow)) arrow(length = unit(5, 'mm')) else NULL,
                  end_cap = circle(12, 'mm')) +
    geom_node_point(aes(color = class), size = node_point_size) +
    geom_node_text(aes(label = name), size = label_size, repel = FALSE) +
    scale_edge_width(range = c(0.5, 3)) +
    scale_edge_color_manual(values = manual_color) +
    scale_color_manual(values = col) +
    theme_void() +
    theme(legend.text = element_text(size = 20)) +
    ggtitle(label = 'TF modules-meta df based network graph')
  
  }
}
                  
 if(layout %in% c('linear')) {   
    
  if(!is.null(arc)) {
      
        # Create the ggraph plot
  p1 <- ggraph(graph, layout = layout, circular = TRUE) + 
    geom_edge_arc(aes(width = weight, color = color),
                  arrow = if (!is.null(arrow)) arrow(length = unit(5, 'mm')) else NULL,
                  start_cap = circle(12, 'mm'),
                  end_cap = circle(12, 'mm')) +
    geom_node_point(aes(color = class), size = node_point_size) +
    geom_node_text(aes(label = name), size = label_size, repel = FALSE) +
    scale_edge_width(range = c(0.5, 3)) +
    scale_edge_color_manual(values = manual_color) +
    scale_color_manual(values = col) +
    theme_void() +
    theme(legend.text = element_text(size = 20)) +
    ggtitle(label = 'TF modules-meta df based network graph')
      
  } else {    
    
  # Create the ggraph plot
  p1 <- ggraph(graph, layout = layout, circular = TRUE) + 
    geom_edge_fan(aes(width = weight, color = color),
                  arrow = if (!is.null(arrow)) arrow(length = unit(5, 'mm')) else NULL,
                  end_cap = circle(12, 'mm')) +
    geom_node_point(aes(color = class), size = node_point_size) +
    geom_node_text(aes(label = name), size = label_size, repel = FALSE) +
    scale_edge_width(range = c(0.5, 3)) +
    scale_edge_color_manual(values = manual_color) +
    scale_color_manual(values = col) +
    theme_void() +
    theme(legend.text = element_text(size = 20)) +
    ggtitle(label = 'TF modules-meta df based network graph')
  
  }
}                       
                  
  return(p1)
}
                  
##########################################################

vis_network_w_thresholded_degrees_w_GGNET2 <- function(meta_DF, 
                                              score_cut = NULL, 
                                              target_genes = NULL,
                                              TF_names = NULL,
                                              node_names = NULL,
                                              seed_value = 123,
                                              label_size = 6,
                                              override_node_size_value = 28,
                                              arrow = NULL,
                                              degree = NULL) {
    

  # To filter the vertices by user-defined degree threshold:  
    
  if (!is.null(degree)) {
      
# Below, I convert the DF to a vector by unlist,
# then, tabulate the frequencies and select the genes that exist more or equal to the input number of degree.  
      
      keep_nodes <- meta_DF %>% unlist %>% table %>% { .[. >= degree] } %>% names
      
      meta_DF <- meta_DF %>% filter(tf %in% keep_nodes | target %in% keep_nodes)
      
    }
    
    
     p1 <- vis_subnetwork_w_GGNET2(meta_DF = meta_DF, 
                                  score_cut = score_cut, 
                                  target_genes = target_genes, 
                                  TF_names = TF_names, 
                                  node_names = node_names,
                                  seed_value = seed_value, 
                                  label_size = label_size, 
                                  override_node_size_value = override_node_size_value,
                                  arrow = arrow)
      
      
    return(p1)
}  

                  
############################################################

# IMPORTANT: I ADDED DEGREE PARAMETER TO ALL FUNCTIONS ABOVE. SO BELOW FUNCTIONS ARE USELESS NOW.
# I STILL KEEP THEM BELOW FOR FUTURE REFERENCE.

############################################################
                  
# Following function contains an alternative code that uses sna package's degree function to calculate
# degrees for each node.
                  
# If the degree is provided, the following function plots the connections that contain those selected nodes,
# where either the target or the tf gene is among those degree-based selected nodes.

# vis_network_w_degree_based_subset_w_GGNET2 <- function(meta_DF, 
#                                               score_cut = NULL, 
#                                               target_genes = NULL,
#                                               TF_names = NULL, 
#                                               seed_value = 123,
#                                               label_size = 3,
#                                               arrow = NULL,
#                                               degree = NULL) {
    
    
# # Create a directed network object from the 'tf' and 'target' columns, interpreted as an edge list.
  
#     net <- network(meta_DF[, c("tf", "target")], directed = TRUE, matrix.type = "edgelist")
  
#   # To filter the vertices by user-defined degree threshold:  
    
#   if (!is.null(degree)) {
      
#      # Find the indices of the vertices with at least the specified number of degree:
      
#      keep_degree <- which(sna::degree(dat = net, cmode = "freeman") >= degree)  # cmode is default. 
      
#      # Filter the network: 
      
#      net <- get.inducedSubgraph(x = net, v = keep_degree)
      
#      vertex_names <-network::get.vertex.attribute(x = net, attrname = 'vertex.names')
     
#      meta_DF <- meta_DF %>% filter(tf %in% vertex_names | target %in% vertex_names)
      
      
      
#      p1 <- vis_subnetwork_w_GGNET2(meta_DF = meta_DF, 
#                                   score_cut = score_cut, 
#                                   target_genes = target_genes, 
#                                   TF_names = TF_names, 
#                                   seed_value = seed_value, 
#                                   label_size = label_size, 
#                                   arrow = arrow)
      
      
#   } else {
      
      
#       p1 <- vis_subnetwork_w_GGNET2(meta_DF = meta_DF, 
#                                   score_cut = score_cut, 
#                                   target_genes = target_genes, 
#                                   TF_names = TF_names, 
#                                   seed_value = seed_value, 
#                                   label_size = label_size, 
#                                   arrow = arrow)
      
      
#       }
    
#     return(p1)
# }  
                 
####################################################

# NOTE: I WILL CHANGE THIS FUNCTION BY INCORPORATING THE LAST IMPROVED VERSION OF GGNET2 FUNCTION INTO THE BELOW FUNCTION !!!!!!!!!!!
                  
# If degree parameter is provided, which specifies the minimum number of connections that a node must have,  
# this function will generate a network plot that includes only these nodes. 

# SO IT ONLY PLOTS THOSE SUBSETTED NODES' CONNECTIONS BETWEEN EACH OTHER, NOT THEIR ALL CONNECTIONS.

vis_network_w_degree_thresholded_vertices_NOT_WORKING_YET <- function(meta_DF, score_cut = NULL, 
                          target_genes = NULL,
                          TF_names = NULL, seed_value = 123,
                          label_size = 3,
                          arrow = NULL,
                          degree = NULL) {
  
  # Filter the network data frame to only include rows where 'target' is in target_genes, if provided.
  if (!is.null(target_genes)) {
    meta_DF <- meta_DF %>% filter(target %in% target_genes)
  }
  
  # Further filter the network data frame to only include rows where 'tf' is in TF_names, if provided.
  if (!is.null(TF_names)) {
    meta_DF <- meta_DF %>% filter(tf %in% TF_names)
  }
  
  # Filter based on the absolute value of the 'estimate' column, if score_cut is provided.
  if (!is.null(score_cut)) {
    meta_DF <- meta_DF %>% filter(abs(estimate) >= score_cut)
    score_cut <- paste0(' >= ', score_cut)  # Record the applied score cut for annotation.
  } else {
    score_cut <- 'unset'  # If no score_cut is provided, note that it's unset.
  }

    # Add check for meta_df having at least one row after the optional filtering steps above:
    
    if (nrow(meta_DF) <= 0) {
    stop("Error: meta_DF should have at least one row.")
  }  
    
  # Identify genes that appear both as transcription factors (TFs) and as targets.
  target_TFs <- intersect(meta_DF$tf, meta_DF$target)

  # Modify the 'tf' and 'target' columns to append ".target-TF" to genes that are both TFs and targets.
  meta_DF$tf <- ifelse(meta_DF$tf %in% target_TFs, paste0(meta_DF$tf, ".target-TF"), meta_DF$tf)
  meta_DF$target <- ifelse(meta_DF$target %in% target_TFs, paste0(meta_DF$target, ".target-TF"), meta_DF$target)
  
  # Rescale the edge weights (based on the 'estimate' column) to a range suitable for visualization (0.5 to 3).
  meta_DF$weight <- scales::rescale(abs(meta_DF$estimate), to = c(0.5, 3))

  # Create a directed network object from the 'tf' and 'target' columns, interpreted as an edge list.
  net <- network(meta_DF[, c("tf", "target")], directed = TRUE, matrix.type = "edgelist")

    
  # To filter the vertices by user-defined degree threshold:  
    
  if (!is.null(degree)) {
      
     # Find the indices of the vertices with at least the specified number of degree:
      
     keep_degree <- which(sna::degree(dat = net, cmode = "freeman") >= degree)  # cmode is default. 
      
     # Filter the network: 
      
     net <- get.inducedSubgraph(x = net, v = keep_degree)
      
  }
    
    

    
  net %>% print()   
 
  # Assign the computed weights to the edges in the network.
  set.edge.value(x = net, attrname = "weight", value = meta_DF$weight)
  
  # Extract the names of the vertices (nodes) in the network.
  vertex_names <- network::get.vertex.attribute(x = net, attrname = 'vertex.names')

  # Determine the identity of each node as 'target', 'TF', or 'target-TF'.
  identity <- ifelse(vertex_names %in% paste0(target_TFs, ".target-TF"), "target-TF", 
                     ifelse(vertex_names %in% meta_DF$target, "target", "TF"))
  
  # Set the class attribute of each vertex to the determined identity ('target', 'TF', or 'target-TF').
  network::set.vertex.attribute(x = net, attrname = 'class', value = identity)
    
  
  # ALTERNATIVE:

# Network object contains a list called 'val'. This list consists of
# list elements. Each of this list elements stores the vertex names 
# as their secon variable.

# For instance the fist vertex name can be accessed via:

# net$val[[1]][[2]]  

# below code uses sapply to access the second element of each nested lists and
# generates a vector consisting of the second elements from each list.

# net %v% "class" <-  
#   ifelse(sapply(net$val,"[[",2) %in% meta_DF$target, "target", "tf")

#

# After we set the classes of vertices (nodes) as target and tf,
# we can now delete the '.' that we previously added to distinguish
# target genes from TFs. Because, there can be target genes that
# are also Transcription factors, we applied this above.
    
    
  # Clean up the vertex names by removing the ".target-TF" suffix.
  network::set.vertex.attribute(net, "vertex.names", gsub(pattern = ".target-TF", replacement = "", x = network::get.vertex.attribute(net, "vertex.names"), fixed = TRUE))

  # Assign edge colors based on whether the correlation ('estimate') is positive (red) or negative (blue).
#  set.edge.attribute(net, "color", ifelse(meta_DF$estimate > 0, "firebrick", "steelblue"))

        # Assign edge colors based on whether the correlation ('estimate') is positive (red) or negative (blue).
set.edge.value(x = net, attrname = "color", value = ifelse(meta_DF$estimate > 0, "firebrick", "steelblue"))
    
    
  # Define colors for the nodes based on their class.
  col <- c("TF" = "gray60", "target" = "sandybrown", "target-TF" = "#AC94F4")

    
# set.seed(seed_value)

# p1 <- 
#   ggnet2(net = net,
#        color.legend = "Class",
#        size.legend = "Degree",
#        label = TRUE,
#        size="class",
#        edge.size = "weight",
#        edge.color = "color",
#        edge.alpha = 0.75,
#        alpha = 0.7,
#        color="class",
#        size.palette = c("target" = 3, "tf" = 1),
#        label.size = 3,
#        palette=col,
#        legend.position = 'right') + 
#   coord_equal() +
#   theme(legend.text = element_text(size = 20))

# p1
    
# I could not add legend for the edge weights therefore I added the following attribute
# to the vertex to use only in the legend. It does not change anything in the structure
# of the network.

    
# I could not add legend for the edge weights therefore I added the following attribute
# to the vertex to use only in the legend. It does not change anything in the structure
# of the network.
    
for_legend <- 
  ifelse(vertex_names %in% meta_DF$target, "activator", "repressor")

network::set.vertex.attribute(x = net, attrname = 'for_legend', value = for_legend)

# Set a random seed for reproducibility in the layout of the network visualization.    
    
set.seed(seed_value)

if (!is.null(arrow)) {
    
    p1 <- 
  ggnet2(net = net,
       color.legend = "Class",
       size.legend = "State",
       label = TRUE,
       size="for_legend",
       edge.size = "weight",
       edge.color = "color",
       edge.alpha = 0.75,
       alpha = 0.7,
       color="class",
       size.palette = c("activator" = 1, "repressor" = 1),
       label.size = label_size,
       palette=col,
       legend.position = 'right',
       arrow.size = 12, # newly added  for edge arrows.
       arrow.gap = 0.025) + # newly added  for edge arrows.
  coord_equal() +
  theme(legend.text = element_text(size = 20)) +
  ggtitle(label = 'TF modules-meta df based network graph')
    
    
  } else {

     p1 <- 
  ggnet2(net = net,
       color.legend = "Class",
       size.legend = "State",
       label = TRUE,
       size="for_legend", # this is just a workaround to add the edge colours to the legend. 
                          # It will be overriden with the two edge  colours after generation of the plot.
       edge.size = "weight",
       edge.color = "color",
       edge.alpha = 0.75,
       alpha = 0.7,
       color="class",
       size.palette = c("activator" = 1, "repressor" = 1),
       label.size = label_size,
       palette=col,
       legend.position = 'right') + 
  coord_equal() +
  theme(legend.text = element_text(size = 20)) +
  ggtitle(label = 'TF modules-meta df based network graph')

}
    

p1 <- p1 + guides(size = guide_legend(override.aes = list(color = c("firebrick", "steelblue")))) # Override the relevant legend for the edges.

p1    
    
}  
                  
###############################################################

message("make sure you installed the GGally, network, igraph and ggraph packages :)")

message("First option: One can directly use only the node_names parameter to specify 
          the gene names to visualize in the graph")

message("Second option: Use TF-names and/or target genes parameters 
         to visualize their connections")