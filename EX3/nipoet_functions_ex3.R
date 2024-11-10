library(plotly)

nipoet_dataset_03 <- function(n=200) {
    cluster1 <- mvrnorm(n, mu = c(2, 2), Sigma = matrix(c(1, 0.5, 0.5, 1), 2))
    cluster2 <- mvrnorm(n, mu = c(7, 7), Sigma = matrix(c(1, -0.5, -0.5, 1), 2))
    cluster3 <- mvrnorm(n, mu = c(12, 2), Sigma = matrix(c(1.5, 0.3, 0.3, 1.5), 2))
    ds <- rbind(cluster1, cluster2, cluster3)
    return (ds)
}

nipoet_generate_half_moon <- function(n, radius = 2.5, noise = 0.2) {
  # Generate x-coordinates
  theta <- runif(n, 0, pi)
  x <- radius * cos(theta) + rnorm(n, 0, noise)
  
  # Generate y-coordinates
  y <- radius * sin(theta) + rnorm(n, 0, noise)
  
  # Shift the data to create a half-moon shape
  y <- y - radius
  
  return(cbind(x, y))
}

#--- Algorithms ---#

create_grid <- function(data, x_num, y_num) {
    x_min <- min(data$X1)
    x_max <- max(data$X1)
    y_min <- min(data$X2)
    y_max <- max(data$X2)
    
    x_grid <- seq(x_min, x_max, length.out = x_num + 1)
    y_grid <- seq(y_min, y_max, length.out = y_num + 1)
    
    regions <- list()
    
    for (i in 1:x_num) {
        for (j in 1:y_num) {
        temp <- list(
            x_min = x_grid[i],
            x_max = x_grid[i + 1],
            y_min = y_grid[j],
            y_max = y_grid[j + 1],
            points = data[data$X1 >= x_grid[i] & data$X1 <= x_grid[i + 1] &
                            data$X2 >= y_grid[j] & data$X2 <= y_grid[j + 1], ],
            cluster = FALSE,
            merged_index = NULL
        )
        regions[[length(regions) + 1]] <- temp
        }
    }
    
    return(regions)
}

identify_clusters <- function(grid, threshold) {
    for (i in seq_along(grid)) {
        grid[[i]]$cluster <- (nrow(grid[[i]]$points) > threshold)
    }
    return(grid)
}

are_neighbours <- function(region1, region2) {
    x_adjacent <- (region1$x_max == region2$x_min) || (region1$x_min == region2$x_max)
    y_adjacent <- (region1$y_max == region2$y_min) || (region1$y_min == region2$y_max)
    return(x_adjacent || y_adjacent)
}

merge_clusters <- function(grid) {
    # Initialize a list to hold clusters
    clusters <- list()
    cluster_id <- 1  # Cluster index
    
    # Iterate over each grid cell and look for unvisited clusters
    for (i in seq_along(grid)) {
        # If the region is not part of a cluster or already merged, skip it
        if (!grid[[i]]$cluster || !is.null(grid[[i]]$merged_index)) next
        
        # Initialize a new cluster and perform a DFS to find all connected regions
        current_cluster <- list()
        stack <- list(i)
        
        while (length(stack) > 0) {
            current_index <- stack[[length(stack)]]
            stack <- stack[-length(stack)]
            
            # Skip if already visited
            if (!is.null(grid[[current_index]]$merged_index)) next
            
            # Mark as visited by setting merged_index to cluster_id
            grid[[current_index]]$merged_index <- cluster_id
            current_cluster <- append(current_cluster, list(grid[[current_index]]))
            
            # Look for neighboring clusters
            for (j in seq_along(grid)) {
                if (is.null(grid[[j]]$merged_index) && grid[[j]]$cluster && are_neighbours(grid[[current_index]], grid[[j]])) {
                    # Add unvisited neighbor to stack
                    stack <- append(stack, j)
                }
            }
        }
        
        # Add the current cluster to the list of clusters
        clusters[[cluster_id]] <- current_cluster
        cluster_id <- cluster_id + 1
    }
    
    return(clusters)
}

visualize_grid <- function(data, regions, merged_regions) {
    plot <- ggplot(data, aes(x = X1, y = X2)) +
        geom_point(size = 2, color = "blue") +
        theme_minimal()

    for (region in regions) {
        plot <- plot + 
        annotate("rect", xmin = region$x_min, xmax = region$x_max, 
                ymin = region$y_min, ymax = region$y_max, 
                alpha = 0.1, fill = ifelse(region$cluster, "gray", "gray"), color = "black")
    }
    
    # THIS DOESN'T WORK PROPERLY????
    colors <- rainbow(length(merged_regions))
    
    for (i in seq_along(merged_regions)) {
        for (region in merged_regions[[i]]) {
            plot <- plot + 
                annotate("rect", xmin = region$x_min, xmax = region$x_max, 
                        ymin = region$y_min, ymax = region$y_max, 
                        alpha = 0.3, fill = colors[i], color = "black")
        }
    }
    
    print(plot)
}

#---- DB-Scan ----#

nipoet_dbscan <- function(data, eps, min_pts) {
  # Initialize variables
  n_points <- nrow(data)
  visited <- rep(FALSE, n_points)
  noise <- rep(FALSE, n_points)
  clusters <- rep(0, n_points)  # 0 means unassigned
  cluster_id <- 0
  
  # function to find neighbors within eps radius
  find_neighbors <- function(point_idx) {
    distances <- sqrt(rowSums((data - matrix(data[point_idx,], 
                                           nrow=n_points, 
                                           ncol=ncol(data), 
                                           byrow=TRUE))^2))
    return(which(distances <= eps))
  }
  
  # main DBSCAN algorithm
  for (point_idx in 1:n_points) {
    if (visited[point_idx]) next
    
    visited[point_idx] <- TRUE
    neighbors <- find_neighbors(point_idx)
    
    # Check if point is noise
    if (length(neighbors) < min_pts) {
      noise[point_idx] <- TRUE
      next
    }
    
    # Start a new cluster
    cluster_id <- cluster_id + 1
    clusters[point_idx] <- cluster_id
    
    # Process neighbors
    seed_set <- neighbors[neighbors != point_idx]
    
    while (length(seed_set) > 0) {
      current_point <- seed_set[1]
      
      if (!visited[current_point]) {
        visited[current_point] <- TRUE
        current_neighbors <- find_neighbors(current_point)
        
        if (length(current_neighbors) >= min_pts) {
          seed_set <- unique(c(seed_set, current_neighbors))
        }
      }
      
      if (clusters[current_point] == 0) {
        clusters[current_point] <- cluster_id
      }
      
      seed_set <- seed_set[-1]
    }
  }
  
  result <- list(
    clusters = clusters,
    noise = noise,
    n_clusters = cluster_id
  )
  
  return(result)
}

nipoet_plot_dbscan_results <- function(data, result) {
  if (ncol(data) != 2) {
    stop("Plotting only works for 2D data")
  }
  
  plot_df <- data.frame(
    x = data[,1],
    y = data[,2],
    cluster = factor(result$clusters),
    point_type = ifelse(result$noise, "Noise", "Cluster Point")
  )
  
  p <- plot_ly() %>%
    # Add cluster points
    add_trace(
      data = subset(plot_df, point_type == "Cluster Point"),
      x = ~x,
      y = ~y,
      color = ~cluster,
      type = 'scatter',
      mode = 'markers',
      marker = list(size = 10),
      name = ~paste("Cluster", cluster),
      showlegend = TRUE
    ) %>%
    # Add noise points
    add_trace(
      data = subset(plot_df, point_type == "Noise"),
      x = ~x,
      y = ~y,
      type = 'scatter',
      mode = 'markers',
      marker = list(
        symbol = 'x',
        size = 10,
        color = 'black'
      ),
      name = 'Noise',
      showlegend = TRUE
    ) %>%
    # Update layout
    layout(
      title = "DBSCAN Clustering Results",
      xaxis = list(title = "X"),
      yaxis = list(title = "Y"),
      hovermode = 'closest'
    )
  
  return(p)
}



