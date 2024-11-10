#setwd() #how am I supposed to do this?
library(FNN)
library(ggplot2)
library(gridExtra) # For arranging plots side-by-side

# distance function from ha1
nipoet_calc_distance <- function(method=NaN, p1, p2) {

  if (length(p1)!=length(p2)) {
    stop("Points have different dimenions")
  }

  if (identical(method, "manhattan")) {
    sum <- 0
    for (i in 1:length(p1)) {
        sum <- sum + abs(p1[i] - p2[i])
    }
    return (sum)
  }
  else if (identical(method, "euclidean")) {
    sum <- 0
    for (i in 1:length(p1)) {
      sum <- sum + (p1[i] - p2[i])^2
    }
    return (sum^0.5)
  }
  else if (identical(method, "chebyshev")) {

    distances<-abs(p1[1]-p2[1])
    for (i in 2:length(p1)) {
        distances[i]<-abs(p1[i]-p2[i])
    }
    sort(distances)
    return (distances[length(p1)])
  }
  
  stop("Invalid method specified.")
}

nipoet_compute_distances <- function(data) {
    size <- nrow(data)
    distances <- rep(Inf, size)

    for (i in 1:size) {
        for (j in 1:size) {
            dist <- nipoet_calc_distance("euclidean", data[i, ], data[j, ])
            if (dist < distances[i] && i!=j) {  # ensure no the same samples
                distances[i] = dist
            }
        }
    }

    return(distances)
}

nipoet_hopkins_statistic <- function(data) {

    n_samples <- nrow(data)

    # get mins and maxs to create Dataset S
    mins <- apply(data, 2, min) # apply min over each column (2)
    maxs <- apply(data, 2, max)

    # create subset of the original dataset (15% of original sample size)
    subset_size <- floor(n_samples*0.15)
    #set.seed(187) # for reproducibility
    sample_idx <- sample(n_samples, subset_size)
    R <- data[sample_idx, ]

    # create dataset S with the same domain as R
    S <-  matrix(nrow = subset_size, ncol = ncol(data))
    for (i in 1:ncol(data)) {
        S[, i] <- runif(subset_size, mins[i], maxs[i])  # for each feature/column draw samples form uniform distribution
    }

    # calculate distances
    beta_distances <- nipoet_compute_distances(S)
    alpha_distances <- nipoet_compute_distances(R)

    hopkins <- sum(beta_distances) / sum(unlist(c(beta_distances, alpha_distances)))

    return(hopkins)
}

nipoet_feature_selection_2d <- function(data) {
  
  combinations <- combn(c(1:ncol(data)), 2)
  
  best_combination <- NULL
  highest_hopkins <- -Inf
  
  # iterate over each combination of features
  for (i in 1:ncol(combinations)) {
    
    hopkins_value <- nipoet_hopkins_statistic(data[, combinations[, i]])
    
    if (hopkins_value > highest_hopkins) {
      highest_hopkins <- hopkins_value
      best_combination <- combinations[, i]
    }

  }
  
  return(list(best_combination = best_combination, highest_hopkins = highest_hopkins))
}

nipoet_plot_3 <- function(data) {
  # Ensure the data has exactly 3 columns
  if (ncol(data) != 3) {
    stop("Data must have exactly 3 columns.")
  }
  
  colnames(data) <- c("Feature1", "Feature2", "Feature3")
  
  # Create each 2D plot for each pairwise combination
  plot1 <- ggplot(data, aes(x = Feature1, y = Feature2)) +
    geom_point() +
    labs(title = "Feature1 vs Feature2") +
    theme_minimal()
  
  plot2 <- ggplot(data, aes(x = Feature1, y = Feature3)) +
    geom_point() +
    labs(title = "Feature1 vs Feature3") +
    theme_minimal()
  
  plot3 <- ggplot(data, aes(x = Feature2, y = Feature3)) +
    geom_point() +
    labs(title = "Feature2 vs Feature3") +
    theme_minimal()
  
  # Arrange the three plots side-by-side
  grid.arrange(plot1, plot2, plot3, ncol = 3)
}