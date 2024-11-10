library(rgl)
library(MASS)
library(here)
source(here::here("./EX1/nipoet_functions_ex1.R"))

get_dataset_03 <- function() {
  n <- 500  # Total number of points
  
  # Generate random samples for each distribution
  x1 <- MASS::mvrnorm(n, mu = c(0, 0, 0), Sigma = diag(3))
  x2 <- MASS::mvrnorm(n, mu = c(5, 5, 5), Sigma = diag(3))
  x3 <- MASS::mvrnorm(n, mu = c(-3, -3, -3), Sigma = diag(3))
  
  # Combine all samples into one dataset
  data <- rbind(x1, x2, x3)
  
  return(data)
}

## Function to plot 3D Gaussians with ellipses
nipoet_plot_3d_gaussians <- function(data, em_result, title = "3D Gaussians") {
  plot3d(data, col = "blue", size = 2, 
        xlab = "X Axis",
        ylab = "Y Axis", 
        zlab = "Z Axis")
  
  for (i in 1:length(em_result$params$mu)) {
    mu <- em_result$params$mu[[i]]
    sigma <- em_result$params$sigma[[i]]
    
    # Ensure sigma is a matrix and draw the 3D ellipse
    ellipse <- ellipse3d(sigma, centre = mu, level = 0.5)
    shade3d(ellipse, col = c("red", "green", "yellow")[i], alpha = 0.3)
  }
  
  title3d(main = title)
}

# Multivariate Gaussian PDF for 3D data
nipoet_gaussian_3d_pdf <- function(x, mu, sigma) {
  k <- length(mu)
  det_sigma <- det(sigma)
  inv_sigma <- solve(sigma)
  centered_x <- t(t(x) - mu)
  exp_term <- exp(-0.5 * rowSums((centered_x %*% inv_sigma) * centered_x))
  return((1 / sqrt((2 * pi)^k * det_sigma)) * exp_term)
}

# E-step: Calculate responsibilities for 3D data
nipoet_e_step <- function(data, params) {
  # Calculate weighted Gaussian densities for each component
  p <- sapply(1:length(params$pi), function(i) {
    params$pi[i] * nipoet_gaussian_3d_pdf(data, params$mu[[i]], params$sigma[[i]])
  })
  
  # Normalize to get responsibilities
  total_p <- rowSums(p)
  gammas <- sweep(p, 1, total_p, "/")
  
  return(gammas)
}

# M-step: Update parameters for 3D data
nipoet_m_step <- function(data, gammas) {
  k <- ncol(gammas)
  
  # Update mu, sigma, and pi for each component
  mu_new <- lapply(1:k, function(i) {
    colSums(gammas[, i] * data) / sum(gammas[, i])
  })
  
  sigma_new <- lapply(1:k, function(i) {
    centered <- sweep(data, 2, mu_new[[i]], "-")
    gamma_diag <- diag(gammas[, i])
    sigma <- (t(centered) %*% gamma_diag %*% centered) / sum(gammas[, i])
    return(sigma)
  })
  
  pi_new <- colMeans(gammas)
  
  updated_params <- list(
    mu = mu_new,
    sigma = sigma_new,
    pi = pi_new
  )
  
  return(updated_params)
}

# EM algorithm for 3D data with plotting, configurable for k clusters
nipoet_em_algorithm <- function(data, k, max_iter = 30, tolerance = 1e-4, plot = TRUE, plotting_delay = 1) {
  n <- nrow(data)
  d <- ncol(data)
  
  # Initial guesses for mu, sigma, and pi
  mu <- lapply(1:k, function(i) data[sample(1:n, 1), ])
  sigma <- lapply(1:k, function(i) diag(d))
  pi <- rep(1 / k, k)

  params <- list(mu = mu, sigma = sigma, pi = pi)
  
  log_likelihood <- numeric(max_iter)
  
  for (i in 1:max_iter) {
    # E-step: Compute responsibilities
    gammas <- nipoet_e_step(data, params)
    
    # M-step: Update parameters
    params <- nipoet_m_step(data, gammas)
    
    # Calculate log-likelihood
    likelihood <- rowSums(sapply(1:k, function(j) {
      params$pi[j] * nipoet_gaussian_3d_pdf(data, params$mu[[j]], params$sigma[[j]])
    }))
    log_likelihood[i] <- sum(log(likelihood))
    
    # Check for convergence by comparing the change in log-likelihood
    if (i > 1 && abs(log_likelihood[i] - log_likelihood[i - 1]) < tolerance) {
      message("Convergence reached at iteration ", i)
      break
    }
    
    # Plot if enabled
    if (plot) {
      title <- paste("Iteration", i)
      nipoet_plot_3d_gaussians(data, list(params = params), title)
      Sys.sleep(plotting_delay)
    }
  }

  # Assign labels based on highest responsibility
  labels <- apply(gammas, 1, which.max)
  
  labeled_data <- cbind(data, cluster = labels)
  
  return(labeled_data)
}




nipoet_kmeans <- function(data, k=NULL, plot=FALSE) {

    if(is.null(k)) {
        stop("No number of clusters provided")
    }

    n_samples <- nrow(data)
    n_dim <- ncol(data)

    # initialize centroids by setting each centroid to random sample
    centroids <- matrix(0, nrow = k, ncol = n_dim)
    r_indx <- sample(1:n_samples, k)    # k unique indices
    for (i in 1:k) {
        centroids[i, ] <- data[r_indx[i], ]
    }

    labels <- rep(-1, n_samples)
    old_labels <- rep(-1, n_samples)

    # repeat until strop criterion is satisfied
    max_iter = 900
    for (iter in 1:max_iter) {
        old_labels <- labels
        for (i in 1:n_samples) {
            best_dist <- Inf
            for (label in 1:k) {
                dist <- nipoet_calc_distance("euclidean", centroids[label, ], data[i, ])
                if (dist < best_dist) {
                    best_dist <- dist
                    labels[i] <- label
                }
            }   
        }
       
        for (label in 1:k) {
            # select all indeces where label = centroid
            indices <- which(labels == label)

            if (length(indices) > 0) {  # Check if there are any points assigned to the cluster
                centroids[label, ] <- apply(data[indices, ], 2, mean)  # get mean for every feature/column of samples with same label
            }
        }

        if (plot) {
          title <- paste("Iteration:", iter)

          p <- plot_ly(x = ~data[,1], y = ~data[,2], z = ~data[,3], type = 'scatter3d', mode = 'markers',
                marker = list(size = 5), color = as.factor(labels)) %>%
              layout(title = title,
                    scene = list(
                      xaxis = list(title = "X1"),
                      yaxis = list(title = "X2"),
                      zaxis = list(title = "X3")
                    ))
          p <- p %>% add_markers(
            x = centroids[, 1],
            y = centroids[, 2], 
            z = centroids[, 3],
            marker = list(size = 20, symbol = 'cross'),
            color = 'balck',
            name = 'Centroids'
          )
          print(p)
        }

        # stopping criterion, if all label assigments remain the same, stop
        if (all(labels == old_labels)) {
            cat("Convergence reached after", iter, "iterations.\n")
            break
        }
    }

    # add labels to dataset
    labeled_data <- cbind(data, labels)

    return(labeled_data)
}