library(here)   # for absolute paths
source(here::here("./EX3/nipoet_functions_ex3.R"))  # get "nipoet_generate_half_moon" function
library(caret)  # for logistic regression model
library(rgl)

nipoet_generate_half_moon_dataset <- function(n, radius = 2.5, noise = 0.2) {
    moon1 <- nipoet_generate_half_moon(n, radius, noise)
    moon1[, 1] <- moon1[, 1] + radius
    moon1[, 2] <- moon1[, 2] + radius # shift the y-coordinates

    moon2 <- nipoet_generate_half_moon(n, radius, noise)
    moon2[, 2] <- -moon2[, 2] - (radius*0.7)   # flip moon and shift

    # Combine both moons and create labels
    data <- rbind(moon1, moon2)

    return(data)
}

nipoet_convert_dataset_to_3d <- function(data_2d) {
    x_val <- data_2d[, 1]
    y_val <- data_2d[, 2]
    z_val <- rep(0, length(x_val))

    for (x in 1:length(x_val)) {
        if (y_val[x] < sin(x_val[x]-1)*1.3) {
            z_val[x] <- 3
        } else {
            z_val[x] <- -3
        }
    }

    data_3d <- cbind(x_val, y_val, z_val)
    
    return(data_3d)
}

nipoet_2d_logistic_regression <- function(data, labels) {

    # convert to dataframe
    df <- data.frame(x = data[,1], y = data[,2], label = factor(labels))

    model <- glm(label ~ x + y, data = df, family = binomial)

    # Confusion matrix output
    predicted_probabilities <- predict(model, type = "response")
    predicted_labels <- ifelse(predicted_probabilities > 0.5, 1, 0)
    conf_matrix <- confusionMatrix(as.factor(predicted_labels), as.factor(labels))

    # USE ACCURACY ONLY IF DATASET IS BALANCED
    message(paste("Accuracy of the model", conf_matrix$overall['Accuracy']))
    print(conf_matrix$table)

    # plotting the classified data with hyperplane
    plot(df$x, df$y, col = ifelse(predicted_labels == 0, "blue", "red"),
     pch = 19, xlab = "x", ylab = "y",
     main = paste("Dataset separation accuracy ", round(conf_matrix$overall['Accuracy'], 4), " in 2D"))

    coef <- coef(model)
    abline(a = -coef[1] / coef[3], b = -coef[2] / coef[3], col = "black", lwd = 2)
}

nipoet_3d_logistic_regression <- function(data, labels) {

    df <- data.frame(x = data[,1] , y = data[,2], z = data[,3], label = factor(labels))

    model <- glm(label ~ x + y + z, data = df, family = binomial)

    # Confusion matrix output
    predicted_probabilities <- predict(model, type = "response")
    predicted_labels <- ifelse(predicted_probabilities > 0.5, 1, 0)
    conf_matrix <- confusionMatrix(as.factor(predicted_labels), as.factor(labels))
    # print(conf_matrix)

    # USE ACCURACY ONLY IF DATASET IT BALANCED
    message(paste("Accuracy of the model", conf_matrix$overall['Accuracy']))
    print(conf_matrix$table)

    # plotting the classified data with hyperplane
    plot3d(df$x, df$y, df$z, col = ifelse(predicted_labels == 0, "blue", "red"), type = "s", radius = 0.1, main="3D Logistic Regression", 
        xlab = "X Axis",
        ylab = "Y Axis", 
        zlab = "Z Axis")
    coeff <- coef(model)
    planes3d(coeff[2], coeff[3], coeff[4], coeff[1], col = 'blue', alpha = 0.3)

    # # Scatter plot of points
    # p <- plot_ly(x = ~df$x, y = ~df$y, z = ~df$z, color = as.factor(predicted_labels), colors = c('blue', 'red')) %>%
    #     add_markers(size = 4) %>%
    #     layout(scene = list(
    #         xaxis = list(title = "X"),
    #         yaxis = list(title = "Y"),
    #         zaxis = list(title = "Z")
    #     ))

    # # Extract coefficients for the plane
    # coeff <- coef(model)
    # grid_size <- 10  # adjust for resolution of the plane
    # x_seq <- seq(min(df$x), max(df$x), length.out = grid_size)
    # y_seq <- seq(min(df$y), max(df$y), length.out = grid_size)
    # z_vals <- outer(x_seq, y_seq, function(x, y) (-coeff[1] - coeff[2] * x - coeff[3] * y) / coeff[4])

    # # Add plane to the plot
    # p <- p %>% add_surface(x = ~x_seq, y = ~y_seq, z = ~z_vals, opacity = 0.5)
    # p
}