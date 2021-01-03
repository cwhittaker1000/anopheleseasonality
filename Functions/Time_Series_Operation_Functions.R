# Entropic measure from Nguyen et al- note, requires time series normalised by total
entropic_measure <- function(timepoints, data) {
  entropy <- 0
  for (i in 1:length(data)) {
    entropy <- entropy + data[i] * log2(data[i]/(1/length(data)))
  }
  return(entropy)
}

# Calculates the Peak Distance from January On a Circle
calculate_peak_distance_from_jan <- function(timepoints, data) {
  peak_index <- which(data == max(data))
  time <- timepoints[peak_index]
  distance_from_jan <- min(abs(time - 0.33), max(timepoints) - time + 0.33) # need to check when it goes with the second option that 0's not implicitly include in the distance there
  return(distance_from_jan)
  # or at the very least, it needs to be consistent with everything else that's going on elsewhere in the code
}

# Calculates the Proportion of Points that Are Above a Specified Threshold (as a Function of the Mean)
points_greater_than_mean_multiple <- function(data, multiple) {
  counts <- na.trim(data)
  mean <- mean(counts)
  greater_than_mean <- sum(counts > (multiple * mean))/length(counts)
  return(greater_than_mean)
} 

# Clustering and Characterising the Time Series Properties 
cluster_characterisation <- function(time_series_properties, number_clusters, time_series) {
  
  set.seed(58)
  par(mfrow = c(1, 1))
  palette(c("#F15025", "#7ACC70", "#00A7E1", "#F2328C", "#E71D36", "#52AD9C", "#7761B5", "#3F220F", "#D6D84F", "#363537", "#0CCE6B", "#2A628F", "#FF9B54"))
  
  # Running PCA on the time series properties and then clustering the outputs
  number_properties <- ncol(time_series_properties)
  property_names <- colnames(time_series_properties)
  PCA <- prcomp(time_series_properties)
  loadings <- PCA$rotation[, 1:number_properties]
  PCA_output <- as.matrix(time_series_properties) %*% loadings
  clustering_results <- kmeans(PCA_output, number_clusters, nstart = 20)
  plot(PCA_output[, 2], PCA_output[, 1], col = clustering_results$cluster, pch = 20, xlab = "PCA Comp 2 (20%)", ylab = "PCA Comp 1 (55%)", cex = 2, las = 1)
  clusters <- as.character(seq(1:number_clusters))
  legend(-2.5, 4, as.character(seq(1:number_clusters)), cex = 1.5, col = palette()[1:number_clusters], pch = 20)

  # Visualising the time series belonging to each cluster
  number_clusters <- max(clustering_results$cluster)
  if (number_clusters %% 3 == 0) {
    number_rows <- number_clusters/3 
  } else {
    number_rows <- ceiling(number_clusters/3)
  }
  par(mfrow = c(number_rows, 3))
  max <- 10 
  for (i in 1:number_clusters) {
    cluster <- time_series[clustering_results$cluster == i, ]
    plot(timepoints, apply(cluster, 2, mean) * 100, type = "l", ylim = c(0, max), lwd = 2, col = palette()[i], las = 1, xaxt = "n")
    for (j in 1:length(cluster[, 1])) {
      lines(timepoints, cluster[j, ] * 100, col = adjustcolor(palette()[i], alpha.f = 0.2))
    }
    number_time_series <- length(cluster[, 1])
    text(1.5, 9.5, paste0("n = ", number_time_series), cex = 1.5, col = "grey20")
  }
  
  # Visualising the distribution of the properties of the time series belonging to each cluster
  par(mfrow = c(number_clusters, number_properties + 1))
  max <- 10 
  mean_operation_values <- matrix(nrow = number_clusters, ncol = number_properties)
  colnames(mean_operation_values) <- property_names
  for (i in 1:number_clusters) {
    if (i < number_clusters) {
      subsetter <- clustering_results$cluster == i
      cluster_time_series <- time_series[subsetter, ]
      cluster_time_series_properties <- time_series_properties[subsetter, ]
      plot(timepoints, apply(cluster_time_series, 2, mean) * 100, type = "l", ylim = c(0, max), col = palette()[i], lwd = 2)
      for (j in 1:number_properties) {
        hist(cluster_time_series_properties[, j], col = palette()[i], main = "", xlim = c(min(cluster_time_series_properties[, j]), max(cluster_time_series_properties[, j])))
      }
      mean_operation_values[i, ] <- apply(cluster_time_series_properties, 2, mean)
    } else {
      subsetter <- clustering_results$cluster == i
      cluster_time_series <- time_series[subsetter, ]
      cluster_time_series_properties <- time_series_properties[subsetter, ]
      plot(timepoints, apply(cluster_time_series, 2, mean) * 100, type = "l", ylim = c(0, max), col = palette()[i], lwd = 2)
      for (j in 1:number_properties) {
        hist(cluster_time_series_properties[, j], col = palette()[i], main = "", xlim = c(min(cluster_time_series_properties[, j]), max(cluster_time_series_properties[, j])))
        mtext(property_names[j], side = 1, outer = FALSE, cex = 1, font = 2, line = 3, col = "grey20")
      }
      mean_operation_values[i, ] <- apply(cluster_time_series_properties, 2, mean)
    }
  }
  
  return(list(PCA_Output = PCA_output,
              Mean_Cluster_Values = as.data.frame(Cluster = seq(1, 8, 1), mean_operation_values),
              Cluster_Numbers = clustering_results$cluster))
}

