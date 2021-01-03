softmax <- function(values) {
  softmax_values <- c()
  for (i in 1:length(values)) {
    softmax_values[i] <- exp(values[i])/sum(exp(values))
  }
  return(softmax_values)
}


# species_coefficient <- as.vector(species_coefficient_matrix[1, ])
# environmental_coefficients <- as.matrix(environmental_coefficients)
# species_raster_value <- ann_distribution[400000]
# environmental_raster_values <- as.matrix(raster_value_storage[400000, ], nrow = 1)
# generate_single_species_predictions(species_coefficient, environmental_coefficients, species_raster_value, environmental_raster_values)

# generates predictions for temporal modality for a single species in a single location
generate_single_species_predictions <- function(species_coefficient, environmental_coefficients, species_raster_value, environmental_raster_values) {
  predicted_probabilities <- c()
  environmental_component <- environmental_raster_values %*% environmental_coefficients
  predicted_prob <- species_coefficient + environmental_component
  predicted_prob <- softmax(as.vector(predicted_prob))
  predicted_probabilities <- species_raster_value * predicted_prob
  return(predicted_probabilities)
}

# generates predictions for temporal modality for all species, for a given location
generate_all_species_predictions <- function(species_coefficient_matrix, environmental_coefficients, species_raster_values, environmental_raster_values) {
  cluster_and_species_probability <- matrix(nrow = 7, ncol = 4)
  for (i in 1:7) {
    cluster_and_species_probability[i, ] <- generate_single_species_predictions(species_coefficient_matrix[i, ], environmental_coefficients, species_raster_values[i], environmental_raster_values)
  }
  return(cluster_and_species_probability)
}


generate_probability_temporal_presence <- function(cluster_and_species_probability) {
  probability_presence <- c()
  for (i in 1:4) {
    prob_absence <- 1
    for (j in 1:7) {
      prob_absence <- prob_absence * (1 - cluster_and_species_probability[j, i])
    } 
    probability_presence[i] <- 1 - prob_absence
  }
  return(probability_presence)
}

## Old
# 
# for(i in 1:4) { 
#   species_component <-  1 * species_coefficient[i] # species contrib to probability of seasonal mode, conditional on presence
#   environmental_component <- environmental_raster_values * environmental_coefficients[, i]
#   predicted_probabilities[i] <- species_component + sum(environmental_component) # note that there are still NAs in the dataset. Might have to do some processing if they occur in India proper
# }
# species_prob_presence <- species_raster_value
# predicted_probabilities <- species_prob_presence * softmax(predicted_probabilities)
