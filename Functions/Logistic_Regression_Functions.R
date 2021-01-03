# Predicted Probabilities for Each Cluster Based on STAN Model Fit
predict_probabilities <- function(Cluster_Identity_Vector, Identifying_String, STAN_Model_Fit, Number_Clusters) {
  pred_probs <- matrix(nrow = length(Cluster_Identity_Vector), ncol = Number_Clusters)
  model_probs <- STAN_Model_Fit[[Identifying_String]]
  for (i in 1:length(Cluster_Identity_Vector)) {
    for (j in 1:Number_Clusters){
      individual_probability <- mean(model_probs[, i, j])
      pred_probs[i, j] <- individual_probability
    }
  }
  return(pred_probs)
}

# Generate Cluster Predictions Based on Predicted Probabilities
predict_cluster <- function(predicted_probabilities) {
  predicted_cluster <- c()
  for (i in 1:length(predicted_probabilities[, 1])) {
    predicted_cluster[i] <- which(predicted_probabilities[i, ] == max(predicted_probabilities[i, ]))
  }
  return(predicted_cluster)
}

# Calculate Whether Predicted Clusters Were Correct Or Not 
calculate_correct_predictions <- function(Actual_Cluster, Predicted_Cluster) {
  counter <- 0
  for (i in 1:length(Actual_Cluster)) {
    if (Actual_Cluster[i] == Predicted_Cluster[i]) {
      counter <- counter + 1
    }
  }
  return(counter)
}

# Calculate Whether Predicted Clusters Were Correct, Treating Bimodal and Flat Clusters As the Same
calculate_correct_predictions_smush <- function(Actual_Cluster, Predicted_Cluster) {
  counter <- 0
  for (i in 1:length(Actual_Cluster)) {
    if ((Actual_Cluster[i] == 1 | Actual_Cluster[i] == 2) & (Predicted_Cluster[i] == 1 | Predicted_Cluster[i] == 2)) {
      if (Actual_Cluster[i] == Predicted_Cluster[i]) {
        counter <- counter + 1
      }
    }
    else {
      if ((Actual_Cluster[i] == 3 | Actual_Cluster[i] == 4) & (Predicted_Cluster[i] == 3 | Predicted_Cluster[i] == 4)) {
        counter <- counter + 1
      }
    }
  }
  return(counter)
}







