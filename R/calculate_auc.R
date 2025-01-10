# Function to calculate AUC given predicted probabilities and actual labels
calculate_auc <- function(pred_probs, actual_labels) {
    # Sort predictions and labels by decreasing probability
    ordered_data <- data.frame(
        pred = pred_probs,
        actual = actual_labels
    )
    ordered_data <- ordered_data[order(-ordered_data$pred), ]
    
    # Get total number of positive and negative cases
    n_pos <- sum(ordered_data$actual == 1)
    n_neg <- sum(ordered_data$actual == 0)
    
    # Calculate true positive rate (TPR) and false positive rate (FPR) at each threshold
    tpr <- cumsum(ordered_data$actual == 1) / n_pos
    fpr <- cumsum(ordered_data$actual == 0) / n_neg
    
    # Calculate AUC using trapezoidal rule
    auc <- sum(diff(fpr) * (tpr[-1] + tpr[-length(tpr)]) / 2)
    
    return(auc)
}

# Example usage:
# pred_probs <- c(0.9, 0.8, 0.7, 0.6, 0.4, 0.3, 0.2, 0.1)
# actual_labels <- c(1, 1, 0, 1, 0, 0, 1, 0)
# auc <- calculate_auc(pred_probs, actual_labels)