# Method 1: Using base R
plot_roc <- function(pred_probs, actual_labels) {
    # Sort by decreasing probability
    ordered_data <- data.frame(
        pred = pred_probs,
        actual = actual_labels
    )
    ordered_data <- ordered_data[order(-ordered_data$pred), ]
    
    # Calculate TPR and FPR
    n_pos <- sum(ordered_data$actual == 1)
    n_neg <- sum(ordered_data$actual == 0)
    tpr <- cumsum(ordered_data$actual == 1) / n_pos
    fpr <- cumsum(ordered_data$actual == 0) / n_neg
    
    # Plot
    plot(fpr, tpr, type = "l", col = "blue",
         xlab = "False Positive Rate", 
         ylab = "True Positive Rate",
         main = "ROC Curve")
    abline(0, 1, lty = 2, col = "gray")  # Add diagonal reference line
    
    # Calculate and display AUC
    auc <- sum(diff(fpr) * (tpr[-1] + tpr[-length(tpr)]) / 2)
    legend("bottomright", 
           legend = sprintf("AUC = %.3f", auc),
           bty = "n")
}

# Method 2: Using pROC package (recommended for practical use)
library(pROC)

plot_roc_pROC <- function(actual_labels, pred_probs, plot_title = "ROC Curve") {
    # Create ROC object
    roc_obj <- roc(actual_labels, pred_probs)
    
    # Plot
    plot(roc_obj, 
         main = plot_title,
         col = "blue",
         lwd = 2)
    
    # Add diagonal reference line
    abline(0, 1, lty = 2, col = "gray")
    
    # Add AUC to plot
    legend("bottomright", 
           legend = sprintf("AUC = %.3f", auc(roc_obj)),
           bty = "n")
}

#' # Example usage:
#' # Generate sample data
#' set.seed(123)
#' pred_probs <- runif(100)
#' actual_labels <- rbinom(100, 1, 0.5)

#' # Create plots using both methods
#' par(mfrow = c(1, 2))  # Side by side plots
#' plot_roc(pred_probs, actual_labels)
#' title(sub = "Base R")
#' plot_roc_pROC(pred_probs, actual_labels)
#' title(sub = "pROC package")