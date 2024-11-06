# Load necessary packages
#if (!require("data.table")) install.packages("data.table")
if (!require("ggplot2")) install.packages("ggplot2")
#library(data.table)
library(ggplot2)

# Set seed for reproducibility
set.seed(42)

# Generate a noisy signal
time <- seq(1, 100, by = 0.1)             # Time points
true_signal <- sin(time)                  # True signal (sine wave)
noise <- rnorm(length(time), mean = 0, sd = 0.5)  # Gaussian noise
noisy_signal <- true_signal + noise       # Noisy signal (true signal + noise)

# Create a DataFrame for the generated data
data <- data.frame(time = time, signal = noisy_signal)

# Save the generated data to a CSV file
write.csv(data, "D:/PROJET R/datasets/generated_signal.csv", row.names = FALSE)

# Load the generated data from the CSV file
selected_data <- read.csv("D:/PROJET R/datasets/generated_signal.csv")

# Parameters for the Kalman filter
A <- matrix(c(1, 0.1, 0, 1), nrow = 2, ncol = 2) # State transition matrix
B <- matrix(c(0, 1), nrow = 2, ncol = 1)         # Input matrix
H <- matrix(c(1, 0), nrow = 1, ncol = 2)         # Measurement matrix
Q <- matrix(c(0.0001, 0, 0, 0.0001), nrow = 2, ncol = 2) # Process noise covariance
R <- matrix(0.01, nrow = 1, ncol = 1)            # Measurement noise covariance
x0 <- matrix(c(selected_data$signal[1], 0), nrow = 2, ncol = 1) # Initial state
P0 <- diag(2)                                    # Initial estimation error covariance

# Check dimensions of the matrices
n <- nrow(A)
if (ncol(A) != n || nrow(B) != n || ncol(B) != 1 || nrow(Q) != n || ncol(Q) != n || nrow(R) != ncol(R) || length(as.vector(x0)) != n) {
  stop("Incorrect matrix dimensions.")
} else {
  cat("Matrix dimensions are correct.\n")
}

# Solve the Riccati equation iteratively
riccati_result <- solve_riccati_iterative(A, B, Q, R)

# Calculate the optimal control gains
K <- solve(R) %*% t(B) %*% riccati_result

# Kalman filter prediction function
kalman_predict <- function(x, P, F, Q) {
  x_pred <- F %*% x
  P_pred <- F %*% P %*% t(F) + Q
  return(list(x = x_pred, P = P_pred))
}

# Kalman filter update function
kalman_update <- function(x_pred, P_pred, H, R, z, K) {
  y <- z - H %*% x_pred                            # Measurement residual
  S <- H %*% P_pred %*% t(H) + R                   # Residual covariance
  Kt <- P_pred %*% t(H) %*% solve(S)               # Kalman gain
  x_upd <- x_pred + t(K) %*% y                        # Updated state estimate
  P_upd <- (diag(nrow(P_pred)) - Kt %*% H) %*% P_pred # Updated estimate covariance
  return(list(x = x_upd, P = P_upd))
}

# Initialize variables for the Kalman filter
x_est <- x0
P_est <- P0
estimations <- matrix(NA, nrow = length(time), ncol = 2)
estimations[1, ] <- as.vector(x_est)

# Kalman filter loop
for (i in 2:length(time)) {
  # Prediction step
  pred <- kalman_predict(x_est, P_est, A, Q)

  # Update step with observation
  z <- matrix(selected_data$signal[i], ncol = 1)
  upd <- kalman_update(pred$x, pred$P, H, R, z, K)
  x_est <- upd$x
  P_est <- upd$P

  # Store estimations
  estimations[i, ] <- as.vector(x_est)
}

# Prepare data for plotting
results <- data.frame(time = time, original_signal = selected_data$signal, filtered_signal = estimations[, 1])

# Plot the results
p1 <- ggplot(results, aes(x = time)) +
  geom_line(aes(y = original_signal, color = "Original Signal"), size = 1.2, alpha = 1) +
  geom_line(aes(y = filtered_signal, color = "Filtered Signal"), size = 1.2) +
  labs(title = "Kalman Filter with Optimal Control Using Riccati Equations", x = "Time", y = "Signal") +
  theme_minimal()

print(p1)

