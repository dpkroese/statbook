# Clear the workspace
rm(list = ls())

# ---------------------------------------
# Define the target density functions
# ---------------------------------------

# A scalar version (used in the MCMC sampling loop)
f_scalar = function(x, y) {
  if (y < 2) {
    return(exp(-4 * (y - x^2)^2 + (y - 1)^2))
  } else {
    return(0)
  }
}

# A vectorized version (used for evaluating the grid for plotting)
f_vec = function(x, y) {
  exp(-4 * (y - x^2)^2 + (y - 1)^2) * (y < 2)
}

# ---------------------------------------
# Part 1: Surface Plot (similar to surfl)
# ---------------------------------------

# Create a grid of x and y values
x_seq = seq(-2, 2, by = 0.05)
y_seq = seq(-2, 2, by = 0.05)

# Compute Z over the grid using the vectorized function
Z = outer(x_seq, y_seq, f_vec)

# Plot a 3D surface using base R's persp.
persp(x_seq, y_seq, Z,
      theta    = 30,       # rotation angle in the x-y plane
      phi      = 30,       # elevation angle
      expand   = 0.5,      # scaling factor for z values
      col      = "gray",   # similar to MATLAB's colormap gray
      shade    = 0.5,      # adds some shading
      ticktype = "detailed",
      main     = "Surface Plot of f(x,y)")

# ---------------------------------------
# Part 2: Random Walk Metropolis Sampling
# ---------------------------------------

N = 10000                          # sample size
xx = matrix(0, nrow = N, ncol = 2)   # initialize matrix with zeros
x = c(0, -1)                       # initial point
xx[1, ] = x

set.seed(123)  # for reproducible results

# Loop to generate MCMC samples using a Random Walk Metropolis algorithm
for (i in 2:N) {
  z = rnorm(2)                 # Draw a standard Gaussian noise vector of length 2
  y = x + z                    # Propose a new state: y = x + z
  fx = f_scalar(x[1], x[2])      # density at current state
  fy = f_scalar(y[1], y[2])      # density at proposed state
  
  # Compute acceptance probability.
  if (fx == 0) {
    alpha = 1
  } else {
    alpha = min(fy / fx, 1)
  }
  
  # Accept the proposal with probability alpha
  if (runif(1) < alpha) {
    x = y
  }
  
  # Store the current (or updated) state
  xx[i, ] = x
}

# ---------------------------------------
# Part 3: Plot the Sampling Results
# ---------------------------------------

# Plot the sample points from the random walk
plot(xx[, 1], xx[, 2],
     pch   = 20,
     col   = "black",
     cex   = 0.3,
     xlab  = "x",
     ylab  = "y",
     main  = "Random Walk Samples with f(x,y) Contours")

# Overlay contour lines of the target density for visual reference
contour(x_seq, y_seq, Z,
        add        = TRUE,
        drawlabels = FALSE,
        col        = "blue",
        lwd        = 2)
