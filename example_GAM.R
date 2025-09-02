library(brms)
set.seed(123)

# Generate some data
dat <- data.frame(
  x = runif(100, 0, 10),
  y = NA
)
dat$y <- 5 + 2 * sin(dat$x) + rnorm(100, 0, 1)

# Fit the brms model
fit_gam <- brm(
  y ~ s(x),
  data = dat,
  family = gaussian(),
  chains = 2,
  cores = 2
)

# Get the full linear predictor (fixed effects + smooth terms)
linpred_all <- posterior_linpred(fit_gam, newdata = dat, summary = FALSE)

# Get the fixed effects part of the linear predictor
# In this simple model, the only fixed effect is the intercept
# The posterior samples for the intercept are in the 'b_Intercept' column
b_intercept_post <- as_draws_df(fit_gam)$b_Intercept
fixed_effects <- matrix(b_intercept_post, nrow = length(b_intercept_post), ncol = nrow(dat))

# The smooth term contribution is the difference
smooth_term_contribution <- linpred_all - fixed_effects

# 'smooth_term_contribution' is a matrix of posterior samples for the smooth term
# Each row is a posterior draw, and each column corresponds to an observation
# You can now calculate the mean or credible intervals for the smooth term
<<<<<<< HEAD
smooth_mean <- apply(smooth_term_contribution, 2, mean)
=======
smooth_mean <- apply(linpred_all, 2, mean)
>>>>>>> 8b5c852852f15ed42e216c4ac9be805a8a9c2300
smooth_lower_ci <- apply(smooth_term_contribution, 2, quantile, probs = 0.025)
smooth_upper_ci <- apply(smooth_term_contribution, 2, quantile, probs = 0.975)

# You can now plot these results
<<<<<<< HEAD
plot(dat$x, smooth_mean, type = "l", 
     main = "Smooth Term Contribution to Linear Predictor", 
     xlab = "x", ylab = "f(x)")
lines(dat$x, smooth_lower_ci, lty = 2)
lines(dat$x, smooth_upper_ci, lty = 2)

=======
plot(dat$x, smooth_mean, type = "p", 
     main = "Smooth Term Contribution to Linear Predictor", 
     xlab = "x", ylab = "f(x)")
points(dat$x, smooth_lower_ci, lty = 2)
points(dat$x, smooth_upper_ci, lty = 2)
dev.off()
>>>>>>> 8b5c852852f15ed42e216c4ac9be805a8a9c2300


# Example for getting design matrix from mgcv
#Step 1: Fit a model with mgcv

#First, we'll fit a model using the gam function from the mgcv package. The structure of the model object from mgcv is much more transparent for this purpose.


library(mgcv)

# Create example data
set.seed(123)
dat <- data.frame(x = runif(100, 0, 10))
dat$y <- 5 + 2 * sin(dat$x) + rnorm(100, 0, 1)

# Fit the GAM model
fit_gam_mgcv <- gam(y ~ s(x), data = dat)

# Step 2: Create new data

#Now, create a new data frame for which you want to generate the design matrix.


# Create new data

new_dat <- data.frame(x = c(1,10))#seq(0, 10, length.out = 50))
# Step 3: Get the design matrix for new data

#Use the predict.gam function with the type = "lpmatrix" argument. This will return the linear predictor design matrix, which includes both the intercept and the smooth term basis functions.


# Generate the design matrix for the new data
X_new <- predict(fit_gam_mgcv, newdata = new_dat, type = "lpmatrix")

# Examine the dimensions of the matrix
dim(X_new)
# This will have dimensions: (number of new data points) x (number of basis functions + 1)
# The first column is for the intercept, and the remaining columns are for the smooth term
#Step 4: Extract the smooth term part

#The first column of the matrix corresponds to the intercept. The rest of the columns are for the smooth term.


# The first column is the intercept
intercept_design <- X_new[, 1]

# The remaining columns are for the smooth term basis functions
smooth_design_matrix <- X_new[, -1]

# You can check the column names to confirm
colnames(smooth_design_matrix)

# Note: In a brms context, the predict function is handled internally. brms extracts the necessary information from the mgcv object it creates during the model fitting process to perform this calculation for you automatically when you use posterior_linpred() or posterior_predict() with new data.

# refitting the brms model cause simulated data is slightly different
fit_gam <- brm(
  y ~ s(x),
  data = dat,
  family = gaussian(),
  chains = 2,
  cores = 2
)

# then extracting the intercept and smooth parameter posteriors which will get multiplied by the smooth design matrix in the linear predictor
post_draws <- as_draws_df(fit_gam)
post_int <- post_draws$b_Intercept[1]

post_smooth_intercept <- post_draws$bs_sx_1[1]
post_smooth_matrix <- as.matrix( post_draws[,grep("s_sx_1\\[", names(post_draws), value = TRUE)])

lin_pred <- post_int + post_smooth_intercept*new_dat + smooth_design_matrix%*% t(post_smooth_matrix[1,])
