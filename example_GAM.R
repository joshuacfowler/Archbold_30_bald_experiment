library(brms)
library(tidyverse)
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

smooth_mean <- apply(linpred_all, 2, mean)
smooth_lower_ci <- apply(linpred_all, 2, quantile, probs = 0.025)
smooth_upper_ci <- apply(linpred_all, 2, quantile, probs = 0.975)

# You can now plot these results

plot(dat$x, smooth_mean, type = "p", 
     main = "Smooth Term Contribution to Linear Predictor", 
     xlab = "x", ylab = "f(x)")
points(dat$x, smooth_lower_ci, lty = 2)
points(dat$x, smooth_upper_ci, lty = 2)
dev.off()


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

new_dat <- data.frame(x = seq(0, 10, length.out = 50))
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



# Another way to generate these similar kinds of matrices which are what is done before fitting model in mgcv
sm <- smoothCon(
  s(x),
  data=new_dat,
  absorb.cons = T,
  diagonal.penalty = T
)

# function that restructures basis in terms of random effects formulation of splines
re <- smooth2random(sm[[1]], "", type=2)

X_new <- re$Xf
Z_new <- re$rand$Xr


brms:::data_sm(brms:::model.frame.brmsfit(fit_gam), newdata = new_dat)

stan_data <- make_standata(fit_gam)




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

# coef <- c(post_draws$b_Intercept[1],
#           post_draws$bs_sx_1[1],
#           post_draws)
post_int <- mean(post_draws$b_Intercept)
post_bX <- mean(post_draws$bs_sx_1)
post_smooth <- unlist(colMeans(post_draws[,6:13])) # "s_sx_1[1]":"s_sx_8[1]"

# post_smooth_matrix <- matrix(data = NA, nrow = 1, ncol = 9)
# post_smooth_matrix[,1] <- post_draws$bs_sx_1[1]
# post_smooth_matrix[,2:9] <- as.matrix(post_draws[,grep("s_sx_1\\[", names(post_draws), value = TRUE)][1,])
# post_smooth_intercept <- post_draws$bs_sx_1[1]
# post_smooth_matrix <- as.matrix( post_draws[,grep("s_sx_1\\[", names(post_draws), value = TRUE)])[1,]
prep <- prepare_predictions(fit_gam, newdata = new_dat)


pred <- tibble(lin_pred = post_int + post_bX*X_new+ Z_new %*% (post_smooth), # 
               lin_pred_prep = mean(prep$dpars$mu$fe$b) + prep$dpars$mu$sm$fe$Xs*mean(prep$dpars$mu$sm$fe$bs) + prep$dpars$mu$sm$re$sx$Zs[[1]] %*% colMeans(prep$dpars$mu$sm$re$sx$s[[1]]),
               lin_rnorm = rnorm(lin_pred_prep, sd = mean(prep$dpars$sigma)),
               x = new_dat$x,
               default = predict(fit_gam, new_dat)[,"Estimate"],
               just_smooth = c(colMeans(brms:::posterior_smooths(fit_gam, newdata = new_dat, smooth = "s(x)"))),
               smooth_plus = post_int + c(colMeans(brms:::posterior_smooths(fit_gam, newdata = new_dat, smooth = "s(x)"))))
                 
          

# plot( pred$x, pred$lin_pred)


ggplot(dat)+
  geom_point(aes(x = x, y = y))+
  geom_line(data = pred, aes(x = x, y = lin_pred) )+
  geom_line(data = pred, aes(x = x, y = default), color = "red")+
  geom_line(data = pred, aes(x = x, y = lin_pred_prep), color = "orange")+
  geom_line(data = pred, aes(x = x, y = lin_pred_prep), color = "blue")
  geom_line(data = pred, aes(x = x, y = just_smooth), color = "blue")+
  geom_line(data = pred, aes(x = x, y = smooth_plus), color = "green")





Xp <- predict(fit_gam_mgcv,new_dat,type="lpmatrix") 

## Xp %*% coef(b) yields vector of predictions
pred <- Xp %*% coef(fit_gam_mgcv)






Xp <-  predict(fit_gam_mgcv, newdata = new_dat, type = "lpmatrix")

xn <- c(.341,.122,.476,.981) ## want prediction at these values
x0 <- 1         ## intercept column
dx <- 1/30      ## covariate spacing in `newd'
for (j in 0:2) { ## loop through smooth terms
  cols <- 1+j*9 +1:9      ## relevant cols of Xp
  i <- floor(xn[j+1]*30)  ## find relevant rows of Xp
  w1 <- (xn[j+1]-i*dx)/dx ## interpolation weights
  ## find approx. predict matrix row portion, by interpolation
  x0 <- c(x0,Xp[i+2,cols]*w1 + Xp[i+1,cols]*(1-w1))
}
dim(x0)<-c(1,28) 
fv <- x0%*%coef(b) + xn[4];fv    ## evaluate and add offset
se <- sqrt(x0%*%b$Vp%*%t(x0));se ## get standard error
## compare to normal prediction
predict(b,newdata=data.frame(x0=xn[1],x1=xn[2],
                             x2=xn[3],x3=xn[4]),se=TRUE)








