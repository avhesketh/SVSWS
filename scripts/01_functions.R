##### FUNCTIONS ######
# Amelia Hesketh #
# June 2024 #

# Function 1: GAM smoother pairwise comparisons
# This code is derived from Gavin Simpson and generates pairwise comparisons
# that show differences in algal cover smoothers from GAM models
# https://fromthebottomoftheheap.net/2017/10/11/difference-splines-i/


smooth_diff <- function(model, newdata, f1, f2, var, alpha = 0.05,
                        unconditional = FALSE) {
  xp <- predict(model, newdata = newdata, type = 'lpmatrix')
  # extract relevant columns (c1 and c2) and rows (r1 and r2) to compare
  c1 <- grepl(f1, colnames(xp))
  c2 <- grepl(f2, colnames(xp))
  r1 <- newdata[[var]] == f1
  r2 <- newdata[[var]] == f2
  ## difference rows of xp matrix for data from comparison
  X <- xp[r1, ] - xp[r2, ]
  ## zero out cols of X related to splines for other treatments
  X[, ! (c1 | c2)] <- 0
  ## zero out the parametric cols
  X[, !grepl('^s\\(', colnames(xp))] <- 0
  # Now that Xp matrix is modified, can obtain predicted values by multiplying matrix by
  # the estimated model coefficients and summing the result row-wise
  dif <- X %*% coef(model)
  # Compute the standard errors of the differences using the variance covariance matrix 
  # of estimated model coefficients
  se <- sqrt(rowSums((X %*% vcov(model, unconditional = unconditional)) * X))
  # compute the 95% confidence interval about the predicted differences
  crit <- qt(alpha/2, df.residual(model), lower.tail = FALSE)
  upr <- dif + (crit * se)
  lwr <- dif - (crit * se)
  data.frame(pair = paste(f1, f2, sep = '-'),
             diff = dif,
             se = se,
             upper = upr,
             lower = lwr)
}


# Function 2: Adding ellipses to NMDS (or capscale) plots 


## adding ellipses to plot

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

