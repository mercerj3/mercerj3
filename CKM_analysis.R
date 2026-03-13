## =============================================================================
## CKM MAIHDA Analysis - Stan Model Post-Processing
## Design-weighted Bayesian multilevel logistic regression
## =============================================================================

library(brms)
library(rstan)
library(dplyr)
library(ggplot2)
library(bayesplot)
library(gridExtra)

## -----------------------------------------------------------------------------
## Helper Functions
## -----------------------------------------------------------------------------

# Rename MCMC array parameters with human-readable group names
renameArray <- function(model, names, start, finish) {
  array <- as.array(model)
  for (i in start:finish) {
    dimnames(array)[[3]][i] <- names[i - start + 1]
  }
  return(array)
}

# VPC calculation (mean-based)
variance <- function(models) {
  vpc <- matrix(nrow = length(models), ncol = 3)
  colnames(vpc) <- c("Mean VPC", "Lower CI", "Upper CI")
  var_sd <- matrix(nrow = length(models), ncol = 2)
  colnames(var_sd) <- c("Variance", "SD")
  for (i in 1:length(models)) {
    model <- models[[i]]
    sd <- rstan::extract(model, pars = "sd_1")
    sd <- sd$sd_1
    sigma2 <- sd^2
    var_sd[i, 1] <- mean(sigma2)
    var_sd[i, 2] <- mean(sd)
    fullPosterior <- sigma2 / (sigma2 + ((pi^2) / 3)) * 100
    vpc[i, 1] <- round(mean(fullPosterior), 2)
    vpc[i, 2:3] <- round(quantile(fullPosterior, prob = c(0.025, 0.975)), 2)
  }
  output <- list(var_sd, vpc)
  names(output) <- c("var_sd", "vpc")
  return(output)
}

# VPC calculation (median-based)
variance_med <- function(models) {
  vpc <- matrix(nrow = length(models), ncol = 3)
  colnames(vpc) <- c("Median VPC", "Lower CI", "Upper CI")
  var_sd <- matrix(nrow = length(models), ncol = 2)
  colnames(var_sd) <- c("Variance", "SD")
  for (i in 1:length(models)) {
    model <- models[[i]]
    sd <- rstan::extract(model, pars = "sd_1")
    sd <- sd$sd_1
    sigma2 <- sd^2
    var_sd[i, 1] <- median(sigma2)
    var_sd[i, 2] <- median(sd)
    fullPosterior <- sigma2 / (sigma2 + ((pi^2) / 3)) * 100
    vpc[i, 1] <- round(median(fullPosterior), 2)
    vpc[i, 2:3] <- round(quantile(fullPosterior, prob = c(0.025, 0.975)), 2)
  }
  output <- list(var_sd, vpc)
  names(output) <- c("var_sd", "vpc")
  return(output)
}

## -----------------------------------------------------------------------------
## 60 Intersectional Strata Names (5 race x 4 gender x 3 age)
## -----------------------------------------------------------------------------

groupName <- c(
  "White cismen 18-29", "White cismen 30-54", "White cismen 55+",
  "White transwomen 18-29", "White transwomen 30-54", "White transwomen 55+",
  "White ciswomen 18-29", "White ciswomen 30-54", "White ciswomen 55+",
  "White transmen 18-29", "White transmen 30-54", "White transmen 55+",
  "Black cismen 18-29", "Black cismen 30-54", "Black cismen 55+",
  "Black transwomen 18-29", "Black transwomen 30-54", "Black transwomen 55+",
  "Black ciswomen 18-29", "Black ciswomen 30-54", "Black ciswomen 55+",
  "Black transmen 18-29", "Black transmen 30-54", "Black transmen 55+",
  "Asian cismen 18-29", "Asian cismen 30-54", "Asian cismen 55+",
  "Asian transwomen 18-29", "Asian transwomen 30-54", "Asian transwomen 55+",
  "Asian ciswomen 18-29", "Asian ciswomen 30-54", "Asian ciswomen 55+",
  "Asian transmen 18-29", "Asian transmen 30-54", "Asian transmen 55+",
  "Other/Multi cismen 18-29", "Other/Multi cismen 30-54", "Other/Multi cismen 55+",
  "Other/Multi transwomen 18-29", "Other/Multi transwomen 30-54", "Other/Multi transwomen 55+",
  "Other/Multi ciswomen 18-29", "Other/Multi ciswomen 30-54", "Other/Multi ciswomen 55+",
  "Other/Multi transmen 18-29", "Other/Multi transmen 30-54", "Other/Multi transmen 55+",
  "Hispanic cismen 18-29", "Hispanic cismen 30-54", "Hispanic cismen 55+",
  "Hispanic transwomen 18-29", "Hispanic transwomen 30-54", "Hispanic transwomen 55+",
  "Hispanic ciswomen 18-29", "Hispanic ciswomen 30-54", "Hispanic ciswomen 55+",
  "Hispanic transmen 18-29", "Hispanic transmen 30-54", "Hispanic transmen 55+"
)

## -----------------------------------------------------------------------------
## M2 Model Formula (for reference)
## brmsformula(
##   as.factor(ckm)|weights(wt_norm) ~
##     1 + year + black + asian + other_multi + hispanic +
##     transwomen + cismen + transmen +
##     age3054 + age55plus + (1 | stratum)
## )
##
## Parameter mapping:
##   b[1]:  year2018       b[2]:  year2019       b[3]:  year2020
##   b[4]:  year2021       b[5]:  year2022       b[6]:  year2023
##   b[7]:  black          b[8]:  asian          b[9]:  other_multi
##   b[10]: hispanic       b[11]: transwomen     b[12]: cismen
##   b[13]: transmen       b[14]: age3054        b[15]: age55plus
##   b_Intercept: Intercept
##
## Reference categories: White, Ciswomen, Age 18-29, Year 2017
## -----------------------------------------------------------------------------

## Corrected parameter labels for M2 fixed effects
## NOTE: b[12] corresponds to "cismen" (Cis Men), NOT "Cis Women".
##       Ciswomen is the REFERENCE category for gender.

m2_parameter_labels <- c(
  "Year 2018", "Year 2019", "Year 2020", "Year 2021",
  "Year 2022", "Year 2023",
  "Black", "Asian", "Other/Multiracial", "Hispanic",
  "Trans Women", "Cis Men", "Trans Men",
  "Age 30-54", "Age 55+", "Intercept"
)

## -----------------------------------------------------------------------------
## M2 Mean Estimates (Odds Ratios with 95% Credible Intervals)
## -----------------------------------------------------------------------------

# Load M2 stanfit object (uncomment when running):
# CKM_m2 <- readRDS("CKM_m2.stan")

# pyMDE_m2_sum <- summary(CKM_m2)$summary

# pyMDE_m2_mean_estimates <- pyMDE_m2_sum %>%
#   as.data.frame() %>%
#   filter(grepl("^b\\[|^b_Intercept", rownames(.))) %>%
#   mutate(
#     parameter = m2_parameter_labels,
#     OR = format(round(exp(mean), digits = 3), nsmall = 2),
#     LB = format(round(exp(`2.5%`), digits = 3), nsmall = 2),
#     UB = format(round(exp(`97.5%`), digits = 3), nsmall = 2),
#     est = paste0(OR, " (", LB, ", ", UB, ")")
#   ) %>%
#   dplyr::select(parameter, est)

## -----------------------------------------------------------------------------
## M2 Median Estimates (Odds Ratios with 95% Credible Intervals)
## -----------------------------------------------------------------------------

# pyMDE_m2_med_estimates <- pyMDE_m2_sum %>%
#   as.data.frame() %>%
#   filter(grepl("^b\\[|^b_Intercept", rownames(.))) %>%
#   mutate(
#     parameter = m2_parameter_labels,
#     OR = format(round(exp(`50%`), digits = 3), nsmall = 2),
#     LB = format(round(exp(`2.5%`), digits = 3), nsmall = 2),
#     UB = format(round(exp(`97.5%`), digits = 3), nsmall = 2),
#     est = paste0(OR, " (", LB, ", ", UB, ")")
#   ) %>%
#   dplyr::select(parameter, est)
