#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Experiments using a low rank Sigma when doing knockoffs for microbiome data
##
## author: sankaran.kris@gmail.com
## date: 09/15/2017

###############################################################################
## Packages and some utility functions
###############################################################################
library("MFKnockoffs")
library("phyloseq")
library("tidyverse")
library("caret")
library("glasso")
library("cate")

scale_colour_discrete <- function(...)
  scale_colour_brewer(..., palette="Set2")
theme_set(theme_bw())
theme_update(
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.ticks = element_blank(),
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 6),
  axis.text = element_text(size = 6),
  axis.title = element_text(size = 8),
  strip.background = element_blank(),
  strip.text = element_text(size = 8),
  legend.key = element_blank()
)

#' Symmetrize a matrix
symm <- function(x) {
  (x + t(x)) / 2
}

#' Construct Knockoffs
#'
#' This is just a wrapper of MFKnockoffs.create.gaussian() that names knockoff
#' variables explicitly
knockoff_wrapper <- function(x, sigma_approx) {
  x_tilde <- MFKnockoffs.create.gaussian(
    x,
    rep(0, ncol(x)),
    symm(sigma_approx),
    method = "asdp"
  )

  colnames(x_tilde) <- paste0("knockoff_", 1:ncol(x_tilde))
  x_tilde
}

#' Fit model on true + knockoffs
fit_wrapper <- function(x, x_tilde, y) {
  train(
    cbind(x, x_tilde),
    y,
    method = "glmnet",
    intercept = FALSE
  )
}

#' Calculate FDP estimates
#'
#' This calculates FDP along the lambda path used by the fit (caret lasso
#' object)
fdp_estimates <- function(fit) {
  beta_hat <- coef(fit$finalModel)
  p <- ncol(beta_hat)

  fdp_hat <- data_frame(
    "lambda" = fit$finalModel$lambda,
    "nonzero" = colSums(beta_hat != 0),
    "fdp" = rep(0, p)
  )

  for (j in seq_len(p)) {
    nonzero_beta <- beta_hat[beta_hat[, j] != 0, j]
    fdp_hat$fdp[j] <- mean(grepl("knockoff", names(nonzero_beta)))
  }

  fdp_hat
}

###############################################################################
## Load and transform the data
###############################################################################
ps <- readRDS("../data/ps.rds")
x <- get_taxa(ps) %>%
  asinh() %>%
  scale(center = TRUE, scale = FALSE)
y <- sample_data(ps) %>%
  .[["age"]]

## factor analysis
fact_res <- factor.analysis(x, 2)
sigma_approx <- fact_res$Gamma %*% t(fact_res$Gamma) + diag(fact_res$Sigma)
x_tilde <- knockoff_wrapper(x, sigma_approx)
fit <- fit_wrapper(x, x_tilde, y)
fdp_hat <- fdp_estimates(fit) %>%
  mutate(method = "factor")

## graphical lasso
glasso_res <- glasso(var(x), rho = 5e-2)
sigma_glasso <- symm(solve(glasso_res$w))
x_tilde <- knockoff_wrapper(x, sigma_approx)
fit <- fit_wrapper(x, x_tilde, y)
fdp_hat <- fdp_estimates(fit) %>%
  mutate(method = "glasso") %>%
  bind_rows(fdp_hat)

ggplot(fdp_hat) +
  geom_point(
    aes(x = nonzero, y = fdp, col = method)
  )
