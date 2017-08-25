#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Experiments using a low rank Sigma when doing knockoffs for microbiome data
##
## author: sankaran.kris@gmail.com
## date: //2017

library("MFKnockoffs")
library("phyloseq")
library("tidyverse")
library("caret")
library("glasso")
library("cate")

symm <- function(x) {
  (x + t(x)) / 2
}

ps <- readRDS("../data/ps.rds")

x <- get_taxa(ps) %>%
  asinh() %>%
  scale(center = TRUE, scale = FALSE)
y <- sample_data(ps) %>%
  .[["age"]]

## factor analysis
fact_res <- factor.analysis(x, 2)
sigma_approx <- fact_res$Gamma %*% t(fact_res$Gamma) + diag(fact_res$Sigma)

x_tilde <- MFKnockoffs.create.gaussian(
  x,
  rep(0, ncol(x)),
  symm(sigma_approx),
  method = "asdp"
)

colnames(x) <- paste0("true_", 1:ncol(x))
colnames(x_tilde) <- paste0("knockoff_", 1:ncol(x_tilde))

## fit model on true + knockoffs
fit <- train(
  cbind(x, x_tilde),
  y,
  method = "glmnet",
  intercept = FALSE
)

## estimate FDP
beta_hat <- coef(fit$finalModel)[, 10]
nonzero_beta <- beta_hat[beta_hat != 0]
fdp_hat <- mean(grepl("knockoff", names(nonzero_beta)))
fdp_hat

## graphical lasso
glasso_res <- glasso(var(x), rho = 5e-2)
sigma_glasso <- symm(solve(glasso_res$w))
x_tilde <- MFKnockoffs.create.gaussian(
  x,
  rep(0, ncol(x)),
  sigma_glasso,
  method = "asdp"
)

## again fit model
fit <- train(
  cbind(x, x_tilde),
  y,
  method = "glmnet",
  intercept = FALSE
)

## estimate FDP
beta_hat <- coef(fit$finalModel)[, 30]
nonzero_beta <- beta_hat[beta_hat != 0]
fdp_hat <- mean(grepl("knockoff", names(nonzero_beta)))
fdp_hat
