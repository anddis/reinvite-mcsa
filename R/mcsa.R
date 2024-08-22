#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# STHLM3MRI Reinvite: Second screening round - Monte Carlo Sensitivity Analysis
# Plot parameters' distributions & Simulate outcome data
# Andrea Discacciati
# 20240816
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#

# Load packages and functions ----
library("tidyr")
library("MASS")
library("RColorBrewer")
source(here::here("R", "functions.R"))

# Set-up ----
n.sim <- 5000

# MCSA1: PSA only at first screening (+75); MCSA2: S3 only at first screening (+31)
data_discordant <- lapply(c(MCSA1 = 75, 
                            MCSA2 = 31), 
                          \(x) data.frame(studieid = paste0("id", sprintf("%02s", seq_len(x)))))
par0 <-  par()

# Prior distributions' parameters ----
prior_parameters <- list(
  p.participation = list(mu = qlogis(0.72), sigma = sqrt(0.01)),
  
  p.mri = list(mu = qlogis(0.90), sigma = sqrt(0.01)),
  
  p.bx = list(mu = qlogis(.999), sigma = sqrt(0.001)),
  
  p.isup2p = list(mu = qlogis(0.80), sigma = sqrt(0.10)),
  
  # ISUP1
  fpf = list(mu = c(psa3p = qlogis(0.30), s3m11p = qlogis(0.30)), 
             sigma = c(psa3p = sqrt(0.15), s3m11p = sqrt(0.15)),
             rho = 0.8),
  
  # ISUP2+
  tpf = list(mu = c(psa3p = qlogis(0.60), s3m11p = qlogis(0.60*1.2)),
             sigma = c(psa3p = sqrt(0.10), s3m11p = sqrt(0.10)),
             rho = 0.8)
)

prior_parameters <- modifyList(prior_parameters,
                               lapply(prior_parameters[c("fpf", "tpf")], 
                                      \(x) list(Sigma = cor2cov(x$rho, x$sigma))))

# Plot distributions ----
layout(matrix(c(1, 2, 3, 0, 4, 5, 6, 7), 4, 2, TRUE))
par(pty = "m", mar = c(4, 1, 3, 1))
plot_univariate_prior_dist(prior_parameters$p.participation$mu, 
                           prior_parameters$p.participation$sigma, 
                           title = "(a) Participation probability")
plot_univariate_prior_dist(prior_parameters$p.mri$mu, 
                           prior_parameters$p.mri$sigma, 
                           title = "(b) MRI compliance probability")
plot_univariate_prior_dist(prior_parameters$p.isup2p$mu, 
                           prior_parameters$p.isup2p$sigma, 
                           title = "(c) ISUP≥2 vs ISUP=1 probability")
plot_univariate_prior_dist(prior_parameters$fpf$mu["psa3p"], 
                           prior_parameters$fpf$sigma["psa3p"], 
                           title = "(d) False Positive Fraction (ISUP=1) PSA≥3 ng/ml")
plot_univariate_prior_dist(prior_parameters$fpf$mu["s3m11p"], 
                           prior_parameters$fpf$sigma["s3m11p"], 
                           title = "(e) False Positive Fraction (ISUP=1) Stockholm3≥0.11")
plot_univariate_prior_dist(prior_parameters$tpf$mu["psa3p"], 
                           prior_parameters$tpf$sigma["psa3p"], 
                           title = "(f) True Positive Fraction (ISUP≥2) PSA≥3 ng/ml")
plot_univariate_prior_dist(prior_parameters$tpf$mu["s3m11p"], 
                           prior_parameters$tpf$sigma["s3m11p"], 
                           title = "(g) True Positive Fraction (ISUP≥2) Stockholm3≥0.11")
par(par0)

layout(matrix(c(1, 2), 1, 2, TRUE))
par(pty = "s")
plot_bivariate_prior_dist(prior_parameters$fpf$mu,
                          prior_parameters$fpf$Sigma,
                          title = "(h) False Positive Fraction (ISUP=1)",
                          lab = "FPF")
plot_bivariate_prior_dist(prior_parameters$tpf$mu,
                          prior_parameters$tpf$Sigma,
                          title = "(i) True Positive Fraction (ISUP≥2)",
                          lab = "TPF")
par(par0)

# Random draws from prior distributions ----
set.seed(1914)
prior_draws <- lapply(data_discordant, \(x) 
                      list(
                        p.participation = do.call("rlogitnormal", c(n = n.sim, prior_parameters$p.participation)),
                        p.mri = do.call("rlogitnormal", c(n = n.sim, prior_parameters$p.mri)),
                        p.isup2p = do.call("rlogitnormal", c(n = n.sim, prior_parameters$p.isup2p)),
                        fpf = do.call("rlogitmvnormal", c(n = n.sim, prior_parameters$fpf[c("mu", "Sigma")])),
                        tpf = do.call("rlogitmvnormal", c(n = n.sim, prior_parameters$tpf[c("mu", "Sigma")])) 
                      )
)

# Simulate outcome data ----
set.seed(1918)
data_simulated <- Map(\(d, n) lapply(seq_len(n.sim), \(i) 
                                     simulate_reinvited_data(d, prior_draws[[n]], i)),
                      data_discordant, names(data_discordant))
