#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# STHLM3MRI Reinvite: Second screening round - Monte Carlo Sensitivity Analysis
# Functions
# Andrea Discacciati
# 20240816
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#

# Functions for logit-normal and logit-bivariate-normal RVs
dlogitnormal <- function(x, mu, sigma) { # x (support): (0, 1), mu: R
  dnorm(qlogis(x), mu, sigma) * (x-x^2)^-1
}

rlogitnormal <- function(n, mu, sigma) { # mu: R
  plogis(rnorm(n, mu, sigma))
}

qlogitnormal <- function(q, mu, sigma){ # mu: R
  plogis(qnorm(q, mu, sigma))
}

rlogitmvnormal <- function(n, mu, Sigma) {
  plogis(MASS::mvrnorm(n, mu, Sigma))
}

dlogitmvnormal <- function(x, y, mu, Sigma) {
  sigma.x <- sqrt(Sigma[1,1])
  sigma.y <- sqrt(Sigma[2,2])
  rho     <- cov2cor(Sigma)[1,2]
  
  mu.x <- mu[1]
  mu.y <- mu[2]
  
  Q <- (sigma.x * sigma.y * (1-rho^2))^(-1) * 
    ((qlogis(x)-mu.x)^2 + (qlogis(y)-mu.y)^2 - 2*rho*(qlogis(x)-mu.x)*(qlogis(y)-mu.y))
  
  f <- exp(-Q/2)/(2 * pi * sigma.x * sigma.y * sqrt(1-rho) * x * y * (1-x) * (1-y))
  f
}

# Correlation -> Covariance matrix
cor2cov <- function(rho, sigma) {
  R <- matrix(c(1, rho, rho, 1), ncol = 2)
  diag(sigma) %*% R %*% diag(sigma)
}

# Plot univariate logit-normal distributions
plot_univariate_prior_dist <- function(mu, sigma, title = "", pl.digits = 2) {
  format_pl <- function(x, digits = 2, text) {
    x <- list(text, x[1], x[2])
    fmt <- gsub("X", digits, "%s %.Xf-%.Xf")
    do.call("sprintf", c(fmt = fmt, x))
  }
  
  myseq <- seq(0.001, 0.999, by = 0.001)
  dln <- dlogitnormal(myseq, mu, sigma)
  dln <- dln/max(dln) # rescale max density to 1 to simplify plotting
  
  # Area between quantiles
  quant <- c(0.005,
             0.995,
             0.025,
             0.975,
             0.100,
             0.900)
  probss <- sapply(quant, \(q) qlogitnormal(q, mu, sigma))
  names(probss) <- quant
  
  # Lower and higher indices on the X-axis
  index <- rep(NA, length(quant))
  index[c(1, 3, 5)] <- sapply(probss[c(1, 3, 5)], \(p) min(which(myseq >= p)))
  index[c(2, 4, 6)] <- sapply(probss[c(2, 4, 6)], \(p) max(which(myseq <= p)))
  names(index) <- quant
  
  fillcolors <- rev(RColorBrewer::brewer.pal(6, "Blues"))
  
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
       ylab = "", xlab = "")
  axis(1, at = seq(0, 1, 0.05))
  for (i in c(1, 3, 5)) {
    polygon(myseq[c(index[i], index[i]:index[i+1], index[i+1])],
            c(0, dln[index[i]:index[i+1]], 0),
            col = fillcolors[i],
            border = NA)
  }
  lines(myseq, dln, lwd = 2, type = "l")
  
  mu.p <- integrate(\(x) x * dlogitnormal(x, mu, sigma), 0, 1)$value # mean mu.p on transformed (0,1) space
  map <- myseq[which.max(dln)] # mode, assuming a unimodal logit-normal distr.
  q50 <- unname(plogis(mu)) # median
  
  #arrows(mu.p, 0.1, mu.p, 0, length = 0.1, angle = 30, lwd = 1.2)
  #arrows(q50, 0.1, q50, 0, length = 0.1, angle = 30, lwd = 2)
  arrows(map, 0.1, map, 0, length = 0.1, angle = 30, lwd = 2)
  
  text(x = 0+1*(plogis(mu)<0.5), y = 0.8, 
       label = paste0(format_pl(c(probss[1], probss[2]), pl.digits, "99% limits: "), "\n",
                      format_pl(c(probss[3], probss[4]), pl.digits, "95% limits: "), "\n",
                      format_pl(c(probss[5], probss[6]), pl.digits, "80% limits: ")), 
       pos = 4-2*(plogis(mu)<0.5))
  title(main = title, cex = 0.95)
  
  probss <- c(probss, `0.5` = q50, mu = mu.p, map = map)
  probss
}

# Plot bivariate logit-normal distributions
plot_bivariate_prior_dist <- function(mu, Sigma, title = "", lab = "") {
  g <- seq(0.001, 0.999, by = 0.005)
  z <- outer(g, g, FUN = "dlogitmvnormal", mu = mu, Sigma = Sigma)
  z <- z/max(z)
  
  levels <- sort(c(10^seq(-12, 0, by = 2), 5e-01, 1e-01))
  
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), 
       axes = FALSE,
       xlab = paste(lab, "for PSA ≥3 ng/ml"), 
       ylab = paste(lab, "for Stockholm3 ≥0.11"),
       xaxs = "i", yaxs = "i")
  .filled.contour(g, g, z,
                  levels = levels,
                  col = hcl.colors(length(levels)-1, "Heat", rev = TRUE))
  contour(g, g, z,
          levels = levels,
          lwd = 0.5,
          drawlabels = FALSE,
          add = TRUE)
  axis(1, at = seq(0, 1, 0.1))
  axis(2, at = seq(0, 1, 0.1), las = 1)
  box(lwd = 1)
  title(main = title, cex = 0.95)
  
  ind.max.density <- which(z == max(z), arr.ind = TRUE)
  points(g[ind.max.density[1]], g[ind.max.density[2]], pch = 16)
  abline(a = 0, b = 1)
}


simulate_reinvited_data <- function(data.discordant, prior.draws, i) {
  #cat(i, "\n")
  df0 <- expand.grid(bl_psa3p = c(0, 1),
                     bl_s3m11p = c(0, 1))
  
  n.discordant <- nrow(data.discordant)
  n.isup2p <- rbinom(1, n.discordant, prior.draws$p.isup2p)
  n.isup1  <- n.discordant - n.isup2p
  n.participated <- rbinom(1, n.discordant, prior.draws$p.participation[i])
  n.notparticipated <- n.discordant - n.participated
  n.mri <- rbinom(1, n.discordant, prior.draws$p.mri[i])
  n.notmri <- n.discordant - n.mri
  
  make.sampling.probs <- function(prior.draws, i) {
    #                 STOCKHOLM3
    # PSA         Negative    Positive 
    # Negative    p_0_0       p_0_1    1-TPF PSA  
    # Positive    p_1_0       TPPF     TPF_PSA
    #             1-TPF_S3    TPF_S3   1
    
    # unconditional probabilities       
    probs.uncond <- rep(0, 4)
    names(probs.uncond) <- c("p.psa3p-.s3m11p-", 
                             "p.psa3p+.s3m11p-",
                             "p.psa3p-.s3m11p+",
                             "p.psa3p+.s3m11p+")
    
    rtbeta <- function(n, a, b, LB, UB) { # random generation from truncated beta
      qbeta(runif(n, pbeta(LB, a, b), pbeta(UB, a, b)), a, b)
    }
    LB <- max(0, prior.draws[[i, "psa3p"]]+prior.draws[[i, "s3m11p"]]-1) # LB TPPR Alonzo, see also Agresti 3.5.1
    UB <- min(prior.draws[i, ]) # UB TPPR Alonzo, see also Agresti 3.5.1
    probs.uncond["p.psa3p+.s3m11p+"] <- rtbeta(1, 8, 4, LB, UB) # alt1: random n from trunc beta [LB, UB]
    probs.uncond["p.psa3p+.s3m11p-"] <- prior.draws[[i, "psa3p"]]-probs.uncond["p.psa3p+.s3m11p+"]
    probs.uncond["p.psa3p-.s3m11p+"] <- prior.draws[[i, "s3m11p"]]-probs.uncond["p.psa3p+.s3m11p+"]
    probs.uncond["p.psa3p-.s3m11p-"] <- -sum(probs.uncond)+1
    
    # conditional probabilities given test+ (either or both)
    probs.test.pos <- probs.uncond[c("p.psa3p+.s3m11p-",
                                     "p.psa3p-.s3m11p+",
                                     "p.psa3p+.s3m11p+")]
    probs.cond <- probs.test.pos / sum(probs.test.pos)
    list(probs.uncond = probs.uncond,
         probs.cond   = probs.cond)
  }
  
  probs.isup1 <- make.sampling.probs(prior.draws$fpf, i)
  counts.isup1 <- rmultinom(1, n.isup1, probs.isup1$probs.cond)
  
  probs.isup2p <- make.sampling.probs(prior.draws$tpf, i)
  counts.isup2p <- rmultinom(1, n.isup2p, probs.isup2p$probs.cond)
  
  uncounts.isup1 <- tidyr::uncount(df0, c(0, counts.isup1))
  if (nrow(uncounts.isup1)>0) { # prevents error if no ISUP1
    uncounts.isup1$cb_pca_isup1  <- 1
    uncounts.isup1$cb_pca_isup2p <- 0
  }
  uncounts.isup2p <- tidyr::uncount(df0, c(0, counts.isup2p))
  uncounts.isup2p$cb_pca_isup1  <- 0
  uncounts.isup2p$cb_pca_isup2p <- 1
  
  uncounts <- cbind(
    data.discordant[sample(n.discordant), "studieid", drop = FALSE],
    rbind(
      uncounts.isup1,
      uncounts.isup2p
    )
  )
  uncounts$participated <- rep(c(0, 1), c(n.notparticipated, n.participated))[sample(n.discordant)]
  uncounts$bl_psa3p <- ifelse(uncounts$participated == 0, 0, uncounts$bl_psa3p)
  uncounts$bl_s3m11p <- ifelse(uncounts$participated == 0, 0, uncounts$bl_s3m11p)
  uncounts$mr_yn <- rep(c(0, 1), c(n.notmri, n.mri))[sample(n.discordant)]
  uncounts$mr_yn <- ifelse(uncounts$participated == 0, 0, uncounts$mr_yn)
  uncounts$cb_yn <- uncounts$mr_yn # 100% biopsy compliance ifelse(uncounts$participated == 0 | uncounts$mr_yn == 0, 0, 1)
  uncounts$cb_pca_isup0  <- 0
  uncounts$cb_pca_isup1 <- ifelse(uncounts$participated == 0 | uncounts$mr_yn == 0, 0, uncounts$cb_pca_isup1)
  uncounts$cb_pca_isup2p <- ifelse(uncounts$participated == 0 | uncounts$mr_yn == 0, 0, uncounts$cb_pca_isup2p)
  uncounts$simulated <- 1
  
  uncounts <- uncounts[, c("studieid", "participated", "bl_psa3p", "bl_s3m11p", 
                           "mr_yn", "cb_yn", "cb_pca_isup0", "cb_pca_isup1", "cb_pca_isup2p",
                           "simulated")]
  
  make.matrix.bloodtest.probs <- function(obj) {
    lapply(obj, \(x) {
      if (length(x) == 3) x <- c(0, x) 
      addmargins(matrix(x, ncol = 2, 
                        dimnames = list(PSA3  = c("-", "+"), 
                                        S3M11 = c("-", "+"))))
    })
  }
  
  # dump stuff
  return(
    list(data = uncounts,
         probs.isup1  = make.matrix.bloodtest.probs(probs.isup1),
         probs.isup2p = make.matrix.bloodtest.probs(probs.isup2p)
    )
  )
}
