# Title: Recursive spatio-temporal analysis for real temperature data  
# Script name: runme_recursiveTemperatureExample.R
# Author: Mario Figueira
# Date: 2025-11-02
# Description: Recursive inference with R-INLA for a real Big data problem.

# Last update: 2025-11-17

remove(list = ls())

# Loading libraries ----

library(Matrix)
library(Rcpp)
library(RcppEigen)
library(parallel)
library(rsm)

library(INLA)
library(INLABMA)
library(fmesher)
library(inlabru)

library(ggplot2)
library(GGally)
library(gridExtra)
library(ggtext)
library(dplyr)

library(sf)
library(geodata)
library(tidyverse)

## Custom functions ----

dens_MVN <- function(theta, mu, Q){
  return((2*pi)^(-nrow(Q)/2) * det(Q)^(1/2) * exp(-1/2 * t(theta - mu) %*% Q %*% (theta - mu)))
}

aug_ccd <- function(ccd_base, aug_coeff){
  require(dplyr)
  aug_coeff <- c(2, aug_coeff) %>% unique(.)
  ccd_res <- lapply(X = 1:length(aug_coeff), FUN = function(i){ccd_base*aug_coeff[i]}) %>% do.call(., what = rbind) %>% unique(.)
  return(ccd_res)
}

colsc <- function(...) {
  fill <- scale_fill_gradientn(
    colours = viridis::turbo(n=11),
    limits = range(..., na.rm = TRUE)
  )
  color <- scale_colour_gradientn(
    colours = viridis::turbo(n=11),
    limits = range(..., na.rm = TRUE)
  )
  return(list(fill=fill,color=color))
}

fix.Q <- function(Q){
  d <- diag(Q)
  Q <- Q + t(Q)
  diag(Q) <- d
  return (Q)
}

theta_muCov <- function(list_output){
  list_Q_hyper <- lapply(X = list_output$Cov, FUN = function(x){solve(x)})
  list_mu_hyper <- lapply(X = list_output$mode, FUN = function(x){x})
  
  Q_hyper <- Reduce(f = "+", x = list_Q_hyper)
  Cov_hyper <- solve(Q_hyper)
  
  mu_hyper <- Cov_hyper %*% (lapply(X = seq_len(length(list_Q_hyper)), FUN = function(idx){list_Q_hyper[[idx]] %*% list_mu_hyper[[idx]]}) %>% Reduce(f = "+", x = .))
  return(list(mu = mu_hyper, Cov = Cov_hyper))
}

theta_ccd <- function(list_output, n_cores = 1){
  list_Q_hyper <- lapply(X = list_output$Cov, FUN = function(x){solve(x)})
  list_mu_hyper <- lapply(X = list_output$mode, FUN = function(x){x})
  
  Q_hyper <- Reduce(f = "+", x = list_Q_hyper)
  Cov_hyper <- solve(Q_hyper)
  
  mu_hyper <- Cov_hyper %*% (lapply(X = seq_len(length(list_Q_hyper)), FUN = function(idx){list_Q_hyper[[idx]] %*% list_mu_hyper[[idx]]}) %>% Reduce(f = "+", x = .))
  
  size_theta <- ncol(Cov_hyper)
  ccd_list <- lapply(X = 1:size_theta, FUN = function(idx){as.formula(paste0("x", idx, "~", "theta", idx))})
  ccd.str <- ccd(basis = length(ccd_list), n0 = 1, alpha = "rotatable",
                 inscribed = TRUE, oneblock = FALSE,
                 coding = ccd_list)
  
  ccd.str <- unique(decode.data(ccd.str)[,1:size_theta+2])
  idx_central <- as.integer(which(apply(X = ccd.str, MARGIN = 1, FUN = function(x){all(x==rep(0,times = size_theta))})))
  ccd.str <- ccd.str[c(idx_central, setdiff(1:nrow(ccd.str), y = idx_central)),]
  
  ccd.ext <- aug_ccd(ccd_base = ccd.str, aug_coeff = c(2))
  # idx.mode <- apply(ccd.ext, MARGIN = 1, FUN = function(x){sum(x == rep(0, length(x))) == length(x)}) %>% which(.)
  Eig_dec <- eigen(Cov_hyper)
  df_theta <- mclapply(X = 1:nrow(ccd.ext), mc.cores = n_cores, FUN = function(i){
    res <- Eig_dec$vectors %*% diag(Eig_dec$values**(1/2)) %*% t(ccd.ext[i,]) + mu_hyper
    return(res)
  }) %>% do.call(., what = cbind) %>% t(.) %>% data.frame(.)
  return(df_theta)
}

## Loading data ----

DFsf_temp <- readRDS(file = "./Temperature_example/DF_temp.RDS")
mesh <- readRDS(file = "./Temperature_example/mesh.RDS")

## Exploring data ----

DF_temp <- data.frame(cbind(DFsf_temp %>% st_drop_geometry(.), DFsf_temp %>% st_coordinates(.)))

plot_func <- function(df, idx, title = NULL, name = "turbo") {
  ggplot(data = df[idx,], aes(x = X, y = Y, fill = max_sst)) +
    geom_tile() + labs(title = title) +
    scale_fill_viridis_c(option = name) + theme_minimal() +
    labs(fill = "Temperature") + xlab(label = "") + ylab(label = "") +
    theme(plot.title = element_text(size = 14, h = 0.5, face = "bold"))
}

list_plots <- lapply(X = 1:36, FUN = function(i){plot_func(df = DF_temp, idx = which(DF_temp$month.id %in% i), title = paste("Month",i), name = "turbo")})

grid.arrange(arrangeGrob(grobs = list_plots, ncol = 4))

DF_temp$month.id %>% max(.)
ntime <- 120 # Number of months to analyse. Change between 480 or 120 to reproduce the results presented in the paper
k_times <- 6 # Number of temporal partitions. Change between 24 and 6 partitions, for the 480 or 120 months respectively
DFsimplified_temp <- DF_temp[which(DF_temp$month.id %in% 1:ntime),]

# Defining the spatial componente for the separable spatio-temporal model (Q_ts = Q_t \otimes Q_s)
spde <- inla.spde2.pcmatern(
  mesh = mesh, alpha = 2,
  prior.range = c(mesh$loc[c(mesh$segm$bnd$idx[,1], mesh$segm$bnd$idx[1,1]),1:2] %>% dist(x = .) %>% as.vector(.) %>% max(.)/5, 0.5),
  prior.sigma = c(1,0.5), constr = FALSE
)

# Complete data analysis (full spatio-temporal model) ----

spde.idx <- inla.spde.make.index(name = "sp.idx", n.spde = spde$n.spde, n.group = ntime)
A_spt <- inla.spde.make.A(mesh = mesh, loc = DFsf_temp[which(DFsf_temp$month.id %in% 1:ntime),] %>% st_coordinates(.), group = DFsimplified_temp$month.id)

stk_total <- inla.stack(data = list(y = DFsimplified_temp$max_sst),
                        A = list(A_spt, 1),
                        effects = list(
                          spde.idx,
                          list(intercept = rep(1, times = nrow(DFsimplified_temp)))
                        ),
                        tag = "stk_inf_total")

formula_inla <- y ~ -1 + intercept + f(sp.idx, model = spde, group = sp.idx.group, control.group = list(model = "ar1"))
full_model_inla <- inla(formula = formula_inla, family = "gaussian", data = inla.stack.data(stk_total), control.predictor = list(A = inla.stack.A(stk_total)), control.compute = list(config = TRUE, waic = TRUE), verbose = FALSE)
t_cpu_full <- as.numeric(full_model_inla$cpu.used[4]/60)

# Recursive inference approach ----

t_cpu <- c()
list_mode_theta <- list()
list_Cov_theta <- list()
theta_mode <- rep(NA, times = 4)
## Compute the Gaussian approximation (proxy) to the actual Gaussian approximation of the whole data ----
for(i in 1:k_times){
  idx <- DFsimplified_temp$month.id %in% ((ntime/k_times*(i-1)+1):(ntime/k_times*i))
  spde.idx <- inla.spde.make.index(name = "sp.idx", n.spde = spde$n.spde, n.group = ntime/k_times)
  ntimes_group <- DFsimplified_temp$month.id[idx] %>% table(.) %>% as.vector(.)
  A_spt <- inla.spde.make.A(mesh = mesh, loc = DFsf_temp[which(idx),] %>% st_coordinates(.), group = rep(1:(ntime/k_times), times = ntimes_group))
  
  stk_inf_i <- inla.stack(data = list(y = DFsimplified_temp$max_sst[idx]),
                          A = list(A_spt, 1),
                          effects = list(
                            spde.idx,
                            list(intercept = rep(1, times = sum(idx)))
                          ),
                          tag = paste0("stk_inf_i"),
                          compress = TRUE,
                          remove.unused = TRUE)
  
  formula_inla_i <- y ~ -1 + intercept + f(sp.idx, model = spde, group = sp.idx.group, control.group = list(model = "ar1"))
  t_cpu[i] <- (system.time(
    model_inla <- inla(formula = formula_inla_i, family = "gaussian",
                       data = inla.stack.data(stk_inf_i),
                       only.hyperparam = TRUE,
                       control.mode = list(theta = theta_mode, fixed = FALSE, restart = TRUE),
                       control.predictor = list(A = inla.stack.A(stk_inf_i)),
                       control.compute = list(config = FALSE, waic = FALSE, return.marginals = FALSE),
                       verbose = FALSE)
  )[3]/60) %>% as.numeric(.) %>% round(., digits = 3)
  list_mode_theta[[i]] <- model_inla$misc$theta.mode
  list_Cov_theta[[i]] <- model_inla$misc$cov.intern
  theta_mode <- model_inla$mode$theta
  cat(paste0("Finalizada iteración ", i, ". Tiempo de ejecución: ", t_cpu[i], ". \n"))
}
hyper_Gapp <- theta_muCov(list_output = list(mode = list_mode_theta, Cov = list_Cov_theta))
hyper_theta_ccd <- theta_ccd(list_output = list(mode = list_mode_theta, Cov = list_Cov_theta), n_cores = 8)

list_gg_hyper <- list()
labs_title <- c("Log prec. likelihood", "Log-range spt.", "Log-sigma spt.", "Intern group rho spt.")
for(i in seq_along(hyper_Gapp$mu)){
  x_seq <- seq(from = hyper_Gapp$mu[i] - 5*sqrt(diag(hyper_Gapp$Cov)[i]), to = hyper_Gapp$mu[i] + 5*sqrt(diag(hyper_Gapp$Cov)[i]), length.out = 1E3)
  den_x <- dnorm(x = x_seq, mean = hyper_Gapp$mu[i], sd = sqrt(diag(hyper_Gapp$Cov)[i]))
  
  list_gg_hyper[[i]] <- ggplot() +
    geom_line(data = full_model_inla$internal.marginals.hyperpar[[i]], mapping = aes(x = x, y = y), color = "red") +
    geom_line(data = data.frame(x = x_seq, y = den_x), mapping = aes(x = x, y = y), color = "blue", linetype = "solid") +
    labs(title = labs_title[i]) + ylab(label = "density") +
    theme_bw() + theme(plot.title = element_text(size = 16, face = "bold", h = 0.5))
}
grid.arrange(arrangeGrob(grobs = list_gg_hyper, ncol = 4))

df_theta_CCD <- data.frame(rbind(as.matrix(hyper_theta_ccd), as.matrix(full_model_inla$joint.hyper[,1:ncol(hyper_theta_ccd)])), group_col = rep(c("blue", "red"), times = c(nrow(hyper_theta_ccd), nrow(full_model_inla$joint.hyper))))
colnames(df_theta_CCD)[1:ncol(hyper_theta_ccd)] <- labs_title
ggparis_hyper <- ggpairs(data = df_theta_CCD,
                         columns = 1:ncol(hyper_theta_ccd),
                         mapping = aes(colour = group_col), upper = list(continuous = "points"), diag = list(continuous = "blankDiag")) +
  scale_colour_manual(values = c('blue','red')) +
  theme_bw() + 
  theme(strip.text.x = element_text(size = 15, face = "bold"),
        strip.text.y = element_text(size = 15, face = "bold"))

diag_idx <- (0:3*4+1:4)
for(i in 1:4){
  ggparis_hyper$plots[[diag_idx[i]]] <- list_gg_hyper[[i]]
}

t0 <- Sys.time()
inf_condpost_ccd <- mclapply(X = seq_len(nrow(hyper_theta_ccd)), mc.cores = 8, FUN = function(idx_ccd){
  # Computing the mu and Q of the conditional Gaussian approximation
  res <- mclapply(X = 1:k_times, mc.cores = floor(k_times/1), FUN = function(i){
    idx <- DFsimplified_temp$month.id %in% ((ntime/k_times*(i-1)+1):(ntime/k_times*i))
    spde.idx <- inla.spde.make.index(name = "sp.idx", n.spde = spde$n.spde, n.group = ntime/k_times)
    ntimes_group <- DFsimplified_temp$month.id[idx] %>% table(.) %>% as.vector(.)
    A_spt <- inla.spde.make.A(mesh = mesh, loc = DFsf_temp[which(idx),] %>% st_coordinates(.), group = rep(1:(ntime/k_times), times = ntimes_group))
    
    stk_inf_i <- inla.stack(data = list(y = DFsimplified_temp$max_sst[idx]),
                            A = list(A_spt, 1),
                            effects = list(
                              spde.idx,
                              list(intercept = rep(1, times = sum(idx)))
                            ),
                            tag = paste0("stk_inf_i"),
                            compress = TRUE,
                            remove.unused = TRUE)
    
    formula_inla_i <- y ~ -1 + intercept + f(sp.idx, model = spde, group = sp.idx.group, control.group = list(model = "ar1"))
    condpost_inla <- inla(formula = formula_inla_i, family = "gaussian",
                          data = inla.stack.data(stk_inf_i),
                          only.hyperparam = FALSE,
                          control.mode = list(theta = as.numeric(hyper_theta_ccd[idx_ccd,]), fixed = TRUE),
                          control.predictor = list(A = inla.stack.A(stk_inf_i)),
                          control.compute = list(config = TRUE, return.marginals = FALSE),
                          verbose = FALSE)
    
    return(list(Q = condpost_inla$misc$configs$config[[1]]$Q, mu = condpost_inla$misc$configs$config[[1]]$improved.mean, mlik = condpost_inla$mlik[1]))
  })

  # Computing the posterior precision matrix and mean of the ensable of partition for the idx_ccd-th support point
  list_Q <- lapply(X = seq_along(res), FUN = function(i){res[[i]]$Q})
  list_mu <- lapply(X = seq_along(res), FUN = function(i){res[[i]]$mu})
  
  list_Qprodmu <- lapply(X = seq_along(res), FUN = function(i){
    Qpm <- drop(fix.Q(res[[i]]$Q) %*% res[[i]]$mu)
    ext_spt <- (length(res[[i]]$mu)-1)
    vec_mu <- rep(0, times = ext_spt*length(res)+1)
    vec_mu[length(vec_mu)] <- last(Qpm)
    vec_mu[((i-1)*ext_spt+1):(i*ext_spt)] <- Qpm[1:ext_spt]
    return(matrix(data = vec_mu, ncol = 1))
  })
  
  list_Qensamb <- lapply(X = list_Q, FUN = function(x){
    x <- as(x, "TsparseMatrix")
    idx.b0 <- (x@j==(ncol(x)-1))
    Qst <- sparseMatrix(i = x@i[!idx.b0]+1, j = x@j[!idx.b0]+1, x = x@x[!idx.b0], dims = c(ncol(x)-1, ncol(x)-1))
    Qb0 <- sparseMatrix(i = x@i[idx.b0]+1, j = rep(1, sum(idx.b0)), x = x@x[idx.b0], dims = c(ncol(x), 1))
    return(list(Qst = Qst, Qb0 = Qb0))
  })
  
  Qst <- Matrix::bdiag(lapply(X = list_Qensamb, FUN = function(x){x$Qst})) %>% as(., "TsparseMatrix")
  Qst <- sparseMatrix(i = Qst@i+1, j = Qst@j+1, x = Qst@x, dims = c(nrow(Qst)+1, ncol(Qst)))
  
  idx_i <- c()
  idx_x <- c()
  diag_b0 <- c()
  max_row <- c()
  for(i in seq_len(length(list_Qensamb))){
    idx_i0 <- c(list_Qensamb[[i]]$Qb0@i)
    
    idx_x <- c(idx_x, list_Qensamb[[i]]$Qb0@x[-which.max(idx_i0)])
    diag_b0 <- c(diag_b0, list_Qensamb[[i]]$Qb0@x[which.max(idx_i0)])
    
    idx_i0 <- idx_i0[-which.max(idx_i0)]
    
    idx_i <- c(idx_i, idx_i0 + sum(max_row))
    
    max_row <- c(max_row, nrow(list_Qensamb[[i]]$Qb0)-1)
  }
  Qb0 <- sparseMatrix(i = c(idx_i+1, nrow(Qst)), j = rep(1, length(c(idx_i+1, nrow(Qst)))), x = c(idx_x, sum(diag_b0)), dims = c(nrow(Qst), 1))
  Qpost <- fix.Q(cbind(Qst, Qb0))
  invQpost <- inla.qinv(Q = Qpost)
  
  Qpmu_post <- Reduce(x = list_Qprodmu, f = "+")
  mu <- inla.qsolve(Qpost, Qpmu_post, method = "solve")
  
  wsi <- lapply(X = seq_along(res), FUN = function(i){log(res[[i]]$mlik)}) %>% Reduce(x = ., f = "+") %>% exp(.)
  
  return(list(mu = mu, invQdiag = sqrt(diag(invQpost)), wsi = wsi))
})
t1 <- Sys.time()
difftime(t1,t0)

t0 <- Sys.time()
f0 <- 2
delta_w <- ((nrow(hyper_theta_ccd)-1)*(f0^2-1)*(1+exp(-ncol(hyper_theta_ccd)*f0^2/2)))**(-1)
delta_w_central <- 1 - (nrow(hyper_theta_ccd)-1)*delta_w  

wsi <- lapply(X = seq_along(inf_condpost_ccd), FUN = function(i){inf_condpost_ccd[[i]]$wsi}) %>% unlist(.) # One approach to compute the weights through the marginal likelihood of the conditional distributions (non-normalized marginal posterior of the hyperparameters)
# wsi <- apply(X = hyper_theta_ccd, MARGIN = 1, FUN = function(x){dens_MVN(theta = x, mu = drop(hyper_Gapp$mu), solve(hyper_Gapp$Cov))}) # Alternative approach way to compute the weights through the marginal likelihood of the conditional distributions (non-normalized marginal posterior of the hyperparameters)
ws <- wsi*rep(c(delta_w_central, delta_w), time = c(1, nrow(hyper_theta_ccd)-1)) # If we want to normalized it we must add (but there is no need for): "/sum(wsi*rep(c(delta_w_central, delta_w), time = c(1, nrow(hyper_theta_ccd)-1)))"

df_post_musd <- lapply(X = seq_along(inf_condpost_ccd), FUN = function(x){
  df <- data.frame(mean = inf_condpost_ccd[[x]]$mu, sd = inf_condpost_ccd[[x]]$invQdiag)
  colnames(df) <- paste0(colnames(df), rep(x, times = 2))
  return(df)
}) %>% do.call(., what = cbind)

mu_post <- apply(X = df_post_musd[,seq_len(ncol(df_post_musd)/2)*2-1], MARGIN = 1, FUN = function(x){sum(ws*x)/sum(ws)})
sum_wmu2 <- apply(X = df_post_musd[,seq_len(ncol(df_post_musd)/2)*2-1], MARGIN = 1, FUN = function(x){sum(ws*x^2)/sum(ws)})
diff_muwmu2 <- sum_wmu2 - mu_post^2
sd_post <- sqrt(((apply(X = df_post_musd[,seq_len(ncol(df_post_musd)/2)*2], MARGIN = 1, FUN = function(x){sum(ws*x^2)/sum(ws)})) + diff_muwmu2))

list_marg <- mclapply(X = seq_along(mu_post), mc.cores = 80, FUN = function(i_lf){
  x_seq <- sort(unique(c(seq(from = mu_post[i_lf] - 5*sd_post[i_lf], to = mu_post[i_lf] + 5*sd_post[i_lf], length.out = 1.5E1), qnorm(p = seq(0.01, 0.99, length.out = 2.5E1), mean = mu_post[i_lf], sd = sd_post[i_lf]))))
  
  log_wpixi_post <- lapply(X = seq_len(ncol(df_post_musd)/2), FUN = function(i_ccd){log(ws[i_ccd]) + dnorm(x_seq, mean = df_post_musd[i_lf,i_ccd*2-1], sd = df_post_musd[i_lf,i_ccd*2], log = TRUE)}) %>% do.call(., what = cbind)
  max_lwpixi <- max(log_wpixi_post)
  log_pix <- max_lwpixi + log(x = apply(X = log_wpixi_post, MARGIN = 1, FUN = function(x){sum(exp(x - max_lwpixi))}))
  
  spline_logf <- splinefun(x_seq, log_pix, method = "natural")
  lw_bound <- min(x_seq); up_bound <- max(x_seq)
  g <- function(xx){exp(spline_logf(xx) - mlog_pix)}
  int_nn_pix <- integrate(f = g, lower = lw_bound, upper = up_bound, rel.tol = 1E-8, abs.tol = 0)$value
  norm_den_pix <- function(xx){
    exp(spline_logf(xx) - mlog_pix) / int_nn_pix
  }
  
  df_pix <- data.frame(x = x_seq, y = norm_den_pix(x_seq))
  return(df_pix)
})

summary_lf <- mclapply(X = list_marg, mc.cores = 80, FUN = function(x){inla.zmarginal(x, silent = TRUE) %>% as.data.frame(.)}) %>% do.call(., what = rbind) %>% as.data.frame(.)
t1 <- Sys.time()
difftime(t1,t0)

xy_grid <- st_sample(x = st_polygon(x = list(mesh$loc[c(mesh$segm$int$idx[,1], mesh$segm$int$idx[1,1]),1:2])), size = 1E4, type = "regular")
A_grid <- fm_basis(x = mesh, loc = xy_grid)

list_ggspt <- list()
for(k in 1:4){
  u_d <- A_grid %*% summary_lf$mean[(1 + (k-1)*mesh$n):(k*mesh$n)]
  u_s <- A_grid %*% full_model_inla$summary.random$sp.idx$mean[(1 + (k-1)*mesh$n):(k*mesh$n)]
  ggs1_mean <- ggplot() + 
    geom_tile(data = data.frame(st_coordinates(xy_grid), mean = drop(u_s), type = "Standard (s)"), mapping = aes(x = X, y = Y, fill = mean)) +
    geom_tile(data = data.frame(st_coordinates(xy_grid), mean = drop(u_d), type = "Distributed (d)"), mapping = aes(x = X, y = Y, fill = mean)) +
    scale_fill_viridis_c(option = "turbo") + 
    facet_wrap(facets = . ~ type, ncol = 2) +
    theme_bw() + theme(strip.text.x = element_text(size = 14, face = "bold"))
  ggs1_meandif <- ggplot() + 
    geom_tile(data = data.frame(st_coordinates(xy_grid), mean = drop(u_d-u_s), type = "Difference (d-s)"), mapping = aes(x = X, y = Y, fill = mean)) +
    scale_fill_viridis_c(option = "turbo") + 
    facet_wrap(facets = . ~ type, ncol = 2) +
    theme_bw() + theme(strip.text.x = element_text(size = 14, face = "bold"))
  
  su_d <- A_grid %*% summary_lf$sd[(1 + (k-1)*mesh$n):(k*mesh$n)]
  su_s <- A_grid %*% full_model_inla$summary.random$sp.idx$sd[(1 + (k-1)*mesh$n):(k*mesh$n)]
  ggs1_sd1 <- ggplot() + 
    geom_tile(data = data.frame(st_coordinates(xy_grid), sd = drop(su_d), type = "Distributed (d)"), mapping = aes(x = X, y = Y, fill = sd)) +
    scale_fill_viridis_c(option = "turbo") + 
    facet_wrap(facets = . ~ type, ncol = 2) +
    theme_bw() + theme(strip.text.x = element_text(size = 14, face = "bold"))
  ggs1_sd2 <- ggplot() + 
    geom_tile(data = data.frame(st_coordinates(xy_grid), sd = drop(su_s), type = "Standard (s)"), mapping = aes(x = X, y = Y, fill = sd)) +
    scale_fill_viridis_c(option = "turbo") + 
    facet_wrap(facets = . ~ type, ncol = 2) +
    theme_bw() + theme(strip.text.x = element_text(size = 14, face = "bold"))
  ggs1_sddif <- ggplot() + 
    geom_tile(data = data.frame(st_coordinates(xy_grid), sd = drop(su_d - su_s), type = "Difference (d-s)"), mapping = aes(x = X, y = Y, fill = sd)) +
    scale_fill_viridis_c(option = "turbo") + 
    facet_wrap(facets = . ~ type, ncol = 2) +
    theme_bw() + theme(strip.text.x = element_text(size = 14, face = "bold"))
  
  list_ggspt[[k]] <- arrangeGrob(grobs = list(ggs1_mean, ggs1_meandif, ggs1_sd1, ggs1_sd2, ggs1_sddif), layout_matrix = matrix(data = c(1,1,2,3,4,5), ncol = 3, byrow = TRUE))
}

grid.arrange(list_ggspt[[4]])

# Graphical representation of the marginal densities ----

gg_int <- ggplot() + 
  geom_line(data = full_model_inla$marginals.fixed$intercept, mapping = aes(x = x, y = y), color = "red") +
  geom_line(data = list_marg[[length(list_marg)]], mapping = aes(x = x, y = y), color = "blue") +
  labs(title = "Intercept") + ylab(label = "density") +
  theme_bw() + theme(plot.title = element_text(face = "bold", size = 16, h = 0.5))

idx_samp <- sort(sample(x = 1:nrow(summary_lf)-1, size = 47, replace = FALSE))
gg_spt <- lapply(X = idx_samp, FUN = function(id){
  idx_mesh <- id %% mesh$n
  idx_temp <- (id-1) %/% mesh$n + 1
  if(!all(idx_mesh != 0)){idx_mesh[which(idx_mesh==0)] <- mesh$n}
  
  gg <- ggplot() + 
    geom_line(data = full_model_inla$marginals.random$sp.idx[[id]], mapping = aes(x = x, y = y), color = "red") +
    geom_line(data = list_marg[[id]], mapping = aes(x = x, y = y), color = "blue") +
    labs(title = paste0("Spt. eff. (s = ", idx_mesh, ", t = ", idx_temp,")")) + ylab(label = "density") +
    theme_bw() + theme(plot.title = element_text(face = "bold", size = 16, h = 0.5))
})

grid.arrange(arrangeGrob(grobs = c(list(gg_int), gg_spt), ncol = 6))

plot(list_marg[[length(list_marg)]], type = "l", col = "blue")
lines(full_model_inla$marginals.fixed$intercept, col = "red")




