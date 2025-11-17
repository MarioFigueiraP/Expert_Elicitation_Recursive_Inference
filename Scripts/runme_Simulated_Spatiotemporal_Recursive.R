# Spatio-temporal simulated scenario to implement a recursive inference approach.

# Loading libraries ----

library(INLA)
library(inlabru)
library(fmesher)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(ggtext)
library(viridis)
library(sf)
library(parallel)

# Some custom functions ----

## Function to armonaised color scales
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

## Funtion to fix the precision matrix
fix.Q <- function(Q){
  d <- diag(Q)
  Q <- Q + t(Q)
  diag(Q) <- d
  return (Q)
}

# First-order random walk 
rw1_sim <- function(seed, n, sd, constr = TRUE){
  if(!missing(seed)){set.seed(seed = seed)}
  x <- rep(0, times = n)
  for(i in 1:n){
    if(i==1){
      x[1] <- rnorm(1, mean = 0, sd = sd)
    } else{
      x[i] <- x[i-1] + rnorm(1, mean = 0, sd = sd) 
    }
  }
  if(constr){
    x <- x - mean(x)
  }
  return(x)
}

# Second-order random walk 
rw2_sim <- function(seed, n, sd, constr = TRUE){
  if(!missing(seed)){set.seed(seed = seed)}
  x <- rep(0, times = n)
  for(i in 1:n){
    if(i==1){
      x[1] <- rnorm(1, mean = 0, sd = sd)
    } else if(i==2){
      x[2] <- x[1] + rnorm(1, mean = 0, sd = sd)
    } else{
      x[i] <- 2*x[i-1] - x[i-2] + rnorm(1, mean = 0, sd = sd) 
    }
  }
  if(constr){
    x <- x - mean(x)
  }
  return(x)
}

# Simulation of the spatio-temporal data ----
global.seed <- 12
set.seed(seed = global.seed)

rho_spatial <- 0.4 # spatial range
sigma_spatial <- 2 # marginal stdev. of the spatial effect
mesh <- fm_mesh_2d_inla(loc.domain=matrix(c(0,0,0,1,1,1,1,0), byrow=TRUE, ncol=2), max.edge=c(0.025,0.15)) # mesh creation
Q <- fm_matern_precision(x=mesh, alpha=2, rho=rho_spatial, sigma=sigma_spatial) # matrix SPDE-FEM
u_s <- fm_sample(n=1, Q=Q, mu=0, constr=NULL) # spatial effect
u_s <- u_s - mean(u_s)

u_t <- rw1_sim(seed = global.seed, n = 60, sd = (20)**(-1/2), constr = TRUE)
plot(1:length(u_t), u_t, type = "l")
u_tt <- rep(u_t, each = 3E2) 

xy <- cbind(runif(60*3E2), runif(60*3E2))
A_xy <- fm_basis(x = mesh, loc = xy)
u_ss <- drop(A_xy %*% u_s)

beta0 <- 2
y_sim <- rnorm(n = length(u_ss), mean = beta0 + u_ss + u_tt, sd = (20)**(-1/2))

# Full data analysis ----
mesh_inf <- fm_mesh_2d_inla(loc.domain=matrix(c(0,0,0,1,1,1,1,0), byrow=TRUE, ncol=2), max.edge=c(0.04,0.15))
spde <- inla.spde2.pcmatern(mesh = mesh_inf, alpha = 2, prior.range = c(0.2,0.5), prior.sigma = c(1,0.5), constr = TRUE)
A_inf <- fm_basis(x = mesh_inf, loc = xy)

inf_stk_full <- inla.stack(data = list(y = y_sim),
                           A = list(A_inf, 1),
                           effects = list(
                             list(sp = 1:mesh_inf$n),
                             list(beta_0 = rep(1, length.out = length(y_sim)),
                                  t = rep(1:60, each = 3E2)
                                  )
                             ),
                           tag = "inf_full")

formula_full <- y ~ -1 + beta_0 + f(t, model = "rw1", constr = TRUE) + f(sp, model = spde)
full_model <- inla(data = inla.stack.data(inf_stk_full), family = "gaussian",
                   formula = formula_full,
                   control.compute = list(config = TRUE),
                   control.predictor = list(A = inla.stack.A(inf_stk_full)))

ggplot() + 
  geom_line(data = data.frame(full_model$summary.random$t), mapping = aes(x = ID, y = mean), color = "red") +
  geom_line(data = data.frame(ID = 1:length(u_t), mean = u_t), mapping = aes(x = ID, y = mean)) + 
  theme_bw()

# Sequential data analysis ----

size_seq <- length(y_sim)/3

k_seq <- 1
y_sim_seq1 <- rep(NA, length(y_sim))
y_sim_seq1[(size_seq*(k_seq-1)+1):(size_seq*k_seq)] <- y_sim[(size_seq*(k_seq-1)+1):(size_seq*k_seq)]
inf_stk_s1 <- inla.stack(data = list(y = y_sim_seq1),
                         A = list(A_inf, 1),
                         effects = list(
                           list(sp = 1:mesh_inf$n),
                           list(beta_0 = rep(1, length.out = length(y_sim)),
                                t = rep(1:60, each = 3E2)
                                )
                           ),
                         compress = TRUE,
                         remove.unused = FALSE,
                         tag = "inf_s1")

formula_seq1 <- y ~ -1 + beta_0 + f(t, model = "rw1", constr = TRUE) + f(sp, model = spde)
seq1_model <- inla(data = inla.stack.data(inf_stk_s1), family = "gaussian",
                   formula = formula_seq1,
                   control.predictor = list(A = inla.stack.A(inf_stk_s1)),
                   control.compute = list(config = TRUE))

mean_gmrf_1 <- seq1_model$misc$configs$config[[1]]$improved.mean
Q_gmrf_1 <- fix.Q(Q = seq1_model$misc$configs$config[[1]]$Q)
Sigma_gmrf_1 <- sqrt(diag(seq1_model$misc$configs$config[[1]]$Qinv))

df_seq1i <- mclapply(X = 1:length(u_t), mc.cores = 4, FUN = function(i){
  data.frame(ID = i, mean = seq1_model$misc$configs$config[[1]]$improved.mean[i], 
             X0.025quant = qnorm(p = 0.025, mean = seq1_model$misc$configs$config[[1]]$improved.mean[i], sd = Sigma_gmrf_1[i]),
             X0.975quant = qnorm(p = 0.975, mean = seq1_model$misc$configs$config[[1]]$improved.mean[i], sd = Sigma_gmrf_1[i])
  )
}) %>% do.call(., what = rbind) %>% as.data.frame(.)

gg_seq1 <- 
  ggplot() + 
  geom_line(data = df_seq1i, mapping = aes(x = ID, y = mean), color = "red") +
  geom_ribbon(data = df_seq1i, mapping = aes(x = ID, ymin = X0.025quant, ymax = X0.975quant), fill = "red", alpha = 0.25) +
  # geom_line(data = data.frame(seq1_model$summary.random$t), mapping = aes(x = ID, y = mean), color = "red") +
  # geom_ribbon(data = data.frame(seq1_model$summary.random$t), mapping = aes(x = ID, ymin = X0.025quant, ymax = X0.975quant), fill = "red", alpha = 0.25) +
  # geom_line(data = data.frame(full_model$summary.random$t), mapping = aes(x = ID, y = mean), color = "blue") + 
  # geom_ribbon(data = data.frame(full_model$summary.random$t), mapping = aes(x = ID, ymin = X0.025quant, ymax = X0.975quant), fill = "blue", alpha = 0.25) +
  geom_line(data = data.frame(ID = 1:length(u_t), mean = u_t), mapping = aes(x = ID, y = mean)) + 
  labs(title = "A. Posterior (first step)") + theme_bw() +
  theme(plot.title = element_markdown(hjust=0, face="bold", size=16, color = "black"))

res_seq <- mclapply(X = 1:nrow(seq1_model$joint.hyper), mc.cores = 4, FUN = function(i){
# res_seq <- mclapply(X = 1, mc.cores = 4, FUN = function(i){
  k_seq <- 2
  y_sim_seq2 <- rep(NA, length(y_sim))
  y_sim_seq2[(size_seq*(k_seq-1)+1):(size_seq*k_seq)] <- y_sim[(size_seq*(k_seq-1)+1):(size_seq*k_seq)]
  
  log_dens <- seq1_model$joint.hyper[i,ncol(seq1_model$joint.hyper)-1]
  
  mean_gmrf_1 <- seq1_model$misc$configs$config[[i]]$improved.mean
  Q_gmrf_1 <- fix.Q(Q = seq1_model$misc$configs$config[[i]]$Q)
  
  A_inf_sp_seq <- fm_basis(x = mesh_inf, loc = xy)
  
  idx_t <- 1:60
  idx_tt <- rep(idx_t, each = 3E2)
  A_inf_t_seq <- lapply(X = 1:length(y_sim), FUN = function(i){x <- rep(0, times = length(idx_t)); x[idx_tt[i]] <- 1; return(x)}) %>% do.call(., what = rbind)
  
  A_seq <- cbind(A_inf_t_seq, A_inf_sp_seq, 1)
  
  formula_seq2 <- y ~ -1 + offset(mean_gmrf_1) + f(gmrf, model = "generic0", Cmatrix = Q_gmrf_1, extraconstr = list(A = seq1_model$misc$configs$constr$A, e = seq1_model$misc$configs$constr$e), hyper = list(prec = list(initial = log(1), fixed = TRUE)))
  seq2_model <- inla(data = list(y = y_sim_seq2, mean_grmf_1 = mean_gmrf_1, gmrf = 1:nrow(Q_gmrf_1)), family = "gaussian",
                     formula = formula_seq2,
                     control.predictor = list(A = A_seq),
                     control.compute = list(config = TRUE),
                     control.family = list(hyper = list(prec = list(initial = seq1_model$joint.hyper[1,1], fixed = TRUE))),
                     verbose = FALSE)
  
  mean_gmrf_2 <- mean_gmrf_1 + seq2_model$misc$configs$config[[1]]$improved.mean
  Q_gmrf_2 <- fix.Q(Q = seq2_model$misc$configs$config[[1]]$Q)
  Sigma_gmrf_2 <- sqrt(diag(seq2_model$misc$configs$config[[1]]$Qinv))
  
  log_dens <- log_dens + seq2_model$mlik[1] + 0.5*as.numeric(determinant(Q_gmrf_2)$modulus)
  
  df_seq2i <- mclapply(X = 1:length(u_t), mc.cores = 4, FUN = function(i){
    data.frame(ID = i, mean = seq2_model$misc$configs$config[[1]]$improved.mean[i], 
               X0.025quant = qnorm(p = 0.025, mean = seq2_model$misc$configs$config[[1]]$improved.mean[i], sd = Sigma_gmrf_2[i]),
               X0.975quant = qnorm(p = 0.975, mean = seq2_model$misc$configs$config[[1]]$improved.mean[i], sd = Sigma_gmrf_2[i])
               )
  }) %>% do.call(., what = rbind) %>% as.data.frame(.)
  
  df_seq2c <- mclapply(X = 1:length(u_t), mc.cores = 4, FUN = function(i){
    data.frame(ID = i, mean = mean_gmrf_2[i],
               X0.025quant = qnorm(p = 0.025, mean = mean_gmrf_2[i], sd = Sigma_gmrf_2[i]),
               X0.975quant = qnorm(p = 0.975, mean = mean_gmrf_2[i], sd = Sigma_gmrf_2[i])
    )
  }) %>% do.call(., what = rbind) %>% as.data.frame(.)
  
  gg_seq2 <- 
    ggplot() + 
    geom_line(data = df_seq2i, mapping = aes(x = ID, y = mean), color = "red") +
    geom_ribbon(data = df_seq2i, mapping = aes(x = ID, ymin = X0.025quant, ymax = X0.975quant), fill = "red", alpha = 0.25) +
    geom_line(data = df_seq2c, mapping = aes(x = ID, y = mean), color = "blue") + 
    geom_ribbon(data = df_seq2c, mapping = aes(x = ID, ymin = X0.025quant, ymax = X0.975quant), fill = "blue", alpha = 0.25) +
    geom_line(data = data.frame(ID = 1:length(u_t), mean = u_t), mapping = aes(x = ID, y = mean)) + 
    labs(title = "B. Posterior (second step)") + theme_bw() +
    theme(plot.title = element_markdown(hjust=0, face="bold", size=16, color = "black"))
  
  k_seq <- 3
  y_sim_seq3 <- rep(NA, length(y_sim))
  y_sim_seq3[(size_seq*(k_seq-1)+1):(size_seq*k_seq)] <- y_sim[(size_seq*(k_seq-1)+1):(size_seq*k_seq)]
  
  formula_seq3 <- y ~ -1 + offset(mean_gmrf_2) + f(gmrf, model = "generic0", Cmatrix = Q_gmrf_2, extraconstr = list(A = seq1_model$misc$configs$constr$A, e = seq1_model$misc$configs$constr$e), hyper = list(prec = list(initial = log(1), fixed = TRUE)))
  seq3_model <- inla(data = list(y = y_sim_seq3, mean_grmf_2 = mean_gmrf_2, gmrf = 1:nrow(Q_gmrf_2)), family = "gaussian",
                     formula = formula_seq3,
                     control.predictor = list(A = A_seq),
                     control.compute = list(config = TRUE),
                     control.family = list(hyper = list(prec = list(initial = seq1_model$joint.hyper[1,1], fixed = TRUE))),
                     verbose = FALSE)
  
  mean_gmrf_3 <- mean_gmrf_2 + seq3_model$misc$configs$config[[1]]$improved.mean
  Q_gmrf_3 <- fix.Q(Q = seq3_model$misc$configs$config[[1]]$Q)
  Sigma_gmrf_3 <- sqrt(diag(seq3_model$misc$configs$config[[1]]$Qinv))
  
  log_dens <- log_dens + seq3_model$mlik[1] + 0.5*as.numeric(determinant(Q_gmrf_3)$modulus)
  
  df_seq3i <- mclapply(X = 1:length(u_t), mc.cores = 4, FUN = function(i){
    data.frame(ID = i, mean = seq3_model$misc$configs$config[[1]]$improved.mean[i], 
               X0.025quant = qnorm(p = 0.025, mean = seq3_model$misc$configs$config[[1]]$improved.mean[i], sd = Sigma_gmrf_3[i]),
               X0.975quant = qnorm(p = 0.975, mean = seq3_model$misc$configs$config[[1]]$improved.mean[i], sd = Sigma_gmrf_3[i])
    )
  }) %>% do.call(., what = rbind) %>% as.data.frame(.)
  
  df_seq3c <- mclapply(X = 1:length(u_t), mc.cores = 4, FUN = function(i){
    data.frame(ID = i, mean = mean_gmrf_3[i],
               X0.025quant = qnorm(p = 0.025, mean = mean_gmrf_3[i], sd = Sigma_gmrf_3[i]),
               X0.975quant = qnorm(p = 0.975, mean = mean_gmrf_3[i], sd = Sigma_gmrf_3[i])
    )
  }) %>% do.call(., what = rbind) %>% as.data.frame(.)
  
  gg_seq3 <- 
    ggplot() + 
    geom_line(data = df_seq3i, mapping = aes(x = ID, y = mean), color = "red") +
    geom_ribbon(data = df_seq3i, mapping = aes(x = ID, ymin = X0.025quant, ymax = X0.975quant), fill = "red", alpha = 0.25) +
    geom_line(data = df_seq3c, mapping = aes(x = ID, y = mean), color = "blue") + 
    geom_ribbon(data = df_seq3c, mapping = aes(x = ID, ymin = X0.025quant, ymax = X0.975quant), fill = "blue", alpha = 0.25) +
    geom_line(data = data.frame(ID = 1:length(u_t), mean = u_t), mapping = aes(x = ID, y = mean)) + 
    labs(title = "C. Posterior (third step)") + theme_bw() +
    theme(plot.title = element_markdown(hjust=0, face="bold", size=16, color = "black"))
  
  gg_seqJ <- 
    ggplot() + 
    geom_line(data = data.frame(full_model$summary.random$t), mapping = aes(x = ID, y = mean), color = "red") +
    geom_ribbon(data = data.frame(full_model$summary.random$t), mapping = aes(x = ID, ymin = X0.025quant, ymax = X0.975quant), fill = "red", alpha = 0.25) +
    geom_line(data = df_seq3c, mapping = aes(x = ID, y = mean), color = "blue") + 
    geom_ribbon(data = df_seq3c, mapping = aes(x = ID, ymin = X0.025quant, ymax = X0.975quant), fill = "blue", alpha = 0.25) +
    geom_line(data = data.frame(ID = 1:length(u_t), mean = u_t), mapping = aes(x = ID, y = mean)) + 
    labs(title = "D. Posterior comparison") + theme_bw() +
    theme(plot.title = element_markdown(hjust=0, face="bold", size=16, color = "black"))
  
  grid.arrange(arrangeGrob(grobs = list(gg_seq1, gg_seq2, gg_seq3, gg_seqJ), ncol = 2))
  
  out <- list(Qpost = Q_gmrf_3, mean_grmf = mean_gmrf_3, log_dens = log_dens)
  return(out)
})

log_post_dens <- c()
for(i in 1:length(res_seq)){
  log_post_dens[i] <- res_seq[[i]]$log_dens
}

# log_post_dens <- log_post_dens - max(log_post_dens)

post_seq_joint_int_hyper <- seq1_model$joint.hyper
post_seq_joint_int_hyper$`Log posterior density` <- log_post_dens

# gg_seq3 <-
#   ggplot() + 
#   geom_line(data = data.frame(ID = 1:length(u_t), mean = res_seq[[which.max(post_seq_joint_int_hyper[,5])]]$mean_grmf[1:length(u_t)]), mapping = aes(x = ID, y = mean), color = "blue") +
#   geom_line(data = data.frame(ID = 1:length(u_t), mean = full_model$misc$configs$config[[1]]$improved.mean[1:length(u_t)]), mapping = aes(x = ID, y = mean), color = "red") +
#   geom_line(data = data.frame(ID = 1:length(u_t), mean = u_t), mapping = aes(x = ID, y = mean)) + 
#   labs(title = "Posterior mean (s3)") + theme_bw() +
#   theme(plot.title = element_markdown(hjust=0.5, face="bold", size=16, color = "black"))

grid.xy <- expand.grid(x = seq(0, 1, length.out = 1E2), y = seq(0, 1, length.out = 1E2))
A_grid <- fm_basis(x = mesh_inf, loc = as.matrix(grid.xy))

seq_sp_grid <- drop(A_grid %*% res_seq[[which.max(post_seq_joint_int_hyper[,5])]]$mean_grmf[length(u_t)+1:mesh_inf$n])
full_sp_grid <- drop(A_grid %*% full_model$misc$configs$config[[1]]$improved.mean[length(u_t)+1:mesh_inf$n])

seq_sp_sd_grid <- drop(A_grid %*% diag(res_seq[[which.max(post_seq_joint_int_hyper[,5])]]$Qpost)[length(u_t)+1:mesh_inf$n]**(-1/2))
full_sp_sd_grid <- drop(A_grid %*% diag(full_model$misc$configs$config[[1]]$Q)[length(u_t)+1:mesh_inf$n]**(-1/2))

list_den_ggplot <- list()
for(i in 1:50){
  xseq <- seq(from = seq_sp_grid[i] - 4*seq_sp_sd_grid[i], to = seq_sp_grid[i] + 4*seq_sp_sd_grid[i], length.out = 1E3)
  
  list_den_ggplot[[i]] <- ggplot() +
    geom_line(data.frame(x = xseq, y = dnorm(x = xseq, mean = seq_sp_grid[i], sd = seq_sp_sd_grid[i])), mapping = aes(x = x, y = y), color = "blue") +
    geom_line(data.frame(x = xseq, y = dnorm(x = xseq, mean = full_sp_grid[i], sd = full_sp_sd_grid[i])), mapping = aes(x = x, y = y), color = "red") +
    theme_bw()
}

grid.arrange(arrangeGrob(grobs = list_den_ggplot))

gg_seq_sp <- ggplot() +
  geom_tile(data = data.frame(grid.xy, sp = seq_sp_grid), mapping = aes(x = x, y = y, fill = sp)) +
  scale_fill_viridis_c(option = "turbo") + 
  labs(title = "Posterior mean (recursive INLA)", fill = "mean") + theme_bw() +
  theme(plot.title = element_markdown(hjust=0.5, face="bold", size=14, color = "black"))

gg_full_sp <- ggplot() +
  geom_tile(data = data.frame(grid.xy, sp = full_sp_grid), mapping = aes(x = x, y = y, fill = sp)) +
  scale_fill_viridis_c(option = "turbo") + 
  labs(title = "Posterior mean (INLA)", fill = "mean") + theme_bw() +
  theme(plot.title = element_markdown(hjust=0.5, face="bold", size=14, color = "black"))

gg_arr1_sp_mean <- ggarrange(plotlist = list(gg_full_sp, gg_seq_sp), ncol = 1, common.legend = TRUE, legend = "right")

gg_seq_sp_sd <- ggplot() +
  geom_tile(data = data.frame(grid.xy, sp = seq_sp_sd_grid), mapping = aes(x = x, y = y, fill = sp)) +
  scale_fill_viridis_c(option = "turbo") + 
  labs(title = "Posterior stdev. (recursive INLA)", fill = "sd") + theme_bw() +
  theme(plot.title = element_markdown(hjust=0.5, face="bold", size=14, color = "black"))

gg_full_sp_sd <- ggplot() +
  geom_tile(data = data.frame(grid.xy, sp = full_sp_sd_grid), mapping = aes(x = x, y = y, fill = sp)) +
  scale_fill_viridis_c(option = "turbo") + 
  labs(title = "Posterior stdev. (INLA)", fill = "sd") + theme_bw() +
  theme(plot.title = element_markdown(hjust=0.5, face="bold", size=14, color = "black"))

gg_arr1_sp_sd <- ggarrange(plotlist = list(gg_full_sp_sd, gg_seq_sp_sd), ncol = 1, common.legend = TRUE, legend = "right")

gg_full_sp_mean_diff <- ggplot() +
  geom_tile(data = data.frame(grid.xy, sp = seq_sp_grid - full_sp_grid), mapping = aes(x = x, y = y, fill = sp)) +
  scale_fill_viridis_c(option = "turbo") + 
  labs(title = "Difference between means", fill = "mean") + theme_bw() +
  theme(plot.title = element_markdown(hjust=0.5, face="bold", size=14, color = "black"))

gg_full_sp_sd_diff <- ggplot() +
  geom_tile(data = data.frame(grid.xy, sp = seq_sp_sd_grid - full_sp_sd_grid), mapping = aes(x = x, y = y, fill = sp)) +
  scale_fill_viridis_c(option = "turbo") + 
  labs(title = "Difference between stdev.", fill = "sd") + theme_bw() +
  theme(plot.title = element_markdown(hjust=0.5, face="bold", size=14, color = "black"))

gg_arr1_sp_diff <- ggarrange(plotlist = list(gg_full_sp_mean_diff, gg_full_sp_sd_diff), ncol = 1, common.legend = FALSE)

grid.arrange(arrangeGrob(grobs = list(gg_arr1_sp_mean, gg_arr1_sp_sd, gg_arr1_sp_diff), ncol = 3)) #layout_matrix = matrix(data = c(1,1,3,3,2,2,3,3), byrow = TRUE)))
