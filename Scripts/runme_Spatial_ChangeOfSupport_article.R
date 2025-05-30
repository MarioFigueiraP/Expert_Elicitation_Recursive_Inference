# Example 1: Expert elicitation with spatial change of support

# Loading libraries ----

library(INLA)
library(inlabru)
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

## Voronoi tesselation given a set of points, used as centroids 
st_voronoi_point <- function(points){
  ## points must be POINT geometry
  # check for point geometry and continue if true
  if(!all(st_geometry_type(points) == "POINT")){
    stop("Input not  POINT geometries")
  }
  g = st_combine(st_geometry(points)) # make multipoint
  v = st_voronoi(g)
  v = st_collection_extract(v)
  return(v[unlist(st_intersects(points, v))])
}

## Creation of polygons with circle shapes
st_circle_polygon <- function(x, r, np){
  theta <- seq(0, 2*pi, length.out = np)
  xy_circle <- cbind(r*cos(theta) + x[1], r*sin(theta) + x[2])
  xy_circle[nrow(xy_circle),] <- xy_circle[1,]
  return(st_polygon(x = list(xy_circle)))
}

## Creation of a sf data frame with multiple disjoint circles as polygons
st_sfc_disjoint_polygon_circles <- function(seed, r, np, n_sim, n_sample){
  if(!missing(seed)){set.seed(seed = seed)}
  
  x <- cbind(runif(n_sim, min = r, max = 1-r), runif(n_sim, min = r, max = 1-r))
  idx_base <- 1:nrow(x)
  list_circles <- mclapply(X = idx_base, mc.cores = detectCores(), FUN = function(idx){st_circle_polygon(x = x[idx,], r = r, np = np)})
  df <- do.call(list_circles, what = st_sfc)
  
  idx_sel <- 1:nrow(x)
  int <- st_intersects(x = df, sparse = TRUE)
  for(i in 1:length(int)){
    if(i %in% idx_sel){
      idx_sel <- setdiff(x = idx_sel, y = setdiff(x = int[[i]], y = i))
    }}
  df <- df[idx_sel,]; df <- df[sample(x = 1:length(df), size = n_sample),]
  return(df)
}

# Setting the global seed ----
globalseed <- 1234
set.seed(globalseed)

# Simulating spatial effect leveraging the SPDE-FEM approach ----
rho_spatial <- 0.4 # spatial range
sigma_spatial <- 2 # marginal stdev. of the spatial effect
mesh <- fm_mesh_2d_inla(loc.domain=matrix(c(0,0,0,1,1,1,1,0), byrow=TRUE, ncol=2), max.edge=c(0.025,0.15)) # mesh creation
Q <- fm_matern_precision(x=mesh, alpha=2, rho=rho_spatial, sigma=sigma_spatial) # matrix SPDE-FEM
u <- fm_sample(n=1, Q=Q, mu=0, constr=NULL) # spatial effect
u <- u - mean(u)

xy.grid <- expand.grid(x=seq(0,1,length.out=5E2),y=seq(0,1,length.out=5E2)) # Grid to plot a continuous map
Asim <- inla.spde.make.A(mesh=mesh, loc=as.matrix(xy.grid)) # Projection matrix from nodes to grid points 
u_sim <- drop(Asim %*% u) # Spatial effect projected to the grid locations

beta0 <- 5 # Intercept
beta1 <- 1 # Beta coefficient 

formula <- function(x,y,s){ # Covariate formula (spatial structure)
  z <- eval(parse(text=s))
}
cov <- formula(x=xy.grid[,1], y=xy.grid[,2], s="2*(x-0.5)**2 + 3*y*log(2+(y-0.5)**2)") # Simulation of the covariate

DFSim <- data.frame( # Simulated data frame
  x=xy.grid$x,
  y=xy.grid$y,
  z=beta0 + beta1*cov + u_sim, # linear predictor 
  cov=cov,
  sp = u_sim
)

gg_sp.sim <- ggplot() + 
  geom_tile(data = data.frame(DFSim[,1:2], sp = u_sim), mapping = aes(x=x, y=y, fill=sp)) + 
  scale_fill_viridis_c(option = "turbo") + 
  theme_bw()

DFSample <- DFSim[sample(x=1:nrow(DFSim),size=2E1),]
prec_obs <- 16
DFSample$ysim <- rnorm(n = nrow(DFSample), mean = DFSample$z, sd = prec_obs**(-1/2))

# First spatial structure (voronoi) for experts ----

envelopment <- st_polygon(x=list(matrix(c(0,0,0,1,1,1,1,0,0,0), ncol=2, byrow=TRUE)))
up_triang <- st_polygon(x=list(matrix(c(0,0,0,1,1,1,0,0), ncol=2, byrow=TRUE)))
low_triang <- st_polygon(x=list(matrix(c(0,0,1,0,1,1,0,0), ncol=2, byrow=TRUE)))
p <-  st_as_sf(data.frame(x = runif(80), y = runif(80)), coords = 1:2)
v <-  st_voronoi_point(p)

v_env <- st_intersection(v, envelopment)
ggplot() + geom_sf(data=v_env, mapping=aes())

DFv_env <- st_sf(data.frame(idx=as.character(1:length(v_env))), v_env)
ggplot() + geom_sf(data=DFv_env, mapping=aes()) + geom_sf_label(data=DFv_env, mapping=aes(label=idx))

indx_SimPoints <- st_intersects(v_env,st_as_sf(DFSim, coords=1:2),sparse=TRUE)
indx_up <- st_intersects(up_triang,v_env,sparse=TRUE)
indx_low <- st_intersects(low_triang,v_env,sparse=TRUE)
indx_low <- setdiff(unlist(indx_low), unlist(indx_up))

Agg <- lapply(X=indx_SimPoints, FUN=function(i){mean(DFSim[i,"z"])})
Agg_cov <- lapply(X=indx_SimPoints, FUN=function(i){mean(DFSim[i,"cov"])}) %>% unlist(.)

Agg.values <- unlist(Agg)

prec_exp1 <- 9
prec_exp2 <- 4
corr_exp12 <- 0.75
Sigma <- Matrix::Matrix(data = rep(0,4), nrow = 2, ncol = 2)
val_sigma <- c(c(prec_exp1,prec_exp2)**(-1),  corr_exp12*c(prec_exp1*prec_exp2,prec_exp1*prec_exp2)**(-1/2))
Sigma[1,1] <- val_sigma[1]
Sigma[2,2] <- val_sigma[2]
Sigma[1,2] <- val_sigma[3]
Sigma[2,1] <- val_sigma[4]

mat_exp12 <- do.call(lapply(X = 1:length(Agg.values), FUN = function(x){MASS::mvrnorm(n = 1, mu = rep(Agg.values[x], times = 2), Sigma = Sigma)}), what = rbind)

exp1 <- exp2 <- rep(NA,times=length(v_env))
exp1[unlist(indx_up)] <- mat_exp12[unlist(indx_up), 1]
exp2[unlist(indx_low)] <- mat_exp12[unlist(indx_low), 2]
DFexp1 <- DFexp1_S1 <- st_sf(data.frame(yexp1=exp1), st_geometry(DFv_env))
DFexp2 <- DFexp2_S1 <- st_sf(data.frame(yexp2=exp2), st_geometry(DFv_env))

csc <- colsc(DFSim$z, unlist(Agg))

gg1 <- 
  ggplot() + geom_point(data=DFSample, mapping=aes(x=x, y=y, color=ysim)) +
  coord_equal() + csc$color +
  labs(title = "<span style = 'color: blue;'>Observational data</span>", color="Y", x="Latitude", y="Longitude") +
  theme_bw() + theme(plot.title = element_markdown(hjust=0.5, face="bold", size=16))
gg2 <- 
  ggplot() + geom_sf(data=DFexp1, mapping=aes(fill=yexp1)) + csc$fill +
  labs(title = "Expert<sub>1</sub> data", color="Y", x="Latitude", y="Longitude") +
  theme_bw() + theme(plot.title = element_markdown(hjust=0.5, face="bold", size=16, color = "blue")) 
gg3 <- 
  ggplot() + geom_sf(data=DFexp2, mapping=aes(fill=yexp2)) + csc$fill + 
  labs(title = "Expert<sub>2</sub> data", color="Y", x="Latitude", y="Longitude") +
  theme_bw() + theme(plot.title = element_markdown(hjust=0.5, face="bold", size=16, color = "blue"))

ggt <- arrangeGrob(grobs=list(gg1,gg2,gg3), ncol=4)
ggarrange(gg1, gg2, gg3, ncol = 3, common.legend = TRUE, legend = "right")

mesh_inf_s1 <- fm_mesh_2d_inla(loc.domain=matrix(c(0,0,0,1,1,1,1,0), byrow=TRUE, ncol=2), max.edge=c(0.04,0.15))
spde <- inla.spde2.pcmatern(mesh = mesh_inf_s1, alpha = 2, prior.range = c(0.2,0.5), prior.sigma = c(1,0.5), constr = FALSE)

Aobs <- inla.spde.make.A(mesh = mesh_inf_s1, loc = as.matrix(DFSample[,1:2]))
Aexp <- inla.spde.make.A(mesh = mesh_inf_s1, loc = as.matrix(DFSim[,1:2]))
indx_blocking <- st_intersects(st_as_sf(DFSim, coords=1:2), v_env, sparse=TRUE)
Aexp_block <- inla.spde.make.block.A(A = Aexp, block = unlist(indx_blocking), rescale = "count")

inf_stk_only_obs <- inla.stack(data = list(y = DFSample$ysim),
                               A = list(Aobs, 1),
                               effects = list(
                                 list(sp = 1:mesh_inf_s1$n),
                                 list(beta0 = rep(1, nrow(DFSample)), beta1 = DFSample$cov)
                               ),
                               tag = "inf_obs")

inf_stk_obs <- inla.stack(data = list(y = cbind(DFSample$ysim, NA)),
                          A = list(Aobs, 1),
                          effects = list(
                            list(sp = 1:mesh_inf_s1$n),
                            list(beta0 = rep(1, nrow(DFSample)), beta1 = DFSample$cov)
                          ),
                          tag = "inf_obs")

inf_stk_exp12 <- inla.stack(data = list(y = cbind(NA, c(exp1, exp2))),
                            A = list(rbind(Aexp_block, Aexp_block), 1),
                            effects = list(
                              list(sp = 1:mesh_inf_s1$n),
                              # list(sp.copy = 1:mesh_inf_s1$n),
                              list(beta0 = rep(1, 2*length(exp1)),
                                   beta1 = c(Agg_cov, Agg_cov), idx_sigma = 1:(2*length(exp1)))
                            ),
                            tag = "inf_exp12")

inf_stk_total <- inla.stack(inf_stk_obs, inf_stk_exp12)

model_formula <- y ~ -1 + beta0 + beta1 + f(sp, model = spde) + f(idx_sigma, model = "iid2d", n = 2*length(exp1))
inf_model_expS1 <- inla(family = c("gaussian", "gaussian"),
                        data = inla.stack.data(inf_stk_total),
                        formula = model_formula,
                        control.predictor = list(A = inla.stack.A(inf_stk_total)),
                        control.family = list(list(),list(hyper=list(prec = list(initial = 14, fixed = TRUE)))),
                        verbose = FALSE)

## Modelling only the observations with the same model components ----

model_formula_only_obs <- y ~ -1 + beta0 + beta1 + f(sp, model = spde)
inf_model_obs <- inla(family = "gaussian",
                      data = inla.stack.data(inf_stk_only_obs),
                      formula = model_formula_only_obs,
                      control.predictor = list(A = inla.stack.A(inf_stk_only_obs)),
                      verbose = FALSE)

# Second spatial structure (circles) for experts ----

## Creating the polygons (circles) for expert data ----
sf_exp1 <- st_sfc_disjoint_polygon_circles(seed = 1234, r = 0.035, np = 30, n_sim = 200, n_sample = 65)
sf_exp2 <- st_sfc_disjoint_polygon_circles(seed = 4321, r = 0.035, np = 30, n_sim = 200, n_sample = 65)

ggplot() + 
  geom_sf(data = sf_exp1, mapping = aes(), fill = "red", alpha = 0.5) +
  geom_sf(data = sf_exp2, mapping = aes(), fill = "blue", alpha = 0.5) +
  theme_bw()

## Aggregating data for experts ----
indx_SimPoints_exp1 <- st_intersects(sf_exp1, st_as_sf(DFSim, coords=1:2), sparse=TRUE)
indx_SimPoints_exp2 <- st_intersects(sf_exp2, st_as_sf(DFSim, coords=1:2), sparse=TRUE)

Agg_exp1 <- lapply(X = indx_SimPoints_exp1, FUN = function(i){mean(DFSim[i,"z"])})
Agg_cov_exp1 <- lapply(X = indx_SimPoints_exp1, FUN = function(i){mean(DFSim[i,"cov"])}) %>% unlist(.)
Agg.values_exp1 <- unlist(Agg_exp1)

Agg_exp2 <- lapply(X = indx_SimPoints_exp2, FUN = function(i){mean(DFSim[i,"z"])})
Agg_cov_exp2 <- lapply(X = indx_SimPoints_exp2, FUN = function(i){mean(DFSim[i,"cov"])}) %>% unlist(.)
Agg.values_exp2 <- unlist(Agg_exp2)

Agg.values <- c(Agg.values_exp1, Agg.values_exp2)

prec_exp1 <- 9
prec_exp2 <- 4
corr_exp12 <- 0.75
Sigma <- Matrix::Matrix(data = rep(0,4), nrow = 2, ncol = 2)
val_sigma <- c(c(prec_exp1,prec_exp2)**(-1),  corr_exp12*c(prec_exp1*prec_exp2,prec_exp1*prec_exp2)**(-1/2))
Sigma[1,1] <- val_sigma[1]
Sigma[2,2] <- val_sigma[2]
Sigma[1,2] <- val_sigma[3]
Sigma[2,1] <- val_sigma[4]

mat_exp12 <- do.call(lapply(X = 1:length(Agg.values), FUN = function(x){MASS::mvrnorm(n = 1, mu = rep(Agg.values[x], times = 2), Sigma = Sigma)}), what = rbind)

exp1 <- exp2 <- rep(NA,times=nrow(mat_exp12))
exp1[1:65] <- mat_exp12[1:65, 1]
exp2[66:130] <- mat_exp12[66:130, 2]
DFexp1 <- DFexp1_S2 <- st_sf(data.frame(yexp1 = exp1), geometry = st_geometry(c(sf_exp1, sf_exp2)))
DFexp2 <- DFexp2_S2 <- st_sf(data.frame(yexp2 = exp2), geometry = st_geometry(c(sf_exp1, sf_exp2)))

csc <- colsc(DFSim$z, c(Agg.values_exp1, Agg.values_exp2))

gg1 <- 
  ggplot() + geom_sf(data = st_as_sf(DFSample, coords = c("x", "y")), mapping = aes(color=ysim)) +
  csc$color +
  labs(title = "<span style = 'color: blue;'>Observational data</span>", color="Y", x="Latitude", y="Longitude") +
  theme_bw() + theme(plot.title = element_markdown(hjust=0.5, face="bold", size=16))
gg2 <- 
  ggplot() + geom_sf(data = DFexp1[1:65,], mapping = aes(fill=yexp1)) + csc$fill +
  labs(title = "Expert<sub>1</sub> data", color = "Y", x = "Latitude", y = "Longitude") +
  theme_bw() + theme(plot.title = element_markdown(hjust=0.5, face="bold", size=16, color = "blue")) 
gg3 <- 
  ggplot() + geom_sf(data = DFexp2[66:130,], mapping = aes(fill=yexp2)) + csc$fill + 
  labs(title = "Expert<sub>2</sub> data", color = "Y", x = "Latitude", y = "Longitude") +
  theme_bw() + theme(plot.title = element_markdown(hjust=0.5, face="bold", size=16, color = "blue")) 

ggt <- arrangeGrob(grobs=list(gg1,gg2,gg3), ncol=4)
ggarrange(gg1, gg2, gg3, ncol = 3, common.legend = TRUE, legend = "right")

mesh_inf_s2 <- fm_mesh_2d_inla(loc.domain=matrix(c(0,0,0,1,1,1,1,0), byrow=TRUE, ncol=2), max.edge=c(0.04,0.15))
spde <- inla.spde2.pcmatern(mesh = mesh_inf_s2, alpha = 2, prior.range = c(0.5,0.5), prior.sigma = c(1,0.5), constr = FALSE)

Aobs <- inla.spde.make.A(mesh = mesh_inf_s2, loc = as.matrix(DFSample[,1:2]))
indx_blocking1 <- st_intersects(st_as_sf(DFSim, coords=1:2), sf_exp1, sparse = TRUE)
Aexp1 <- inla.spde.make.A(mesh = mesh_inf_s2, loc = as.matrix(DFSim[(lapply(X = indx_blocking1, FUN = length) %>% unlist(.))!=0,1:2]))
Aexp_block1 <- inla.spde.make.block.A(A = Aexp1, block = unlist(indx_blocking1), rescale = "count")
indx_blocking2 <- st_intersects(st_as_sf(DFSim, coords=1:2), sf_exp2, sparse = TRUE)
Aexp2 <- inla.spde.make.A(mesh = mesh_inf_s2, loc = as.matrix(DFSim[(lapply(X = indx_blocking2, FUN = length) %>% unlist(.))!=0,1:2]))
Aexp_block2 <- inla.spde.make.block.A(A = Aexp2, block = unlist(indx_blocking2), rescale = "count")

inf_stk_obs <- inla.stack(data = list(y = cbind(DFSample$ysim, NA)),
                          A = list(Aobs, 1),
                          effects = list(
                            list(sp = 1:mesh_inf_s2$n),
                            list(beta0 = rep(1, nrow(DFSample)), beta1 = DFSample$cov)
                          ),
                          tag = "inf_obs")

inf_stk_exp12 <- inla.stack(data = list(y = cbind(NA, c(exp1, exp2))),
                            A = list(rbind(Aexp_block1, Aexp_block1, Aexp_block2, Aexp_block2), 1),
                            effects = list(
                              list(sp = 1:mesh_inf_s2$n),
                              list(beta0 = rep(1, 2*length(exp1)),
                                   beta1 = c(c(Agg_cov_exp1, rep(NA, 65)), c(Agg_cov_exp2, rep(NA, 65))), idx_sigma = 1:(2*length(exp1)))
                            ),
                            tag = "inf_exp12")

inf_stk_total <- inla.stack(inf_stk_obs, inf_stk_exp12)

model_formula <- y ~ -1 + beta0 + beta1 + f(sp, model = spde) + f(idx_sigma, model = "iid2d", n = 2*length(exp1))
inf_model_expS2 <- inla(family = c("gaussian", "gaussian"),
                        data = inla.stack.data(inf_stk_total),
                        formula = model_formula,
                        control.predictor = list(A = inla.stack.A(inf_stk_total)),
                        control.family = list(list(),list(hyper=list(prec = list(initial = 14, fixed = TRUE)))),
                        verbose = FALSE)


# Graphical results ----

## Graphics for the expert information ----

csc_obs_exp <- colsc(DFSample$z, DFexp1_S1$yexp1, DFexp1_S2$yexp1, DFexp2_S1$yexp2, DFexp2_S2$yexp2)

gg_obs <- 
  ggplot() + geom_point(data=DFSample, mapping=aes(x=x, y=y, color=ysim), size = 3) +
  coord_equal() + csc_obs_exp$color +
  labs(title = "<span style = 'color: black;'>Observational data</span>", color="Values", x="Latitude", y="Longitude") +
  theme_bw() + theme(plot.title = element_markdown(hjust=0.5, face="bold", size=16))
gg_exp1_s1 <- 
  ggplot() + geom_sf(data = DFexp1_S1, mapping = aes(fill=yexp1)) + csc_obs_exp$fill +
  labs(title = "Expert<sub>1</sub> data (S<sub>1</sub>)", color = "Y", x = "Latitude", y = "Longitude") +
  theme_bw() + theme(plot.title = element_markdown(hjust=0.5, face="bold", size=16, color = "black")) +
  theme(legend.position="none")
gg_exp2_s1 <- 
  ggplot() + geom_sf(data = DFexp2_S1, mapping = aes(fill=yexp2)) + csc_obs_exp$fill + 
  labs(title = "Expert<sub>2</sub> data (S<sub>1</sub>)", color = "Y", x = "Latitude", y = "Longitude") +
  theme_bw() + theme(plot.title = element_markdown(hjust=0.5, face="bold", size=16, color = "black")) +
  theme(legend.position="none")
gg_exp1_s2 <- 
  ggplot() + geom_sf(data = DFexp1_S2[1:65,], mapping = aes(fill=yexp1)) + csc_obs_exp$fill +
  labs(title = "Expert<sub>1</sub> data (S<sub>2</sub>)", color = "Y", x = "Latitude", y = "Longitude") +
  theme_bw() + theme(plot.title = element_markdown(hjust=0.5, face="bold", size=16, color = "black")) +
  theme(legend.position="none")
gg_exp2_s2 <- 
  ggplot() + geom_sf(data = DFexp2_S2[66:130,], mapping = aes(fill=yexp2)) + csc_obs_exp$fill + 
  labs(title = "Expert<sub>2</sub> data (S<sub>2</sub>)", color = "Y", x = "Latitude", y = "Longitude") +
  theme_bw() + theme(plot.title = element_markdown(hjust=0.5, face="bold", size=16, color = "black")) + 
  theme(legend.position="none")

gg_obs_exp <- arrangeGrob(grobs = list(gg_obs, gg_exp1_s1, gg_exp1_s2, gg_exp2_s1, gg_exp2_s2), layout_matrix = matrix(data = c(2,3,4,5,1,1), ncol = 3))
grid.arrange(gg_obs_exp)

## Graphical results for the spatial field ----

sp_obs_mean <- drop(inla.spde.make.A(mesh = mesh_inf_s1, loc = as.matrix(DFSim[,1:2])) %*% inf_model_obs$summary.random$sp$mean)
sp_expS1_mean <- drop(inla.spde.make.A(mesh = mesh_inf_s1, loc = as.matrix(DFSim[,1:2])) %*% inf_model_expS1$summary.random$sp$mean)
sp_expS2_mean <- drop(inla.spde.make.A(mesh = mesh_inf_s2, loc = as.matrix(DFSim[,1:2])) %*% inf_model_expS2$summary.random$sp$mean)

sp_obs_sd <- drop(inla.spde.make.A(mesh = mesh_inf_s1, loc = as.matrix(DFSim[,1:2])) %*% inf_model_obs$summary.random$sp$sd)
sp_expS1_sd <- drop(inla.spde.make.A(mesh = mesh_inf_s1, loc = as.matrix(DFSim[,1:2])) %*% inf_model_expS1$summary.random$sp$sd)
sp_expS2_sd <- drop(inla.spde.make.A(mesh = mesh_inf_s2, loc = as.matrix(DFSim[,1:2])) %*% inf_model_expS2$summary.random$sp$sd)

csc_sp_mean_models <- colsc(sp_obs_mean, sp_expS1_mean, sp_expS2_mean)
csc_sp_sd_models <- colsc(sp_obs_sd, sp_expS1_sd, sp_expS2_sd)

gg_sp_sim <- 
  ggplot() + 
  geom_tile(data = DFSim, mapping = aes(x=x, y=y, fill=sp)) + scale_fill_viridis_c(option = "turbo") +
  scale_color_viridis_c(option = "turbo") + coord_equal() + labs(fill = "Values") +
  labs(title = "Simulated", color="value", x="Latitude", y="Longitude") +
  theme_bw() + theme(plot.title = element_markdown(hjust=0.5, face="bold", size=16, color = "black"))

gg_sp_obs_mean <- 
  ggplot() + 
  geom_tile(data = data.frame(DFSim[,1:2], sp = sp_obs_mean), mapping = aes(x=x, y=y, fill=sp)) + scale_fill_viridis_c(option = "turbo") +
  scale_color_viridis_c(option = "turbo") + coord_equal() + csc_sp_mean_models$fill +
  labs(title = "Observational (mean)", fill="mean", x="Latitude", y="Longitude") +
  theme_bw() + theme(plot.title = element_markdown(hjust=0.5, face="bold", size=16, color = "black"))

gg_sp_obs_sd <- 
  ggplot() + 
  geom_tile(data = data.frame(DFSim[,1:2], sp = sp_obs_sd), mapping = aes(x=x, y=y, fill=sp)) + scale_fill_viridis_c(option = "turbo") +
  scale_color_viridis_c(option = "turbo") + coord_equal() + csc_sp_sd_models$fill +
  labs(title = "Observational (sd)", fill="sd", x="Latitude", y="Longitude") +
  theme_bw() + theme(plot.title = element_markdown(hjust=0.5, face="bold", size=16, color = "black"))

gg_sp_expS1_mean <- 
  ggplot() + 
  geom_tile(data = data.frame(DFSim[,1:2], sp = sp_expS1_mean), mapping = aes(x=x, y=y, fill=sp)) + scale_fill_viridis_c(option = "turbo") +
  scale_color_viridis_c(option = "turbo") + coord_equal() + csc_sp_mean_models$fill +
  labs(title = "Expert S<sub>1</sub> (mean)", color="value", x="Latitude", y="Longitude") +
  theme_bw() + theme(plot.title = element_markdown(hjust=0.5, face="bold", size=16, color = "black"), legend.position="none")

gg_sp_expS1_sd <- 
  ggplot() + 
  geom_tile(data = data.frame(DFSim[,1:2], sp = sp_expS1_sd), mapping = aes(x=x, y=y, fill=sp)) + scale_fill_viridis_c(option = "turbo") +
  scale_color_viridis_c(option = "turbo") + coord_equal() + csc_sp_sd_models$fill +
  labs(title = "Expert S<sub>1</sub> (sd)", color="value", x="Latitude", y="Longitude") +
  theme_bw() + theme(plot.title = element_markdown(hjust=0.5, face="bold", size=16, color = "black"), legend.position="none")

gg_sp_expS2_mean <- 
  ggplot() + 
  geom_tile(data = data.frame(DFSim[,1:2], sp = sp_expS2_mean), mapping = aes(x=x, y=y, fill=sp)) + scale_fill_viridis_c(option = "turbo") +
  scale_color_viridis_c(option = "turbo") + coord_equal() + csc_sp_mean_models$fill +
  labs(title = "Expert S<sub>2</sub> (mean)", color="value", x="Latitude", y="Longitude") +
  theme_bw() + theme(plot.title = element_markdown(hjust=0.5, face="bold", size=16, color = "black"), legend.position="none")

gg_sp_expS2_sd <- 
  ggplot() + 
  geom_tile(data = data.frame(DFSim[,1:2], sp = sp_expS2_sd), mapping = aes(x=x, y=y, fill=sp)) + scale_fill_viridis_c(option = "turbo") +
  scale_color_viridis_c(option = "turbo") + coord_equal() + csc_sp_sd_models$fill +
  labs(title = "Expert S<sub>2</sub> (sd)", color="value", x="Latitude", y="Longitude") +
  theme_bw() + theme(plot.title = element_markdown(hjust=0.5, face="bold", size=16, color = "black"), legend.position="none")

gg_sp_obs_exp <- arrangeGrob(grobs = list(gg_sp_sim, gg_sp_obs_mean, gg_sp_expS1_mean, gg_sp_expS2_mean, gg_sp_obs_sd, gg_sp_expS1_sd, gg_sp_expS2_sd), 
                             ncol = 3, layout_matrix = matrix(data = c(NA, 1, NA, 2, 3, 4, 5, 6, 7), ncol = 3, byrow = TRUE))
grid.arrange(gg_sp_obs_exp)
