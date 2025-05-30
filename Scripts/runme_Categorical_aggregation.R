# Categorical aggregation. Example 1. Full model inference and sequential Bayesian inference

library(INLA)
library(inlabru)
library(ggplot2)
library(ggtext)
library(gridExtra)

A_cat <- function(idx, idx.agg){
  if(missing(idx.agg)){
    return(lapply(X = idx %>% unique(.), FUN = function(i){as.numeric(idx==i)}) %>% do.call(what = cbind, .))
  } else{
    A_agg <- lapply(X = idx.agg %>% unique(.), FUN = function(i){as.numeric(idx.agg==i)}) %>% do.call(what = cbind, .)
    A_idx <- matrix(0, nrow = nrow(A_agg), ncol = idx %>% unique(.) %>% length(.))
    for(i in unique(idx.agg)){
      A_idx[,idx[idx.agg==i] %>% unique(.)] <- A_agg[,i]
    }
    return(A_idx)
  }
}

# Examples are given for three different types of data in "categorical crossing scales" ----

## Simulation of non-aggregated categories ----

cat_values <- rnorm(1:5, sd = 1) %>% matrix(data = ., nrow = 1) %>% apply(X = ., MARGIN = 1, FUN = function(x){x-mean(x)})

idx.cat <- rep(1:5, each = 10)
lin_pred <- -1 + cat_values %>% rep(x = ., each = 10)

### Gaussian data ----
prec_gauss <- 36
y_gauss <- rnorm(n = length(lin_pred), mean = lin_pred, sd = prec_gauss**(-1/2))

## Simulation of aggregated categories ----

cat.agg_values <- c(cat_values[1:3] %>% sum(.), cat_values[-c(1:3)])
idx.cat.agg <- rep(x = 1:3, times = c(30, 10, 10))
lin_pred.agg <- -1 + cat.agg_values %>% rep(x = ., times = c(30, 10, 10))

### Gaussian data ----
prec_gauss <- 36
y_gauss.agg <- rnorm(n = length(lin_pred.agg), mean = lin_pred.agg, sd = prec_gauss**(-1/2))

## Inference ----

### Gaussian data analysis ----
model_gaussian <- 
  inla(formula = y ~ -1 + beta0 + f(idx.cat, model = "iid", constr = TRUE),
       data = list(y = y_gauss, idx.cat = idx.cat, beta0 = rep(1, length(lin_pred))),
       control.compute = list(config = TRUE),
       family=  "gaussian", verbose = FALSE)

model_gaussian.agg <- 
  inla(formula = y ~ -1 + beta0 + f(idx.cat, model = "iid", constr = TRUE),
       data = list(y = y_gauss.agg, idx.cat = idx.cat.agg, beta0 = rep(1, length(lin_pred.agg))),
       family=  "gaussian", verbose = FALSE)

ggplot() +
  geom_line(data = model_gaussian$marginals.fixed$beta0, mapping = aes(x = x, y = y), color = "red") + 
  geom_line(data = model_gaussian.agg$marginals.fixed$beta0, mapping = aes(x = x, y = y), color = "blue") +
  theme_bw()

ggplot() +
  geom_line(data = model_gaussian$marginals.random$idx.cat$index.4, mapping = aes(x = x, y = y), color = "red") + 
  geom_line(data = model_gaussian.agg$marginals.random$idx.cat$index.2, mapping = aes(x = x, y = y), color = "blue") +
  theme_bw()

ggplot() +
  geom_line(data = model_gaussian$marginals.random$idx.cat$index.5, mapping = aes(x = x, y = y), color = "red") + 
  geom_line(data = model_gaussian.agg$marginals.random$idx.cat$index.3, mapping = aes(x = x, y = y), color = "blue") +
  theme_bw()

#### Joint model ----

inf_stk <- inla.stack(data = list(y = c(y_gauss, y_gauss.agg)),
                      A = list(1, rbind(A_cat(idx = idx.cat), A_cat(idx = idx.cat, idx.agg = idx.cat.agg))),
                      effects = list(
                        list(beta0 = rep(1, times = length(c(y_gauss, y_gauss.agg)))),
                        list(idx.cat = idx.cat %>% unique(.))
                      ),
                      tag = "inf_stk")

model_gaussian_joint <- 
  inla(formula = y ~ -1 + beta0 + f(idx.cat, model = "iid", constr = TRUE),
       data = inla.stack.data(inf_stk),
       control.predictor = list(A = inla.stack.A(inf_stk)),
       control.compute = list(config = TRUE),
       family=  "gaussian", verbose = FALSE)

gg_beta0 <- ggplot() +
  geom_line(data = model_gaussian_joint$marginals.fixed$beta0, mapping = aes(x = x, y = y), color = "red") +
  geom_line(data = model_gaussian$marginals.fixed$beta0, mapping = aes(x = x, y = y), color = "blue") + 
  geom_line(data = model_gaussian.agg$marginals.fixed$beta0, mapping = aes(x = x, y = y), color = "green") +
  geom_vline(xintercept = -1) + theme_bw() + labs(title = "Intercept") +
  theme(plot.title = element_markdown(hjust=0.5, face="bold", size=16, color = "black"))

gg_cat1 <- ggplot() +
  geom_line(data = model_gaussian_joint$marginals.random$idx.cat$index.1, mapping = aes(x = x, y = y), color = "red") +
  geom_line(data = model_gaussian$marginals.random$idx.cat$index.1, mapping = aes(x = x, y = y), color = "blue") + 
  geom_vline(xintercept = cat_values[1]) + theme_bw() + labs(title = "Level 1 (category)") +
  theme(plot.title = element_markdown(hjust=0.5, face="bold", size=16, color = "black"))

gg_cat2 <- ggplot() +
  geom_line(data = model_gaussian_joint$marginals.random$idx.cat$index.2, mapping = aes(x = x, y = y), color = "red") +
  geom_line(data = model_gaussian$marginals.random$idx.cat$index.2, mapping = aes(x = x, y = y), color = "blue") + 
  geom_vline(xintercept = cat_values[2]) + theme_bw() + labs(title = "Level 2 (category)") +
  theme(plot.title = element_markdown(hjust=0.5, face="bold", size=16, color = "black"))

gg_cat3 <- ggplot() +
  geom_line(data = model_gaussian_joint$marginals.random$idx.cat$index.3, mapping = aes(x = x, y = y), color = "red") +
  geom_line(data = model_gaussian$marginals.random$idx.cat$index.3, mapping = aes(x = x, y = y), color = "blue") + 
  geom_vline(xintercept = cat_values[3]) + theme_bw() + labs(title = "Level 3 (category)") +
  theme(plot.title = element_markdown(hjust=0.5, face="bold", size=16, color = "black"))

gg_cat4 <- ggplot() +
  geom_line(data = model_gaussian_joint$marginals.random$idx.cat$index.4, mapping = aes(x = x, y = y), color = "red") +
  geom_line(data = model_gaussian$marginals.random$idx.cat$index.4, mapping = aes(x = x, y = y), color = "blue") +
  geom_line(data = model_gaussian.agg$marginals.random$idx.cat$index.2, mapping = aes(x = x, y = y), color = "green") +
  geom_vline(xintercept = cat_values[4]) + theme_bw() + labs(title = "Level 4 (category)") +
  theme(plot.title = element_markdown(hjust=0.5, face="bold", size=16, color = "black"))

gg_cat5 <- ggplot() +
  geom_line(data = model_gaussian_joint$marginals.random$idx.cat$index.5, mapping = aes(x = x, y = y), color = "red") +
  geom_line(data = model_gaussian$marginals.random$idx.cat$index.5, mapping = aes(x = x, y = y), color = "blue") + 
  geom_line(data = model_gaussian.agg$marginals.random$idx.cat$index.3, mapping = aes(x = x, y = y), color = "green") +
  geom_vline(xintercept = cat_values[5]) + theme_bw() + labs(title = "Level 5 (category)") +
  theme(plot.title = element_markdown(hjust=0.5, face="bold", size=16, color = "black"))

gg_post <- arrangeGrob(grobs = list(gg_beta0, gg_cat1, gg_cat2, gg_cat3, gg_cat4, gg_cat5), ncol = 3)
grid.arrange(gg_post)
