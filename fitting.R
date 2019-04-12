# Notes -------------------------------------------------------------------

# 

# Initialisation ----------------------------------------------------------

rm(list = ls()) # Clear Workspace

source("functions.R")

seed <- 154815468 # seed also used for stan # 462528635
set.seed(seed)

run <- FALSE

stan_code <- "" # Path to stan code
res_file <- "" # Path to results file (load or save)

n_it <- 2000
n_chains <- 5

# Data --------------------------------------------------------------------

dfA <- read.csv("Data/Gause-yeist-competition-exp-1.csv", comment.char="#")
dfB <- read.csv("Data/Gause-yeist-competition-exp-2.csv", comment.char="#")

df <- process_data(dfA, dfB)

plot_data(df)
# ggsave("data.jpg", width = 28, height = 21, units = "cm", dpi = 300)

# Model -------------------------------------------------------------------

format_data <- function(df){
  df1_A <- subset(df, Condition == "Single" & Species == "Saccharomyces" & Experiment == "A")
  df1_B <- subset(df, Condition == "Single" & Species == "Saccharomyces" & Experiment == "B")
  
  df2_A <- subset(df, Condition == "Single" & Species == "Schixosachararomyces" & Experiment == "A")
  df2_B <- subset(df, Condition == "Single" & Species == "Schixosachararomyces" & Experiment == "B")
  
  df12_A <- subset(df, Condition == "Mixed" & Experiment == "A")
  df12_B <- subset(df, Condition == "Mixed" & Experiment == "B")
  
  list(
    D = 2,
    
    N1_A = nrow(df1_A),
    t1_A = df1_A$Age,
    y1_A = df1_A$Volume,
    
    N1_B = nrow(df1_B),
    t1_B = df1_B$Age,
    y1_B = df1_B$Volume,
    
    N2_A = nrow(df2_A),
    t2_A = df2_A$Age,
    y2_A = df2_A$Volume,
    
    N2_B = nrow(df2_B),
    t2_B = df2_B$Age,
    y2_B = df2_B$Volume,
    
    N12_A = nrow(df12_A),
    t12_A = df12_A$Age,
    y12_A = cbind(subset(df12_A, Species == "Saccharomyces", Volume),
                  subset(df12_A, Species == "Schixosachararomyces", Volume)),
    
    N12_B = nrow(df12_B),
    t12_B = df12_B$Age,
    y12_B = cbind(subset(df12_B, Species == "Saccharomyces", Volume),
                  subset(df12_B, Species == "Schixosachararomyces", Volume)),
    
    N_rep = 150,
    t_rep = 1:150
  )
}

data_stan <- format_data(df)

param <- c("r", "alpha", "sigma",
           "sigma_f0", "sigma_alpha",
           "k_uni", "k_multi",
           "f0",
           "f1_A", "f1_B", "f2_A", "f2_B", "f1_A", "f12_B",
           "y1_rep", "y2_rep", "y12_rep",
           "log_lik")

library(rstan)
# Parallel computing
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

if (run){
  fit <- stan(file = stan_code, data = data_stan, iter = n_it, chains = n_chains, pars = param, seed = seed,
              control = list(adapt_delta = 0.99, max_treedepth = 12))
  
  saveRDS(fit, res_file)
}else{
  fit <- readRDS(res_file)
}

# Explore solution -----------------------------------------------------------------

key_param <- c("r", "alpha", "sigma")

pairs(fit, pars = key_param)
plot(fit, pars = key_param, plotfun = "trace")
plot(fit, pars = c(key_param, "sigma_f0", "sigma_alpha"))

print(fit)

# shinystan::launch_shinystan(fit)

# Posterior Predictive Checks -------------------------------------------------------------------

# rep <- process_replications_density(data_stan, fit, maxVolume = 15)
rep <- process_replications_spaghetti(data_stan, fit, draws = 100)
plot_PPC(rep, df)

# Compare models ----------------------------------------------------------

library(loo)

log_lik1 <- extract_log_lik(fit, merge_chains = FALSE)
r_eff1 <- relative_eff(exp(log_lik1), chain_id = rep(1:dim(log_lik1)[2], each = dim(log_lik1)[1]))
(loo1 <- loo(log_lik1, r_eff = r_eff1))
# plot(loo1)

fit2 <- readRDS("Results/fit_LV12_exp2_v3.rds")
log_lik2 <- extract_log_lik(fit2, merge_chains = FALSE)
r_eff2 <- relative_eff(exp(log_lik2), chain_id = rep(1:dim(log_lik2)[2], each = dim(log_lik2)[1]))
(loo2 <- loo(log_lik2, r_eff = r_eff2))
comp <- compare(loo1, loo2)
print(comp) # negative elpd favors first model, to compared with SE
