## check h2 on observed scale
library(tidyverse)
setwd("/project/nmancuso_8/elezhang/projects/jaxqtl/code/sim/params")

calc_h2obs <- function(V_a, V_re, V_disp=0, baseline_mu=0){
  tot_var = V_a + V_re + V_disp
  lamb = exp(baseline_mu + tot_var / 2.0)
  h2g_obs = lamb * V_a / (lamb * (exp(tot_var) - 1) + 1)
  return(h2g_obs)
}

num_cells <- c(10, 50, 100, 200, 500, 1000) # fix 100
maf <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5)

beta0 <- c(-3, -2, -1, 0) # fix 0
Va <- c(seq(0, 0.1, 0.02), 0.2, 0.5, 0.7, 1) # fix 0.1
Vre <- c(0, 0.05, 0.1, 0.2, 0.5, 0.7, 1, 4) # fix 0.2

expand.grid(Va = Va, Vre = Vre) %>% 
  mutate(h2obs = calc_h2obs(Va, Vre)) %>% 
  ggplot(aes(x = Va, y = h2obs)) +
  geom_point(size=0.7) + geom_line(linewidth=0.5)+ geom_hline(yintercept = 0.1, color = "blue", linetype="dashed")+
  facet_grid(. ~ Vre) +
  axis_theme + theme(axis.text.x = element_text(angle=30, size = 6))

ggsave2(paste0(plots_dir, "h2obs_Va_Vre.png"),
        width = onecol_width*2, height = onecol_height, units = 'cm', dpi = 300)

libsize <- 1 # can use empirical library size from onek1k


####### examine effect of params_Va ########
params_file <- "params_Va"
num_cells <- 100
maf <- 0.2
libsize <- 1
beta0 <- 0
Va <- c(seq(0, 0.1, 0.02), 0.2, 0.5, 0.7, 1)
Vre <- c(0.2, 1, 4)
rep <- 200

params <- expand.grid(num_cells = num_cells,
                      libsize = libsize,
                      maf = maf, 
                      beta0 = beta0,
                      Va = Va, 
                      Vre = Vre) %>% 
  mutate(h2obs = calc_h2obs(Va, Vre))

n_scenario <- nrow(params)

params %>% 
  write_tsv(params_file, col_names=F)

# write params for sqtl
params %>% 
  slice(rep(1:n(), each=rep)) %>%
  mutate(sim_idx = rep(1:n_scenario, each=rep),
         idx = rep(1:rep, n_scenario)) %>%
  write_tsv(paste0(params_file, "_sqtl"), col_names=F)

####### examine effect of params_Vre ########
params_file <- "params_Vre"
num_cells <- 100
maf <- 0.2
libsize <- 1
beta0 <- 0
Va <- c(0.05, 0.1)
Vre <- c(0, 0.05, 0.1, 0.2, 0.5, 0.7, 1, 4)
rep <- 200

params <- expand.grid(num_cells = num_cells,
            libsize = libsize,
            maf = maf, 
            beta0 = beta0,
            Va = Va, 
            Vre = Vre) %>% 
  mutate(h2obs = calc_h2obs(Va, Vre))

n_scenario <- nrow(params)

params %>% 
  write_tsv(params_file, col_names=F)

# write params for sqtl
params %>% 
  slice(rep(1:n(), each=rep)) %>%
  mutate(sim_idx = rep(1:n_scenario, each=rep),
         idx = rep(1:rep, n_scenario)) %>%
  write_tsv(paste0(params_file, "_sqtl"), col_names=F)
  

####### examine type I error of params_type1 ########
params_file <- "params_type1"
num_cells <- 100
maf <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5) # c(0.05, 0.2, 0.5)
libsize <- 1
beta0 <- 0
Va <- 0
Vre <- c(0, 0.05, 0.1, 0.2, 0.5, 0.7, 1, 4)
rep <- 200

params <- expand.grid(num_cells = num_cells,
                      libsize = libsize,
                      maf = maf, 
                      beta0 = beta0,
                      Va = Va, 
                      Vre = Vre) %>% 
  mutate(h2obs = calc_h2obs(Va, Vre))

n_scenario <- nrow(params)

params %>% 
  write_tsv(params_file, col_names=F)

# write params for sqtl
params %>% 
  slice(rep(1:n(), each=rep)) %>%
  mutate(sim_idx = rep(1:n_scenario, each=rep),
         idx = rep(1:rep, n_scenario)) %>%
  write_tsv(paste0(params_file, "_sqtl"), col_names=F)


####### examine number of cells: params_num_cells ########
params_file <- "params_num_cells"
num_cells <- c(10, 50, 100, 200, 500, 1000)
maf <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5)  # c(0.2)
libsize <- 1
beta0 <- 0
Va <- 0.05
Vre <- c(0, 0.05, 0.1, 0.2, 0.5, 0.7, 1, 4)
rep <- 200

params <- expand.grid(num_cells = num_cells,
                      libsize = libsize,
                      maf = maf, 
                      beta0 = beta0,
                      Va = Va, 
                      Vre = Vre) %>% 
  mutate(h2obs = calc_h2obs(Va, Vre))

n_scenario <- nrow(params); n_scenario

params %>% 
  write_tsv(params_file, col_names=F)

# write params for sqtl
params %>% 
  slice(rep(1:n(), each=rep)) %>%
  mutate(sim_idx = rep(1:n_scenario, each=rep),
         idx = rep(1:rep, n_scenario)) %>%
  write_tsv(paste0(params_file, "_sqtl"), col_names=F)


####### examine baseline mean: params_beta0 ########
params_file <- "params_beta0"
num_cells <- 100
maf <- 0.2
libsize <- 1
beta0 <- c(-3, -2, -1, 0)
Va <- 0.05
Vre <- c(0, 0.05, 0.1, 0.2, 0.5, 0.7, 1, 4)
rep <- 200

params <- expand.grid(num_cells = num_cells,
                      libsize = libsize,
                      maf = maf, 
                      beta0 = beta0,
                      Va = Va, 
                      Vre = Vre) %>% 
  mutate(h2obs = calc_h2obs(Va, Vre))

n_scenario <- nrow(params); n_scenario

params %>% 
  write_tsv(params_file, col_names=F)

# write params for sqtl
params %>% 
  slice(rep(1:n(), each=rep)) %>%
  mutate(sim_idx = rep(1:n_scenario, each=rep),
         idx = rep(1:rep, n_scenario)) %>%
  write_tsv(paste0(params_file, "_sqtl"), col_names=F)


####### Round2: examine baseline mean: params_beta0 ########
params_file <- "params_Va_Vre_cells_beta0_more"
num_cells <- c(10, 50, 100, 200, 500)
maf <- 0.2
libsize <- 1
beta0 <- c(-2, -1, 0) # c(-20, -15, -10, -5)
Va <- c(seq(0, 0.1, 0.02), 0.2, 0.5, 1)
Vre <- c(0, 0.1, 0.2, 0.7, 1, 2)
rep <- 200

params <- expand.grid(num_cells = num_cells,
                      libsize = libsize,
                      maf = maf, 
                      beta0 = beta0,
                      Va = Va, 
                      Vre = Vre) %>% 
  mutate(h2obs = calc_h2obs(Va, Vre))

n_scenario <- nrow(params); n_scenario

params %>% distinct(Va, Vre, h2obs)

params %>% 
  write_tsv(params_file, col_names=F)

# write params for sqtl
params %>% 
  slice(rep(1:n(), each=rep)) %>%
  mutate(sim_idx = rep(1:n_scenario, each=rep),
         idx = rep(1:rep, n_scenario)) %>%
  write_tsv(paste0(params_file, "_sqtl"), col_names=F)

####### Round3: bulk ########
params_file <- "params_bulk_Va_beta0_alpha_libsizefix"
num_cells <- 100
maf <- c(0.2)
libsize <- 1
rep <- 300

maf <- c(0.05, 0.2, 0.5)
libsize <- 1
beta0 <- c(-3, -2, -1, 0, 0.5, 1, 3, 5)
Va <- c(0, 0.01, 0.05, 0.1, 0.2, 0.5, 1)
alpha <- c(0.01, 0.1, 0.2, 0.3, 0.6, 1)
N <- c(50, 100, 200, 500, 1000)

params <- expand.grid(libsize = libsize,
                      maf = maf,
                      beta0 = beta0,
                      Va = Va,
                      alpha = alpha,
                      N = N)

n_scenario <- nrow(params); n_scenario

# fix library size (1)
params %>% write_tsv(params_file, col_names=F)


# different library size (most recent)
params_file <- "params_bulk_celltype_libsize"
libsize <- "./CD4_ET.libsize.tsv"
libsize <- c("CD4_NC", "B_IN", "Plasma")
rep <- 300

maf <- 0.2
# CD4_ET, 0.07, 0.16, 0.96
beta0 <- c(-16, -15, -14, -13, -12, -11, -10)
# Va <- c(0, 0.01, 0.05, 0.1, 0.2, 0.5, 1)
Va <- c(0, 0.01, 0.05, 0.1, 0.5)

alpha <- c(0.01, 0.1, 0.2, 0.3, 0.6, 1)
N <- c(50, 100, 200, 500, 1000)

params_large_libsize <- expand.grid(libsize = c("CD4_NC"),
                                    maf = maf,
                                    beta0 = c(-18, -16, -15, -14, -13, -12, -11),
                                    Va = Va,
                                    alpha = alpha,
                                    N = N)

params_medium_libsize <- expand.grid(libsize = c("B_IN"),
                                    maf = maf,
                                    beta0 = c(-16, -15, -14, -13, -12, -11, -10),
                                    Va = Va,
                                    alpha = alpha,
                                    N = N)

# params_large_libsize <- expand.grid(libsize = c("CD4_NC", "B_IN"),
#                                     maf = maf,
#                                     beta0 = beta0,
#                                     Va = Va,
#                                     alpha = alpha,
#                                     N = N)

params_small_libsize <- expand.grid(libsize = c("Plasma"),
                                    maf = maf,
                                    beta0 = c(-16, -14, -12, -11, -10, -9, -8),
                                    Va = Va,
                                    alpha = alpha,
                                    N = N)

params %>% write_tsv(params_file, col_names=F)

# diff MAF
params_file <- "params_bulk_celltype_libsize_maf"
maf <- c(0.2)
libsize <- c("B_IN")
rep <- 300

maf <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
beta0 <- c(-16, -15, -14, -13, -12, -11, -10)
Va <- c(0, 0.05, 0.1)
alpha <- 0.3
N <- 1000

params_maf <- expand.grid(libsize = c("B_IN"),
                          maf = maf,
                          beta0 = beta0,
                          Va = Va,
                          alpha = alpha,
                          N = N)


####### Round3: bulk ########
# params_file <- "params_bulk_Va_beta0_alpha"
# num_cells <- 100
# maf <- 0.2
# libsize <- 1
# beta0 <- c(-5, -2, -1, 0, 1, 2, 3, 4, 5)
# Va <- c(0, 0.02, 0.05, 0.06, 0.08, 0.1, 0.2, 0.5, 1)
# Vre <- 0
# alpha <- c(0.01, 0.1, 0.2, 0.6, 1)
# rep <- 200
# 
# params <- expand.grid(num_cells = num_cells,
#                       libsize = libsize,
#                       maf = maf, 
#                       beta0 = beta0,
#                       Va = Va, 
#                       Vre = Vre,
#                       alpha = alpha) %>% 
#   mutate(h2obs = calc_h2obs(Va, Vre))
# 
# n_scenario <- nrow(params); n_scenario
# 
# params %>% distinct(Va, Vre, h2obs)
# 
# params %>% 
#   write_tsv(params_file, col_names=F)

##### collect unfinished sqtl params file ###
# gather: 
files <- c("params_Vre", "params_num_cells", "params_type1", "params_beta0")

allparams <- data.frame()

for (i in files){
  params <- read_tsv(paste0(i, "_sqtl"), F)
  colnames(params) <- c("num_cells", "libsize", "maf", "beta0", "Va", "Vre", "h2obs",
                        "sim_idx", "idx")
  allfiles <- list.files(paste0("../output/", i))
  sqtl_files <- allfiles[grepl("*_sqtl.rda", allfiles)]
  
  newparams <- tibble(params_name = gsub("_sqtl.rda", "", sqtl_files)) %>% 
    mutate(params_name = gsub("sim|pheno", "", params_name)) %>% 
    separate(params_name, into=c("sim_idx", "idx"), sep="\\.") %>% 
    mutate(sim_idx = as.integer(sim_idx),
           idx = as.integer(idx)) 
  
  df <- params %>% 
    anti_join(newparams, by = c("sim_idx", "idx")) %>% 
    mutate(name = i)
  
  allparams <- bind_rows(allparams, df)
}
