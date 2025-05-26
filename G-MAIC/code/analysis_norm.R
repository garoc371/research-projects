library(dplyr)
library(rsimsum)
library(ggplot2)
library(patchwork)
library(knitr)
library(forcats)
library(viridis)
library(patchwork)

N_sim = 2000   # 2000 or 5000, 2000 MC replications for the main analysis
N_AC <- c(100,200, 600)
N_BC <- 600
b_trt <- log(0.25)
b_EM <- -log(0.6)
b_1 <- -log(0.8)
sdX = 0.4  # sd for both trials in the Normal scenario
meanX_BC = 0.6
rhoX = 0.25
meanX_AC = c(0.45,0.375,0.25)
event_rate = 0.4

param.combinations <- expand.grid(N_AC=N_AC, meanX_AC=meanX_AC)
pc <- param.combinations

scenarios <- nrow(pc) # number of simulation scenarios


file_names_norm <- list.files("results/normal3")
scenarios_names <- vector(mode = "character", length = scenarios)
methods <- c("bucher", "maic", "bayes_gcomp", "ce","mis_bayes_gcomp")
n_methods <- length(methods)

for (i in 1:scenarios){
  file.id <- paste0("N_AC", pc$N_AC[i], "meanX_AC", pc$meanX_AC[i])
  scenarios_names[i] <- file.id
}



mk_dataset <- function(index){
  # index - the index of chosen scenaio
  file_idx <- grep(scenarios_names[index], file_names_norm)
  list_of_file_names <- file_names_norm[file_idx]
  obj_names <- lapply(paste0("results/normal3/", list_of_file_names), load, .GlobalEnv)
  res1 <- bind_rows(bucher.results, maic.results, bayes_gcomp.results, CE.results, mis_bayes_gcomp.results) %>% 
      mutate(se = sqrt(var),
             dataset = rep(1:N_sim, n_methods),
             method = c(rep("Bucher", N_sim), rep("MAIC", N_sim),
                        rep("Bayes_Gcomp", N_sim), rep("G-MAIC", N_sim),
                        rep("(Mis)Bayes_Gcomp", N_sim)),
             N_AC = rep(pc$N_AC[index], n_methods*N_sim),
             meanX_AC = rep(pc$meanX_AC[index], n_methods*N_sim))
    
}

results_all_norm <- vector(mode = "list", length = scenarios)

for(i in 1:scenarios){
  results_all_norm[[i]]<- mk_dataset(i)
}

sim_all_norm <- bind_rows(results_all_norm)
s1 <- simsum(sim_all_norm, estvarname = "mean", se = "se", methodvar = "method",
             by = c("N_AC", "meanX_AC"), x = T, true = 0, ref = "Bucher")

s2 <- summary(s1) %>% tidy(stats = "bias")


norm_bias <- autoplot(summary(s1), type = "lolly", stats = "bias") + 
  xlab("Estimation Bias") +
  ylab("") +
  ggtitle("Bias across scenarios\n Multivariate Normal covariate structure") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size = 14, hjust = 0.5,colour = "#003D4C", face = "bold"),
    axis.title.x = element_text(size = 12,colour = "#003D4C", face = "bold"),
    plot.title.position = "plot"
  ) +
  scale_color_viridis(discrete = TRUE)

norm_bias
# ggsave("plots/normal/bias.png", width = 5, height = 8)

norm_cover <-autoplot(summary(s1), type = "lolly", stats = "cover") + 
  xlab("Coverage") +
  ylab("Method")+
  theme_bw()+
  ggtitle("Coverage across scenarios\n under multivariate Normal covariate structure") +
  theme(
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 14, hjust = 0.5,colour = "#8DB9CA", face = "bold"),
    axis.title = element_text(size = 12, colour = "#8DB9CA", face = "bold")
  ) +
  scale_color_viridis(discrete = TRUE)

norm_cover


# ggsave("plots/normal/coverage.png", width = 5, height = 8)


norm_empse <-autoplot(summary(s1), type = "lolly", stats = "empse") + 
  xlab("Empirical Std. Errors") +
  ylab("") +
  ggtitle("Empirical Std. Errors across scenarios\n Multivariate Normal covariate structure") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size = 14, hjust = 0,colour = "#8DB9CA", face = "bold"),
    axis.title.x = element_text(size = 12,colour = "#8DB9CA", face = "bold"),
    plot.title.position = "plot"
  ) +
  scale_color_viridis(discrete = TRUE)

norm_empse

# ggsave("plots/normal/norm_empse.png")





norm_se_sum <- tidy(summary(s1, stats = c("empse", "modelse")))
library(kableExtra)
norm_se_sum$stat <- rep(c("Empirical SE", "Model SE"),9)
kable(norm_se_sum, format = "latex",
      col.names = c("Statistics", "Sample size", "IPD location", "Bucher", "Bayes Gcomp", "G-MAIC",
                    "MAIC", "Bayes Gcomp(mis-specified)")) %>% 
  kable_styling(bootstrap_options = c("striped", "hover"))


