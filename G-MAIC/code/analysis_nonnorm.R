library(dplyr)
library(rsimsum)
library(ggplot2)
library(viridis)
library(patchwork)


N_sim = 2000 # 2000 or 5000, 2000 MC replications for the main analysis
N_AC <- c(100,200,600)
N_BC <- 600
b_trt <- log(0.25)
b_EM <- -log(0.6)
b_1 <- -log(0.8)
X1_AC = 0
var1_AC = 0.2
var1_BC = 0.4
X1_BC = c(0.275,0.4,0.6)
event_rate = 0.4

param.combinations <- expand.grid(N_AC=N_AC, X1_BC=X1_BC)
pc <- param.combinations

scenarios <- nrow(pc) # number of simulation scenarios

file_names_nonnorm <- list.files("results/non_norm3")
scenarios_names <- vector(mode = "character", length = scenarios)
methods <- c("bucher", "maic", "bayes_gcomp", "ce", "mis_bayes_gcomp")
n_methods <- length(methods)

for (i in 1:scenarios){
  file.id <- paste0("N_AC", pc$N_AC[i], "X1_BC", pc$X1_BC[i])
  scenarios_names[i] <- file.id
}


mk_dataset <- function(index){
  # index - the index of chosen scenaio
  file_idx <- grep(scenarios_names[index], file_names_nonnorm)
  list_of_file_names <- file_names_nonnorm[file_idx]
  obj_names <- lapply(paste0("results/non_norm3/", list_of_file_names), load, .GlobalEnv)
  res1 <- bind_rows(bucher.results, maic.results, bayes_gcomp.results, CE.results, mis_bayes_gcomp.results) %>% 
    mutate(se = sqrt(var),
           dataset = rep(1:N_sim, n_methods),
           method = c(rep("Bucher", N_sim), rep("MAIC", N_sim),
                      rep("Bayes_Gcomp", N_sim), rep("G-MAIC", N_sim),
                      rep("(Mis)Bayes_Gcomp", N_sim)),
           N_AC = rep(pc$N_AC[index], n_methods*N_sim),
           X1_BC = rep(pc$X1_BC[index], n_methods*N_sim))
  
}



results_all_nonnorm <- vector(mode = "list", length = scenarios)
for(i in 1:scenarios){
  results_all_nonnorm[[i]]<- mk_dataset(i)
}

sim_all_nonnorm <- bind_rows(results_all_nonnorm)
s_nonnorm <- simsum(sim_all_nonnorm, estvarname = "mean", se = "se", methodvar = "method",
             by = c("N_AC", "X1_BC"), x = T, true = 0, ref = "Bucher", dropbig = T, control = list(dropbig.max = 20))




non_norm_bias <-autoplot(summary(s_nonnorm), type = "lolly", stats = "bias") +
  ggtitle("Estimation Bias across scenarios, Non-normal covariate structure")+
  xlab("Estimation Bias") +
  ylab("")+
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 14, hjust = 0,colour = "#4B384C", face = "bold"),
    axis.title.x = element_text(size = 12,colour = "#4B384C", face = "bold"),
    plot.title.position = "plot"
  ) +
  scale_color_viridis(discrete = TRUE)
non_norm_bias
# ggsave("plots/non_norm/bias.png")  
  


non_norm_cover <-autoplot(summary(s_nonnorm), type = "lolly", stats = "cover") + 
  ggtitle("Coverage across scenarios under non-normal covariate structure")+
  xlab("Coverage") +
  ylab("Method")+
  theme_bw()+
  theme(
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 14, hjust = 0.5),
    axis.title = element_text(size = 12)
  ) +
  scale_color_viridis(discrete = TRUE)
non_norm_cover
# ggsave("plots/non_norm/coverage.png") 

non_norm_empse <-autoplot(summary(s_nonnorm), type = "lolly", stats = "empse") + 
  ggtitle("Empirical Std. Error across scenarios,Non-normal covariate structure")+
  xlab("Empirical Std Error") +
  ylab("")+
  theme_bw()+
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 14, hjust = 0, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    plot.title.position = "plot"
  ) +
  scale_color_viridis(discrete = TRUE)
non_norm_empse
# ggsave("plots/non_norm/empse.png") 






se_nonnorm_summary <- tidy(summary(s_nonnorm, stats = c("empse", "modelse")))
library(kableExtra)
se_nonnorm_summary$stat <- rep(c("Empirical SE", "Model SE"),9)
kable(se_nonnorm_summary, format = "latex",
      col.names = c("Statistics", "Sample size", "IPD location", "Bucher", "Bayes Gcomp", "G-MAIC",
                    "MAIC", "Bayes Gcomp(mis-specified)")) %>% 
  kable_styling(bootstrap_options = c("striped", "hover"))


