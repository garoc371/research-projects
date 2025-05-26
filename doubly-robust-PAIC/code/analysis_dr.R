library(tidyverse)
library(rsimsum)
library(here)

setwd(here("results2"))

N_sim = 2000
N_sample <- c(100,200,600)
b_trt <- log(0.25)
b_EM <- -log(0.6)
b_1 <- -log(0.8)
sdX_AC <- 0.6
X1_var_target <- 0.3
X1_BC <- c(0.2, 0.5) 
event_rate <- 0.4
meanX_trial <- 0

pc <- expand.grid(N_sample = N_sample, X1_BC = X1_BC)



file.names <- list.files()

scenario.names <- vector(mode = "character", length = nrow(pc))

for(i in 1:nrow(pc)){
  file.id <- paste0("N_sample_", pc$N_sample[i], "_X1_BC_", pc$X1_BC[i])
  scenario.names[i] <- file.id
}

# base case analysis assuming marginalizing over correct marginals
# function to extract estimate based on type of target population from
# regression based methods
# margin - correct marginals
# oracle - true ALD population
# mis-specified - mis-specifying one marginal (using normal instead of gamma)


extract_target_results <- function(results.all, target = "margin"){
  
  est_name <- paste0("est_", target)
  var_name <- paste0("var_", target)
  
  results <- list(mean = results.all[[est_name]], 
                  var = results.all[[var_name]])
  
  return(results)
  
}

methods_to_extract <- c("aipw_ml", "bayes_gcomp", "tmle", "weighted_dr")


process_scenario <- function(scenario_name, file_names, target = "margin") {
  scenario_results <- list()
  
  # Find files for this scenario
  scenario_files <- file.names[grep(scenario_name, file.names, ignore.case = TRUE)]
  
  # Process methods that need extraction
  for (method in methods_to_extract) {
    method_file <- scenario_files[grep(method, scenario_files, ignore.case = TRUE)]
      # Load the file
    obj_name <- load(method_file)
      
      # Extract base results
    results_object <- get(obj_name)
    base_results <- lapply(results_object, function(x) extract_target_results(x, target = target))
      
      # Store in the list
    scenario_results[[method]] <- base_results
    
  }
  
  # Handle 'maic' separately
  maic_file <- scenario_files[grep("maic", scenario_files, ignore.case = TRUE)]
  obj_name <- load(maic_file)
  scenario_results[["maic"]] <- get(obj_name)  # Store maic results as is
  
  
  # Convert results to a data frame
  results_df <- bind_rows(scenario_results)
  method <- names(scenario_results)
  
  # Extract sample size and BC location from scenario name
  sample_size <- as.numeric(gsub(".*N_sample_([0-9]+).*", "\\1", scenario_name))
  bc_location <- as.numeric(gsub(".*BC_([0-9.]+).*", "\\1", scenario_name))
  
  # Add extra columns
  results_df <- results_df %>%
    mutate(
      method = rep(method, each = 2000),
      se = sqrt(var),
      dataset = rep(1:N_sim, length(unique(method))),
      N_sample = sample_size,
      BC_location = bc_location
    )
  
  return(results_df)
}

library(kableExtra)
library(dplyr)

create_performance_table <- function(data, format = "html",
                                     metric_full_name,
                                     column_names = c("N_sample", "BC_location", "MAIC", "AIPW", 
                                                      "OM", "TMLE", "Weighted DR"),
                                     methods = c("MAIC", "AIPW", "Bayes_gcomp", "TMLE", "Weighted DR"),
                                     method_full_names = c(
                                       "MAIC: Matching-Adjusted Indirect Comparison",
                                       "AIPW: Augmented Inverse Probability Weighting",
                                       "OM: Outcome modelling using Bayesian G-computation",
                                       "TMLE: Targeted Maximum Likelihood Estimation",
                                       "Weighted DR: Doubly Robust estimation with weighted regressions"
                                     )) {
  
  # Prepare header vectors
  header1 <- c(" " = 2, "Method" = length(methods))
  
  header2 <- setNames(length(column_names), 
                      paste(metric_full_name, "Comparison Across Methods"))
  
  data %>%
    kbl(format = format, 
        col.names = column_names,
        caption = NULL) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                  full_width = FALSE) %>%
    add_header_above(header1) %>%
    add_header_above(header2, 
                     bold = TRUE, 
                     color = "black", 
                     background = "white", 
                     font_size = 20) %>%
    pack_rows("BC_location = 0.2", 1, 3) %>%
    pack_rows("BC_location = 0.5", 4, 6) %>%
    footnote(
      general = sprintf("<span style='font-size: 0.8em;'>%s values are reported for each method and scenario combination.</span>", metric_full_name),
      number = sapply(method_full_names, function(x) sprintf("<span style='font-size: 0.8em;'>%s</span>", x)),
      general_title = "",
      number_title = "",
      escape = FALSE
    )
}


sim_all <- list()
for (scenario in scenario.names){
  sim_all[[scenario]] <- process_scenario(scenario, file.names, target = "margin")
}


results_all <- bind_rows(sim_all)

s1 <- simsum(results_all, estvarname = "mean", se = "se", methodvar = "method",
             by = c("N_sample", "BC_location"), x = T, true = 0, ref = "maic",dropbig = T, control = list(dropbig.max = 20))


autoplot(summary(s1), stats = "cover", type = "lolly")
autoplot(summary(s1), stats = "empse", type = "lolly")
autoplot(summary(s1), stats = "modelse", type = "lolly")


sim_oracle <- list()
for (scenario in scenario.names){
  sim_oracle[[scenario]] <- process_scenario(scenario, file.names, target = "oracle")
}

results_oracle <- bind_rows(sim_oracle)
s_oracle <- simsum(results_oracle, estvarname = "mean", se = "se", methodvar = "method",
             by = c("N_sample", "BC_location"), x = T, true = 0, ref = "maic",dropbig = T, control = list(dropbig.max = 20))

# sim_misspecified <- list()
# for (scenario in scenario.names){
#   sim_misspecified[[scenario]] <- process_scenario(scenario, file.names, target = "misspecified")
# }
# 
# results_misspecified <- bind_rows(sim_misspecified)
# s_misspecified <- simsum(results_misspecified, estvarname = "mean", se = "se", methodvar = "method",
#                    by = c("N_sample", "BC_location"), x = T, true = 0, ref = "maic",dropbig = T, control = list(dropbig.max = 20))
# 

autoplot(summary(s_oracle), stats = "cover", type = "lolly")
autoplot(summary(s_oracle), stats = "empse", type = "lolly")
autoplot(summary(s_oracle), stats = "modelse", type = "lolly")


library(patchwork)

bias_oracle <- tidy(s_oracle) %>% filter(stat == "bias") %>% pull(est)
bias_margin <- tidy(s1) %>% filter(stat == "bias") %>% pull(est)
# bias_misspecified <- tidy(s_misspecified) %>% filter(stat == "bias") %>% pull(est)

# Determine the overall range of bias values
bias_range <- range(c(bias_oracle, bias_margin), na.rm = TRUE)
# bias_range2 <- range(c(bias_misspecified, bias_margin), na.rm = TRUE)
# Add a small buffer to the range for better visualization
bias_range <- bias_range + c(-1, 1) * 0.05 * diff(bias_range)
# bias_range2 <- bias_range2 + c(-1, 1) * 0.05 * diff(bias_range2)

bias_oracle <- autoplot(summary(s_oracle), stats = "bias", type = "lolly")+
  coord_cartesian(xlim = bias_range) +
  labs(title = "Bias w/ true ALD population")
bias_margin <- autoplot(summary(s1), stats = "bias", type = "lolly") +
  coord_cartesian(xlim = bias_range) +
  labs(title = "Bias w/ correct marginals")
# bias_misspecified <- autoplot(summary(s_misspecified), stats = "bias", type = "lolly")+
#   coord_cartesian(xlim = bias_range2) +
#   labs(title = "Misspecified Bias")

bias_oracle /bias_margin
# bias_margin / bias_misspecified


empse_oracle <- tidy(summary(s_oracle), stats = "empse")
empse_oracle <- empse_oracle %>%
  select(N_sample, BC_location, method, est) %>%
  pivot_wider(
    names_from = method,
    values_from = est,
    names_prefix = "empse_"
  ) %>%
  mutate(across(starts_with("empse_"), ~round(., 3))) %>%
  arrange(BC_location, N_sample) %>% 
  
  create_performance_table(., format = "latex",
                           metric_full_name = "Empirical Standard Errors")

empse_oracle


empse_margin <- tidy(summary(s1), stats = "empse")
empse_margin <- empse_margin %>%
  select(N_sample, BC_location, method, est) %>%
  pivot_wider(
    names_from = method,
    values_from = est,
    names_prefix = "empse_"
  ) %>%
  mutate(across(starts_with("empse_"), ~round(., 3))) %>%
  arrange(BC_location, N_sample) %>% 
  
  create_performance_table(., format = "latex",
                           metric_full_name = "Empirical Standard Errors")

empse_margin

modelse_margin <- tidy(summary(s1), stats = "modelse")
modelse_margin <- modelse_margin %>%
  select(N_sample, BC_location, method, est) %>%
  pivot_wider(
    names_from = method,
    values_from = est,
    names_prefix = "empse_"
  ) %>%
  mutate(across(starts_with("empse_"), ~round(., 3))) %>%
  arrange(BC_location, N_sample) %>% 
  
  create_performance_table(., format = "html",
                           metric_full_name = "Model Standard Errors")

modelse_margin

relerror_margin <- tidy(summary(s1), stats = "relerror")
relerror_margin <- relerror_margin %>%
  select(N_sample, BC_location, method, est) %>%
  pivot_wider(
    names_from = method,
    values_from = est,
    names_prefix = "relerror_"
  ) %>%
  mutate(across(starts_with("relerror_"), ~round(., 3))) %>%
  arrange(BC_location, N_sample) %>% 
  
  create_performance_table(., format = "latex",
                           metric_full_name = "Relative Errors in SE")

relerror_margin


mse_oracle <- tidy(summary(s_oracle), stats = "mse")
mse_oracle <- mse_oracle %>%
  select(N_sample, BC_location, method, est) %>%
  pivot_wider(
    names_from = method,
    values_from = est,
    names_prefix = "mse_"
  ) %>%
  mutate(across(starts_with("mse_"), ~round(., 3))) %>%
  arrange(BC_location, N_sample) %>% 
  
  create_performance_table(., format = "html",
                           metric_full_name = "Mean Squared Errors")

mse_oracle

empse_margin <- tidy(summary(s1), stats = "empse")
empse_margin <- empse_margin %>%
  select(N_sample, BC_location, method, est) %>%
  pivot_wider(
    names_from = method,
    values_from = est,
    names_prefix = "empse_"
  ) %>%
  mutate(across(starts_with("empse_"), ~round(., 3))) %>%
  arrange(BC_location, N_sample) %>% 
  
  create_performance_table(., format = "html",
                           metric_full_name = "Empirical Standard Errors")

empse_margin


bec_margin <- tidy(summary(s1), stats = "becover")
bec_margin <- bec_margin %>%
  dplyr::select(N_sample, BC_location, method, est) %>%
  pivot_wider(
    names_from = method,
    values_from = est,
    names_prefix = "bec_"
  ) %>%
  mutate(across(starts_with("bec_"), ~round(., 3))) %>%
  arrange(BC_location, N_sample) %>% 
  
  create_performance_table(., format = "latex",
                           metric_full_name = "Bias eliminated coverage")

bec_margin

cover_margin <- tidy(summary(s1), stats = "cover")
cover_margin <- cover_margin %>%
  dplyr::select(N_sample, BC_location, method, est) %>%
  pivot_wider(
    names_from = method,
    values_from = est,
    names_prefix = "cov_"
  ) %>%
  mutate(across(starts_with("cov_"), ~round(., 3))) %>%
  arrange(BC_location, N_sample) %>% 
  
  create_performance_table(., format = "latex",
                           metric_full_name = "Coverage")

cover_margin


mse_margin <- tidy(summary(s1), stats = "mse")
mse_margin <- mse_margin %>%
  dplyr::select(N_sample, BC_location, method, est) %>%
  pivot_wider(
    names_from = method,
    values_from = est,
    names_prefix = "mse_"
  ) %>%
  mutate(across(starts_with("mse_"), ~round(., 3))) %>%
  arrange(BC_location, N_sample) %>% 
  
  create_performance_table(., format = "html",
                           metric_full_name = "Mean squared errors")

mse_margin

#############################################################################

#zip plot

z1 <- results_all %>% 
  filter(BC_location==0.2 & method != "maic") %>% 
  simsum(., estvarname = "mean", se = "se", methodvar = "method",
             by = c("N_sample", "BC_location"), x = T, true = 0, 
         dropbig = T, control = list(dropbig.max = 20)) %>% 
  autoplot(., type = "zip", zoom = 0.2)


z2 <- results_all %>% 
  filter(BC_location==0.5 & method != "maic") %>% 
  simsum(., estvarname = "mean", se = "se", methodvar = "method",
         by = c("N_sample", "BC_location"), x = T, true = 0, 
         dropbig = T, control = list(dropbig.max = 10)) %>% 
  autoplot(., type = "zip", zoom = 0.2)


z1 + z2
