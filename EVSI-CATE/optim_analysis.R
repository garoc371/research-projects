library(tidyverse)
library(data.table)
library(patchwork)

source("analysis_functions.R")
source("simulation_functions.R")
source("config.R")

annual_population <- 2e4
discount_rate <- 0.03
time_horizon <- 10

ttl_population <- sum(annual_population/(1 + discount_rate)^(0:(time_horizon-1)))
ttl_population <- round(ttl_population,0)

sample_sizes <- seq(100, 1000, by = 100)
results4_list_g8 <- vector("list", length(sample_sizes))
names(results4_list_g8) <- sample_sizes

results5_list_g8 <- vector("list", length(sample_sizes))
names(results5_list_g8) <- sample_sizes

costs_g5 <- c(5, 6, 8, 10, 20)*10
costs_g8 <- c(5, 5, 6, 6, 8, 8, 10, 20)*100
fixed_costs <- 1e6

for (n in sample_sizes) {
    file_res4 <- sprintf("data/results4/logs_subgroups_8/lhs_log_%d.csv", n)
    file_res5 <- sprintf("data/results5/logs_subgroups_8/lhs_log_%d.csv", n)

    results4_list_g8[[as.character(n)]] <- read_optim_results(
        file_res4,
        n,
        costs_g8,
        ttl_population,
        fixed_cost = fixed_costs
    )

    results5_list_g8[[as.character(n)]] <- read_optim_results(
        file_res5,
        n,
        costs_g8,
        ttl_population,
        fixed_cost = fixed_costs
    )
}

# Combine all results
final_results4_g8 <- bind_rows(results4_list_g8)
final_results5_g8 <- bind_rows(results5_list_g8)


final_results4_g8 <- final_results4_g8 %>%
    group_by(allocation, allocation_num, sample_size) %>%
    summarise(
        evsi = median(evsi),
        enbs = median(enbs),
        evsi_sd = median(evsi_sd)
    ) %>%
    ungroup()


final_results5_g8 <- final_results5_g8 %>%
    group_by(allocation, allocation_num, sample_size) %>%
    summarise(
        evsi = median(evsi),
        enbs = median(enbs),
        evsi_sd = median(evsi_sd)
    ) %>%
    ungroup()

max4_evsi_cate_g8 <- final_results4_g8 %>%
    group_by(sample_size) %>%
    slice_max(order_by = evsi, n = 1) %>%
    ungroup()

max4_evsi_cate_g8

max5_evsi_cate_g8 <- final_results5_g8 %>%
    group_by(sample_size) %>%
    slice_max(order_by = evsi, n = 1) %>%
    ungroup()

max5_evsi_cate_g8

max4_enbs_cate_g8 <- final_results4_g8 %>%
    group_by(sample_size) %>%
    slice_max(order_by = enbs, n = 1) %>%
    ungroup()

max4_enbs_cate_g8

max5_enbs_cate_g8 <- final_results5_g8 %>%
    group_by(sample_size) %>%
    slice_max(order_by = enbs, n = 1) %>%
    ungroup()

max5_enbs_cate_g8
#
#

# EVPPI for each subgroup-specific treatment effects
# read in the posterior treatment effects and INB

# utilize build-in functions from voi pacakge

library(voi)

res4_input_param_g8 <- read_rds(
    "data/results4/logs_subgroups_8/treatment_effects.RDs"
)
res5_input_param_g8 <- read_rds(
    "data/results5/logs_subgroups_8/treatment_effects.RDs"
)

# cor_eff_g5 <- plot_treatment_effects_correlation(input_param_g5)
# cor_eff_g8 <- plot_treatment_effects_correlation(input_param_g8)

# control_probs_g5 <- read_rds("data/results5/logs_subgroups_5/ctrl_probs.RDs")
# control_probs_g8 <- read_rds("data/results5/logs_subgroups_8/ctrl_probs.RDs")

res4_age_groups_g8 <- colnames(res4_input_param_g8)
res5_age_groups_g8 <- colnames(res5_input_param_g8)
INB_post_g8_res4 <- read_rds("data/results4/logs_subgroups_8/INB_post.RDs")
INB_post_g8_res5 <- read_rds("data/results5/logs_subgroups_8/INB_post.RDs")

evppi_g8_res4 <- compute_evppi(res4_input_param_g8, INB_post_g8_res4)
evppi_g8_res5 <- compute_evppi(res5_input_param_g8, INB_post_g8_res5)

p4_evppi_g8 <- plot_evppi(evppi_g8_res4, age_groups = res4_age_groups_g8)
p5_evppi_g8 <- plot_evppi(evppi_g8_res5, age_groups = res5_age_groups_g8)

# standard can be computed by re-fitting the bart model given inputs and
# outcomes. This could introduce additional Monte Carlo errors, and a more
# direct approach is to average over the evsi values given by corner candidate
# during the search

res4_evsi_corner_g8 <- final_results4_g8 %>%
    filter(allocation == "c(0, 0, 0, 0, 0, 0, 0, 1)") %>%
    dplyr::select(evsi, sample_size, enbs, evsi_sd) %>%
    mutate(measure = "Standard EVSI")

res5_evsi_corner_g8 <- final_results5_g8 %>%
    filter(allocation == "c(0, 0, 0, 0, 0, 0, 0, 1)") %>%
    dplyr::select(evsi, sample_size, evsi_sd, enbs) %>%
    mutate(measure = "Standard EVSI")


#############################################################################


# 
# comp_g5 <- max_evsi_cate_g5 %>%
#     dplyr::select(evsi, sample_size,enbs, evsi_sd) %>%
#     mutate(measure = rep("EVSI-CATE", length(sample_sizes))) %>%
#     bind_rows(df_evsi_corner_g5)
# 
# p_evsi_comp_g5 <- plot_evsi_comparison(comp_g5)
# 
# p_evppi_g5 + p_evsi_comp_g5
# 
# ggsave(filename = "data/results5/comb_voi_g5.png", width = 7, height = 5, units = "in", dpi = 300)


res4_comp_g8 <- max4_evsi_cate_g8 %>%
  dplyr::select(evsi, sample_size,enbs, evsi_sd) %>%
  mutate(measure = rep("EVSI-CATE", length(sample_sizes))) %>%
  bind_rows(res4_evsi_corner_g8)

p4_evsi_comp_g8 <- plot_evsi_comparison(res4_comp_g8)

res5_comp_g8 <- max5_evsi_cate_g8 %>%
    dplyr::select(evsi, sample_size,enbs, evsi_sd) %>%
    mutate(measure = rep("EVSI-CATE", length(sample_sizes))) %>%
    bind_rows(res5_evsi_corner_g8)

p5_evsi_comp_g8 <- plot_evsi_comparison(res5_comp_g8)

(p4_evppi_g8 + p4_evsi_comp_g8)/(p5_evppi_g8 + p5_evsi_comp_g8)



ggsave(filename = "data/comp_voi.png", width = 7, height = 7, units = "in", dpi = 300)


enbs_res4_g8 <- max4_enbs_cate_g8 %>%
  dplyr::select(evsi, sample_size, enbs) %>%
  mutate(measure = rep("EVSI-CATE", length(sample_sizes))) %>%
  bind_rows(res4_evsi_corner_g8)


p4_enbs_g8 <- plot_enbs_comparison(enbs_res4_g8)


enbs_res5_g8 <- max5_enbs_cate_g8 %>% 
  dplyr::select(evsi, sample_size, enbs) %>%
  mutate(measure = rep("EVSI-CATE", length(sample_sizes))) %>%
  bind_rows(res5_evsi_corner_g8)


p5_enbs_g8 <-plot_enbs_comparison(enbs_res5_g8)

(p4_evppi_g8 + p4_enbs_g8) / (p5_evppi_g8 + p5_enbs_g8)
ggsave(filename = "data/comp_enbs.png", width = 7, height = 7, units = "in", dpi = 300)




top_res_g8 <- final_results_g8 |>
    filter(sample_size == 1000) |>
    slice_max(order_by = evsi, n = 200) %>% 
  arrange(desc(evsi)) %>% 
  mutate(group = ((row_number() -1) %/% 50) +1)

top_res_list <- group_split(top_res_g8, group)

plot_list <- imap(top_res_list, ~ {
  
  start_rank <- (.y -1) * 50 + 1
  end_rank <- .y * 50
  
  plot_parallel_simplex(
    results_df = .x,
    age_groups = age_groups_g8,
    order_axes = c(8,7,6,5,1,3,4,2),
    line_scaling = 1.5
  ) + labs(title = sprintf("Allocation Strategies: top %d–%d", 
                           start_rank, end_rank))
})

wrap_plots(plot_list)
ggsave(filename = "data/results5/parcorrd_top200.png", width = 8, height = 6, units = "in", dpi = 300)
