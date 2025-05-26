parse_allocations <- function(s_strings) {
  allocations <- lapply(s_strings, function(s) {
    as.numeric(stringr::str_extract_all(s, "\\d*\\.?\\d+")[[1]])
  })
  do.call(rbind, allocations)
}


read_optim_results <- function(file_path, n, costs, n_target, fixed_cost) {
  dt <- fread(file_path)
  allocation_matrix <- parse_allocations(dt$S)

  tibble(
    allocation = dt$S,
    # Use the already parsed numeric allocations
    allocation_num = split(allocation_matrix, seq(nrow(allocation_matrix))),
    evsi = dt$evsi,
    evsi_sd = dt$evsi_sd,
    sample_size = dt$n_trial,
  ) %>%
    mutate(
      total_cost = map_dbl(allocation_num, ~ sum(.x * n * costs)),
      enbs = evsi * n_target - total_cost - fixed_costs
    )
}


plot_treatment_effects_correlation <- function(treatment_effects) {
  corr_df <- treatment_effects %>%
    cor() %>%
    round(2) %>%
    as.data.frame() %>%
    rownames_to_column("age_group_row") %>%
    pivot_longer(
      -age_group_row,
      names_to = "age_group_col",
      values_to = "correlation"
    )

  ggplot(
    corr_df,
    aes(x = age_group_row, y = age_group_col, fill = correlation)
  ) +
    geom_tile() +
    scale_fill_gradient2(
      low = "white",
      high = "#002855",
      mid = "white",
      midpoint = 0,
      limit = c(0, 1),
      name = "Correlation",
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = 10,
        barheight = 0.5,
        frame.colour = "white",
        ticks.colour = "white"
      )
    ) +
    geom_text(aes(label = sprintf("%.2f", correlation)), size = 3) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.spacing.x = unit(0.5, "cm"),
      axis.title = element_blank(),
      plot.margin = margin(t = 1, r = 1, b = 1, l = 1)
    ) +
    coord_fixed()
}


compute_evppi <- function(treatment_effects, INB_post) {
  input_param <- as.data.frame(treatment_effects)
  # Rename columns (e.g. for 8 subgroups) so that they are p1, p2, …, p8.
  input_param <- setNames(input_param, paste0("p", seq_len(ncol(input_param))))
  nb_output <- data.frame(ctrl = 0, trt = INB_post)
  evppi_res <- evppi(
    nb_output,
    input_param,
    pars = as.list(names(input_param)),
    method = "bart",
    k = 4,
    se = TRUE
  )
  evppi_res
}

plot_evppi <- function(evppi_results, age_groups) {
  evppi_results %>%
    mutate(age_group = age_groups) %>%
    ggplot(aes(x = age_groups, y = evppi)) +
    geom_col(fill = "#3B5873", width = 0.6, alpha = 0.6) +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = NULL, y = "EVPPI (£)", title = "") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "black"),
      plot.title = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 11)
    )
}

compute_evsi_standard <- function(
  sample_sizes,
  treatment_effects,
  control_prob,
  INB_post,
  subgroup_index,
  bart_k = 3,
  age_groups
) {
  evsi_values <- numeric(length(sample_sizes))
  for (i in seq_along(sample_sizes)) {
    S_vector <- rep(0, ncol(treatment_effects))
    S_vector[subgroup_index] <- 1
    trial_data <- lapply(seq_len(nrow(treatment_effects)), function(j) {
      gen_future_trial(
        N = sample_sizes[i],
        S = S_vector,
        theta_s = treatment_effects[j, ],
        p_control = control_prob[j, ],
        age_groups = age_groups
      )
    }) %>%
      bind_rows() %>%
      dplyr::select(ends_with(paste0("log_OR_g", subgroup_index)))

    bart_mod <- bart(x.train = trial_data, y.train = INB_post, k = bart_k)
    bart_fitted <- as.numeric(fitted(bart_mod))
    evsi_values[i] <- mean(pmax(bart_fitted, 0)) - pmax(mean(bart_fitted), 0)
  }
  evsi_values
}

make_comparison <- function(evsi_cate_df, standard_evsi_values, sample_sizes) {
  df_standard <- tibble(
    evsi = round(standard_evsi_values, 0),
    sample_size = sample_sizes,
    measure = "Standard EVSI"
  )
  # It is assumed that evsi_cate_df already contains a measure column (e.g. "EVSI-CATE").
  bind_rows(evsi_cate_df, df_standard)
}

plot_evsi_comparison <- function(combined_df) {
  combined_df %>%
    ggplot(aes(x = sample_size, y = evsi, colour = measure)) +
    geom_smooth(
      data = combined_df %>% filter(measure == "EVSI-CATE"),
      linewidth = 1,
      se = F
    ) +
    geom_smooth(
      data = combined_df %>% filter(measure == "Standard EVSI"),
      linewidth = 1,
      se = F
    ) +
    scale_colour_manual(
      values = c("EVSI-CATE" = "#9B4B9D", "Standard EVSI" = "#002855")
    ) +
    scale_x_continuous(
      breaks = seq(0, max(combined_df$sample_size), by = 100),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(
      labels = scales::comma,
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    labs(
      x = "Sample Size (n)",
      y = "Expected Value of Sample Information (£)",
      colour = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom")
}



plot_enbs_comparison <- function(enbs_df) {
  enbs_df %>%
    ggplot(aes(x = sample_size, y = enbs, colour = measure)) +
    geom_smooth(
      data = enbs_df %>% filter(measure == "EVSI-CATE"),
      linewidth = 1
      ) +
    geom_smooth(
      data = enbs_df %>% filter(measure == "Standard EVSI"),
      linewidth = 1
    ) +
    scale_colour_manual(
      values = c("EVSI-CATE" = "#9B4B9D", "Standard EVSI" = "#002855")
    ) +
    scale_x_continuous(
      breaks = seq(0, max(enbs_df$sample_size), by = 100),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(
      labels = scales::comma,
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    labs(
      x = "Sample Size (n)",
      y = "Expected Net Benefits of Sampling (£)",
      colour = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom")
}
