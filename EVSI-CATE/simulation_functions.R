# simulation_functions.R

# Generate subgroup proportions
# gen_group_proportions <- function(n) {
#   if (n < 1) stop("n must be a positive integer")
#   if (n == 1) return(1)
#   weights <- (n - 1):1
#   remaining_prob <- 1 - 1 / (2 * n)
#   scaling_factor <- remaining_prob / sum(weights)
#   c(weights * scaling_factor, 1 / (2 * n))
# }

gen_group_proportions <- function(n, p = 0.5) {
  if(n < 1) stop("n must be a positive integer")
  weights <- (n:1)^p            # Highest weight for youngest age group
  sample_props <- weights / sum(weights)
  return(sample_props)
}


# Create integer age breaks given an age grid and number of subgroups
create_age_breaks <- function(age_grid, n_subgroups) {
  age_min <- min(age_grid)
  age_max <- max(age_grid) + 1  # so the last bin is exclusive
  raw_breaks <- seq(age_min, age_max, length.out = n_subgroups + 1)
  int_breaks <- unique(round(raw_breaks))
  if (int_breaks[length(int_breaks)] < age_max) {
    int_breaks[length(int_breaks)] <- age_max
  }
  return(int_breaks)
}

# Generate a dataset based on subgroup proportions and treatment effects
generate_dataset <- function(N, S, theta_s, trt, group_labels) {
  raw_allocations <- N * S
  integer_allocations <- floor(raw_allocations)
  fractional_allocations <- raw_allocations - integer_allocations
  allocated_samples <- integer_allocations
  remaining_samples <- N - sum(allocated_samples)
  if (remaining_samples > 0) {
    ranked_indices <- order(fractional_allocations, decreasing = TRUE)
    allocated_samples[ranked_indices[1:remaining_samples]] <-
      allocated_samples[ranked_indices[1:remaining_samples]] + 1
  }
  data_list <- mapply(
    function(n_i, label, theta) {
      if (n_i > 0) {
        outcomes <- rbinom(n_i, 1, theta)
        data.frame(
          age_group = label,
          outcome   = outcomes,
          trt       = trt
        )
      } else {
        data.frame(
          age_group = label,
          outcome   = 0,
          trt       = trt
        )
      }
    },
    n_i    = allocated_samples,
    label  = group_labels,
    theta  = theta_s,
    SIMPLIFY = FALSE
  )
  dataset <- dplyr::bind_rows(data_list)
  return(dataset)
}

# Generate future trial data given parameters and age group labels
gen_future_trial <- function(N, S, theta_s, p_control, age_groups) {
  N_trt <- N_ctrl <- N / 2
  df_ctrl <- generate_dataset(N_ctrl, S, theta_s = p_control,
                              trt = 0, group_labels = age_groups)
  odds_control   <- p_control / (1 - p_control)
  odds_treatment <- odds_control * exp(theta_s)
  p_treatment    <- odds_treatment / (1 + odds_treatment)
  df_trt <- generate_dataset(N_trt, S, theta_s = p_treatment,
                             trt = 1, group_labels = age_groups)
  df <- dplyr::bind_rows(df_ctrl, df_trt)
  
  # Map original age groups to short labels (e.g. "g1", "g2", …)
  unique_groups <- sort(unique(df$age_group))
  age_group_map <- setNames(paste0("g", seq_along(unique_groups)), unique_groups)
  df <- df %>% dplyr::mutate(age_group_short = age_group_map[as.character(age_group)])
  
  df_summary <- df %>%
    dplyr::group_by(age_group_short) %>%
    dplyr::summarise(
      y_trt  = sum(outcome[trt == 1]),
      y_ctrl = sum(outcome[trt == 0]),
      n_trt  = sum(trt == 1),
      n_ctrl = sum(trt == 0),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      y_trt_adj  = y_trt + 1,
      y_ctrl_adj = y_ctrl + 1,
      n_trt_adj  = n_trt + 2,
      n_ctrl_adj = n_ctrl + 2,
      log_OR = log((y_trt_adj / (n_trt_adj - y_trt_adj)) /
                     (y_ctrl_adj / (n_ctrl_adj - y_ctrl_adj)))
    ) %>%
    dplyr::select(age_group_short, log_OR)
  
  df_wide <- tidyr::pivot_wider(df_summary, names_from = age_group_short,
                                values_from = log_OR,
                                names_prefix = "log_OR_")
  return(df_wide)
}

# Multi-cohort state transition model (STM)
run_multi_cohort_sequential <- function(start_ages, final_tp_ctrl, final_tp_trt, n_cycles,
                                        utilities, cost_ctrl, cost_trt, wtp,
                                        extended_age_grid, n_target, age_groups) {
  multi_cohort_stm <- lapply(start_ages, function(start_age) {
    cohort_list <- list()
    age_idx_range <- which(extended_age_grid >= start_age & extended_age_grid < (start_age + 40))
    for (i in 1:nrow(final_tp_ctrl)) {
      ctrl_probs <- final_tp_ctrl[i, ]
      trt_probs  <- final_tp_trt[i, ]
      cohort_list[[i]] <- cohort_state_transition(ctrl_probs, trt_probs, utilities,
                                                  cost_ctrl, cost_trt, n_cycles,
                                                  age_idx_range, wtp, n_target, extended_age_grid)
    }
    return(cohort_list)
  })
  return(multi_cohort_stm)
}

cohort_state_transition <- function(ctrl_probs, trt_probs, utilities, cost_ctrl, cost_trt,
                                    time_horizon, age_idx_range, wtp, n_target, extended_age_grid) {
  if (length(ctrl_probs) == 1) {
    ctrl_probs <- rep(ctrl_probs, length(extended_age_grid))
  }
  if (length(trt_probs) == 1) {
    trt_probs <- rep(trt_probs, length(extended_age_grid))
  }
  ctrl_probs <- matrix(ctrl_probs, nrow = 1)
  trt_probs  <- matrix(trt_probs, nrow = 1)
  ctrl_probs <- ctrl_probs[, age_idx_range]
  trt_probs  <- trt_probs[, age_idx_range]
  
  trans_array_trt <- trans_array_ctrl <- array(0, dim = c(3, 3, time_horizon))
  for (t in 1:time_horizon) {
    p_ctrl <- ctrl_probs[t]
    p_trt  <- trt_probs[t]
    trans_mat_ctrl <- matrix(c(1 - 0.5 * p_ctrl, 0.5 * p_ctrl, 0,
                               0, 1 - p_ctrl, p_ctrl,
                               0, 0, 1), nrow = 3, byrow = TRUE)
    trans_mat_trt <- matrix(c(1 - 0.5 * p_trt, 0.5 * p_trt, 0,
                              0, 1 - p_trt, p_trt,
                              0, 0, 1), nrow = 3, byrow = TRUE)
    trans_array_trt[,,t] <- trans_mat_trt
    trans_array_ctrl[,,t] <- trans_mat_ctrl
  }
  trace_mat_ctrl <- trace_mat_trt <- matrix(NA, nrow = 3, ncol = time_horizon + 1)
  trace_mat_ctrl[, 1] <- trace_mat_trt[, 1] <- c(n_target, 0, 0)
  for (s in 1:time_horizon) {
    trace_mat_ctrl[, s + 1] <- trace_mat_ctrl[, s] %*% trans_array_ctrl[,,s]
    trace_mat_trt[, s + 1]  <- trace_mat_trt[, s] %*% trans_array_trt[,,s]
  }
  discount_rate <- 0.03
  discount_vector <- (1 / (1 + discount_rate)) ^ (1:time_horizon)
  discounted_trace_ctrl <- sweep(trace_mat_ctrl[,-1], 2, discount_vector, FUN = "*")
  discounted_trace_trt  <- sweep(trace_mat_trt[,-1], 2, discount_vector, FUN = "*")
  cum_qalys_ctrl <- sum(utilities %*% discounted_trace_ctrl)
  cum_costs_ctrl <- sum(cost_ctrl %*% discounted_trace_ctrl)
  cum_qalys_trt  <- sum(utilities %*% discounted_trace_trt)
  cum_costs_trt  <- sum(cost_trt %*% discounted_trace_trt)
  
  NB_ctrl <- wtp * cum_qalys_ctrl - cum_costs_ctrl
  NB_trt  <- wtp * cum_qalys_trt - cum_costs_trt
  
  return(list(NB_ctrl = NB_ctrl, NB_trt = NB_trt))
}

# Average across multiple cohorts
multi_cohort_average <- function(result, weights) {
  weighted_NB_ctrl <- Reduce('+',
                             lapply(seq_along(result), function(i) {
                               sapply(result[[i]], function(x) x$NB_ctrl) * weights[i]
                             }))
  weighted_NB_trt <- Reduce('+',
                            lapply(seq_along(result), function(i) {
                              sapply(result[[i]], function(x) x$NB_trt) * weights[i]
                            }))
  list(NB_ctrl = weighted_NB_ctrl, NB_trt = weighted_NB_trt)
}
