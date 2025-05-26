library(tidyverse)
library(GGally)
library(voi)

# Source required functions
source("analysis_functions.R")
source("simulation_functions.R")

#' Create a parallel coordinate plot of the simplex space
#' @param results_df A tibble from read_optim_results for a specific sample size
#' @param age_groups Vector of age group names
#' @param order_axes Optional vector of indices to reorder the axes (default: NULL for original order)
#' @param line_scaling Factor to control line thickness scaling (higher = more contrast between high/low EVSI)
#' @return A ggplot object
plot_parallel_simplex <- function(
  results_df,
  age_groups,
  order_axes = NULL,
  line_scaling = 1.5
) {
  n_groups <- length(age_groups)

  # Convert list column to regular columns
  plot_df <- results_df %>%
    dplyr::mutate(
      # Name the elements in each vector using age group names
      allocation_num = purrr::map(allocation_num, function(vec) {
        # Ensure we only take the first n_groups elements
        vec <- vec[seq_len(n_groups)]
        names(vec) <- age_groups
        return(vec)
      })
    ) %>%
    tidyr::unnest_wider(allocation_num) %>%
    dplyr::select(dplyr::all_of(age_groups), evsi)

  # If order_axes is provided, reorder the columns directly
  if (
    !is.null(order_axes) &&
      is.numeric(order_axes) &&
      length(order_axes) == n_groups
  ) {
    ordered_groups <- age_groups[order_axes]
    # Reorder the columns directly in the dataframe
    plot_df <- plot_df %>%
      dplyr::select(dplyr::all_of(ordered_groups), evsi)
    cols_to_use <- names(plot_df)[1:n_groups] # Get actual column names after reordering
  } else {
    # Use original order
    cols_to_use <- age_groups
  }

  # Add linewidth based on EVSI
  max_evsi <- max(plot_df$evsi)
  min_evsi <- min(plot_df$evsi)

  plot_df <- plot_df %>%
    dplyr::mutate(
      # Normalize EVSI to 0-1 range
      norm_evsi = (evsi - min_evsi) / (max_evsi - min_evsi),
      # Scale linewidth - higher values get thicker lines
      linewidth = 0.1 + norm_evsi * line_scaling,
      # More aggressive alpha scaling - make low values very transparent
      alpha_value = norm_evsi^2 * 0.9,  # Square the normalized value to make low values more transparent
      # Add row ID for grouping
      id = dplyr::row_number()
    ) %>%
    # Sort by EVSI so higher values are drawn first (underneath)
    dplyr::arrange(desc(evsi))

  # Prepare data in long format for ggplot
  plot_long <- plot_df %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(cols_to_use),
      names_to = "variable",
      values_to = "value"
    ) %>%
    # Convert variable to factor to maintain order
    dplyr::mutate(variable = factor(variable, levels = cols_to_use))

  # Create parallel coordinate plot manually with ggplot
  ggplot(plot_long, aes(x = variable, y = value, group = id, color = evsi)) +
    # Add lines with varying width and transparency
    geom_line(aes(linewidth = linewidth, alpha = alpha_value)) +
    # Add points with matching transparency
    geom_point(aes(alpha = alpha_value), size = 1) +
    # Set color scale - multiple grey shades for better differentiation
    scale_color_gradientn(
      name = "EVSI (£)",
      colors = c("#FFFFFF", "#CCCCCC", "#999999", "#666666", "#333333", "#000000"),
      values = c(0, 0.4, 0.6, 0.8, 0.9, 1),  # More aggressive progression to emphasize high values
      guide = "colorbar"
    ) +
    # Use the calculated linewidth directly
    scale_linewidth_identity() +
    # Use the calculated alpha directly - no need for scale_alpha when we handle it manually
    scale_alpha_identity() +
    # Hide the linewidth and alpha legends
    guides(linewidth = "none", alpha = "none") +
    # Set y-axis limits
    scale_y_continuous(limits = c(0, 1)) +
    # Add labels
    labs(
      x = "Age Group",
      y = "Allocation Proportion",
      title = "Allocation Strategy Comparison",
      subtitle = "Line thickness indicates EVSI value (thicker = higher value)"
    ) +
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 9)
    )
}

#' Create ternary plots showing slices of the simplex space
#' @param results_df A tibble from read_optim_results for a specific sample size
#' @param dims Vector of three indices indicating which dimensions to plot
#' @param age_groups Vector of age group names
#' @param other_sums Sequence of values for sum of other dimensions (default: seq(0, 0.8, by = 0.2))
#' @return A list of ggplot objects, one for each slice
plot_ternary_simplex <- function(results_df, dims, age_groups, 
                                other_sums = seq(0, 0.2, by = 0.05)) {
  if (length(dims) != 3) {
    stop("dims must be a vector of exactly three indices")
  }
  
  n_groups <- length(age_groups)
  if (max(dims) > n_groups) {
    stop("dims contains indices larger than the number of age groups")
  }
  
  # Get indices of dimensions not in dims
  other_dims <- setdiff(seq_len(n_groups), dims)
  
  # Extract shortened age group labels for readability
  labels <- stringr::str_extract(age_groups[dims], "\\d+[-–]\\d+|\\d+\\+")
  if (any(is.na(labels))) {
    labels <- age_groups[dims]  # Fall back to original labels if extraction fails
  }
  
  # Function to create a single ternary plot for a specific slice
  create_slice_plot <- function(df, other_sum) {
    # Calculate the remaining allocation that must be split among our three dimensions
    remaining_alloc <- 1 - other_sum
    
    # Create breaks from 0 to remaining_alloc
    break_seq <- seq(0, remaining_alloc, by = remaining_alloc*0.2)

    # Normalize EVSI for point sizing
    max_evsi <- max(df$evsi)
    min_evsi <- min(df$evsi)
    
    # Transform EVSI values to better separate high values
    df <- df %>%
      dplyr::mutate(
        norm_evsi = (evsi - min_evsi) / (max_evsi - min_evsi),
        # More subtle point size scaling
        point_size = 1 + norm_evsi * 2,
        # Adjust transparency to emphasize high values
        alpha_value = 0.1 + norm_evsi * 0.8
      ) %>%
      dplyr::arrange(evsi)  # Sort by EVSI so higher values are on top
    
    # Create the ternary plot
    p <- ggtern::ggtern(
      data = df,
      mapping = aes(
        x = dim1,
        y = dim2,
        z = dim3,
        color = evsi,
        size = point_size,
        alpha = alpha_value
      )
    ) +
      geom_point() +
      # Use refined blue gradient for better separation of high values
      scale_color_gradientn(
        name = "EVSI (£)",
        colors = c("#FFFFFF", "#E6E6FF", "#B3B3FF", "#6666FF", "#0000FF", "#000066"),
        values = c(0, 0.2, 0.4, 0.6, 0.8, 1),
        guide = "colorbar"
      ) +
      # Set breaks
      ggtern::scale_T_continuous(breaks = seq(0,1, by = 0.2),
                               labels = break_seq) +
      ggtern::scale_L_continuous(breaks = seq(0,1, by = 0.2),
                               labels = break_seq) +
      ggtern::scale_R_continuous(breaks = seq(0,1, by = 0.2),
                               labels = break_seq) +
      scale_size_identity() +
      scale_alpha_identity() +
      labs(
        title = sprintf("Other dimensions sum = %.2f", other_sum),
        subtitle = sprintf(
          "Remaining allocation (%.2f) split between %s",
          remaining_alloc,
          paste(labels, collapse = ", ")
        )
      ) +
      ggtern::Llab(labels[1]) +
      ggtern::Tlab(labels[2]) +
      ggtern::Rlab(labels[3]) +
      theme_bw() +
      theme(tern.axis.title.show = F) +
      ggtern::theme_arrowsmall()
    
    return(p)
  }
  
  # Create plots for each slice
  plot_list <- list()
  
  # Convert list column to regular columns for all dimensions
  full_df <- results_df %>%
    dplyr::mutate(
      allocation_num = purrr::map(allocation_num, function(vec) {
        result <- vec[seq_len(n_groups)]
        names(result) <- age_groups
        return(result)
      })
    ) %>%
    tidyr::unnest_wider(allocation_num)
  
  for (other_sum in other_sums) {
    # Calculate sum of other dimensions
    other_dims_sum <- rowSums(full_df[, age_groups[other_dims], drop = FALSE])
    
    # Filter rows where other dimensions sum exactly to other_sum
    slice_df <- full_df %>%
      dplyr::filter(abs(other_dims_sum - other_sum) < 1e-10)
    
    if (nrow(slice_df) >= 5) {  # Only create plot if we have enough points
      # Extract the three dimensions we want to plot
      plot_data <- slice_df %>%
        dplyr::mutate(
          dim1 = !!sym(age_groups[dims[1]]),
          dim2 = !!sym(age_groups[dims[2]]),
          dim3 = !!sym(age_groups[dims[3]])
        ) %>%
        dplyr::select(dim1, dim2, dim3, evsi)
      
      # Verify that these three dimensions sum to 1 - other_sum
      three_sum <- rowSums(plot_data[, c("dim1", "dim2", "dim3")])
      if (all(abs(three_sum - (1 - other_sum)) < 1e-10)) {
        plot_list[[sprintf("slice_%.2f", other_sum)]] <- create_slice_plot(plot_data, other_sum)
      }
    }
  }
  
  # Add class attribute to make it clear this is a list of grobs
  class(plot_list) <- c("tern_plot_list", "list")
  return(plot_list)
}
