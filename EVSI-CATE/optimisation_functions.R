library(data.table)



adapt_rdn_search <- function(niter, n_trial, nb_samples, eff_samples, 
                             aux_samples, workers, log_file, age_groups,
                             method = "LHS", 
                             batch_size = 200, 
                             top_percent = 10, 
                             convergence_threshold = 5, 
                             patience = 4,
                             corner_init = T,
                             increment = 0.01,
                             scaling_factor = .1) {
  
  
  
  
  # Set up parallel plan
  plan(multisession, workers = workers)
  
  # Initialize log file with headers if it doesn't exist
    if (!file.exists(log_file)) {
    fwrite(
      data.frame(timestamp = character(), 
                 S = character(), 
                 evsi = numeric(), 
                 evsi_sd = numeric(),
                 n_trial = numeric()), 
      file = log_file, 
      sep = ",", 
      col.names = TRUE
    )
  }
  
  # If a previous log exists, load previous results to continue search
  all_allocations <- vector("list", niter)
  all_evsi <- numeric(niter)
  all_evsi_sd <- numeric(niter)
  alloc_counter <- 1  # default if no previous results
  top_allocations <- NULL
  best_evsi_history <- c()
  
  if (file.exists(log_file)) {
    previous_log <- fread(log_file)
    if (nrow(previous_log) > 0) {
      # Convert stored allocation strings back to numeric vectors
      previous_allocs <- lapply(previous_log$S, function(x) {
        x_clean <- sub("^c\\(", "", x)
        x_clean <- sub("\\)$", "", x_clean)
        as.numeric(unlist(strsplit(x_clean, ", ")))
      })
      previous_evsi <- previous_log$evsi
      n_prev <- length(previous_evsi)
      if(n_prev > 0) {
        all_allocations[1:n_prev] <- previous_allocs
        all_evsi[1:n_prev] <- previous_evsi
        alloc_counter <- n_prev + 1
        best_evsi_history <- max(previous_evsi)
        top_n <- ceiling((top_percent / 100) * n_prev)
        prev_alloc_matrix <- do.call(rbind, previous_allocs)
        top_indices <- order(previous_evsi, decreasing = TRUE)[1:top_n]
        top_allocations <- prev_alloc_matrix[top_indices, , drop = FALSE]
      }
    }
  }
  
  # Calculate remaining iterations based on already computed allocations
  remaining_iter <- niter - (alloc_counter - 1)
  num_batches <- ceiling(remaining_iter / batch_size)
  
  # Function to generate simplex samples based on the chosen method
  generate_simplex <- function(n, dim, method) {
    if (method == "LHS") {
      # Latin Hypercube Sampling
      sample_matrix <- randomLHS(n, dim - 1)
      
    } else if (method == "sobol") {
      # Sobol Sequence Sampling
      sobol_points <- sobol(n = n, dim = dim - 1, scrambling = 1, seed = 123)
      
      # Ensure all points are within [0,1]
      sobol_points[sobol_points < 0] <- 0
      sobol_points[sobol_points > 1] <- 1
      
      sample_matrix <- sobol_points
    } else {
      stop("Unsupported sampling method. Choose 'LHS' or 'Sobol'.")
    }
    
    # Sort each row
    sorted_matrix <- t(apply(sample_matrix, 1, sort))
    
    # Add 0 and 1 to each row
    augmented_matrix <- cbind(0, sorted_matrix, 1)
    
    # Compute simplex samples
    simplex_samples <- augmented_matrix[, 2:(dim + 1)] - augmented_matrix[, 1:dim]      
    
    return(simplex_samples)
  }
  
  # Function to adjust allocations to specified increment
  adjust_allocations <- function(allocations, increment = 0.01) {
    if(!is.matrix(allocations)) allocations <- matrix(allocations, nrow = 1)
    
    alloc_rounded <- round(allocations / increment) * increment
    residuals <- 1 - rowSums(alloc_rounded)
    
    # For rows with residuals, distribute them based on original allocation values
    for(i in which(abs(residuals) >= (increment / 2))) {  # Adjust threshold as needed
      # Explicitly round n_steps to avoid floating-point issues
      n_steps <- round(abs(residuals[i]) / increment)
      if(n_steps == 0) next  # No adjustment needed
      
      if(residuals[i] > 0) {
        idx <- order(allocations[i,])[1:n_steps]  # smallest first
      } else {
        idx <- order(allocations[i,], decreasing = TRUE)[1:n_steps]  # largest first
      }
      alloc_rounded[i, idx] <- alloc_rounded[i, idx] + sign(residuals[i]) * increment
    }
    
    return(alloc_rounded)
  }
  
  # Function to generate allocations.
  # If center_alloc is provided, use a Dirichlet distribution with parameters based on the centroid.
  # Function to generate allocations.
  # If center_alloc is provided, use a Dirichlet distribution with parameters based on the centroid.
  generate_allocations <- function(n,
                                   dim,
                                   method = "LHS",
                                   increment = 0.01,
                                   center_alloc = NULL,
                                   scaling_factor = .5,
                                   overshoot_factor = 1.2) {
    
    # 1. Determine how many points to generate (overshoot factor).
    total_count <- as.integer(ceiling(overshoot_factor * n))
    
    # 2. Generate raw samples.
    if (is.null(center_alloc)) {
      # All points come from LHS (or Sobol).
      sample_matrix <- generate_simplex(total_count, dim, method = "LHS")
    } else {
      # Half from LHS, half from Dirichlet around centroid.
      half_count <- as.integer(ceiling(total_count / 2))
      
      # LHS (or Sobol).
      lhs_samples <- generate_simplex(half_count, dim, method = "LHS")
      
      # Dirichlet-based sampling around the centroid.
      centroid <- colMeans(center_alloc)
      alpha_vec <- centroid * scaling_factor
      dirichlet_samples <- MCMCpack::rdirichlet(half_count, alpha_vec)
      
      sample_matrix <- rbind(lhs_samples, dirichlet_samples)
    }
    
    # 3. Round and remove duplicates.
    adjusted_allocations <- adjust_allocations(sample_matrix, increment)
    unique_allocations   <- unique(adjusted_allocations)
    
    # 4. Check the number of unique points left after rounding.
    #    If too few, sample *with replacement* to avoid aborting the job.
    num_unique <- nrow(unique_allocations)
    
    if (num_unique < n) {
      warning(paste(
        "Not enough unique allocation points available after rounding:",
        "only", num_unique, "unique points for", n, "requested.",
        "Sampling with replacement to avoid aborting the job."
      ))
      final_alloc <- unique_allocations[sample(num_unique, n, replace = TRUE), ]
    } else {
      final_alloc <- unique_allocations[sample(num_unique, n, replace = FALSE), ]
    }
    
    return(final_alloc)
  }
  
  # Define single_iter function outside the loop for efficiency
    single_iter <- function(S, n_trial, eff_samples, aux_samples, nb_samples) {
    # Generate future trial data
    dat <- lapply(seq_len(nrow(eff_samples)), function(j) {
      gen_future_trial(
        N = n_trial,
        S = S,
        theta_s = eff_samples[j, ],
        p_control = aux_samples[j, ],
        age_groups = age_groups
      )
    }) %>% bind_rows()
    
    # Fit BART model
    bart_mod <- bart(x.train = dat, y.train = nb_samples, k = 3)
    # Assume bart_post is a B x nsim matrix, where each row is a posterior draw.
    bart_post <- dbarts::extract(bart_mod)
    
    # Compute a replicate EVSI for each posterior draw:
    evsi_rep <- apply(bart_post, 1, function(preds) {
      # preds is a vector (length = nsim) of predicted incremental net benefits.
      # Compute the "within-draw" mean of the positive parts:
      Z <- pmax(preds, 0)
      replicate_evsi <- mean(Z) - pmax(mean(preds), 0)
      replicate_evsi
    })
    
    # Calculate the overall EVSI mean and variance across posterior replicates
    evsi_mean <- mean(evsi_rep)
    evsi_sd <- sd(evsi_rep)
    
    # Return the EVSI mean and variance (logging and subsequent steps remain unchanged)
    return(c(evsi_mean = evsi_mean, evsi_sd = evsi_sd))
  }
  
  # Determine the dimension of the simplex
  dim_simplex <- ncol(eff_samples)
  
  
  
  for (batch in 1:num_batches) {
    cat("Processing batch", batch, "of", num_batches, "\n")
    
    remaining_iter <- niter - ((batch - 1) * batch_size)
    current_batch_size <- ifelse(remaining_iter >= batch_size, batch_size, remaining_iter)
    
    if (is.null(top_allocations) && corner_init) {
      # e.g. corner is [0,...,0,1]
      corner_alloc <- numeric(dim_simplex)
      corner_alloc[dim_simplex] <- 1
      
      # Number of near-corner candidates needed
      num_near_corner <- current_batch_size - 1
      
      # Use a Dirichlet concentration that concentrates points near the corner.
      # Here, for dim_simplex = d, we use: c(rep(1, d - 1), 10)
      alpha_vec <- c(rep(1, dim_simplex - 1), 2)
      
      # Generate near-corner candidates
      near_corner_candidates <- MCMCpack::rdirichlet(num_near_corner, alpha_vec)
      
      # Optionally, round candidates to the desired increment
      near_corner_candidates <- adjust_allocations(near_corner_candidates, increment)
      
      # Combine the exact corner candidate with the near-corner candidates.
      allocation_matrix <- rbind(corner_alloc, near_corner_candidates)
      
    } else if (is.null(top_allocations)) {
      # Original random approach
      allocation_matrix <- generate_allocations(current_batch_size, dim_simplex, method)
    } else {
      # Subsequent batches: half random, half around top
      half_size <- floor(current_batch_size / 2)
      random_alloc <- generate_allocations(half_size, dim_simplex, method)
      sampled_indices <- sample(nrow(top_allocations), half_size, replace = TRUE)
      sampled_top <- top_allocations[sampled_indices, , drop = FALSE]
      perturbed_alloc <- generate_allocations(half_size, dim_simplex,
                                              method = "LHS",
                                              center_alloc = sampled_top)
      allocation_matrix <- rbind(random_alloc, perturbed_alloc)
    }
    
    
    # Execute iterations in parallel
    results_batch <- future_sapply(1:nrow(allocation_matrix), function(i) {
      single_iter(allocation_matrix[i, ], n_trial, eff_samples, aux_samples, nb_samples)
    }, future.seed = TRUE)
    
    # Ensure results_batch is a matrix with 2 rows: first row = evsi_mean, second row = evsi_sd.
    if(nrow(results_batch) != 2) results_batch <- t(results_batch)
    evsi_means_batch <- results_batch[1, ]
    evsi_sds_batch <- results_batch[2, ]
    
    # Store results
    allocations_list <- asplit(allocation_matrix, 1)
    end_index <- alloc_counter + length(allocations_list) - 1
    all_allocations[alloc_counter:end_index] <- allocations_list
    all_evsi[alloc_counter:end_index] <- evsi_means_batch
    if (length(all_evsi_sd) == 0) all_evsi_sd <- numeric(niter)
    all_evsi_sd[alloc_counter:end_index] <- evsi_sds_batch
    alloc_counter <- end_index + 1
    
    # Log batch results
        batch_log <- data.frame(
      timestamp = Sys.time(),
      S = apply(allocation_matrix, 1, function(s) paste0("c(", paste(round(s, 2), collapse = ", "), ")")),
      evsi = evsi_means_batch,
      evsi_sd = evsi_sds_batch,
      n_trial = n_trial
    )
    
    fwrite(batch_log, file = log_file, append = TRUE, sep = ",")
    
    # Identify top allocations from current batch
    top_n <- ceiling((top_percent / 100) * length(evsi_means_batch))
    top_indices <- order(evsi_means_batch, decreasing = TRUE)[1:top_n]
    top_allocations <- allocation_matrix[top_indices, , drop = FALSE]
    current_max_evsi <- max(evsi_means_batch, na.rm = TRUE)
    
    # Update best EVSI history
    best_evsi_history <- c(best_evsi_history, current_max_evsi)
    
    # Check for convergence
    if (length(best_evsi_history) > patience) {
      recent_history <- tail(best_evsi_history, patience)
      improvement <- recent_history[length(recent_history)] - recent_history[1]
      
      if (improvement < convergence_threshold) {
        cat("Convergence achieved at batch", batch, ". Stopping search.\n")
        break
      }
    }
    
    # Optional: Print batch summary
    cat("Batch", batch, "completed. Top EVSI:", current_max_evsi, "\n")
  }
  
  # Trim the pre-allocated lists to actual size
  all_allocations <- all_allocations[1:(alloc_counter - 1)]
  all_evsi <- all_evsi[1:(alloc_counter - 1)]
  all_evsi_sd <- all_evsi_sd[1:(alloc_counter - 1)]
  
  # Return the results
  return(list(allocation = all_allocations, evsi = all_evsi, evsi_sd = all_evsi_sd))
}
