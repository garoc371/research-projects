# Load required libraries
library(tidyverse)

# Define the target directory
target_directory <- "MC_sim/results/unrestricted_spline/"

# List the files in the target directory
existing_files <- list.files(target_directory, pattern = "\\.RDs$", full.names = TRUE)

# Extract the base names of the existing files and remove the suffix
existing_scenarios <- basename(existing_files) %>%
  str_remove("(_extended|_limited)\\.RDs$") %>%
  unique()

# Define the control and treatment types
control_types <- c("non_linear_mon_increase", "non_linear_mon_decrease", "non_linear_non_mon")
treatment_types <- c("constant", "non_linear_mon_increase", "non_linear_mon_decrease", "non_linear_non_mon")

# Create the combinations
combinations <- expand.grid(control_type = control_types, treatment_type = treatment_types)


# Add short versions to the combinations data frame
combinations <- combinations %>%
  mutate(
    control_short = gsub("non_linear", "non", gsub("linear", "lin", as.character(control_type))),
    treatment_short = gsub("non_linear", "non", gsub("linear", "lin", as.character(treatment_type))),
    scenario = paste0(control_short, "_ctrl_", treatment_short, "_trt")
  )

# Filter out the existing scenarios
missing_combinations <- combinations %>%
  filter(!scenario %in% existing_scenarios)
