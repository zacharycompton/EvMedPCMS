# Load necessary libraries
library(dplyr)
library(ggplot2)
library(cowplot)

# Read in CSV files for each evolutionary model
resultsOU <- read.csv("./pglsLambdaSimResultsOU123.csv")
resultsBrown <- read.csv("./pglsLambdaSimResultsBrown123.csv")
resultsWhite <- read.csv("pglsLambdaSimResultsWhite123.csv")
resultsPagel <- read.csv("./pglsLambdaSimResultsPagel123.csv")

# Combine datasets and add evolutionary model information as a column
resultsOU$evo_model <- "OU"
resultsBrown$evo_model <- "Brownian"
resultsWhite$evo_model <- "Random Noise"
resultsPagel$evo_model <- "Pagel's Lambda"

# Combine all data into one dataframe
all_results <- bind_rows(resultsOU, resultsBrown, resultsWhite, resultsPagel)

# Adding error type category based on `simslope` and `pval` conditions
all_results <- all_results %>%
  mutate(error_type = case_when(
    simslope != 0 & pval > 0.05 ~ "False Negative",
    simslope == 0 & pval < 0.05 ~ "False Positive",
    TRUE ~ "Correct"
  ))

# Calculate error rates by `min`, `model`, and `evo_model`
error_rates <- all_results %>%
  group_by(min, model, evo_model) %>%
  mutate(total = n()) %>%
  group_by(min, model, evo_model, error_type) %>%
  summarise(count = n(), total = first(total), .groups = 'drop') %>%
  ungroup() %>%
  mutate(rate = count / total)

error_rate_false<-filter(error_rates, error_type != "Correct")

write.csv(error_rate_false, "error_rate_false150.csv")

# Filter for False Negative and False Positive rates
filtered_rates <- error_rates %>%
  filter(error_type %in% c("False Negative", "False Positive")) %>%
  mutate(percent = rate * 100)

# Filter for False Negative and False Positive rates and calculate average error rate by model
average_error_rate_by_model <- error_rates %>%
  filter(error_type %in% c("False Negative", "False Positive")) %>%
  group_by(model) %>%
  summarise(avg_rate = sum(rate) / 16, .groups = 'drop') %>%  # Divide by 4 to account for the four model types
  mutate(percent = avg_rate * 100)

# Expand the data to include all combinations of min, model, evo_model, and error_type
expanded_rates <- expand.grid(
  min = unique(filtered_rates$min),
  model = unique(filtered_rates$model),
  evo_model = unique(filtered_rates$evo_model),
  error_type = c("False Negative", "False Positive")
)

# Filter expanded_rates to keep only valid combinations
expanded_rates <- expanded_rates %>%
  filter(
    (evo_model %in% c("Brownian", "OU", "Pagel's Lambda") & error_type == "False Negative") |
      (evo_model == "Random Noise" & error_type == "False Positive")
  )


brownian_zero_error <- expand.grid(
  min = unique(filtered_rates$min),
  model = unique(filtered_rates$model),
  evo_model = "Brownian",
  error_type = "False Negative"
) %>%
  mutate(percent = 0)

# Merge with the existing filtered_rates data
filtered_rates_with_gaps <- expanded_rates %>%
  left_join(filtered_rates, by = c("min", "model", "evo_model", "error_type")) %>%
  mutate(percent = ifelse(is.na(percent), 0, percent))

filtered_rates_with_gaps <- filtered_rates_with_gaps %>%
  bind_rows(brownian_zero_error)

# Plot with gaps and specific error type rules
ggplot(filtered_rates_with_gaps, aes(x = as.factor(min), y = percent, fill = evo_model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), aes(group = interaction(evo_model, error_type))) +
  facet_wrap(~model, scales = "free_y") +
  labs(
    title = "Error Rates by Evolutionary Model for Each Statistical Model 150 Species",
    x = "Minimum Necropsies",
    y = "Error Rate (%)",
    fill = "Evolutionary Model"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    panel.grid = element_blank(),
    plot.title = element_text(size = 25),         # Increase title text size
    axis.title.x = element_text(size = 20),       # Increase x-axis label text size
    axis.title.y = element_text(size = 20),       # Increase y-axis label text size
    axis.text = element_text(size = 20),          # Increase tick text size
    strip.text = element_text(size = 20),         # Increase facet label text size
    legend.text = element_text(size = 20),        # Increase legend text size
    legend.title = element_text(size = 20)        # Increase legend title text size
  ) +
  scale_y_continuous(limits = c(0, 60), breaks = seq(0, 60, by = 5))

ggsave(file = "error150Species.pdf", width = 13, height = 10)