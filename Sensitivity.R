library(dplyr)
library(stringr)
library(tidyr)
library(lubridate)
library(purrr)
library(ggplot2)
library(patchwork)
library(readxl)

# Read data
tumours_raw <- read_excel("E:/data/RARECAREnet_list_of_rare_cancers.xlsx")
tumours_final <- read.csv("E:/data/Splitted2/tumours_final.csv")
research_scores_df <- read.csv("E:/data/Matched/Research_scores.csv")
research_scores_s1_df <- read.csv("E:/data/Matched/Research_scores_s1.csv")
research_scores_s2_df <- read.csv("E:/data/Matched/Research_scores_s2.csv")
research_scores_s3_df <- read.csv("E:/data/Matched/Research_scores_s3.csv")
results_df <- read.csv("E:/data/Matched/Results2.csv", sep = ",") %>%
  mutate(Tumours.ID = as.character(Tumours.ID))

# Get RARECARE tumour name
tumour_names <- tumours_final %>%
  mutate(Tumour.ID = as.character(Tumour.ID)) %>%
  select(Tumour.ID, Synonym_main, Crude.incidence.rate.per.100.000) %>%
  rename(Tumours.ID = Tumour.ID) %>%
  distinct()

# Original RARECARE names
tumours_raw_filtered <- tumours_raw %>%
  filter(`R=rare` == "R") %>%
  filter(Tier == 2) %>%
  mutate(Tumour = str_remove(Tumour, "\\*")) %>%
  select(Tumour, `Crude incidence rate per 100,000`)

# Split scores into tertiles and put tumours in groups
big_df <- tumour_names %>%
  mutate(Tumours.ID = as.character(Tumours.ID)) %>%
  full_join(research_scores_df %>%mutate(Tumours.ID = as.character(Tumours.ID)), by = "Tumours.ID") %>%
  distinct() %>%
  select(Tumour, Tumours.ID, Crude.incidence.rate.per.100.000, Research_score, Normalised_research_score) %>%
  arrange(desc(Research_score)) %>%
  mutate(
    research_group = case_when(
      Research_score == 0 ~ "Non-researched",
      ntile(if_else(Research_score > 0,
                    Research_score,
                    NA_real_), 3) == 1 ~ "Under-researched",
      ntile(if_else(Research_score > 0,
                    Research_score,
                    NA_real_), 3) == 2 ~ "Moderately-researched",
      ntile(if_else(Research_score > 0,
                    Research_score,
                    NA_real_), 3) == 3 ~ "Well-researched"
    ),
    rank = rank(-Research_score, ties.method = "first"),
    # Include top and bottom 10 tumour information
    top_bottom_group = case_when(
      rank <= 10 ~ "Top 10",
      rank > n() - 10 ~ "Bottom 10",
      TRUE ~ NA_character_
    )
  )

# Retrieve top and bottom 10 tumorus by (normalised) research score
top_10_RS <- big_df %>%
  slice_max(order_by = as.numeric(Research_score), n = 10, with_ties = FALSE)
top_10_nRS <- big_df %>%
  slice_max(order_by = as.numeric(Normalised_research_score), n = 10, with_ties = FALSE)
bottom_10_RS <- big_df %>%
  slice_min(order_by = as.numeric(Research_score), n = 10, with_ties = FALSE)
bottom_10_nRS <- big_df %>%
  slice_min(order_by = as.numeric(Normalised_research_score), n = 10, with_ties = FALSE)

# Plot
p1 <- ggplot(top_10_RS, aes(x = reorder(Tumour, Research_score), y = Research_score)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) +
  labs(
    x = "Tumour",
    y = "Research score",
    title = "Baseline Analysis"
  ) +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 22),
    plot.title = element_text(size = 26)
  )

p2 <- ggplot(top_10_nRS, aes(x = reorder(Tumour, Normalised_research_score), y = Normalised_research_score)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) +
  labs(
    x = "Tumour",
    y = "Normalised research score",
    title = "Baseline Analysis"
  ) +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 22),
    plot.title = element_text(size = 26)
  )

p3 <- ggplot(bottom_10_RS, aes(x = reorder(Tumour, Research_score), y = Research_score)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) +
  labs(
    x = "Tumour",
    y = "Research score",
    title = "Baseline Analysis"
  ) +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 22),
    plot.title = element_text(size = 26)
  )

p4 <- ggplot(bottom_10_nRS, aes(x = reorder(Tumour, Normalised_research_score), y = Normalised_research_score)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) +
  labs(
    x = "Tumour",
    y = "Normalised research score",
    title = "Baseline Analysis"
  ) +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 22),
    plot.title = element_text(size = 26)
  )

# Split scores into tertiles and put tumours in groups
big_df1 <- tumour_names %>%
  full_join(research_scores_s1_df %>% mutate(Tumours.ID = as.character(Tumours.ID)), by = "Tumours.ID") %>%
  distinct() %>%
  mutate(Tumour = as.character(Tumour)) %>%
  select(Tumour, Tumours.ID, Crude.incidence.rate.per.100.000, Research_score, Normalised_research_score) %>%
  arrange(desc(Research_score)) %>%
  mutate(
    research_group = case_when(
      Research_score == 0 ~ "Non-researched",
      ntile(if_else(Research_score > 0,
                    Research_score,
                    NA_real_), 3) == 1 ~ "Under-researched",
      ntile(if_else(Research_score > 0,
                    Research_score,
                    NA_real_), 3) == 2 ~ "Moderately-researched",
      ntile(if_else(Research_score > 0,
                    Research_score,
                    NA_real_), 3) == 3 ~ "Well-researched"
    ),
    rank = rank(-Research_score, ties.method = "first"),
    # Include top and bottom 10 tumour information
    top_bottom_group = case_when(
      rank <= 10 ~ "Top 10",
      rank > n() - 10 ~ "Bottom 10",
      TRUE ~ NA_character_
    )
  )

# Retrieve top and bottom 10 tumorus by (normalised) research score
top_10_RS1 <- big_df1 %>%
  slice_max(order_by = as.numeric(Research_score), n = 10, with_ties = FALSE)
top_10_nRS1 <- big_df1 %>%
  slice_max(order_by = as.numeric(Normalised_research_score), n = 10, with_ties = FALSE)
bottom_10_RS1 <- big_df1 %>%
  slice_min(order_by = as.numeric(Research_score), n = 10, with_ties = FALSE)
bottom_10_nRS1 <- big_df1 %>%
  slice_min(order_by = as.numeric(Normalised_research_score), n = 10, with_ties = FALSE)

# Plot
p1s1 <- ggplot(top_10_RS1, aes(x = reorder(Tumour, Research_score), y = Research_score)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) +
  labs(
    x = "Tumour",
    y = "Research score",
    title = "MA More Important"
  ) +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 22),
    plot.title = element_text(size = 26)
  )

p2s1 <- ggplot(top_10_nRS1, aes(x = reorder(Tumour, Normalised_research_score), y = Normalised_research_score)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) +
  labs(
    x = "Tumour",
    y = "Normalised research score",
    title = "MA More Important"
  ) +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 22),
    plot.title = element_text(size = 26)
  )

p3s1 <- ggplot(bottom_10_RS1, aes(x = reorder(Tumour, Research_score), y = Research_score)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) +
  labs(
    x = "Tumour",
    y = "Research score",
    title = "MA More Important"
  ) +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 22),
    plot.title = element_text(size = 26)
  )

p4s1 <- ggplot(bottom_10_nRS1, aes(x = reorder(Tumour, Normalised_research_score), y = Normalised_research_score)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) +
  labs(
    x = "Tumour",
    y = "Normalised research score",
    title = "MA More Important"
  ) +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 22),
    plot.title = element_text(size = 26)
  )

# Split scores into tertiles and put tumours in groups
big_df2 <- tumour_names %>%
  full_join(research_scores_s2_df %>% mutate(Tumours.ID = as.character(Tumours.ID)), by = "Tumours.ID") %>%
  distinct() %>%
  mutate(Tumour = as.character(Tumour)) %>%
  select(Tumour, Tumours.ID, Crude.incidence.rate.per.100.000, Research_score, Normalised_research_score) %>%
  arrange(desc(Research_score)) %>%
  mutate(
    research_group = case_when(
      Research_score == 0 ~ "Non-researched",
      ntile(if_else(Research_score > 0,
                    Research_score,
                    NA_real_), 3) == 1 ~ "Under-researched",
      ntile(if_else(Research_score > 0,
                    Research_score,
                    NA_real_), 3) == 2 ~ "Moderately-researched",
      ntile(if_else(Research_score > 0,
                    Research_score,
                    NA_real_), 3) == 3 ~ "Well-researched"
    ),
    rank = rank(-Research_score, ties.method = "first"),
    # Include top and bottom 10 tumour information
    top_bottom_group = case_when(
      rank <= 10 ~ "Top 10",
      rank > n() - 10 ~ "Bottom 10",
      TRUE ~ NA_character_
    )
  )

# Retrieve top and bottom 10 tumorus by (normalised) research score
top_10_RS2 <- big_df2 %>%
  slice_max(order_by = as.numeric(Research_score), n = 10, with_ties = FALSE)
top_10_nRS2 <- big_df2 %>%
  slice_max(order_by = as.numeric(Normalised_research_score), n = 10, with_ties = FALSE)
bottom_10_RS2 <- big_df2 %>%
  slice_min(order_by = as.numeric(Research_score), n = 10, with_ties = FALSE)
bottom_10_nRS2 <- big_df2 %>%
  slice_min(order_by = as.numeric(Normalised_research_score), n = 10, with_ties = FALSE)

# Plot
p1s2 <- ggplot(top_10_RS2, aes(x = reorder(Tumour, Research_score), y = Research_score)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) +
  labs(
    x = "Tumour",
    y = "Research score",
    title = "MA Removed"
  ) +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 22),
    plot.title = element_text(size = 26)
  )

p2s2 <- ggplot(top_10_nRS2, aes(x = reorder(Tumour, Normalised_research_score), y = Normalised_research_score)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) +
  labs(
    x = "Tumour",
    y = "Normalised research score",
    title = "MA Removed"
  ) +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 22),
    plot.title = element_text(size = 26)
  )

p3s2 <- ggplot(bottom_10_RS2, aes(x = reorder(Tumour, Research_score), y = Research_score)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) +
  labs(
    x = "Tumour",
    y = "Research score",
    title = "MA Removed"
  ) +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 22),
    plot.title = element_text(size = 26)
  )

p4s2 <- ggplot(bottom_10_nRS2, aes(x = reorder(Tumour, Normalised_research_score), y = Normalised_research_score)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) +
  labs(
    x = "Tumour",
    y = "Normalised research score",
    title = "MA Removed"
  ) +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 22),
    plot.title = element_text(size = 26)
  )

# Split scores into tertiles and put tumours in groups
big_df3 <- tumour_names %>%
  full_join(research_scores_s3_df %>% mutate(Tumours.ID = as.character(Tumours.ID)), by = "Tumours.ID") %>%
  distinct() %>%
  mutate(Tumour = as.character(Tumour)) %>%
  select(Tumour, Tumours.ID, Crude.incidence.rate.per.100.000, Research_score, Normalised_research_score) %>%
  arrange(desc(Research_score)) %>%
  mutate(
    research_group = case_when(
      Research_score == 0 ~ "Non-researched",
      ntile(if_else(Research_score > 0,
                    Research_score,
                    NA_real_), 3) == 1 ~ "Under-researched",
      ntile(if_else(Research_score > 0,
                    Research_score,
                    NA_real_), 3) == 2 ~ "Moderately-researched",
      ntile(if_else(Research_score > 0,
                    Research_score,
                    NA_real_), 3) == 3 ~ "Well-researched"
    ),
    rank = rank(-Research_score, ties.method = "first"),
    # Include top and bottom 10 tumour information
    top_bottom_group = case_when(
      rank <= 10 ~ "Top 10",
      rank > n() - 10 ~ "Bottom 10",
      TRUE ~ NA_character_
    )
  )

# Retrieve top and bottom 10 tumorus by (normalised) research score
top_10_RS3 <- big_df3 %>%
  slice_max(order_by = as.numeric(Research_score), n = 10, with_ties = FALSE)
top_10_nRS3 <- big_df3 %>%
  slice_max(order_by = as.numeric(Normalised_research_score), n = 10, with_ties = FALSE)
bottom_10_RS3 <- big_df3 %>%
  slice_min(order_by = as.numeric(Research_score), n = 10, with_ties = FALSE)
bottom_10_nRS3 <- big_df3 %>%
  slice_min(order_by = as.numeric(Normalised_research_score), n = 10, with_ties = FALSE)

# Plot
p1s3 <- ggplot(top_10_RS3, aes(x = reorder(Tumour, Research_score), y = Research_score)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) +
  labs(
    x = "Tumour",
    y = "Research score",
    title = "MA Less Important"
  ) +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 22),
    plot.title = element_text(size = 26)
  )

p2s3 <- ggplot(top_10_nRS3, aes(x = reorder(Tumour, Normalised_research_score), y = Normalised_research_score)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) +
  labs(
    x = "Tumour",
    y = "Normalised research score",
    title = "MA Less Important"
  ) +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 22),
    plot.title = element_text(size = 26)
  )

p3s3 <- ggplot(bottom_10_RS3, aes(x = reorder(Tumour, Research_score), y = Research_score)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) +
  labs(
    x = "Tumour",
    y = "Research score",
    title = "MA Less Important"
  ) +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 22),
    plot.title = element_text(size = 26)
  )

p4s3 <- ggplot(bottom_10_nRS3, aes(x = reorder(Tumour, Normalised_research_score), y = Normalised_research_score)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) +
  labs(
    x = "Tumour",
    y = "Normalised research score",
    title = "MA Less Important"
  ) +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 22),
    plot.title = element_text(size = 26)
  )

# Function that changes text width
wrapper <- function(x, ...) {
  paste(strwrap(x, width = 40), collapse = "\n")
}

# Add plots together 
plot1 <- (p1 / p1s1) | (p1s2 / p1s3) 
plot1 +
  plot_annotation(
  title = "Sensitivity Analysis: Top 10 Rare Tumours by Research Score",
  tag_levels = "A",
  theme = theme(
    plot.title = element_text(size = 32, face = "bold", hjust = 0.5, margin = margin(b = 20))
  )
)

plot2 <- (p2/p2s1) | (p2s2/p2s3) 
plot2 +
  plot_annotation(
    title = "Sensitivity Analysis: Top 10 Rare Tumours by Normalised Research Score",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 32, face = "bold", hjust = 0.5, margin = margin(b = 20))
    )
  )

plot3 <- (p3/p3s1) | (p3s2/p3s3) 
plot3 +
  plot_annotation(
    title = "Sensitivity Analysis: Bottom 10 Rare Tumours by Research Score",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 32, face = "bold", hjust = 0.5, margin = margin(b = 20))
    )
  )

plot4 <- (p4/p4s1) | (p4s2/p4s3) 
plot4 +
  plot_annotation(
    title = "Sensitivity Analysis: Bottom 10 Rare Tumours by Normalised Research Score",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 32, face = "bold", hjust = 0.5, margin = margin(b = 20))
    )
  )

# Function that enables multiple (normalised) research scores and group to be added to one dataframe
prep_df <- function(df, i) {
  df %>%
    select(Tumour, Research_score, Normalised_research_score, research_group, Tumours.ID) %>%
    rename(
      !!paste0("Research_score_", i) := Research_score,
      !!paste0("Normalised_research_score_", i) := Normalised_research_score,
      !!paste0("research_group_", i) := research_group
    )
}

# Add all sensitivity analysis data together
big_df_all <- list(
  prep_df(big_df, 1),
  prep_df(big_df1, 2),
  prep_df(big_df2, 3),
  prep_df(big_df3, 4)
) %>%
  reduce(full_join, by = c("Tumour", "Tumours.ID")) %>%
  arrange(desc(Research_score_1)) %>%
  left_join(results_df %>% select(Tumours.ID, group), by = "Tumours.ID")

# Save data
write.csv2(big_df_all, "Results_sens.csv", row.names = FALSE, quote = TRUE)

