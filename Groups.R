library(dplyr)
library(stringr)
library(tidyr)
library(lubridate)
library(purrr)
library(ggplot2)
library(patchwork)
library(readxl)

# Read data
base_df <- read.csv("E:/data/Matched/Results.csv", sep = ";") 

# Determine which kind of activity is present in tumours
base_binary <- base_df %>%
  mutate(
    pubmed = if_else(total_publications > 0, 1, 0),
    pubmed_drug = if_else(total_publications_od > 0, 1, 0),
    trials = if_else(n_trials > 0, 1, 0),
    orphan_trials = if_else(n_orphan_true > 0, 1, 0),
    designations = if_else(n_orphan_designations > 0, 1, 0),
    approvals = if_else(n_orphan_approvals > 0, 1, 0),
    
  )

# Put tumours into orphan drug development stage groups
base_groups <- base_binary %>%
  mutate(
    group = case_when(
      pubmed == 0 & trials == 0 & designations == 0 & approvals == 0 ~ "None",
      pubmed == 1 & trials == 0 & designations == 0 & approvals == 0 ~ "Only published literature",
      pubmed == 0 & trials == 1 & designations == 0 & approvals == 0 ~ "Only clinical trials",
      pubmed == 0 & trials == 0 & designations == 1 & approvals == 0 ~ "Only orphan designations",
      pubmed == 1 & trials == 1 & designations == 0 & approvals == 0 ~ "Published literature and clinical trials",
      pubmed == 1 & trials == 0 & designations == 1 & approvals == 0 ~ "Published literature and orphan designations",
      pubmed == 0 & trials== 1 & designations == 1 & approvals == 0 ~ "Clinical trials and orphan designations",
      pubmed == 1 & trials == 1 & designations == 1 & approvals == 0 ~ "Published literature, clinical trials and orphan designations",
      pubmed == 1 & trials == 1 & designations == 1 & approvals == 1 ~ "Published literature, clinical trials, orphan designations and marketing authorisations",
      pubmed == 0 & trials == 0 & designations == 1 & approvals == 1 ~ "Full research"
    )
  ) %>%
  select(-pubmed, -pubmed_drug, -trials, -orphan_trials, -designations, -approvals)

# Add tumour families
base_groups2 <- base_groups %>%
  mutate(tumor_family = case_when(
    str_detect(Tumour, "Squamous") ~ "Squamous cell carcinoma",
    str_detect(Tumour, "Adeno") ~ "Adenocarcinoma",
    str_detect(Tumour, "Undifferentiated") ~ "Undifferentiated carcinoma",
    str_detect(Tumour, "leukemia|lymphoma|Myeloma|Hodgkin|Myeloproliferative|Myelodysplastic|mast cell|histiocytic") ~ "Haematological malignancies",
    str_detect(Tumour, "sarcoma|Ewing|Osteogenic|Chondrogenic|rhabdomyo|vascular") ~ "Sarcomas",
    str_detect(Tumour, "germ cell|seminoma|teratoma") ~ "Germ cell tumours",
    str_detect(Tumour, "melanoma") ~ "Melanoma",
    TRUE ~ "Other epithelial tumours"
  ))

# Save data
write.csv(base_groups, "Results2.csv", row.names = FALSE)

# Calculate how often certain groups appear in the data
plot_groups <- base_groups %>%
  count(group, research_group) %>%
  mutate(percent = n / sum(n) * 100)

plot_groups2 <- base_groups2 %>%
  count(tumor_family, group) 

# # Ranks
# family_order <- c(
#   "Haematological malignancies",
#   "CNS tumours",
#   "Sarcomas",
#   "Germ cell tumours",
#   "Neuroendocrine tumours",
#   "Squamous cell carcinoma",
#   "Adenocarcinoma",
#   "Undifferentiated carcinoma",
#   "Melanoma",
#   "Other epithelial tumours")

# Rank from lowest to highest orphan drug development stage
phase_order <- c(
  "None",
  "Only published literature",
  "Published literature and orphan designations",
  "Only orphan designations",
  "Only clinical trials",
  "Published literature and clinical trials",
  "Clinical trials and orphan designations",
  "Published literature, clinical trials and orphan designations",
  "Published literature, clinical trials, orphan designations and marketing authorisations"
)

# Create tumour families
base_groups2 <- base_groups %>%
  mutate(tumor_family = case_when(
    str_detect(str_to_lower(Tumour), "squamous") ~ "Squamous cell carcinomas",
    str_detect(str_to_lower(Tumour), "adeno") ~ "Adenocarcinomas",
    str_detect(str_to_lower(Tumour), "undifferentiated") ~ "Undifferentiated carcinomas",
    str_detect(str_to_lower(Tumour), "leuk|lymph|aml|myeloma|hodgkin|myeloproliferative|myelodysplastic|mast cell|histiocytic") ~ "Haematological malignancies",
    str_detect(str_to_lower(Tumour), "sarcoma") ~ "Sarcomas",
    str_detect(str_to_lower(Tumour), "germ cell|teratoma") ~ "Germ cell tumours",
    str_detect(str_to_lower(Tumour), "mullerian") ~ "Müllerian mixed tumours",
    str_detect(Tumour, "Semino|seminoma") ~ "Germ cell tumours",
    str_detect(str_to_lower(Tumour), "blastoma") ~ "Blastomas",
    str_detect(str_to_lower(Tumour), "salivary") ~ "Salivary gland type tumours",
    TRUE ~ "Other tumours"
  ))
  #filter(tumor_family != "Other tumours")

write.csv(base_groups2, "Results3.csv", row.names = FALSE)

# Calculate per tumour family and group
plot_groups_ <- base_groups2 %>%
  count(tumor_family, group)

# Calculate incidence per group
incidence_labels <- base_groups2 %>%
  mutate(Crude.incidence.rate.per.100.000 =
           as.numeric(gsub(",", ".", Crude.incidence.rate.per.100.000))) %>%
  group_by(tumor_family) %>%
  summarise(
    total_incidence = sum(Crude.incidence.rate.per.100.000, na.rm = TRUE),
    .groups = "drop"
  )

# Combine
plot_groups2 <- plot_groups_ %>%
  left_join(incidence_labels, by = "tumor_family")

# Compute tumour family ordering based on stage score
family_order <- plot_groups2 %>%
  mutate(score = case_when(
    group == "Published literature, clinical trials, orphan designations and marketing authorisations" ~ 125,
    group == "Published literature, clinical trials and orphan designations" ~ 20,
    group == "Clinical trials and orphan designations" ~ 6,
    group == "Published literature and clinical trials" ~ 4,
    group == "Published literature and orphan designations" ~ 4,
    group == "Only orphan designations" ~ 3,
    group == "Only clinical trials" ~ 3,
    group == "Only published literature" ~ 1,
    group == "None" ~ 0,
    TRUE ~ 0
  )) %>%
  group_by(tumor_family) %>%
  summarise(total_score = sum(n * score, na.rm = TRUE)) %>%
  arrange(desc(total_score)) %>%
  pull(tumor_family)

# Labels
label_df <- plot_groups2 %>%
  group_by(tumor_family) %>%
  summarise(
    y_pos = sum(n),
    total_incidence = first(total_incidence),
    .groups = "drop"
  )

# Plot the results
ggplot(plot_groups2,
       aes(x = factor(tumor_family, levels = family_order),
           y = n,
           fill = factor(group, levels = phase_order))) +
  geom_col() +
  geom_text(
    data = label_df,
    aes(
      x = factor(tumor_family, levels = family_order),
      y = y_pos + 2,  
      label = round(total_incidence, 1)
    ),
    inherit.aes = FALSE,
    size = 5
  ) +
  scale_fill_discrete(
    labels = \(x) stringr::str_wrap(x, width = 35)
    ) +
  labs(
    x = "Tumour group",
    y = "Number of tumours",
    fill = "Orphan drug development stage\n(least to most advanced)",
    title = "Orphan Drug Development Stage Progression Across Several Identified Tumour Groups"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 12),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 18),
        strip.text = element_text(size = 18),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 18),
        plot.title = element_text(size = 26, face = "bold")
  )

# group_order <- plot_groups %>%
#   group_by(group) %>%
#   summarise(total = sum(percent)) %>%
#   arrange(desc(total)) %>%
#   pull(group)
# 
# wrapper <- function(x, ...) {
#   paste(strwrap(x, width = 20), collapse = "\n")
# }
# 
# ggplot(plot_groups, aes(x = factor(group, levels = group_order),
#                         y = percent,
#                         fill = research_group)) +
#   geom_col() +
#   scale_x_discrete(labels = function(x) sapply(x, wrapper)) +
#   labs(x = "Research phase group", y = "Percentage", title = "Distribution of Research Phase and Research Activity Groups", fill = "Research activity group") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         text = element_text(size = 12),
#         axis.title = element_text(size = 22),
#         axis.text = element_text(size = 18),
#         strip.text = element_text(size = 18),
#         legend.title = element_text(size = 22),
#         legend.text = element_text(size = 18),
#         plot.title = element_text(size = 26, face = "bold")
#   )

# variables <- c(
#   "n_orphan_designations",
#   "n_orphan_approvals",
#   "approval_rate",
#   "total_publications",
#   "Crude.incidence.rate.per.100.000",
#   "n_trials",
#   "n_orphan_true",
#   "total_publications_od",
#   "total_publications",
#   "Research_score"
# )

# # Calculate mean scores per orphan drug development stage
# summary_research_group <- base_groups %>%
#   mutate(
#     Crude.incidence.rate.per.100.000 = as.numeric(gsub(",", ".", Crude.incidence.rate.per.100.000)),
#     Research_score = as.numeric(gsub(",", ".", Research_score)),
#     approval_rate = as.numeric(gsub(",", ".", approval_rate))
#   ) %>%
#   group_by(research_group) %>%
#   summarise(
#     across(
#       all_of(variables),
#       ~ mean(.x, na.rm = TRUE)
#     ),
#     n = n(),
#     .groups = "drop"
#   )
# 
# # Calculate median scores per orphan drug development stage
# summary_group <- base_groups %>%
#   mutate(
#     Crude.incidence.rate.per.100.000 = as.numeric(gsub(",", ".", Crude.incidence.rate.per.100.000)),
#     Research_score = as.numeric(gsub(",", ".", Research_score)),
#     approval_rate = as.numeric(gsub(",", ".", approval_rate))
#   ) %>%
#   group_by(group) %>%
#   summarise(
#     across(
#       all_of(variables),
#       ~ mean(.x, na.rm = TRUE)
#     ),
#     n = n(),
#     .groups = "drop"
#   )

# Calculate median scores and IQR per orphan drug development stage
df_plot <- base_groups %>%
  filter(!group %in% c("Only clinical trials", "None")) %>%
  group_by(group) %>%
  summarise(
    median_publications = median(total_publications, na.rm = TRUE),
    q1 = quantile(total_publications, 0.25, na.rm = TRUE),
    q3 = quantile(total_publications, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_publications)) %>%
  mutate(group = factor(group, levels = group))

# Plot the results
ggplot(df_plot, aes(x = group, y = median_publications)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(
    aes(ymin = q1, ymax = q3),
    width = 0.2
  ) +
  scale_x_discrete(
    labels = \(x) stringr::str_wrap(x, width = 25)
  ) +
  labs(
    x = "Orphan drug development stage",
    y = "Median number of publications",
    title = "Published Literature Density Across Orphan Drug Development Stages"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 12),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 18),
        strip.text = element_text(size = 18),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 18),
        plot.title = element_text(size = 26, face = "bold")
  )
