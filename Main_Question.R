library(dplyr)
library(stringr)
library(tidyr)
library(lubridate)
library(purrr)
library(ggplot2)
library(patchwork)
library(readxl)

# Read data and only include trials that are included in products
tumours_raw <- read_excel("E:/data/RARECAREnet_list_of_rare_cancers.xlsx")
tumours_final <- read.csv("E:/data/Splitted2/tumours_final.csv")
products_df <- read.csv("E:/data/Matched/products.csv")
indications_df <- read.csv("E:/data/Matched/indications.csv") %>%
  semi_join(products_df, by = c("CT.EDU.ID", "NCT.ID"))
orphan_indications_df <- read.csv("E:/data/Matched/orphan_designation_better.csv")
sponsor_df <- read.csv("E:/data/Splitted2/Sponsors.csv") %>%
  semi_join(products_df, by = c("CT.EDU.ID", "NCT.ID"))
trial_information_df <- read.csv("E:/data/Splitted2/Trial_information.csv") %>%
  semi_join(products_df, by = c("CT.EDU.ID", "NCT.ID"))
PubMed_df <- read.csv("E:/data/Matched/PubMed.csv") %>%
  select(Tumours.ID, PubMed.Count.General) %>%
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

# Add together
tumour_names <- bind_cols(tumours_raw_filtered, tumour_names) %>%
  select(Tumours.ID, Crude.incidence.rate.per.100.000, Tumour, Synonym_main)

# Retrieve trial IDs
nct_ids <- indications_df$NCT.ID
edu_ids <- indications_df$CT.EDU.ID

# Filter other datasets on rare tumour trials based on trial IDs
sponsor_df_filtered <- sponsor_df %>%
  filter((!is.na(NCT.ID) & NCT.ID %in% nct_ids) |
           (!is.na(CT.EDU.ID) & CT.EDU.ID %in% edu_ids))

trial_information_df_filtered <- trial_information_df %>%
  filter((!is.na(NCT.ID) & NCT.ID %in% nct_ids) |
           (!is.na(CT.EDU.ID) & CT.EDU.ID %in% edu_ids))

products_df_filtered <- products_df %>%
  filter((!is.na(NCT.ID) & NCT.ID %in% nct_ids) |
           (!is.na(CT.EDU.ID) & CT.EDU.ID %in% edu_ids))

# Add missing trial information status data
trial_information_update <- tibble(
  CT.EDU.ID = c(
    "2005-004599-19",
    "2012-001477-82",
    "2015-002916-34",
    "2017-001804-31",
    "2017-003434-87",
    "2021-000677-89"
  ),
  Member.State = c(
    "france",
    NA,
    "germany",
    NA,
    NA,
    "netherlands"
  ),
  Trial.Status = c(
    "completed",
    "completed",
    "completed",
    "completed",
    "unknown",
    "ongoing"
  ),
  End.of.Trial = c(
    "01 Jan 2023",
    NA,
    NA,
    "20 Jun 2007",
    NA,
    NA
  )
) %>%
  mutate(End.of.Trial = dmy(End.of.Trial))

# Join missing trial information data with old trial information data
trial_information_df_filtered <- trial_information_df_filtered %>%
  left_join(trial_information_update, by = "CT.EDU.ID", suffix = c("", ".new")) %>%
  mutate(
    End.of.Trial = ymd(End.of.Trial),
    Trial.Status = coalesce(Trial.Status.new, Trial.Status),
    Trial.ID = coalesce(NCT.ID, CT.EDU.ID),
    End.of.Trial = coalesce(End.of.Trial.new, End.of.Trial)
  ) %>%
  select(-Trial.Status.new, -End.of.Trial.new)

# Separate tumours and add to clinical trial data
indications_long <- indications_df %>%
  separate_rows(tumours_id, sep = ",") %>%
  rename(Tumours.ID = tumours_id) %>%
  mutate(Tumours.ID = str_trim(Tumours.ID),
         Tumours.ID = str_squish(Tumours.ID),
         Trial.ID = coalesce(NCT.ID, CT.EDU.ID)) %>%
  left_join(trial_information_df_filtered, by = c("NCT.ID", "CT.EDU.ID"), relationship = "many-to-many") %>%
  rename(Trial.ID = Trial.ID.x) %>%
  group_by(Tumours.ID, Trial.ID) %>%
  # Get last trial start date
  slice_max(Start.of.Trial, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(Trial.ID, NCT.ID, CT.EDU.ID, Indication, Tumours.ID, Start.of.Trial)

# Separate tumours
orphan_long <- orphan_indications_df %>%
  separate_rows(Tumours.ID, sep = ",") %>%
  mutate(Tumours.ID = str_squish(Tumours.ID),
         Orphan.ID = coalesce(as.character(FDA.ID), EMA.ID))
  
# PM_scores <- orphan_long %>%
#   group_by(Tumours.ID) %>%
#   summarise(
#     PubMed_Count = sum(PubMed.Count, na.rm = TRUE)
#   ) %>%
#   rename(PM_score = PubMed_Count)

# Main research question:
# Determine PubMed score:
PM_scores <- PubMed_df %>%
  rename(PM_score = PubMed.Count.General) %>%
  mutate(
    PM_score = if_else(Tumours.ID == 197,
                                   14365,
                       PM_score)
  )
  
# Determine OD score:
OD_scores <- orphan_long %>%
  filter(Orphan.Designation.Status == "designated") %>%
  group_by(Tumours.ID) %>%
  summarise(
    OD_score = n(),
    .groups = "drop"
  )

# Determine MA score:
MA_scores <- orphan_long %>%
  filter(Orphan.Designation.Status == "approved") %>%
  group_by(Tumours.ID) %>%
  summarise(
    MA_score = n(),
    .groups = "drop"
  )

# Determine CT score:
CT_scores <- indications_long %>%
  # Combine trial information with tumour information
  #left_join(
  #  trial_information_df_filtered, by = c("NCT.ID", "CT.EDU.ID"), relationship = "many-to-many"
  #) %>%
  # Prevent double trials
  #distinct(Tumours.ID, NCT.ID, CT.EDU.ID, Trial.Status, Trial.Phase) %>%
  #mutate(
    # Phase score
  #  phase_score = case_when(
  #    str_detect(Trial.Phase, "IV") ~ 1.0,
  #    str_detect(Trial.Phase, "III") ~ 0.75,
  #    str_detect(Trial.Phase, "II") ~ 0.5,
  #    str_detect(Trial.Phase, "I") ~ 0.25,
  #  TRUE ~ 0.25                                                               # No trial phase
  #  ),
    
  #  # Status score
  #  status_score = case_when(
  #    Trial.Status %in% c("completed") ~ 1.0,
  #    Trial.Status %in% c("ongoing", "recruiting", "authorised", "restarted", "trial now transitioned", "eea") ~ 0.6,
  #    Trial.Status %in% c("temporarily halted", "unknown") ~ 0.3,
  #    TRUE ~ 0.1                                                                # Prematurely ended and prohibited
  #  ),
    
    # Calculate trial score
  #  trial_score = phase_score * status_score
  #) %>%
  group_by(Tumours.ID) %>%
  # Sum the trial scores per tumour
  summarise(
    #CT_score = sum(trial_score, na.rm = TRUE),
    CT_score = n(),
    .groups = "drop"
  )

# Min-max scaling function
minmax_scale <- function(x, log = FALSE) {
  if (log) x <- log1p(x)
  if (max(x, na.rm = TRUE) == min(x, na.rm = TRUE)) {
    return(rep(0, length(x)))
  }
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# Calculate total research score:
Research_score <- OD_scores %>%
  full_join(CT_scores, by = "Tumours.ID") %>%
  full_join(MA_scores, by = "Tumours.ID") %>%
  full_join(PM_scores, by = "Tumours.ID") %>%
  full_join(tumour_names, by = "Tumours.ID") %>%
  # Change NA scores to 0
  mutate(
    OD_score = ifelse(is.na(OD_score), 0, OD_score),
    MA_score = ifelse(is.na(MA_score), 0, MA_score),
    PM_score = ifelse(is.na(PM_score), 0, PM_score),
    CT_score = ifelse(is.na(CT_score), 0, CT_score)
  ) %>%
  # Scale scores with min-max scaling
  mutate(
    OD_score_scaled = minmax_scale(OD_score),
    MA_score_scaled = minmax_scale(MA_score),
    PM_score_scaled = minmax_scale(PM_score, log = TRUE),
    CT_score_scaled = minmax_scale(CT_score)
  ) %>%
  # Calculate research score
  mutate(
    Research_score = round((0.3 * OD_score_scaled + 0.4 * MA_score_scaled + 0.2 * CT_score_scaled + 0.1 * PM_score_scaled), 2),
    Research_score = ifelse(
      OD_score == 0 & MA_score == 0 & CT_score == 0 & PM_score == 0,
      0,
      Research_score)
  ) %>%
  mutate(
    Normalised_research_score =
      round(Research_score / sqrt(Crude.incidence.rate.per.100.000), 2)
  )
  
# Correlation test
cor.test(Research_score$Research_score,
         Research_score$Crude.incidence.rate.per.100.000,
         method = "pearson")

# Plot the linear relationship between incidence rate and research score
ggplot(Research_score, aes(x = Crude.incidence.rate.per.100.000, y = Research_score)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Crude incidence per 100k", y = "Research score",
       title = "Relationship between tumour incidence and research score")

# Only include specific columns for saving
Research_score_select <- Research_score %>%
  select(Tumours.ID, OD_score_scaled, MA_score_scaled, PM_score_scaled, CT_score_scaled, Research_score, Normalised_research_score)

# Save data
write.csv(Research_score_select, "Research_scores.csv", row.names = FALSE)

# Sensitivity analysis:
# MA more important:
Research_score <- OD_scores %>%
  full_join(CT_scores, by = "Tumours.ID") %>%
  full_join(MA_scores, by = "Tumours.ID") %>%
  full_join(PM_scores, by = "Tumours.ID") %>%
  full_join(tumour_names, by = "Tumours.ID") %>%
  # Change NA scores to 0
  mutate(
    OD_score = ifelse(is.na(OD_score), 0, OD_score),
    MA_score = ifelse(is.na(MA_score), 0, MA_score),
    PM_score = ifelse(is.na(PM_score), 0, PM_score),
    CT_score = ifelse(is.na(CT_score), 0, CT_score)
  ) %>%
  # Scale scores with min-max scaling
  mutate(
    OD_score_scaled = minmax_scale(OD_score),
    MA_score_scaled = minmax_scale(MA_score),
    PM_score_scaled = minmax_scale(PM_score, log = TRUE),
    CT_score_scaled = minmax_scale(CT_score)
  ) %>%
  # Calculate research score
  mutate(
    Research_score = round((0.2 * OD_score_scaled + 0.6 * MA_score_scaled + 0.1 * CT_score_scaled + 0.1 * PM_score_scaled), 2),
    Research_score = ifelse(
      OD_score == 0 & MA_score == 0 & CT_score == 0 & PM_score == 0,
      0,
      Research_score)
  ) %>%
  mutate(
    Normalised_research_score =
      round(Research_score / sqrt(Crude.incidence.rate.per.100.000), 2)
)

# Only include specific columns for saving
Research_score_select <- Research_score %>%
  select(Tumours.ID, OD_score_scaled, MA_score_scaled, PM_score_scaled, CT_score_scaled, Research_score, Normalised_research_score)

# Save data
write.csv(Research_score_select, "Research_scores_s1.csv", row.names = FALSE)

# MA removed
Research_score <- OD_scores %>%
  full_join(CT_scores, by = "Tumours.ID") %>%
  full_join(MA_scores, by = "Tumours.ID") %>%
  full_join(PM_scores, by = "Tumours.ID") %>%
  full_join(tumour_names, by = "Tumours.ID") %>%
  # Change NA scores to 0
  mutate(
    OD_score = ifelse(is.na(OD_score), 0, OD_score),
    MA_score = ifelse(is.na(MA_score), 0, MA_score),
    PM_score = ifelse(is.na(PM_score), 0, PM_score),
    CT_score = ifelse(is.na(CT_score), 0, CT_score)
  ) %>%
  # Scale scores with min-max scaling
  mutate(
    OD_score_scaled = minmax_scale(OD_score),
    MA_score_scaled = minmax_scale(MA_score),
    PM_score_scaled = minmax_scale(PM_score, log = TRUE),
    CT_score_scaled = minmax_scale(CT_score)
  ) %>%
  # Calculate research score
  mutate(
    Research_score = round((0.5 * OD_score_scaled + 0.0 * MA_score_scaled + 0.3 * CT_score_scaled + 0.2 * PM_score_scaled), 2),
    Research_score = ifelse(
      OD_score == 0 & MA_score == 0 & CT_score == 0 & PM_score == 0,
      0,
      Research_score)
  ) %>%
  mutate(
    Normalised_research_score =
      round(Research_score / sqrt(Crude.incidence.rate.per.100.000), 2)
  )

# Only include specific columns for saving
Research_score_select <- Research_score %>%
  select(Tumours.ID, OD_score_scaled, MA_score_scaled, PM_score_scaled, CT_score_scaled, Research_score, Normalised_research_score)

# Save data
write.csv(Research_score_select, "Research_scores_s2.csv", row.names = FALSE)

# MA less important
Research_score <- OD_scores %>%
  full_join(CT_scores, by = "Tumours.ID") %>%
  full_join(MA_scores, by = "Tumours.ID") %>%
  full_join(PM_scores, by = "Tumours.ID") %>%
  full_join(tumour_names, by = "Tumours.ID") %>%
  # Change NA scores to 0
  mutate(
    OD_score = ifelse(is.na(OD_score), 0, OD_score),
    MA_score = ifelse(is.na(MA_score), 0, MA_score),
    PM_score = ifelse(is.na(PM_score), 0, PM_score),
    CT_score = ifelse(is.na(CT_score), 0, CT_score)
  ) %>%
  # Scale scores with min-max scaling
  mutate(
    OD_score_scaled = minmax_scale(OD_score),
    MA_score_scaled = minmax_scale(MA_score),
    PM_score_scaled = minmax_scale(PM_score, log = TRUE),
    CT_score_scaled = minmax_scale(CT_score)
  ) %>%
  # Calculate research score
  mutate(
    Research_score = round((0.4 * OD_score_scaled + 0.3 * MA_score_scaled + 0.2 * CT_score_scaled + 0.1 * PM_score_scaled), 2),
    Research_score = ifelse(
      OD_score == 0 & MA_score == 0 & CT_score == 0 & PM_score == 0,
      0,
      Research_score)
  ) %>%
  mutate(
    Normalised_research_score =
      round(Research_score / sqrt(Crude.incidence.rate.per.100.000), 2)
  )

# Only include specific columns for saving
Research_score_select <- Research_score %>%
  select(Tumours.ID, OD_score_scaled, MA_score_scaled, PM_score_scaled, CT_score_scaled, Research_score, Normalised_research_score)

# Save data
write.csv(Research_score_select, "Research_scores_s3.csv", row.names = FALSE)

