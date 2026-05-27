library(readxl)
library(dplyr)
library(openxlsx)
library(stringr)
library(tidyr)
library(lubridate)
library(purrr)
library(jsonlite)

# Cleaning function
clean_text <- function(x){
  x <- tolower(x)
  for(i in names(abbr)){
    x <- str_replace_all(x, paste0("\\b", i, "\\b"), abbr[i])
  }
  x <- str_replace_all(x, "\\btumours\\b", "tumor")
  x <- str_replace_all(x, "\\btumors\\b", "tumor")
  x <- str_replace_all(x, "\\bleukaemia\\b", "leukemia")
  x <- str_replace_all(x, "\\bhaematology\\b", "hematology")
  x <- str_replace_all(x, "\\bhaematological\\b", "hematological")
  x <- str_replace_all(x, "\\bhaematologic\\b", "hematologic")
  x <- str_replace_all(x, "\\bpaediatric\\b", "pediatric")
  x <- str_remove(x, regex("^diagnosis of ", ignore_case = TRUE))
  x <- str_remove(x, regex("^prevention of ", ignore_case = TRUE))
  x <- str_remove(x, "(?i)^.*?treatment of\\s*")
  x <- str_replace_all(x, "[0-9]+", " ")
  x <- str_replace_all(x, "[^a-z0-9 ]+", " ")
  x <- str_replace_all(x, "\\s+", " ")
  x <- str_remove_all(x, "\\bwith variants( of)?\\b") 
  x <- str_remove_all(x, "\\bof\\b")
  x <- str_remove_all(x, "\\bwith\\b")   
  x <- str_remove_all(x, "\\band\\b")                         
  x <- str_remove_all(x, "\\bthe\\b")
  x <- str_remove_all(x, "\\bin\\b")
  x <- str_remove_all(x, "\\bfrom\\b")
  x <- str_remove_all(x, "\\bfor\\b")
  x <- str_remove_all(x, "\\bdiagnosis\\b")
  x <- str_remove_all(x, "\\bor\\b")
  x <- str_remove_all(x, "\\bpatients\\b")
  x <- str_remove_all(x, "\\b[a-z]{1,2}\\b")
  x <- str_remove_all(x, "(?i)\\b(as a )?diagnostic for( the| clinical)? management( of)?\\b")
  x <- str_remove_all(x, "(?i)for use as hematopoietic support\\s*")
  x <- str_remove_all(x, "(?i)diagnostic agent for the clinical management\\s*")
  x <- str_split(x, " ") %>%
    lapply(unique) %>%
    lapply(paste, collapse = " ") %>%
    unlist()
  x <- na_if(x, "NA")
  x <- na_if(x, "")
  x <- str_squish(x)
  return(x)
}

# Sorting function
sort_words <- function(x){
  sapply(x, function(txt){
    if(is.na(txt)) return(NA)
    paste(sort(unlist(strsplit(txt, " "))), collapse = " ")
  })
}

# Plural removal function
normalise_plural <- function(text, terms) {
  for(term in terms) {
    pattern <- paste0("\\b", term, "(s|es)?\\b")
    text <- str_replace_all(text, regex(pattern, ignore_case = TRUE), term)
  }
  text
}

# Read data
rarecare_df <- read_excel("E:/data/RARECAREnet_list_of_rare_cancers.xlsx")
EMA_df <- read_excel("E:/data/EMA_drugs.xlsx")
colnames(EMA_df) <- EMA_df[2, ] %>% as.character()
EMA_df <- EMA_df[-c(1,2), ]
file_path <- "E:/data/Medicines_output_orphan_designations_en.json"
data_list <- fromJSON(file_path)
EMA_2_df <- data_list[[2]]
FDA_df <- read_excel("E:/data/FDA_drugs.xlsx")
CT_EU_old_df <- read.csv("E:/data/clinical_trials_eu_old.csv")
CT_EU_new_df <- read.csv("E:/data/clinical_trials_eu_new.csv")
CT_US_df <- read.csv("E:/data/clinical_trials_us.csv")
CT_EU_new_extra_df <- read.csv("E:/data/clinical_trials_eu_new_extra.csv")
ORPHA_df <- read_excel("E:/data/First Filtering/ORPHAnomenclature_MasterFile_en_2025.xlsx")
name_differenes_df <- read_excel("E:/data/First Filtering/Tumour name differences.xlsx")

# Terms from RARECARE list and https://www.cancer.gov/types
cancer_terms <- c("tumor", "tumour", "tumors", "tumours", "glioma", "cancer", "blastoma", "carcinoma", "paget", "thymoma", "teratoma", "seminoma", "malignancy", "malignancies", "mesothelioma", "melanoma", "blastoma", "sarcoma", "adamantinoma", "mycosis fungoides", "mf", "neoplasia", "sézary syndrome", "carcinoid", "ganglioma", "cytoma", "meningioma", "lymphoma", "leukaemia", "leukemia", "myeloma", "aml", "neoplasm", "myelodysplastic syndrome", "mds", "cml")
abbr <- c(
  "aml" = "acute myeloid leukemia",
  "cml" = "chronic myeloid leukemia",
  "mds" = "myelodysplastic syndrome",
  "mf"  = "mycosis fungoides",
  "amkl" = "acute megakaryoblastic leukemia",
  "cll" = "chronic leukemia lymphocytic"
)

# Filter on only rare, tier 2 tumours
rarecare_filtered_df <- rarecare_df %>%
  filter(`R=rare` == "R") %>%
  filter(Tier == 2) %>%
  select(Tumour, `Crude incidence rate per 100,000`)

# Match with the manually created name differences file
rarecare_orphanet_df <- rarecare_filtered_df %>%
  left_join(name_differenes_df,
            by = c("Tumour" = "RARECARE name"))

# Edit data fame so that only one row per preferred term remains
ORPHA_synonyms <- ORPHA_df %>%
  filter(!is.na(Synonyms)) %>%
  distinct(PreferredTerm, Synonyms) %>%   
  group_by(PreferredTerm) %>%
  mutate(syn_id = row_number()) %>%       
  ungroup() %>%
  pivot_wider(
    names_from = syn_id,
    values_from = Synonyms,
    names_prefix = "Synonym_"
  )

# Combine the RARECARE list with the Orphanet synonyms
rarecare_orphanet_synonyms_df <- rarecare_orphanet_df %>%
  separate_rows(`Orphanet name`, sep = " \\+ ") %>%
  left_join(ORPHA_synonyms,
            by = c("Orphanet name" = "PreferredTerm"),
            relationship = "many-to-many")

# Clean data frame
rare_care_filtered_tumours <- rarecare_orphanet_synonyms_df %>%
  mutate(across(-`Crude incidence rate per 100,000`,
                ~ sort_words(normalise_plural(clean_text(.), cancer_terms)))) %>%
  mutate(Tumour.ID = as.numeric(factor(Tumour))) %>%
  select(where(~ !all(is.na(.)))) %>%
  mutate(across(starts_with("Synonym_"),
                ~ na_if(., "variants"))) %>%
  distinct()

rare_care_names <- rare_care_filtered_tumours %>%
  select(Tumour.ID, Tumour, `Orphanet name`, `Crude incidence rate per 100,000`) %>%
  distinct()
  

# Rename and separate columns, normalise dates and select only relevant columns
EMA_filtered_df <- EMA_df %>%
  rename(
    EMA.ID = `EU #`,
    Generic.Name = Product,
    Designation.Date = `Designation date`) %>%
  separate(`Tradename - EU product # - Implemented on`, into = c("Brandname", "EU.Product.ID", "Implemented.on"), sep = "-", fill = "right") %>%
  mutate(Designation.Date = dmy(str_to_lower(Designation.Date)),
         Implemented.on = dmy(str_to_lower(Implemented.on)),
         str_squish(Implemented.on),
         EU.Product.ID = na_if(str_squish(EU.Product.ID), "")) %>%
  select(EMA.ID, Indication, Sponsor, Generic.Name, Designation.Date, Brandname, EU.Product.ID, Implemented.on)

# Join with normal EMA database, rename columns and edit the orphan designation status
EMA_filtered_df <- EMA_filtered_df %>%
  left_join(EMA_2_df %>%
              select(eu_designation_number, status),
            by = c("EMA.ID" = "eu_designation_number")) %>%
  distinct() %>%
  rename(Orphan.Designation.Status = status) %>%
  mutate(Orphan.Designation.Status = case_when(
    !is.na(Implemented.on) ~ "implemented",
    str_to_lower(str_trim(Orphan.Designation.Status)) == "positive" ~ "designated",
    TRUE ~ Orphan.Designation.Status
  ))

# Filter on cancer terms and clean all text
EMA_filtered_df <- EMA_filtered_df %>%
  filter(str_detect(
    str_to_lower(Indication),
    paste(cancer_terms, collapse = "|"))) %>%
  mutate(
    Indication = sort_words(clean_text(Indication)),
    Indication = normalise_plural(Indication, cancer_terms),
    Sponsor = clean_text(Sponsor),
    Generic.Name = clean_text(Generic.Name),
    Brandname = clean_text(Brandname)
  )

# Rename and separate columns, normalise dates and select only relevant columns
FDA_filtered_df <- FDA_df %>%
  rename(
    FDA.ID = `CF Grid Key`,
    Brandname = `Trade Name`,
    Designation.Date = `Date Designated`,
    Indication = `Orphan Designation`,
    Sponsor = `Sponsor Company`,
    Generic.Name = `Generic Name`,
    Orphan.Designation.Status = `Orphan Designation Status`
  ) %>%
  mutate(
    Designation.Date = case_when(
      str_detect(Designation.Date, "^\\d{1,2}/\\d{1,2}/\\d{2,4}$") ~ mdy(Designation.Date),
      str_detect(Designation.Date, "^\\d+$") ~ as.Date(as.numeric(Designation.Date), origin = "1899-12-30"),
      TRUE ~ as.Date(NA)
    ),
    `Marketing Approval Date` = case_when(
      str_detect(`Marketing Approval Date`, "^\\d{1,2}/\\d{1,2}/\\d{2,4}$") ~ mdy(`Marketing Approval Date`),
      str_detect(`Marketing Approval Date`, "^\\d+$") ~ as.Date(as.numeric(`Marketing Approval Date`), origin = "1899-12-30"),
      TRUE ~ as.Date(NA)
    ),
    `Date Designation Withdrawn or Revoked` = case_when(
      str_detect(`Date Designation Withdrawn or Revoked`, "^\\d{1,2}/\\d{1,2}/\\d{2,4}$") ~ 
        mdy(`Date Designation Withdrawn or Revoked`),
      str_detect(`Date Designation Withdrawn or Revoked`, "^\\d+$") ~ 
        as.Date(as.numeric(`Date Designation Withdrawn or Revoked`), origin = "1899-12-30"),
      TRUE ~ as.Date(NA)
    )
  ) %>%
  select(FDA.ID, Generic.Name, Brandname, Designation.Date, Indication, Orphan.Designation.Status, `Date Designation Withdrawn or Revoked`, `Marketing Approval Date`, Sponsor)

# Filter on cancer terms and clean all text
FDA_filtered_df <- FDA_filtered_df %>%
  filter(str_detect(
    str_to_lower(`Indication`),
    paste(cancer_terms, collapse = "|")
  )) %>%
  mutate(
    Indication = sort_words(clean_text(Indication)),
    Indication = normalise_plural(Indication, cancer_terms),
    Sponsor = clean_text(Sponsor),
    Generic.Name = clean_text(Generic.Name),
    Brandname = clean_text(Brandname),
    # Generalise orphan designation status
    Orphan.Designation.Status = case_when(
      str_detect(str_to_lower(Orphan.Designation.Status), "withdrawn|revoked") ~ "withdrawn",
      str_detect(str_to_lower(Orphan.Designation.Status), "approved") ~ "approved",
      TRUE ~ str_to_lower(Orphan.Designation.Status)
      )
    )

# Join the two orphan drug designation databases based on their shared column names
orphan_designation_df <- full_join(
  EMA_filtered_df,
  FDA_filtered_df,
  by = c("Indication", "Sponsor", "Generic.Name", "Designation.Date", "Brandname", "Orphan.Designation.Status")
) %>%
  select(EMA.ID, FDA.ID, Indication, Generic.Name, Brandname, Orphan.Designation.Status, Designation.Date, `Date Designation Withdrawn or Revoked`, Implemented.on, `Marketing Approval Date`, EU.Product.ID, Sponsor)

# Rename and combine trial phase columns, normalise dates and select only relevant columns
CT_EU_old_filtered_df <- CT_EU_old_df %>%
  filter(e12_meddra_classification.e13_condition_being_studied_is_a_rare_disease == "Yes") %>%
  rename(`CT.EDU.ID` = a2_eudract_number) %>%
  rename(`NCT.ID` = a52_us_nct_clinicaltrialsgov_registry_number) %>%
  rename(`Member.State` = a1_member_state_concerned) %>%
  rename(`Product.Name` = dimp.d31_product_name) %>%
  rename(`Indication` = e11_medical_conditions_being_investigated) %>%
  rename(`Minors.Included` = f11_trial_has_subjects_under_18) %>%
  rename(`Sponsor` = b1_sponsor.b11_name_of_sponsor) %>%
  rename(`Controlled` = e81_controlled) %>%
  rename(`Randomised` = e811_randomised) %>%
  rename(`Open` = e812_open) %>%
  rename(`Paediatric.Plan` = a7_trial_is_part_of_a_paediatric_investigation_plan) %>%
  rename(`Parallel.Group` = e815_parallel_group) %>%
  rename(`Multiple.Member.States` = e85_the_trial_involves_multiple_member_states) %>%
  rename(`Trial.Status` = p_end_of_trial_status) %>%
  rename(`End.of.Trial` = p_date_of_the_global_end_of_the_trial) %>%
  rename(`Start.of.Trial` = .startDate) %>%
  mutate(across(c(`End.of.Trial`, `Start.of.Trial`),
                ~ parse_date_time(.,
                                  orders = c("dmy", "d B Y", "Ymd", "mdY"),
                                  quiet = TRUE))) %>%
  rowwise() %>%
  mutate(`Trial.Phase` = paste(
      c("I","II","III","IV")[which(c(e71_human_pharmacology_phase_i, e72_therapeutic_exploratory_phase_ii, e73_therapeutic_confirmatory_phase_iii, e74_therapeutic_use_phase_iv))],
      collapse = " and ")
  ) %>%
  ungroup() %>%
  mutate(Trial.Pase = na_if(Trial.Phase, ""),
         # Remove mistakes in NCT IDs
         NCT.ID = if_else(
           NCT.ID %in% c("NCT00000000", "NCT12345678"),
           NA_character_,
           NCT.ID
         )
  ) %>%
  select(`CT.EDU.ID`, `NCT.ID`, `Member.State`, `Product.Name`, Indication, `Minors.Included`, Sponsor, Randomised, Open, `Paediatric.Plan`, `Parallel.Group`, `Multiple.Member.States`, `Trial.Status`, `Trial.Phase`, `End.of.Trial`, `Start.of.Trial`)

# Filter on cancer terms
CT_EU_old_filtered_df <- CT_EU_old_filtered_df %>%
  filter(str_detect(
    str_to_lower(Indication),
    paste(cancer_terms, collapse = "|")))

# Merge both new EU CT databases
CT_EU_new_merged_df <- CT_EU_new_df %>%
  left_join(
    CT_EU_new_extra_df,
    by = c("Trial.number" = "X_id")
  )

# Rename columns, normalise dates, clean trial status, phases, regions, and age, and select only relevant columns
CT_EU_new_filtered_df <- CT_EU_new_merged_df %>%
  rename(Rare.Disease = authorizedApplication.authorizedPartI.medicalConditions.isConditionRareDisease) %>%
  rename(CT.EDU.ID = Trial.number) %>%
  rename(Member.State = Location.s..and.recruitment.status) %>%
  rename(Trial.Status = Overall.trial.status) %>%
  rename(Minors.Included = Age.group) %>%
  rename(Number.of.Participants = Number.of.participants.enrolled) %>%
  rename(Trial.Region = Trial.region) %>%
  rename(Indication = Medical.conditions) %>%
  rename(Trial.Phase = Trial.phase) %>%
  rename(Product.Name = Product) %>%
  rename(End.of.Trial = End.date) %>%
  rename(Start.of.Trial = Start.date) %>%
  rename(Sponsor = Sponsor.Co.Sponsors) %>%
  rename(Sponsor.Type = Sponsor.type) %>%
  rename(NCT.ID = authorizedApplication.authorizedPartI.trialDetails.clinicalTrialIdentifiers.secondaryIdentifyingNumbers.nctNumber.number) %>%
  filter(if_any(Rare.Disease, ~ str_detect(str_to_lower(.), "true"))) %>%
  mutate(across(c(End.of.Trial, Start.of.Trial),
                ~ parse_date_time(.,
                                  orders = c("dmy", "d B Y", "Ymd", "mdY"),
                                  quiet = TRUE))) %>%
  mutate(
    Trial.Status = str_to_lower(Trial.Status),
    Trial.Status = str_squish(Trial.Status),
    Trial.Status = case_when(
        str_detect(Trial.Status, "authorised") ~ "authorised",
        str_detect(Trial.Status, "ongoing") ~ "ongoing",
        str_detect(Trial.Status, "recruit") ~ "recruitment",
        str_detect(Trial.Status, "ended") ~ "completed",
        TRUE ~ Trial.Status
      )
  ) %>%
  mutate(
    Minors.Included = case_when(
      str_detect(Minors.Included, "0-17") ~ "TRUE",
      TRUE ~ "FALSE"
    )
  ) %>%
  mutate(
    Trial.Phase = str_to_lower(Trial.Phase),
    Trial.Phase = str_squish(Trial.Phase),
    Trial.Phase = case_when(
      str_detect(Trial.Phase, "phase i and phase ii") ~ "I and II",
      str_detect(Trial.Phase, "phase ii and phase iii") ~ "II and III",
      str_detect(Trial.Phase, "phase iii and phase iv") ~ "III and IV",
      str_detect(Trial.Phase, "phase i") ~ "I",
      str_detect(Trial.Phase, "phase ii") ~ "II",
      str_detect(Trial.Phase, "phase iii") ~ "III",
      str_detect(Trial.Phase, "phase iv") ~ "IV",
      TRUE ~ Trial.Phase
    )
  ) %>%
  mutate(
    Trial.Region = str_to_lower(Trial.Region),
    Trial.Region = str_squish(Trial.Region),
    Trial.Region = case_when(
      str_detect(Trial.Region, "only") ~ "eea",
      TRUE ~ "eea and non eea")) %>%
  mutate(Product.Name = str_to_lower(Product.Name) %>%
           str_squish() %>%
           na_if("n/a") %>%
           na_if("")
  ) %>%
  # Remove NCT ID mistakes
  mutate(
    NCT.ID = if_else(
      NCT.ID %in% c("NCT00000000", "NCT12345678"),
      NA_character_,
      NCT.ID
    )
  ) %>%
  select(CT.EDU.ID, Rare.Disease, Member.State, Trial.Status, Gender, Minors.Included, Number.of.Participants, Trial.Region, Indication, Trial.Phase, Product.Name, Start.of.Trial, End.of.Trial, Sponsor, Sponsor.Type, NCT.ID)

# Filter on cancer terms
CT_EU_new_filtered_df <- CT_EU_new_filtered_df %>%
  filter(str_detect(
    str_to_lower(Indication),
    paste(cancer_terms, collapse = "|")))

# Rename columns, normalise dates, clean trial status, trial phases, and age, and select only relevant columns
CT_US_filtered_df <- CT_US_df %>%
  rename(NCT.ID = NCT.Number) %>%
  rename(Trial.Status = Study.Status) %>%
  rename(Product.Name = Interventions) %>%
  rename(Indication = Conditions) %>%
  rename(Number.of.Participants = Enrollment) %>%
  rename(Gender = Sex) %>%
  rename(Trial.Phase = Phases) %>%
  rename(Sponsor.Type = Funder.Type) %>%
  rename(End.of.Trial = Completion.Date) %>%
  rename(Start.of.Trial = Start.Date) %>%
  rename(Study.Design = Study.Design) %>%
  rename(Minors.Included = Age) %>%
  mutate(
    Start.of.Trial = str_trim(as.character(Start.of.Trial)),
    End.of.Trial   = str_trim(as.character(End.of.Trial))
  ) %>%
  mutate(
    Start.of.Trial = na_if(Start.of.Trial, ""),
    End.of.Trial   = na_if(End.of.Trial, "")
  ) %>%
  mutate(across(c(End.of.Trial, Start.of.Trial),
                ~ parse_date_time(.,
                                  orders = c("dmy", "d B Y", "Ymd", "mdY"),
                                  quiet = TRUE))) %>%
  mutate(
    Trial.Status = str_to_lower(Trial.Status),
    Trial.Status = str_squish(Trial.Status),
    Trial.Status = case_when(
      str_detect(Trial.Status, "recruiting") ~ "recruitment",
      str_detect(Trial.Status, "active") ~ "ongoing",
      str_detect(Trial.Status, "enrolling") ~ "ongoing",
      str_detect(Trial.Status, "terminated") ~ "prematurely ended",
      str_detect(Trial.Status, "withdrawn") ~ "prematurely ended",
      str_detect(Trial.Status, "suspended") ~ "temporarily halted",
      TRUE ~ Trial.Status)
    ) %>%
  mutate(
    Trial.Phase = str_to_lower(Trial.Phase),
    Trial.Phase = str_squish(Trial.Phase),
    Trial.Phase = case_when(
      str_detect(Trial.Phase, "phase1\\|phase2|phase2\\|phase1") ~ "I and II",
      str_detect(Trial.Phase, "phase2\\|phase3|phase3\\|phase2") ~ "II and III",
      str_detect(Trial.Phase, "phase3\\|phase4|phase4\\|phase3") ~ "III and IV",
      str_detect(Trial.Phase, "phase1") ~ "I",
      str_detect(Trial.Phase, "phase2") ~ "II",
      str_detect(Trial.Phase, "phase3") ~ "III",
      str_detect(Trial.Phase, "phase4") ~ "IV",
      TRUE ~ Trial.Phase
    )
  ) %>%
  mutate(
    Minors.Included = case_when(
      str_detect(Minors.Included, "CHILD") ~ "TRUE",
      TRUE ~ "FALSE"
    )
  ) %>%
  select(NCT.ID, Indication, Trial.Status, Number.of.Participants, Gender, Product.Name, Sponsor, Minors.Included, Trial.Phase, Number.of.Participants, Sponsor.Type, Study.Design, Start.of.Trial, End.of.Trial)

# Filter on cancer terms
CT_US_filtered_df <- CT_US_filtered_df %>%
  filter(str_detect(
    str_to_lower(Indication),
    paste(cancer_terms, collapse = "|")
  )) 

# Split the study design column
CT_US_filtered_df <- CT_US_filtered_df %>%
  mutate(Study.Design = str_to_lower(Study.Design)) %>%
  separate_rows(Study.Design, sep = "\\|") %>%
  separate(Study.Design,
           into = c("key", "value"),
           sep = ":",
           fill = "right") %>%
  mutate(
    key = str_squish(key),
    value = str_squish(value)
  ) %>%
  pivot_wider(
    names_from = key,
    values_from = value,
    values_fill = NA
  ) %>%
  mutate(
    randomised = case_when(
      str_detect(allocation, "non_random") ~ FALSE,
      str_detect(allocation, "random") ~ TRUE,
      is.na(allocation) ~ NA
    )
  ) %>%
  select(-masking) %>%
  rename(Model = `intervention model`,
         Purpose = `primary purpose`,
         Randomised = randomised)

# Save all data to files
write.csv(CT_EU_new_filtered_df, "new_CT_EU.csv", row.names = FALSE)
write.csv(CT_EU_old_filtered_df, "old_CT_EU.csv", row.names = FALSE)
write.csv(CT_US_filtered_df, "CT_US.csv", row.names = FALSE)
write.csv(orphan_designation_df, "orphan_designation.csv", row.names = FALSE)
write.csv(rare_care_filtered_tumours, "tumours_synonyms.csv", row.names = FALSE)

# Prioritise Orphanet name over RARECARE name, if available, and save data to a file
write.csv(rare_care_filtered_tumours %>% mutate(Tumour = coalesce(`Orphanet name`, Tumour)), "tumours_names.csv", row.names = FALSE)
