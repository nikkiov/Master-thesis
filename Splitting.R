library(readxl)
library(dplyr)
library(openxlsx)
library(stringr)
library(tidyr)
library(lubridate)
library(dplyr)
library(stringdist)

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
  x <- str_remove(x, regex("^treatment of ", ignore_case = TRUE))
  x <- str_remove(x, regex("^diagnosis of ", ignore_case = TRUE))
  x <- str_remove(x, regex("^prevention of ", ignore_case = TRUE))
  x <- str_remove(x, "(?i)^.*?treatment of\\s*")
  x <- str_remove(x, "(?i)^.*\\bin\\b\\s*")
  x <- str_replace_all(x, "[0-9]+", " ")
  x <- str_replace_all(x, "[^a-z0-9 ]+", " ")
  x <- str_replace_all(x, "\\s+", " ")
  x <- str_remove_all(x, "\\bwith variants of\\b") 
  x <- str_remove_all(x, "\\bwith variants\\b")   
  x <- str_remove_all(x, "\\bof\\b")
  x <- str_remove_all(x, "\\bwith\\b")   
  x <- str_remove_all(x, "\\band\\b")                         
  x <- str_remove_all(x, "\\bthe\\b")
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
new_CT_EU_df <- read.csv("E:/data/First Filtering2/new_CT_EU.csv")
old_CT_EU_df <- read.csv("E:/data/First Filtering2/old_CT_EU.csv")
CT_US_df <- read.csv("E:/data/First Filtering2/CT_US.csv")

# Cancer terms from RARECARE list and https://www.cancer.gov/types
cancer_terms <- c("tumor", "tumour", "tumors", "tumours", "glioma", "cancer", "blastoma", "carcinoma", "paget", "thymoma", "teratoma", "seminoma", "malignancy", "malignancies", "mesothelioma", "melanoma", "blastoma", "sarcoma", "adamantinoma", "mycosis fungoides", "mf", "neoplasia", "sézary syndrome", "carcinoid", "ganglioma", "cytoma", "meningioma", "lymphoma", "leukaemia", "leukemia", "myeloma", "aml", "neoplasm", "myelodysplastic syndrome", "mds", "cml")
abbr <- c(
  "aml" = "acute myeloid leukemia",
  "cml" = "chronic myeloid leukemia",
  "mds" = "myelodysplastic syndrome",
  "mf"  = "mycosis fungoides",
  "amkl" = "acute megakaryoblastic leukemia",
  "cll" = "chronic leukemia lymphocytic"
)

# Separate product fields (including product type and product name) and clean it
CT_US_product <- CT_US_df %>%
  select(NCT.ID, Product.Name) %>%
  separate_rows(Product.Name, sep = "\\|") %>%
  separate(Product.Name,
           into = c("Product.Type", "Product.Name"),
           sep = ":",
           fill = "right",
           extra = "merge") %>%
  mutate(
    Product.Type = str_squish(Product.Type),
    Product.Name = str_squish(Product.Name)
  ) %>%
  distinct(NCT.ID, Product.Type, Product.Name, .keep_all = TRUE) %>%
  mutate(Product.Type = clean_text(Product.Type),
         Product.Name = clean_text(Product.Name)) %>%
  # Only include biological products
  filter(Product.Type %in% c("drug", "biological"),
         !is.na(Product.Name))

# Separate indication fields, first filter and clean the indication field
CT_US_indication <- CT_US_df %>%
  select(NCT.ID, Indication) %>%
  separate_rows(Indication, sep = "\\|") %>%
  filter(str_detect(
    str_to_lower(Indication),
    paste(cancer_terms, collapse = "|"))) %>%
  mutate(Indication = sort_words(clean_text(Indication)),
         Indication = normalise_plural(Indication, cancer_terms)) %>%
  # Only include IDs from drug trials
  semi_join(CT_US_product, by = "NCT.ID")

# Clean sponsor fields
CT_US_sponsor <- CT_US_df %>%
  select(NCT.ID, Sponsor, Sponsor.Type) %>%
  mutate(
    Sponsor = clean_text(Sponsor),
    Sponsor.Type = str_to_lower(Sponsor.Type)
  ) %>%
  mutate(across(everything(), str_squish)) %>%
  distinct(NCT.ID, Sponsor, Sponsor.Type, .keep_all = TRUE) %>%
  # Only include IDs from drug trials
  semi_join(CT_US_product, by = "NCT.ID")

# Clean trial information fields
CT_US_trial_information <- CT_US_df %>%
  select(NCT.ID, Trial.Status, Number.of.Participants, Gender, Minors.Included, Trial.Phase, Start.of.Trial, End.of.Trial, Model, Purpose, Randomised) %>%
  mutate(across(c(Trial.Status, Minors.Included, Gender, Model, Purpose, Randomised), clean_text)) %>%
  # Only include IDs from drug trials
  semi_join(CT_US_product, by = "NCT.ID")

# Separate product fields and clean it
CT_EU_old_product <- old_CT_EU_df %>%
  select(CT.EDU.ID, NCT.ID, Product.Name) %>%
  mutate(
    Product.Name = str_to_lower(Product.Name),
    Product.Name = str_split(Product.Name, " / "),
  ) %>%
  mutate(Product.Name = clean_text(Product.Name)) %>%
  distinct(CT.EDU.ID, NCT.ID, Product.Name, .keep_all = TRUE) %>%
  # Only include drugs
  filter(!is.na(Product.Name))

# Separate indication fields, first filter and clean the indication field
CT_EU_old_indication <- old_CT_EU_df %>%
  select(CT.EDU.ID, NCT.ID, Indication) %>%
  separate_rows(Indication, sep = "\\|") %>%
  filter(str_detect(
    str_to_lower(Indication),
    paste(cancer_terms, collapse = "|"))) %>%
  mutate(Indication = sort_words(clean_text(Indication)),
         Indication = normalise_plural(Indication, cancer_terms)) %>%
  distinct(CT.EDU.ID, Indication, .keep_all = TRUE) %>%
  # Only include IDs from drug trials
  semi_join(CT_EU_old_product, by = "CT.EDU.ID")

# Clean sponsor fields
CT_EU_old_sponsor <- old_CT_EU_df %>%
  select(CT.EDU.ID, NCT.ID, Sponsor) %>%
  mutate(
    Sponsor = clean_text(Sponsor)
  ) %>%
  mutate(across(everything(), str_squish)) %>%
  distinct(CT.EDU.ID, Sponsor, .keep_all = TRUE) %>%
  # Only include IDs from drug trials
  semi_join(CT_EU_old_product, by = "CT.EDU.ID")

# Clean trial information fields
CT_EU_old_trial_information <- old_CT_EU_df %>%
  select(CT.EDU.ID, NCT.ID, Member.State, Trial.Status, Minors.Included, Trial.Phase, Start.of.Trial, End.of.Trial, Open, Parallel.Group, Randomised, Multiple.Member.States) %>%
  rename(Model = Parallel.Group) %>%
  mutate(Model = if_else(Model == TRUE, "parallel", NA_character_),
         Member.State = Member.State %>% 
           str_split(" - ", simplify = TRUE) %>%
           .[,1] %>%                                                            # Only retrieve information in front of the -
           str_to_lower(),
         Member.State = if_else(
           str_detect(Member.State, "czech republic"),
           "czechia",
           Member.State)
         ) %>%
  mutate(across(c(Trial.Status, Minors.Included, Open, Randomised, Multiple.Member.States), clean_text)) %>%
  # Only include IDs from drug trials
  semi_join(CT_EU_old_product, by = "CT.EDU.ID")

# Separate product fields and clean it
CT_EU_new_product <- new_CT_EU_df %>%
  select(CT.EDU.ID, NCT.ID, Product.Name) %>%
  mutate(
    Product.Name = str_to_lower(Product.Name)
  ) %>%
  distinct(CT.EDU.ID, NCT.ID, Product.Name, .keep_all = TRUE) %>%
  mutate(Product.Name = clean_text(Product.Name)) %>%
  # Only include drugs
  filter(!is.na(Product.Name))

# Separate indication fields, first filter and clean the indication field
CT_EU_new_indication <- new_CT_EU_df %>%
  select(CT.EDU.ID, NCT.ID, Indication) %>%
  separate_rows(Indication, sep = "\\,") %>%
  filter(str_detect(
    str_to_lower(Indication),
    paste(cancer_terms, collapse = "|"))) %>%
  mutate(Indication = sort_words(clean_text(Indication)),
         Indication = normalise_plural(Indication, cancer_terms)) %>%
  # Only include IDs from drug trials
  semi_join(CT_EU_new_product, by = "CT.EDU.ID")

# Separate sponsor fields and clean the sponsor fields
CT_EU_new_sponsor <- new_CT_EU_df %>%
  select(CT.EDU.ID, NCT.ID, Sponsor, Sponsor.Type) %>%
  mutate(
    Sponsor = str_to_lower(Sponsor),
    Sponsor.Type = str_to_lower(Sponsor.Type)
  ) %>%
  mutate(
    Sponsor = str_split(Sponsor, ","),
    Sponsor.Type = str_split(Sponsor.Type, ","),
  ) %>%
  unnest_longer(Sponsor.Type) %>%
  unnest_longer(Sponsor) %>%
  mutate(across(everything(), str_squish)) %>%
  distinct(CT.EDU.ID, NCT.ID, Sponsor, Sponsor.Type, .keep_all = TRUE) %>%
  mutate(Sponsor.Type = clean_text(Sponsor.Type),
         Sponsor = clean_text(Sponsor)) %>%
  # Only include IDs from drug trials
  semi_join(CT_EU_new_product, by = "CT.EDU.ID")

# Clean trial information fields
CT_EU_new_trial_information <- new_CT_EU_df %>%
  select(CT.EDU.ID, NCT.ID, Member.State, Trial.Status, Gender, Minors.Included, Number.of.Participants, Trial.Region, Trial.Phase, Start.of.Trial, End.of.Trial) %>%
  mutate(Gender = if_else(Gender == "Female, Male", "all", str_to_lower(Gender)),
         Member.State = str_extract_all(Member.State, "[^:,]+(?=:)") %>%        # Function that removes everything after the : 
             sapply(function(x) paste(str_to_lower(trimws(x)), collapse = ", "))
         ) %>%
  mutate(across(c(Trial.Status, Trial.Region, Gender, Minors.Included), clean_text)) %>%
  # Only include IDs from drug trials
  semi_join(CT_EU_new_product, by = "CT.EDU.ID")

# Combine EU trial data
EU_combined_trial <- bind_rows(CT_EU_new_trial_information, CT_EU_old_trial_information)

# Join all trial information together, take EU information over US information
Trial_information_df <- CT_US_trial_information %>%
  full_join(
    EU_combined_trial,
    by = "NCT.ID",
    suffix = c("_US", "_EU")
  ) %>%
  mutate(
    Trial.Status = coalesce(Trial.Status_EU, Trial.Status_US),
    Minors.Included = coalesce(Minors.Included_EU, Minors.Included_US),
    Trial.Phase = coalesce(Trial.Phase_EU, Trial.Phase_US),
    Trial.Phase = na_if(Trial.Phase, ""),
    Gender = coalesce(Gender_EU, Gender_US),
    Start.of.Trial = coalesce(Start.of.Trial_EU, Start.of.Trial_US),
    End.of.Trial = coalesce(End.of.Trial_EU, End.of.Trial_US),
    Number.of.Participants = coalesce(Number.of.Participants_EU, Number.of.Participants_US),
    Model = coalesce(Model_EU, Model_US),
    Randomised = coalesce(Randomised_EU, Randomised_US),
    End.of.Trial = if_else(as.Date(End.of.Trial) < as.Date(Start.of.Trial), as.Date(NA), as.Date(End.of.Trial))     # If start date is later than end date, remove trial end date
  ) %>%
  distinct(NCT.ID, CT.EDU.ID, Member.State, Trial.Status, Gender, Minors.Included,
           Number.of.Participants, Trial.Phase, Start.of.Trial, End.of.Trial,
           Model, Randomised, Purpose, Open, Trial.Region, Multiple.Member.States, .keep_all = TRUE) %>%
  select(NCT.ID, CT.EDU.ID, Member.State, Trial.Status, Gender, Minors.Included,
         Number.of.Participants, Trial.Phase, Start.of.Trial, End.of.Trial,
         Model, Randomised, Purpose, Open, Trial.Region, Multiple.Member.States)

# Combine EU sponsor data
EU_combined_sponsor <- bind_rows(CT_EU_new_sponsor, CT_EU_old_sponsor)

# Join all sponsor information together, take EU over US
Sponsor_df <- CT_US_sponsor %>%
  full_join(
    EU_combined_sponsor,
    by = "NCT.ID",
    suffix = c("_US", "_EU")
  ) %>%
  mutate(
    Sponsor = coalesce(Sponsor_EU, Sponsor_US),
    Sponsor.Type = if_else(!is.na(Sponsor.Type_EU), Sponsor.Type_EU, Sponsor.Type_US)
  ) %>%
  distinct(NCT.ID, CT.EDU.ID, Sponsor.Type, Sponsor, .keep_all = TRUE) %>%
  select(NCT.ID, CT.EDU.ID, Sponsor.Type, Sponsor)

# Combine EU indication data
EU_combined_indication <- bind_rows(CT_EU_new_indication, CT_EU_old_indication)

# Join all indication information together, take US information over EU information (less noise)
Indication_df <- CT_US_indication %>%
  full_join(
    EU_combined_indication,
    by = "NCT.ID",
    suffix = c("_US", "_EU"),
    relationship = "many-to-many"
  ) %>%
  mutate(
    Indication = coalesce(Indication_US, Indication_EU)
  ) %>%
  distinct(NCT.ID, CT.EDU.ID, Indication, .keep_all = TRUE) %>%
  select(NCT.ID, CT.EDU.ID, Indication) %>%
  filter(str_detect(str_to_lower(Indication), paste(cancer_terms, collapse = "|")))

# Combine EU product data
EU_combined_product <- bind_rows(CT_EU_new_product, CT_EU_old_product)

# Join all product information together, take EU information over US information
Product_df <- CT_US_product %>%
  full_join(
    EU_combined_product,
    by = "NCT.ID",
    suffix = c("_US", "_EU"),
    relationship ="many-to-many"
  ) %>%
  mutate(
    Product.Name = coalesce(Product.Name_US, Product.Name_EU)
  ) %>%
  filter(!is.na(Product.Name)) %>%
  distinct(NCT.ID, CT.EDU.ID, Product.Type, Product.Name, .keep_all = TRUE) %>%
  select(NCT.ID, CT.EDU.ID, Product.Type, Product.Name)

# Save data
write.csv(Product_df, "Products.csv", row.names = FALSE)
write.csv(Indication_df, "Indications.csv", row.names = FALSE)
write.csv(Sponsor_df, "Sponsors.csv", row.names = FALSE)
write.csv(Trial_information_df, "Trial_information.csv", row.names = FALSE)

