library(stringr)
library(stringdist)
library(purrr)
library(lubridate)
library(dplyr)
library(purrr)
library(tidyr)
library(fuzzyjoin)
library(text2vec)

# Function that cleans data
clean_text <- function(x){
  x <- tolower(x)
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
  x <- str_remove_all(x, "\\bin\\b")
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

# Function that sorts words
sort_words <- function(x){
  sapply(x, function(txt){
    if(is.na(txt)) return(NA)
    paste(sort(unlist(strsplit(txt, " "))), collapse = " ")
  })
}

# Read files
Products_df <- read.csv("data/Splitted2/Products.csv")
Orphan_designations_df <- read.csv("data/Splitted2/orphan_designation.csv") %>%
  mutate(EU.Product.ID = str_to_lower(EU.Product.ID))
drugbank_df <- read.csv("Drugbank.csv")

# Clean up the drugbank database
drugbank_clean <- drugbank_df %>%
  mutate(
    name = str_to_lower(str_squish(name)),
    synonyms = str_to_lower(str_squish(synonyms)),
    products = str_to_lower(str_squish(products)),
    ema_ma_number = str_to_lower(str_squish(ema_ma_number)),
    categories = str_remove(categories, "D0.*")
  ) %>%
  rowwise() %>%
  # Split all columns
  mutate(
    synonyms_list = list(str_split(synonyms, ";\\s*")[[1]]),
    products_list = list(str_split(products, ";\\s*")[[1]]), 
    ema_list = list(str_split(ema_ma_number, ";\\s*")[[1]]),
    # Include all possible names for one drug in one column
    all_names = list(unique(c(name, unlist(synonyms_list), unlist(products_list))))
  ) %>%
  ungroup() %>%
  unnest(all_names) %>%
  group_by(drugbank_id) %>%
  # Remove duplicate EMA numbers
  mutate(ema_ma_number = list(unique(unlist(ema_list)))) %>%
  ungroup() %>%
  unnest(ema_ma_number) %>%
  select(drugbank_id, all_names, classification, groups, ema_ma_number) %>%
  # Remove drugs with no name(s)
  filter(
    !is.na(all_names) & all_names != ""
  ) %>%
  # Remove duplicates
  distinct() %>%
  # Remove last numbers
  mutate(ema_ma_number = str_remove(ema_ma_number, "/\\d{3}$")) %>%
  # Remove latter row(s) where one has an EMA number while the other(s) do(es) not
  distinct(drugbank_id, all_names, classification, groups, ema_ma_number, .keep_all = FALSE) %>%
  group_by(drugbank_id, all_names) %>%
  filter(!(is.na(ema_ma_number) & any(!is.na(ema_ma_number)))) %>%
  ungroup()

# List all EMA numbers and names
drugbank_wide <- drugbank_clean %>%
  group_by(drugbank_id) %>%
  summarise(
    all_ema_numbers = list(unique(ema_ma_number)),
    all_names = list(unique(all_names))
  ) %>%
  ungroup()

# Find max length of EMA numbers and names
max_ema <- max(lengths(drugbank_wide$all_ema_numbers))
max_names <- max(lengths(drugbank_wide$all_names))

# Give every unique EMA number / drug name its own field
drugbank_final <- drugbank_wide %>%
  mutate(
    ema = map(all_ema_numbers, ~ {
      x <- .x
      length(x) <- max_ema
      return(x)
    })
  ) %>%
  unnest_wider(ema, names_sep = "_") %>%
  rename_with(~ paste0("ema_", 1:max_ema), starts_with("ema_")) %>%
  mutate(
    names = map(all_names, ~ {
      x <- .x
      length(x) <- max_names
      return(x)
    })
  ) %>%
  unnest_wider(names, names_sep = "_") %>%
  rename_with(~ paste0("name_", 1:max_names), starts_with("names_")) %>%
  select(-all_ema_numbers, -all_names)

# Add all names to one column
drugbank_all_names <- drugbank_final %>%
  pivot_longer(
    cols = starts_with("name_"),
    names_to = "name_column",
    values_to = "name",
    values_drop_na = TRUE
  ) %>%
  select(drugbank_id, name) %>%
  distinct() %>%
  # Clean and sort drug names
  mutate(
    name_clean = clean_text(name)
  ) %>%
  # Remove empty fields
  filter(!is.na(name_clean) & name_clean != "") %>%
  select(drugbank_id, name = name_clean) %>%
  distinct()

# Take first row of the cleaned Drugbank, select fields and clean the text
drugbank_unique <- drugbank_clean %>%
  group_by(drugbank_id) %>%
  slice(1) %>%
  ungroup() %>%
  select(drugbank_id, 
         Drugbank_official_name = all_names,
         classification, 
         groups) %>%
  mutate(Drugbank_official_name = clean_text(Drugbank_official_name),
         classification = str_to_lower(classification))

# Add all EMA numbers to one column
ema_mapping <- drugbank_final %>%
  select(drugbank_id, starts_with("ema_")) %>%
  pivot_longer(
    cols = starts_with("ema_"),
    names_to = "ema_column",
    values_to = "ema_number",
    values_drop_na = TRUE
  ) %>%
  select(drugbank_id, ema_number) %>%
  distinct() %>%
  filter(!is.na(ema_number) & ema_number != "")

# Use batches of size 100 to increase performance
batch_size <- 100
batches <- split(1:nrow(Orphan_designations_df), ceiling(seq_len(nrow(Orphan_designations_df)) / batch_size))
match_list_orphan <- list()

for(i in seq_along(batches)) {
  message("Batch ", i, " of ", length(batches))
  batch_idx <- batches[[i]]
  
  # Go over every orphan designation drug name within the batch
  batch_result <- Orphan_designations_df[batch_idx, ] %>%
    rowwise() %>%
    mutate(
      # Create a list with match scores
      match = list({
        # If there is an EU Product ID, try to match on that first
        if(!is.na(EU.Product.ID) && EU.Product.ID != "") {
          ema_match <- ema_mapping %>%
            filter(ema_number == EU.Product.ID)
          
          # If there is a match
          if(nrow(ema_match) > 0) {
            tibble(
              fuzzy_drugbank_id = ema_match$drugbank_id[1],
              # Get first name = Generic name
              fuzzy_name = drugbank_all_names %>%
                filter(drugbank_id == ema_match$drugbank_id[1]) %>%
                pull(name) %>%
                first(),
              fuzzy_similarity = 1.0,
              match_method = "EMA_number"
            )
          # If there is no match, go to distance matching
          } else {
            names_to_try <- c(Generic.Name, Brandname)                          # Try multiple drug names
            names_to_try <- names_to_try[!is.na(names_to_try) & names_to_try != ""]
            
            if(length(names_to_try) > 0) {
              best_dist <- Inf
              best_index <- NA
              best_source <- NA_character_
              
              # Go over names
              for(nm in names_to_try) {
                dist <- stringdist::stringdist(nm, drugbank_all_names$name, method = "jw")
                best <- which.min(dist)                                         # Shortest distance
                
                if(length(best) > 0 && dist[best] < best_dist) {                # Update best match
                  best_dist <- dist[best]
                  best_index <- best
                  best_source <- nm
                }
              }
              
              if(!is.na(best_index) && best_dist < 0.25) {                      # Thresholds for acceptable match
                tibble(
                  fuzzy_drugbank_id = drugbank_all_names$drugbank_id[best_index],
                  fuzzy_name = drugbank_all_names$name[best_index],
                  fuzzy_similarity = 1 - best_dist
                )
                # If no acceptable distance --> No match
              } else {
                tibble(
                  fuzzy_drugbank_id = NA_character_,
                  fuzzy_name = NA_character_,
                  fuzzy_similarity = NA_real_
                )
              }
              # If there is no name available --> No match
            } else {
              tibble(
                fuzzy_drugbank_id = NA_character_,
                fuzzy_name = NA_character_,
                fuzzy_similarity = NA_real_,
                match_method = "no_name"
              )
            }
          }
          
          # If no EU Product ID --> Fuzzy matching
        } else {
          names_to_try <- c(Generic.Name, Brandname)                            # Try multiple drug names
          names_to_try <- names_to_try[!is.na(names_to_try) & names_to_try != ""]
          
          if(length(names_to_try) > 0) {
            best_dist <- Inf
            best_index <- NA
            
            # Go over names
            for(nm in names_to_try) {
              dist <- stringdist::stringdist(nm, drugbank_all_names$name, method = "jw")
              best <- which.min(dist)                                           # Shortest distance
              
              if(length(best) > 0 && dist[best] < best_dist) {                  # Update best match
                best_dist <- dist[best]                                         
                best_index <- best
              }
            }
            
            if(!is.na(best_index) && best_dist < 0.25) {                        # Thresholds for acceptable match
              tibble(
                fuzzy_drugbank_id = drugbank_all_names$drugbank_id[best_index],
                fuzzy_name = drugbank_all_names$name[best_index],
                fuzzy_similarity = 1 - best_dist
              )
              # If no acceptable distance --> No match
            } else {
              tibble(
                fuzzy_drugbank_id = NA_character_,
                fuzzy_name = NA_character_,
                fuzzy_similarity = NA_real_
              )
            }
            # If there is no name available --> No match
          } else {
            tibble(
              fuzzy_drugbank_id = NA_character_,
              fuzzy_name = NA_character_,
              fuzzy_similarity = NA_real_
            )
          }
        }
      })
    ) %>%
    unnest(match) %>%
    select(EMA.ID, FDA.ID, Generic.Name, Brandname, EU.Product.ID,
           fuzzy_drugbank_id, fuzzy_name, fuzzy_similarity, match_method)
  
  # Save batch to list
  match_list_orphan[[i]] <- batch_result
}

# Combine all rows
matches_orphan <- bind_rows(match_list_orphan)

# Filter matches
matches_orphan <- matches_orphan %>%
  filter(fuzzy_similarity > 0.96)

# Match the matches with Drugbank data based on Drugbank ID
matches_enhanced_orphan <- matches_orphan %>%
  left_join(drugbank_unique, by = c("fuzzy_drugbank_id" = "drugbank_id")) %>%
  # Take data from Drugbank if there is a match
  mutate(
    Drugbank_official_name = ifelse(is.na(fuzzy_drugbank_id), NA, Drugbank_official_name),
    classification = ifelse(is.na(fuzzy_drugbank_id), NA, classification),
    groups = ifelse(is.na(fuzzy_drugbank_id), NA, groups),
    # Take generic name from orphan drug designation database if there is no Drugbank match
    Drugbank_official_name = ifelse(
      is.na(Drugbank_official_name) | Drugbank_official_name == "",
      Generic.Name,
      Drugbank_official_name
    )
  ) %>%
  select(EMA.ID, FDA.ID, Drugbank_official_name, fuzzy_drugbank_id, classification, groups) %>%
  rename(Generic.Name = Drugbank_official_name,
         Drugbank.ID = fuzzy_drugbank_id,
         Classification = classification,
         Groups = groups)

# Use batches of size 100 to increase performance
batch_size <- 100
batches <- split(1:nrow(Products_df), ceiling(seq_len(nrow(Products_df)) / batch_size))
match_list_products <- list()

for(i in seq_along(batches)) {
  message("Batch ", i, " of ", length(batches))
  batch_idx <- batches[[i]]
  batch <- Products_df[batch_idx, ]
  
  # Go over every clinical trial drug name within the batch
  batch_result <- batch %>%
    rowwise() %>%
    mutate(
      # Create a list with match scores
      match = list({
        # Distance matching between Drugbank names and clinical trial name
        if(!is.na(Product.Name) && Product.Name != "") {
          dist <- stringdist::stringdist(Product.Name, drugbank_all_names$name, method = "jw")
          best <- which.min(dist)                                               # Shortest distance
          
          if(length(best) > 0 && dist[best] < 0.25) {                           # Thresholds for acceptable match
            # Find Drugbank ID, generic name and distance
            tibble(
              fuzzy_drugbank_id = drugbank_all_names$drugbank_id[best],
              fuzzy_name = drugbank_all_names$name[best],
              fuzzy_similarity = 1 - dist[best]
            )
            # If no acceptable distance --> No match
          } else {
            tibble(
              fuzzy_drugbank_id = NA_character_,
              fuzzy_name = NA_character_,
              fuzzy_similarity = NA_real_
            )
          }
          # If there is no generic name
        } else {
          tibble(
            fuzzy_drugbank_id = NA_character_,
            fuzzy_name = NA_character_,
            fuzzy_similarity = NA_real_
          )
        }
      })
    ) %>%
    # Change matches list to normal rows
    unnest(match) %>%
    select(NCT.ID, CT.EDU.ID, Product.Type, Product.Name,
           fuzzy_drugbank_id, fuzzy_name, fuzzy_similarity)
  
  # Save batch results to list
  match_list_products[[i]] <- batch_result
}

# Combine all rows from the list
matches_products <- bind_rows(match_list_products)

# Filter matches
matches_products <- matches_products %>%
  filter(fuzzy_similarity > 0.96)

# Match the matches with Drugbank data based on Drugbank ID
matches_enhanced_products <- matches_products %>%
  left_join(drugbank_unique, by = c("fuzzy_drugbank_id" = "drugbank_id")) %>%
  # Take data from Drugbank if there is a match
  mutate(
    Drugbank_official_name = ifelse(is.na(fuzzy_drugbank_id), NA, Drugbank_official_name),
    classification = ifelse(is.na(fuzzy_drugbank_id), NA, classification),
    groups = ifelse(is.na(fuzzy_drugbank_id), NA, groups),
    # Take generic name from clinical trial database if there is no Drugbank match
    Drugbank_official_name = ifelse(
      is.na(Drugbank_official_name) | Drugbank_official_name == "",
      Product.Name,
      Drugbank_official_name
    )
  ) %>%
  select(NCT.ID, CT.EDU.ID, Product.Type, Drugbank_official_name, fuzzy_drugbank_id, classification, groups) %>%
  rename(Generic.Name = Drugbank_official_name,
         Drugbank.ID = fuzzy_drugbank_id,
         Classification = classification,
         Groups = groups)

# Save files
write.csv(matches_enhanced_orphan, "matches_orphan_products.csv", row.names = FALSE)
write.csv(matches_enhanced_products, "matches_products.csv", row.names = FALSE)
