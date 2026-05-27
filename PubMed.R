library(pubmedR)
library(dplyr)
library(tidyr)
library(stringr)

# Get unique Tumour-Drug combinations
orphan_indications_df <- read.csv("E:/data/Matched/orphan_designation.csv") 
orphan_indications_better_df <- read.csv("E:/data/Matched/orphan_designation_better.csv") 
tumours_final <- read.csv("E:/data/Splitted2/tumours_final.csv") %>%
  select(Tumour.ID, Synonym_main, Synonym_orpha) %>%
  # Prefer Orpha synonym over RARECARE name
  mutate(Tumour = coalesce(Synonym_orpha, Synonym_main)) %>%
  select(-Synonym_main, -Synonym_orpha)%>%
  rename(Tumours.ID = Tumour.ID)

# Split AML and APML
split_aml <- orphan_indications_df %>%
  filter(Tumours.ID == "5, 4") %>%
  mutate(Tumours.ID = str_split(Tumours.ID, ",\\s*")) %>%
  unnest(Tumours.ID) %>%
  mutate(Tumours.ID = str_trim(Tumours.ID),
         Indication = if_else(
           Tumours.ID == "4",
           "acute leukemia myeloid",
           Indication
           )
         )

# Bind with other data
orphan_indications <- orphan_indications_df %>%
  filter(Tumours.ID != "5, 4") %>%
  bind_rows(split_aml)
  
# Retrieve unique tumour - drug matches
orphan_indications_unique_combinations <- orphan_indications %>%
  filter(!is.na(Generic.Name),
         !is.na(Drugbank.ID)) %>%                                               # Only those with Drugbank matches
  distinct(Indication, Generic.Name)

tumours_unique <- orphan_indications %>%
  distinct(
    Tumour = Indication,
    Tumours.ID = as.character(Tumours.ID)
  ) %>%
  bind_rows(
    tumours_final %>%
      distinct(
        Tumour,
        Tumours.ID = as.character(Tumours.ID)
      )
  ) %>%
  mutate(Tumours.ID = str_split(Tumours.ID, ",\\s*")) %>%
  unnest(Tumours.ID) %>%
  distinct(Tumours.ID, Tumour)

# Function to make query, include both tumour and drug
make_query <- function(tumour, drug) {
  paste0(
    tumour,
    " AND \"", drug, "\"[tiab]"
  )
}
# Initialise list
results_list <- list()

# Go over all combinations
for (i in 1:nrow(orphan_indications_unique_combinations)) {
  
  tumour <- orphan_indications_unique_combinations$Indication[i]
  drug   <- orphan_indications_unique_combinations$Generic.Name[i]
  
  # Make query
  query <- make_query(tumour, drug)
  
  # Return number of papers
  count <- tryCatch({
    pmQueryTotalCount(query = make_query(tumour, drug))
    
  }, error = function(e1) {
    return(0)   
  })
    
  # Bind results to data frame
  results_list[[i]] <- data.frame(
    Tumour = tumour,
    Generic.Name = drug,
    PubMed_Count = count,
    Query = query
  )
  
  # Progress
  print(i)  
  
  # Prevent rate limiting
  Sys.sleep(0.5)
}

# Bind results
results <- bind_rows(results_list)

# Only select important columns
results_select <- results %>%
  select(Tumour, Generic.Name, PubMed_Count.total_count) %>%
  rename(PubMed.Count = PubMed_Count.total_count) %>%
  # Join with original database
  right_join(orphan_indications, by = c("Tumour" = "Indication", "Generic.Name")) %>%
  mutate(PubMed.Count = replace_na(PubMed.Count, 0))

# Function to make query, only tumour
make_query <- function(tumour) {
  query <- tumour
  return(query)
}

# Initialise list
results_list2 <- list()

# Go over all combinations
for (i in 1:nrow(tumours_unique)) {
  
  Tumours_id <- tumours_unique$Tumours.ID[i]
  tumour <- tumours_unique$Tumour[i]
  
  # Make query
  query <- make_query(tumour)
  
  # Return number of papers
  count <- tryCatch({
    pmQueryTotalCount(query = make_query(tumour))
    
  }, error = function(e1) {
    return(0)   
  })
  
  # Bind results to data frame
  results_list2[[i]] <- data.frame(
    Tumours.ID = Tumours_id,
    Tumour = tumour,
    PubMed_Count = count,
    Query = query
  )
  
  # Progress
  print(i)  
  
  # Prevent rate limiting
  Sys.sleep(0.5)
}

# Bind results
results2 <- bind_rows(results_list2)

# Only select important columns
results_select2 <- results2 %>%
  select(Tumours.ID, Tumour, PubMed_Count.total_count) %>%
  rename(PubMed.Count.General = PubMed_Count.total_count) %>%
  mutate(PubMed.Count.General = replace_na(PubMed.Count.General, 0)) %>%
  group_by(Tumours.ID) %>%
  summarise(
    PubMed.Count.General = sum(PubMed.Count.General),
    .groups = "drop"
  )

# Save data
write.csv(results_select, "orphan_designation_better.csv", row.names = FALSE)
write.csv(results_select2, "PubMed.csv", row.names = FALSE)
