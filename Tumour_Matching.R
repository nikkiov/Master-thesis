library(stringr)
library(stringdist)
library(purrr)
library(lubridate)
library(dplyr)
library(purrr)
library(tidyr)
library(fuzzyjoin)
library(text2vec)

# Read files
Indications_df <- read.csv("data/Splitted2/Indications.csv")
Orphan_designations_df <- read.csv("data/Splitted2/orphan_designation.csv")
tumours_df <- read.csv("data/Splitted2/tumours_synonyms.csv")

# Rename columns
tumours_df <- tumours_df %>%
  mutate(Synonym_main = Tumour) %>%
  rename(Synonym_orpha = Orphanet.name)

# Put all tumour names in one column
tumours_all_names <- tumours_df %>%
  pivot_longer(
    cols = starts_with("Synonym_"),
    names_to = "Synonym_column",
    values_to = "Synonym",
    values_drop_na = TRUE
  ) %>%
  select(Tumour.ID, Synonym) %>%
  filter(!is.na(Synonym) & Synonym != "") %>%
  # Remove duplicates
  distinct() %>%
  mutate(Synonym = str_squish(Synonym))

# Select tumours from orphan designation databases
Orphan_designations_names <- Orphan_designations_df %>%
  select(EMA.ID, FDA.ID, Indication) %>%
  mutate(row_id = row_number())

# Make one large vector with all text that should be compared and splits these sentences into words (tokens)
all_terms <- c(tumours_all_names$Synonym, Orphan_designations_names$Indication)
it <- itoken(all_terms, 
             tokenizer = word_tokenizer, 
             progressbar = FALSE)

# Create a vocabulary and a mapping for it
vocab <- create_vocabulary(it)
vectorizer <- vocab_vectorizer(vocab)

# Create a Document-Term Matrix that shows how often words appear in certain sentences
dtm <- create_dtm(it, vectorizer)

# Use Term Frequency-Inverse Document Frequency to determine the importance of words --> common words are less important
tfidf <- TfIdf$new()
dtm_tfidf <- fit_transform(dtm, tfidf)

# Split matrices in RARECARE tumours and orphan designation tumours
n_tumours <- nrow(tumours_all_names)
tumour_vectors <- dtm_tfidf[1:n_tumours, ]
indication_vectors <- dtm_tfidf[(n_tumours + 1):nrow(dtm_tfidf), ]

# Calculate cosine similarity between indication vector and tumour vector using a pairwise similarity matrix
cosine_sim <- as.matrix(sim2(indication_vectors, tumour_vectors, method = "cosine"))

# Use batches of size 100 to increase performance
batch_size <- 100
batches <- split(1:nrow(Orphan_designations_names), ceiling(seq_len(nrow(Orphan_designations_names)) / batch_size))
match_list_orphan <- list()                                                     # Save batch results in a list

for(i in seq_along(batches)) {
  message("Batch ", i, " of ", length(batches))
  batch_idx <- batches[[i]]
  
  # Go over every orphan designation tumour name within the batch
  batch_result <- Orphan_designations_names[batch_idx, ] %>%
    rowwise() %>%
    mutate(
      # Create a list with match scores
      matches = list({
        # Cosine scores
        cos_scores <- cosine_sim[row_id, ]
        # Split the tumour name (orphan designation) in words
        ind_words <- str_split(Indication, "\\s+")[[1]]
        # Put scores per tumour row
        scores <- tibble(
          idx = 1:length(cos_scores),
          tumour_id = tumours_all_names$Tumour.ID,
          tumour_name = tumours_all_names$Synonym,
          cosine = cos_scores
        ) %>%
          rowwise() %>%
          mutate(
            # Split the tumour name (RARECARE) in words
            tumour_words = list(str_split(tumour_name, "\\s+")[[1]]),
            # Calculate Jaccard similarity
            jaccard = length(intersect(ind_words, tumour_words)) / 
              length(union(ind_words, tumour_words)),
            # Combine scores and gives it a weight
            combined = 0.4 * cosine + 0.6 * jaccard
          ) %>%
          ungroup() %>%
          # Keep only relatively good matches
          filter(combined > 0.2) %>%
          arrange(desc(combined))
        
        if(nrow(scores) > 0) {
          # Group on tumour names, finds all tumour IDs for this tumour and adds the best scores
          scores %>%
            group_by(tumour_name) %>%
            summarise(
              all_tumour_ids = list(unique(tumour_id)),
              cosine = max(cosine),
              jaccard = max(jaccard),
              combined = max(combined),
              .groups = 'drop'
            ) %>%
            arrange(desc(combined)) %>%
            
            # Take top 3 matches
            slice(1:min(3, n())) %>%
            mutate(rank = row_number()) %>%
            # Change tumour ID list to a string 
            mutate(
              fuzzy_tumours_id = sapply(all_tumour_ids, function(x) {
                paste(x, collapse = ", ")
              }),
              fuzzy_tumours_id_list = all_tumour_ids
            ) %>%
            select(rank, fuzzy_tumours_id, fuzzy_tumours_id_list,
                   fuzzy_name = tumour_name,
                   similarity = combined,
                   cosine_score = cosine,
                   jaccard_score = jaccard)
          
          # In case of no matches
        } else {
          tibble(
            rank = NA_integer_,
            fuzzy_tumours_id = NA_character_,
            fuzzy_tumours_id_list = list(NA_character_),
            fuzzy_name = NA_character_,
            similarity = NA_real_,
            cosine_score = NA_real_,
            jaccard_score = NA_real_
          )
        }
      })
    ) %>%
    # Change matches list to normal rows
    unnest(matches) %>%
    select(EMA.ID, FDA.ID, Indication, rank, fuzzy_name, 
           fuzzy_tumours_id, fuzzy_tumours_id_list,
           similarity, cosine_score, jaccard_score)
  
  # Save batch result to list
  match_list_orphan[[i]] <- batch_result
}

# Combine all rows from the list
matches_orphan <- bind_rows(match_list_orphan) %>%
  select(-fuzzy_tumours_id_list)

# Rename columns
Indications_names <- Indications_df %>%
  select(NCT.ID, CT.EDU.ID, Indication) %>%
  mutate(Indication = clean_text(Indication)) %>%
  mutate(row_id = row_number())

# Make one large vector with all text that should be compared and splits these sentences into words (tokens)
all_terms <- c(tumours_all_names$Synonym, Indications_names$Indication)
it <- itoken(all_terms, 
             tokenizer = word_tokenizer, 
             progressbar = FALSE)

# Create a vocabulary and a mapping for it
vocab <- create_vocabulary(it)
vectorizer <- vocab_vectorizer(vocab)

# Create a Document-Term Matrix that shows how often words appear in certain sentences
dtm <- create_dtm(it, vectorizer)

# Use Term Frequency-Inverse Document Frequency to determine the importance of words --> common words are less important
tfidf <- TfIdf$new()
dtm_tfidf <- fit_transform(dtm, tfidf)

# Split matrices in RARECARE tumours and clinical trial tumours
n_tumours <- nrow(tumours_all_names)
tumour_vectors <- dtm_tfidf[1:n_tumours, ]
indication_vectors <- dtm_tfidf[(n_tumours + 1):nrow(dtm_tfidf), ]

# Calculate cosine similarity between indication vector and tumour vector using a pairwise similarity matrix
cosine_sim <- as.matrix(sim2(indication_vectors, tumour_vectors, method = "cosine"))

# Use batches of size 100 to increase performance
batch_size <- 100
batches <- split(1:nrow(Indications_names), ceiling(seq_len(nrow(Indications_names)) / batch_size))
match_list_Indications <- list()

for(i in seq_along(batches)) {
  message("Batch ", i, " of ", length(batches))
  batch_idx <- batches[[i]]
  
  # Go over every clinical trial tumour name within the batch
  batch_result <- Indications_names[batch_idx, ] %>%
    rowwise() %>%
    mutate(
      # Create a list with match scores
      matches = list({
        # Cosine scores
        cos_scores <- cosine_sim[row_id, ]
        # Split the tumour name (clinical trials) in words
        ind_words <- str_split(Indication, "\\s+")[[1]]
        # Put scores per tumour row
        scores <- tibble(
          idx = 1:length(cos_scores),
          tumour_id = tumours_all_names$Tumour.ID,
          tumour_name = tumours_all_names$Synonym,
          cosine = cos_scores
        ) %>%
          rowwise() %>%
          mutate(
            # Split the tumour name (RARECARE) in words
            tumour_words = list(str_split(tumour_name, "\\s+")[[1]]),
            # Calculate Jaccard similarity
            jaccard = length(intersect(ind_words, tumour_words)) / 
              length(union(ind_words, tumour_words)),
            # Combine scores and give it a weight
            combined = 0.4 * cosine + 0.6 * jaccard
          ) %>%
          ungroup() %>%
          # Keep only relatively good matches
          filter(combined > 0.2) %>%
          arrange(desc(combined))
        
        if(nrow(scores) > 0) {
          # Group on tumour names, finds all tumour IDs for this tumour and adds the best scores
          scores %>%
            group_by(tumour_name) %>%
            summarise(
              all_tumour_ids = list(unique(tumour_id)),
              cosine = max(cosine),
              jaccard = max(jaccard),
              combined = max(combined),
              .groups = 'drop'
            ) %>%
            arrange(desc(combined)) %>%
            
            # Take top 3 matches
            slice(1:min(3, n())) %>%
            mutate(rank = row_number()) %>%
            # Change tumour ID list to a string
            mutate(
              fuzzy_tumours_id = sapply(all_tumour_ids, function(x) {
                paste(x, collapse = ", ")
              }),
              fuzzy_tumours_id_list = all_tumour_ids
            ) %>%
            select(rank, fuzzy_tumours_id, fuzzy_tumours_id_list,
                   fuzzy_name = tumour_name,
                   similarity = combined,
                   cosine_score = cosine,
                   jaccard_score = jaccard)
          
          # In case of no matches
        } else {
          tibble(
            rank = NA_integer_,
            fuzzy_tumours_id = NA_character_,
            fuzzy_tumours_id_list = list(NA_character_),
            fuzzy_name = NA_character_,
            similarity = NA_real_,
            cosine_score = NA_real_,
            jaccard_score = NA_real_
          )
        }
      })
    ) %>%
    # Change matches list to normal rows
    unnest(matches) %>%
    select(NCT.ID, CT.EDU.ID, Indication, rank, fuzzy_name, 
           fuzzy_tumours_id, fuzzy_tumours_id_list,
           similarity, cosine_score, jaccard_score)
  
  # Save batch results to list
  match_list_Indications[[i]] <- batch_result
}

# Combine all rows from the list
matches_Indications <- bind_rows(match_list_Indications) %>%
  select(-fuzzy_tumours_id_list)

# Save files
write.csv(tumours_final, "tumours_final.csv", row.names = FALSE)
write.csv(matches_orphan, "matches_orphan_indications.csv", row.names = FALSE)
write.csv(matches_Indications, "matches_indications.csv", row.names = FALSE)
