library(dplyr)
library(stringr)
library(tidyr)
library(lubridate)
library(purrr)

# Read files
matches_orphan_indications <- read.csv("E:/data/Splitted2/matches_orphan_indications.csv")
matches_orphan_products <- read.csv("E:/data/Splitted2/matches_orphan_products.csv")
matches_products <- read.csv("E:/data/Splitted2/matches_products.csv")
matches_indications <- read.csv("E:/data/Splitted2/matches_indications.csv")
tumours_final <- read.csv("E:/data/Splitted2/tumours_final.csv")
products_df <- read.csv("E:/data/Splitted2/Products.csv")
indications_df <- read.csv("E:/data/Splitted2/Indications.csv")
sponsor_df <- read.csv("E:/data/Splitted2/Sponsors.csv")
trial_information_df <- read.csv("E:/data/Splitted2/Trial_information.csv")
orphan_designation_df <- read.csv("E:/data/Splitted2/orphan_designation.csv")

# Change type of Tumour ID for matching
tumours_final <- tumours_final %>%
  mutate(Tumour.ID = as.character(Tumour.ID))

# Only take the best ranked match
best_matches <- matches_orphan_indications %>%
  filter(rank == 1)

# Change scores of badly scored matches
good_matches_orphan_indications <- best_matches %>%
  mutate(
    similarity = if_else(str_to_lower(Indication) == "carcinoma thymic thymoma", 0, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "carcinoma thymic thymoma", 0, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "carcinoma thymic thymoma", 0, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "adenocarcinoma advanced gastroesophageal her including junction overexpressing stomach", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "adenocarcinoma advanced gastroesophageal her including junction overexpressing stomach", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "adenocarcinoma advanced gastroesophageal her including junction overexpressing stomach", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "adequate alone are cancer differentiated doses intolerant levo thyroid thyroxine well who", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "adequate alone are cancer differentiated doses intolerant levo thyroid thyroxine well who", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "adequate alone are cancer differentiated doses intolerant levo thyroid thyroxine well who", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "arid cancer containing deficient domain fallopian interactive ovarian peritoneal primary protein rich tube", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "arid cancer containing deficient domain fallopian interactive ovarian peritoneal primary protein rich tube", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "arid cancer containing deficient domain fallopian interactive ovarian peritoneal primary protein rich tube", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "acute antibiotic consolidation duration fever following hospitalization induction leukemia myeloid neutropenia use
e", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "acute antibiotic consolidation duration fever following hospitalization induction leukemia myeloid neutropenia use
", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "acute antibiotic consolidation duration fever following hospitalization induction leukemia myeloid neutropenia use
", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "myelodysplastic requiring syndromes therapy", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "myelodysplastic requiring syndromes therapy", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "myelodysplastic requiring syndromes therapy", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "bcl cell diffuse dlbcl grade hgbcl high large lymphoma myc rearrangements", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "bcl cell diffuse dlbcl grade hgbcl high large lymphoma myc rearrangements", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "bcl cell diffuse dlbcl grade hgbcl high large lymphoma myc rearrangements", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "detection localization neuroblastoma scintigraphic staging", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "detection localization neuroblastoma scintigraphic staging", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "detection localization neuroblastoma scintigraphic staging", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "anemia are associated blood cell dependent myelofibrosis myeloproliferative neoplasm persons red tranfusion who", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "anemia are associated blood cell dependent myelofibrosis myeloproliferative neoplasm persons red tranfusion who", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "anemia are associated blood cell dependent myelofibrosis myeloproliferative neoplasm persons red tranfusion who", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "carcinoma cell lung small", 0, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "carcinoma cell lung small", 0, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "carcinoma cell lung small", 0, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "neuroblastoma pediatric", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "neuroblastoma pediatric", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "neuroblastoma pediatric", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "hodgkin lymphoma non", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "hodgkin lymphoma non", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "hodgkin lymphoma non", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "cell enktl extra killer lymphoma natural nodal", 0, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "cell enktl extra killer lymphoma natural nodal", 0, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "cell enktl extra killer lymphoma natural nodal", 0, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "diagnostic management multiple myeloma", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "diagnostic management multiple myeloma", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "diagnostic management multiple myeloma", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "ewing family sarcoma tumor", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "ewing family sarcoma tumor", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "ewing family sarcoma tumor", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "adenocarcinoma cancer gastric gastroesophageal including junction", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "adenocarcinoma cancer gastric gastroesophageal including junction", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "adenocarcinoma cancer gastric gastroesophageal including junction", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "advanced anaplastic carcinoma follicular locally medullary metastatic papillary thyroid", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "advanced anaplastic carcinoma follicular locally medullary metastatic papillary thyroid", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "advanced anaplastic carcinoma follicular locally medullary metastatic papillary thyroid", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "cell hodgkin lymphoma non", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "cell hodgkin lymphoma non", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "cell hodgkin lymphoma non", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "arising cell chronic leukemia lymphocytic prolymphocytic teatment", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "arising cell chronic leukemia lymphocytic prolymphocytic teatment", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "arising cell chronic leukemia lymphocytic prolymphocytic teatment", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "sarcoma soft tissues", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "sarcoma soft tissues", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "sarcoma soft tissues", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "bcl cell diffuse grade high large lymphoma myc rearrangements", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "bcl cell diffuse grade high large lymphoma myc rearrangements", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "bcl cell diffuse grade high large lymphoma myc rearrangements", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "chromosome chronic leukemia myeloid philadelphia positive", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "chromosome chronic leukemia myeloid philadelphia positive", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "chromosome chronic leukemia myeloid philadelphia positive", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "cell lymphoma peripheral ptcl", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "cell lymphoma peripheral ptcl", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "cell lymphoma peripheral ptcl", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "adult cell leukemia lymphoma", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "adult cell leukemia lymphoma", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "adult cell leukemia lymphoma", 1, jaccard_score),
    cosine_score = if_else(str_to_lower(Indication) == "acute leukemia myelocytic subtypes", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "acute leukemia myelocytic subtypes", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "acute leukemia myelocytic subtypes", 1, similarity),
    similarity = if_else(str_to_lower(Indication) == "acute apl leukemia promyelocytic", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "acute apl leukemia promyelocytic", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "acute apl leukemia promyelocytic", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "cell lymphoma mantle mcl", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "cell lymphoma mantle mcl", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "cell lymphoma mantle mcl", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == " gastro intestinal malignant stromal tumor", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == " gastro intestinal malignant stromal tumor", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == " gastro intestinal malignant stromal tumor", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "ewings sarcoma", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "ewings sarcoma", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "ewings sarcoma", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "cell lymphoma not otherwise peripheral specified", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "cell lymphoma not otherwise peripheral specified", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "cell lymphoma not otherwise peripheral specified", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "bone cancer osteosarcoma", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "bone cancer osteosarcoma", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "bone cancer osteosarcoma", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "hodgkin lymphoma non refractory relapsed", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "hodgkin lymphoma non refractory relapsed", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "hodgkin lymphoma non refractory relapsed", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "diagnostic management neuroblastoma", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "diagnostic management neuroblastoma", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "diagnostic management neuroblastoma", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "inappropriate multiple myeloma oral therapy whom", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "inappropriate multiple myeloma oral therapy whom", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "inappropriate multiple myeloma oral therapy whom", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "mutliple myeloma", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "mutliple myeloma", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "mutliple myeloma", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "astrocytic tumor", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "astrocytic tumor", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "astrocytic tumor", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "astrocytic glioma", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "astrocytic glioma", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "astrocytic glioma", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "clinical diagnostic ewing management sarcoma", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "clinical diagnostic ewing management sarcoma", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == " clinical diagnostic ewing management sarcoma", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "advanced associated hiv kaposi sarcoma", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "advanced associated hiv kaposi sarcoma", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "advanced associated hiv kaposi sarcoma", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "cell follicular including lymphoma peripheral", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "cell follicular including lymphoma peripheral", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "cell follicular including lymphoma peripheral", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "disease hodgkin lymphoma non", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "disease hodgkin lymphoma non", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "disease hodgkin lymphoma non", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "adult atll cell leukemia lymphoma", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "adult atll cell leukemia lymphoma", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "adult atll cell leukemia lymphoma", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "associated fungoides mycosis pruritus", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "associated fungoides mycosis pruritus", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "associated fungoides mycosis pruritus", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "hodgkin lymphoma malignant non", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "hodgkin lymphoma malignant non", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "hodgkin lymphoma malignant non", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "acute all leukemia lymphoblastic", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "acute all leukemia lymphoblastic", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "acute all leukemia lymphoblastic", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "acute childhood leukemia lymphocytic refractory", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "acute childhood leukemia lymphocytic refractory", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "acute childhood leukemia lymphocytic refractory", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "acute induction leukemia myeloid therapy undergoing", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "acute induction leukemia myeloid therapy undergoing", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "acute induction leukemia myeloid therapy undergoing", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "acute leukemia lymphoblastic philadelphia positive", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "acute leukemia lymphoblastic philadelphia positive", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "acute leukemia lymphoblastic philadelphia positive", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "acute leukemia lymphocytic positive tdt", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "acute leukemia lymphocytic positive tdt", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "acute leukemia lymphocytic positive tdt", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "acute leukemia myelocytic", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "acute leukemia myelocytic", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "acute leukemia myelocytic", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "cell lymphoma peripheral ptcl", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "cell lymphoma peripheral ptcl", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "cell lymphoma peripheral ptcl", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "accessible administration anatomical carcinoma cavity cell local oral oropharynx sites squamous topical", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "accessible administration anatomical carcinoma cavity cell local oral oropharynx sites squamous topical", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "accessible administration anatomical carcinoma cavity cell local oral oropharynx sites squamous topical", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "gastro intestinal malignant stromal tumor", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "gastro intestinal malignant stromal tumor", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "gastro intestinal malignant stromal tumor", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "adenocarcinoma advanced expressing gastroesophageal her junction stomach", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "adenocarcinoma advanced expressing gastroesophageal her junction stomach", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "adenocarcinoma advanced expressing gastroesophageal her junction stomach", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "clinical diagnostic management neuroblastoma", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "clinical diagnostic management neuroblastoma", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "clinical diagnostic management neuroblastoma", 1, jaccard_score)) %>%
  # Filter on scores
  filter(similarity > 0.5385,
         cosine_score > 0.4678,
         jaccard_score > 0.42)

# Change scores of too high-ranked matches
bad_matches_orphan_indications <- best_matches %>%
  mutate(
    similarity = if_else(str_to_lower(Indication) == "carcinoma thymic thymoma", 0, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "carcinoma thymic thymoma", 0, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "carcinoma thymic thymoma", 0, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "adenocarcinoma advanced gastroesophageal her including junction overexpressing stomach", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "adenocarcinoma advanced gastroesophageal her including junction overexpressing stomach", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "adenocarcinoma advanced gastroesophageal her including junction overexpressing stomach", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "adequate alone are cancer differentiated doses intolerant levo thyroid thyroxine well who", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "adequate alone are cancer differentiated doses intolerant levo thyroid thyroxine well who", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "adequate alone are cancer differentiated doses intolerant levo thyroid thyroxine well who", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "arid cancer containing deficient domain fallopian interactive ovarian peritoneal primary protein rich tube", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "arid cancer containing deficient domain fallopian interactive ovarian peritoneal primary protein rich tube", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "arid cancer containing deficient domain fallopian interactive ovarian peritoneal primary protein rich tube", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "acute antibiotic consolidation duration fever following hospitalization induction leukemia myeloid neutropenia use
e", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "acute antibiotic consolidation duration fever following hospitalization induction leukemia myeloid neutropenia use
", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "acute antibiotic consolidation duration fever following hospitalization induction leukemia myeloid neutropenia use
", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "myelodysplastic requiring syndromes therapy", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "myelodysplastic requiring syndromes therapy", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "myelodysplastic requiring syndromes therapy", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "bcl cell diffuse dlbcl grade hgbcl high large lymphoma myc rearrangements", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "bcl cell diffuse dlbcl grade hgbcl high large lymphoma myc rearrangements", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "bcl cell diffuse dlbcl grade hgbcl high large lymphoma myc rearrangements", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "detection localization neuroblastoma scintigraphic staging", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "detection localization neuroblastoma scintigraphic staging", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "detection localization neuroblastoma scintigraphic staging", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "anemia are associated blood cell dependent myelofibrosis myeloproliferative neoplasm persons red tranfusion who", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "anemia are associated blood cell dependent myelofibrosis myeloproliferative neoplasm persons red tranfusion who", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "anemia are associated blood cell dependent myelofibrosis myeloproliferative neoplasm persons red tranfusion who", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "carcinoma cell lung small", 0, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "carcinoma cell lung small", 0, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "carcinoma cell lung small", 0, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "neuroblastoma pediatric", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "neuroblastoma pediatric", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "neuroblastoma pediatric", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "hodgkin lymphoma non", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "hodgkin lymphoma non", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "hodgkin lymphoma non", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "cell enktl extra killer lymphoma natural nodal", 0, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "cell enktl extra killer lymphoma natural nodal", 0, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "cell enktl extra killer lymphoma natural nodal", 0, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "diagnostic management multiple myeloma", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "diagnostic management multiple myeloma", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "diagnostic management multiple myeloma", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "ewing family sarcoma tumor", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "ewing family sarcoma tumor", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "ewing family sarcoma tumor", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "adenocarcinoma cancer gastric gastroesophageal including junction", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "adenocarcinoma cancer gastric gastroesophageal including junction", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "adenocarcinoma cancer gastric gastroesophageal including junction", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "advanced anaplastic carcinoma follicular locally medullary metastatic papillary thyroid", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "advanced anaplastic carcinoma follicular locally medullary metastatic papillary thyroid", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "advanced anaplastic carcinoma follicular locally medullary metastatic papillary thyroid", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "cell hodgkin lymphoma non", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "cell hodgkin lymphoma non", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "cell hodgkin lymphoma non", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "arising cell chronic leukemia lymphocytic prolymphocytic teatment", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "arising cell chronic leukemia lymphocytic prolymphocytic teatment", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "arising cell chronic leukemia lymphocytic prolymphocytic teatment", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "sarcoma soft tissues", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "sarcoma soft tissues", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "sarcoma soft tissues", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "bcl cell diffuse grade high large lymphoma myc rearrangements", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "bcl cell diffuse grade high large lymphoma myc rearrangements", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "bcl cell diffuse grade high large lymphoma myc rearrangements", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "chromosome chronic leukemia myeloid philadelphia positive", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "chromosome chronic leukemia myeloid philadelphia positive", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "chromosome chronic leukemia myeloid philadelphia positive", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "cell lymphoma peripheral ptcl", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "cell lymphoma peripheral ptcl", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "cell lymphoma peripheral ptcl", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "adult cell leukemia lymphoma", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "adult cell leukemia lymphoma", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "adult cell leukemia lymphoma", 1, jaccard_score),
    cosine_score = if_else(str_to_lower(Indication) == "acute leukemia myelocytic subtypes", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "acute leukemia myelocytic subtypes", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "acute leukemia myelocytic subtypes", 1, similarity),
    similarity = if_else(str_to_lower(Indication) == "acute apl leukemia promyelocytic", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "acute apl leukemia promyelocytic", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "acute apl leukemia promyelocytic", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "cell lymphoma mantle mcl", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "cell lymphoma mantle mcl", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "cell lymphoma mantle mcl", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == " gastro intestinal malignant stromal tumor", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == " gastro intestinal malignant stromal tumor", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == " gastro intestinal malignant stromal tumor", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "ewings sarcoma", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "ewings sarcoma", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "ewings sarcoma", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "cell lymphoma not otherwise peripheral specified", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "cell lymphoma not otherwise peripheral specified", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "cell lymphoma not otherwise peripheral specified", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "bone cancer osteosarcoma", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "bone cancer osteosarcoma", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "bone cancer osteosarcoma", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "hodgkin lymphoma non refractory relapsed", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "hodgkin lymphoma non refractory relapsed", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "hodgkin lymphoma non refractory relapsed", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "diagnostic management neuroblastoma", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "diagnostic management neuroblastoma", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "diagnostic management neuroblastoma", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "inappropriate multiple myeloma oral therapy whom", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "inappropriate multiple myeloma oral therapy whom", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "inappropriate multiple myeloma oral therapy whom", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "mutliple myeloma", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "mutliple myeloma", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "mutliple myeloma", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "astrocytic tumor", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "astrocytic tumor", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "astrocytic tumor", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "astrocytic glioma", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "astrocytic glioma", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "astrocytic glioma", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "clinical diagnostic ewing management sarcoma", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "clinical diagnostic ewing management sarcoma", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == " clinical diagnostic ewing management sarcoma", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "advanced associated hiv kaposi sarcoma", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "advanced associated hiv kaposi sarcoma", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "advanced associated hiv kaposi sarcoma", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "cell follicular including lymphoma peripheral", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "cell follicular including lymphoma peripheral", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "cell follicular including lymphoma peripheral", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "disease hodgkin lymphoma non", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "disease hodgkin lymphoma non", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "disease hodgkin lymphoma non", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "adult atll cell leukemia lymphoma", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "adult atll cell leukemia lymphoma", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "adult atll cell leukemia lymphoma", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "associated fungoides mycosis pruritus", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "associated fungoides mycosis pruritus", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "associated fungoides mycosis pruritus", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "hodgkin lymphoma malignant non", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "hodgkin lymphoma malignant non", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "hodgkin lymphoma malignant non", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "acute all leukemia lymphoblastic", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "acute all leukemia lymphoblastic", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "acute all leukemia lymphoblastic", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "acute childhood leukemia lymphocytic refractory", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "acute childhood leukemia lymphocytic refractory", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "acute childhood leukemia lymphocytic refractory", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "acute induction leukemia myeloid therapy undergoing", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "acute induction leukemia myeloid therapy undergoing", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "acute induction leukemia myeloid therapy undergoing", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "acute leukemia lymphoblastic philadelphia positive", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "acute leukemia lymphoblastic philadelphia positive", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "acute leukemia lymphoblastic philadelphia positive", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "acute leukemia lymphocytic positive tdt", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "acute leukemia lymphocytic positive tdt", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "acute leukemia lymphocytic positive tdt", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "acute leukemia myelocytic", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "acute leukemia myelocytic", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "acute leukemia myelocytic", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "cell lymphoma peripheral ptcl", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "cell lymphoma peripheral ptcl", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "cell lymphoma peripheral ptcl", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "accessible administration anatomical carcinoma cavity cell local oral oropharynx sites squamous topical", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "accessible administration anatomical carcinoma cavity cell local oral oropharynx sites squamous topical", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "accessible administration anatomical carcinoma cavity cell local oral oropharynx sites squamous topical", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "gastro intestinal malignant stromal tumor", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "gastro intestinal malignant stromal tumor", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "gastro intestinal malignant stromal tumor", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "adenocarcinoma advanced expressing gastroesophageal her junction stomach", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "adenocarcinoma advanced expressing gastroesophageal her junction stomach", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "adenocarcinoma advanced expressing gastroesophageal her junction stomach", 1, jaccard_score),
    similarity = if_else(str_to_lower(Indication) == "clinical diagnostic management neuroblastoma", 1, similarity),
    cosine_score = if_else(str_to_lower(Indication) == "clinical diagnostic management neuroblastoma", 1, cosine_score),
    jaccard_score = if_else(str_to_lower(Indication) == "clinical diagnostic management neuroblastoma", 1, jaccard_score)) %>%
  # Filter on scores
  filter(similarity <= 0.5385,
         cosine_score <= 0.4678,
         jaccard_score <= 0.42)

# Join with the RARECARE-Orphanet mapping table
good_matches_orphan_indications <- good_matches_orphan_indications %>%
  select(EMA.ID, FDA.ID, Indication, fuzzy_tumours_id) %>%
  separate_rows(fuzzy_tumours_id, sep = ",\\s*") %>%
  rename(tumours_id = fuzzy_tumours_id) %>%
  mutate(tumours_id = as.character(tumours_id)) %>%
  left_join(
    tumours_final %>% select(Tumour.ID, Synonym_main, Synonym_orpha),
    by = c("tumours_id" = "Tumour.ID"),
    relationship = "many-to-many"
  ) %>%
  # Prefer Orpha synonym over RARECARE name
  mutate(Indication = coalesce(Synonym_orpha, Synonym_main)) %>%
  select(-Synonym_main, -Synonym_orpha)

# Join the tumour matches with the old orphan designation data
orphan_designation_1_df <- good_matches_orphan_indications %>%
  left_join(
    orphan_designation_df %>% select(-Indication), by = c("EMA.ID" = "EMA.ID", "FDA.ID" = "FDA.ID"), relationship = "many-to-many") %>%
  filter(!is.na(tumours_id)) %>%
  group_by(EMA.ID, FDA.ID) %>%
  summarise(
    tumours_id = paste(unique(tumours_id), collapse = ", "),
    across(everything(), first),
    .groups = "drop"
  )

# Join with product matches data
orphan_designation_2_df <- matches_orphan_products %>%
  right_join(
    orphan_designation_1_df, by = c("EMA.ID", "FDA.ID"), relationship = "many-to-many") %>%
  mutate(
    # Take generic name from the product matches if available
    Generic.Name = coalesce(Generic.Name.x, Generic.Name.y),
    # Fix authorised/implemented drug information manually
    Generic.Name = case_when(
      Generic.Name == "fludarabine phosphate oral tablets" ~ "fludarabine",
      Generic.Name == "tisagenlecleucel autologous cells transduced lentiviral vector containing chimeric antigen receptor directed against" ~ "tisagenlecleucel",
      TRUE ~ Generic.Name
    ),
    Drugbank.ID = case_when(
      Generic.Name == "fludarabine" ~ "DB01073",
      Generic.Name == "tisagenlecleucel" ~ "DB13881",
      TRUE ~ Drugbank.ID
    ),
    # Change all NAs in designation status to designated
    Orphan.Designation.Status = ifelse(
      is.na(Orphan.Designation.Status),
      "designated",
      Orphan.Designation.Status),
    # Change all mentions of implemented to approved
    Orphan.Designation.Status = ifelse(
      Orphan.Designation.Status == "implemented",
      "approved",
      Orphan.Designation.Status),
    # Combine marketing approval date and implemented on columns
    Marketing.Authorisation.Date = coalesce(Marketing.Approval.Date, Implemented.on)
  ) %>%
  # Remove data fields which are not in the RARECARE list
  filter(!is.na(tumours_id)) %>% 
  rename(Tumours.ID = tumours_id) %>%
  distinct() %>%
  select(EMA.ID, FDA.ID, Drugbank.ID, Generic.Name, Classification, Brandname, EU.Product.ID, Tumours.ID, Indication, Groups, Orphan.Designation.Status, Designation.Date, Date.Designation.Withdrawn.or.Revoked, Marketing.Authorisation.Date, Sponsor)

# Only take the best ranked match
best_matches <- matches_indications %>%
  filter(rank == 1)

# Filter on scores
good_matches_indications <- best_matches %>%
  filter(jaccard_score > 0.74 | similarity > 0.74 | (cosine_score > 0.95 & similarity > 0.676))

# Join with the RARECARE-Orphanet mapping table
good_matches_indications <- good_matches_indications %>%
  select(NCT.ID, CT.EDU.ID, Indication, fuzzy_tumours_id) %>%
  separate_rows(fuzzy_tumours_id, sep = ",\\s*") %>%
  rename(tumours_id = fuzzy_tumours_id) %>%
  mutate(tumours_id = as.character(tumours_id)) %>%
  left_join(
    tumours_final %>% select(Tumour.ID, Synonym_main, Synonym_orpha),
    by = c("tumours_id" = "Tumour.ID"),
    relationship = "many-to-many"
  ) %>%
  # Prefer Orpha synonym over RARECARE name
  mutate(Indication = coalesce(Synonym_orpha, Synonym_main)) %>%
  select(-Synonym_main, -Synonym_orpha)

# Join the tumour matches with the old tumour data
indications_1_df <- indications_df %>%
  select(-Indication) %>%
  left_join(
    good_matches_indications,
    by = c("NCT.ID", "CT.EDU.ID"),
    relationship = "many-to-many"
  ) %>%
  filter(!is.na(tumours_id)) %>%
  group_by(NCT.ID, CT.EDU.ID) %>%
  summarise(
    tumours_id = paste(unique(tumours_id), collapse = ", "),
    across(everything(), first),
    .groups = "drop"
  )

# Join clinical trial product data with product matches data
products_2_df <- products_df %>%
  left_join(
    matches_products, by = c("NCT.ID" = "NCT.ID", "CT.EDU.ID" = "CT.EDU.ID", "Product.Type" = "Product.Type"), relationship = "many-to-many") %>%
  # Prefer matched names over already existing ones
  mutate(
    Generic.Name = coalesce(Generic.Name, Product.Name)
  ) %>%
  select(-Product.Name, -EMA.ID) %>%
  # Remove NCT ID error
  mutate(NCT.ID = na_if(NCT.ID, "NCT00000000")) %>%
  # Remove all mentions of "cohort"
  mutate(Generic.Name = str_remove_all(Generic.Name, regex("cohort", ignore_case = TRUE))) %>%
  mutate(Generic.Name = str_trim(Generic.Name)) %>%
  filter(Generic.Name != "" & !is.na(Generic.Name))
  filter(!str_detect(Generic.Name, "placebo|not yet assigned|not available|cohort not previously treated"),
         !Generic.Name == "iron") %>%
  distinct() 

# Save files
write.csv(orphan_designation_2_df, "orphan_designation.csv", row.names = FALSE)
write.csv(indications_1_df, "indications.csv", row.names = FALSE)
write.csv(products_2_df, "products.csv", row.names = FALSE)

