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
tumour_name <- tumours_final %>%
  mutate(Tumour.ID = as.character(Tumour.ID)) %>%
  select(Tumour.ID, Synonym_main, Crude.incidence.rate.per.100.000) %>%
  rename(Tumours.ID = Tumour.ID) %>%
  distinct()

# Original RARECARE names
tumours_raw_filtered <- tumours_raw %>%
  filter(`R=rare` == "R") %>%
  filter(Tier == 2) %>%
  mutate(Tumour = str_remove(Tumour, "\\*")) %>%
  select(Tumour, `Crude incidence rate per 100,000`) %>%
  mutate(tumor = c(
    "Squamous cell carcinoma with variants of nasal cavity and sinuses",
    "Lymphoepithelial carcinoma of nasal cavity and sinuses",
    "Undifferentiated carcinoma of nasal cavity and sinuses",
    "Intestinal type adenocarcinoma of nasal cavity and sinuses",
    "Squamous cell carcinoma with variants of nasopharynx",
    "Papillary adenocarcinoma of nasopharynx",
    "Epithelial tumours of major salivary glands",
    "Salivary gland type tumours of head and neck",
    "Squamous cell carcinoma with variants of hypopharynx",
    "Squamous cell carcinoma with variants of larynx",
    "Squamous cell carcinoma with variants of oropharynx",
    "Squamous cell carcinoma with variants of oral cavity",
    "Squamous cell carcinoma with variants of lip",
    "Squamous cell carcinoma with variants of oesophagus",
    "Adenocarcinoma with variants of oesophagus",
    "Salivary gland type tumours of oesophagus",
    "Undifferentiated carcinoma of oesophagus",
    "Squamous cell carcinoma with variants of stomach",
    "Salivary gland type tumours of stomach",
    "Undifferentiated carcinoma of stomach",
    "Adenocarcinoma with variants of small intestine",
    "Squamous cell carcinoma with variants of small intestine",
    "Squamous cell carcinoma with variants of colon",
    "Fibromyxoma and low-grade mucinous adenocarcinoma (pseudomyxoma peritonei) of the appendix",
    "Squamous cell carcinoma with variants of rectum",
    "Squamous cell carcinoma with variants of anal canal",
    "Adenocarcinoma with variants of anal canal",
    "Paget's disease of anal canal",
    "Squamous cell carcinoma with variants of pancreas",
    "Acinar cell carcinoma of pancreas",
    "Mucinous cystadenocarcinoma of pancreas",
    "Intraductal papillary mucinous carcinoma (invasive) of pancreas",
    "Solid pseudopapillary carcinoma of pancreas",
    "Serous cystadenocarcinoma of pancreas",
    "Carcinoma with osteoclast-like giant cells of pancreas",
    "Hepatocellular carcinoma of liver and intrahepatic bile tract",
    "Fibrolamellar hepatocellular carcinoma",
    "Cholangiocarcinoma of intrahepatic bile tract",
    "Adenocarcinoma with variants of liver and intrahepatic bile tract",
    "Undifferentiated carcinoma of liver and intrahepatic bile tract",
    "Squamous cell carcinoma with variants of liver and intrahepatic bile tract",
    "Bile duct cystadenocarcinoma of intrahepatic bile tract",
    "Adenocarcinoma with variants of gallbladder",
    "Adenocarcinoma with variants of extrahepatic bile tract",
    "Squamous cell carcinoma of gallbladder and extrahepatic bile tract",
    "Squamous cell carcinoma with variants of trachea",
    "Adenocarcinoma with variants of trachea",
    "Salivary gland type tumours of trachea",
    "Adenosquamous carcinoma of lung",
    "Large cell carcinoma of lung",
    "Salivary gland type tumours of lung",
    "Sarcomatoid carcinoma of lung",
    "Malignant thymoma",
    "Squamous cell carcinoma of thymus",
    "Undifferentiated carcinoma of thymus",
    "Lymphoepithelial carcinoma of thymus",
    "Adenocarcinoma with variants of thymus",
    "Mammary Paget’s disease of breast",
    "Special types of adenocarcinoma of breast",
    "Metaplastic carcinoma of breast",
    "Salivary gland type tumours of breast",
    "Epithelial tumour of male breast",
    "Squamous cell carcinoma with variants of corpus uteri",
    "Adenoid cystic carcinoma of corpus uteri",
    "Clear cell adenocarcinoma, NOS",
    "Serous (papillary) carcinoma",
    "Müllerian mixed tumour",
    "Squamous cell carcinoma with variants of cervix uteri",
    "Adenocarcinoma with variants of cervix uteri",
    "Undifferentiated carcinoma of cervix uteri",
    "Müllerian mixed tumour of cervix uteri",
    "Adenocarcinoma with variants of ovary",
    "Mucinous adenocarcinoma of ovary",
    "Clear cell adenocarcinoma of ovary",
    "Primary peritoneal serous/papillary carcinoma",
    "Müllerian mixed tumour of ovary",
    "Adenocarcinoma with variants of fallopian tube",
    "Sex cord tumours of ovary",
    "Malignant/immature teratomas of ovary",
    "Germ cell tumour of ovary",
    "Squamous cell carcinoma with variants of vulva and vagina",
    "Adenocarcinoma with variants of vulva and vagina",
    "Paget's disease of vulva and vagina",
    "Undifferentiated carcinoma of vulva and vagina",
    "Choriocarcinoma of placenta",
    "Squamous cell carcinoma with variants of prostate",
    "Infiltrating duct carcinoma of prostate",
    "Transitional cell carcinoma of prostate",
    "Basal cell adenocarcinoma of prostate",
    "Paratesticular adenocarcinoma with variants",
    "Non-seminomatous testicular cancer",
    "Seminomatous testicular cancer",
    "Spermatocytic seminoma",
    "Teratoma with malignant transformation",
    "Testicular sex cord tumour",
    "Squamous cell carcinoma with variants of penis",
    "Adenocarcinoma with variants of penis",
    "Squamous cell carcinoma, spindle cell type, of kidney",
    "Squamous cell carcinoma with variants of kidney",
    "Transitional cell carcinoma of renal pelvis and ureter",
    "Squamous cell carcinoma with variants of renal pelvis and ureter",
    "Adenocarcinoma with variants of renal pelvis and ureter",
    "Transitional cell carcinoma of urethra",
    "Squamous cell carcinoma with variants of urethra",
    "Adenocarcinoma with variants of urethra",
    "Squamous cell carcinoma with variants of bladder",
    "Adenocarcinoma with variants of bladder",
    "Salivary gland type tumours of bladder",
    "Squamous cell carcinoma with variants of eye and adnexa",
    "Adenocarcinoma with variants of eye and adnexa",
    "Squamous cell carcinoma with variants of middle ear",
    "Adenocarcinoma with variants of middle ear",
    "Mesothelioma of pleura and pericardium",
    "Mesothelioma of peritoneum and tunica vaginalis",
    "Malignant melanoma of mucosa and extracutaneous sites",
    "Malignant melanoma of uvea",
    "Adnexal carcinoma of skin",
    "Neuroblastoma and ganglioneuroblastoma",
    "Nephroblastoma",
    "Retinoblastoma",
    "Hepatoblastoma",
    "Pleuropulmonary blastoma",
    "Pancreatoblastoma",
    "Olfactory neuroblastoma",
    "Odontogenic malignant tumours",
    "Non-seminomatous germ cell tumours",
    "Seminomatous germ cell tumours",
    "Germ cell tumours of central nervous system (CNS)",
    "Soft tissue sarcoma of head and neck",
    "Soft tissue sarcoma of limbs",
    "Soft tissue sarcoma of superficial trunk",
    "Soft tissue sarcoma of mediastinum",
    "Soft tissue sarcoma of heart",
    "Soft tissue sarcoma of breast",
    "Soft tissue sarcoma of uterus",
    "Soft tissue sarcoma of paratestis",
    "Soft tissue sarcomas of other genitourinary sites",
    "Soft tissue sarcoma of viscera",
    "Soft tissue sarcoma of retroperitoneum and peritoneum",
    "Soft tissue sarcoma of pelvis",
    "Soft tissue sarcoma of skin",
    "Soft tissue sarcoma of paraorbital region",
    "Soft tissue sarcoma of brain and nervous system",
    "Embryonal rhabdomyosarcoma",
    "Alveolar rhabdomyosarcoma",
    "Ewing sarcoma of soft tissue",
    "Osteosarcoma",
    "Chondrosarcoma",
    "Chordoma",
    "Vascular sarcoma",
    "Ewing sarcoma",
    "Adamantinoma",
    "High-grade sarcoma (fibrosarcoma, undifferentiated pleomorphic sarcoma)",
    "Gastrointestinal stromal tumour (GIST)",
    "Kaposi sarcoma",
    "Well-differentiated neuroendocrine tumour (NET) of pancreas and digestive tract",
    "Well-differentiated functioning neuroendocrine tumour of pancreas and digestive tract",
    "Poorly differentiated neuroendocrine carcinoma of pancreas and digestive tract",
    "Mixed neuroendocrine-non-neuroendocrine neoplasm (MiNEN) of pancreas and digestive tract",
    "Endocrine carcinoma of thyroid gland",
    "Neuroendocrine carcinoma of skin",
    "Typical and atypical carcinoid of the lung",
    "Neuroendocrine carcinoma of other sites",
    "Malignant pheochromocytoma",
    "Paraganglioma",
    "Pituitary carcinoma",
    "Carcinoma of thyroid gland",
    "Parathyroid carcinoma",
    "Adrenal cortical carcinoma",
    "Astrocytic tumours of CNS",
    "Oligodendroglial tumours of CNS",
    "Ependymal tumours of CNS",
    "Neuronal and mixed neuronal-glial tumours of CNS",
    "Choroid plexus carcinoma of CNS",
    "Malignant meningioma",
    "Embryonal tumours of CNS",
    "Classical Hodgkin lymphoma",
    "Nodular lymphocyte predominant Hodgkin lymphoma",
    "Precursor B/T lymphoblastic leukaemia/lymphoma (including Burkitt leukaemia/lymphoma)",
    "Cutaneous T-cell lymphoma (Sézary syndrome, mycosis fungoides)",
    "Other T-cell lymphomas and NK-cell neoplasms",
    "Diffuse large B-cell lymphoma",
    "Follicular lymphoma",
    "Hairy cell leukaemia",
    "Multiple myeloma / Plasmacytoma",
    "Mantle cell lymphoma",
    "B-cell prolymphocytic leukaemia",
    "Acute promyelocytic leukaemia (APL)",
    "Acute myeloid leukaemia (AML)",
    "Chronic myeloid leukaemia (CML)",
    "Other myeloproliferative neoplasms",
    "Mast cell tumour",
    "Myelodysplastic syndrome with isolated del(5q)",
    "Other myelodysplastic syndromes",
    "Chronic myelomonocytic leukaemia",
    "Atypical chronic myeloid leukaemia, BCR-ABL1 negative",
    "Histiocytic malignancies",
    "Lymph node accessory cell tumours"
  )
  ) %>%
  mutate(abbreviation = c(
    "SCC-Nasal", "LEC-Nasal", "UC-Nasal", "ITAC-Nasal", "SCC-NPC", "PAC-NPC",
    "Salivary-Major", "Salivary-HN", "SCC-Hypopharynx", "SCC-Larynx", "SCC-OPC",
    "SCC-Oral", "SCC-Lip", "SCC-Oesophagus", "ADC-Oesophagus", "Salivary-Oesophagus",
    "UC-Oesophagus", "SCC-Stomach", "Salivary-Stomach", "UC-Stomach", "ADC-SI",
    "SCC-SI", "SCC-Colon", "LAMN-Appendix", "SCC-Rectum", "SCC-Anal", "ADC-Anal",
    "Paget-Anal", "SCC-Pancreas", "Acinar-Pancreas", "MCAC-Pancreas", "IPMC-Pancreas",
    "SPPC-Pancreas", "SCAC-Pancreas", "OGC-Pancreas", "HCC", "FL-HCC", "CCA",
    "ADC-Liver", "UC-Liver", "SCC-Liver", "BDCAC", "ADC-GB", "ADC-EBT", "SCC-GB/EBT",
    "SCC-Trachea", "ADC-Trachea", "Salivary-Trachea", "ADSCC-Lung", "LCLC", 
    "Salivary-Lung", "Sarcomatoid-Lung", "Malignant-Thymoma", "SCC-Thymus", 
    "UC-Thymus", "LEC-Thymus", "ADC-Thymus", "MPD-Breast", "Special-ADC-Breast",
    "Metaplastic-Breast", "Salivary-Breast", "Male-Breast-Ca", "SCC-Endometrium",
    "ACC-Endometrium", "CCC", "Serous-Ca", "MMMT", "SCC-Cervix", "ADC-Cervix",
    "UC-Cervix", "MMMT-Cervix", "ADC-Ovary", "MAC-Ovary", "CCC-Ovary", 
    "PPSC", "MMMT-Ovary", "ADC-Fallopian", "SCST-Ovary", "Immature-Teratoma-Ovary",
    "GCT-Ovary", "SCC-Vulva/Vagina", "ADC-Vulva/Vagina", "Paget-Vulva", 
    "UC-Vulva/Vagina", "Choriocarcinoma", "SCC-Prostate", "IDC-Prostate", 
    "TCC-Prostate", "BCAC-Prostate", "Paratesticular-ADC", "NSGCT", "SGCT", 
    "Spermatocytic-Seminoma", "TMT", "Testicular-SCT", "SCC-Penis", "ADC-Penis",
    "Spindle-SCC-Kidney", "SCC-Kidney", "TCC-RenalPelvis", "SCC-RenalPelvis",
    "ADC-RenalPelvis", "TCC-Urethra", "SCC-Urethra", "ADC-Urethra", "SCC-Bladder",
    "ADC-Bladder", "Salivary-Bladder", "SCC-Eye", "ADC-Eye", "SCC-MiddleEar",
    "ADC-MiddleEar", "Mesothelioma-Pleura", "Mesothelioma-Peritoneum", 
    "Mucosal-Melanoma", "Uveal-Melanoma", "Adnexal-Ca-Skin", "Neuroblastoma",
    "Nephroblastoma", "Retinoblastoma", "Hepatoblastoma", "PPB", "Pancreatoblastoma",
    "ONB", "Odontogenic-Ca", "NSGCT", "SGCT", "CNS-GCT", "STS-HN", "STS-Limbs",
    "STS-Trunk", "STS-Mediastinum", "STS-Heart", "STS-Breast", "STS-Uterus",
    "STS-Paratestis", "STS-GU-Other", "STS-Viscera", "STS-Retroperitoneum",
    "STS-Pelvis", "STS-Skin", "STS-Paraorbital", "STS-CNS", "ERMS", "ARMS",
    "EWS-STS", "OS", "CS", "Chordoma", "Angiosarcoma", "EWS", "Adamantinoma",
    "HighGrade-Sarcoma", "GIST", "KS", "NET-Pancreas", "fNET-Pancreas", 
    "PDNEC", "MiNEN", "Thyroid-Endocrine-Ca", "NEC-Skin", "Carcinoid-Lung",
    "NEC-Other", "Malignant-Pheo", "Paraganglioma", "Pituitary-Ca", "Thyroid-Ca",
    "Parathyroid-Ca", "ACC-Adrenal", "Astrocytoma-CNS", "Oligodendroglioma-CNS",
    "Ependymoma-CNS", "Neuronal-Glial-CNS", "CPC", "Malignant-Meningioma",
    "Embryonal-CNS", "cHL", "NLPHL", "LBL-Burkitt", "CTCL", "T/NK-Cell-Lymphoma",
    "DLBCL", "FL", "HCL", "MM", "MCL", "B-PLL", "APL", "AML", "CML", 
    "MPN-Other", "MastCell", "MDS-5q", "MDS", "CMML", "aCML", "Histiocytic-Sarcoma",
    "AccessoryCell-Tumour"
  ))

# Bind the two data frames
tumour_names <- bind_cols(tumours_raw_filtered, tumour_name) %>%
  select(Tumours.ID, Crude.incidence.rate.per.100.000, tumor, abbreviation) %>%
  rename(Tumour = tumor) 

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
  select(Tumour, Tumours.ID, abbreviation, Crude.incidence.rate.per.100.000, Research_score, Normalised_research_score) %>%
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

# Retrieve top and bottom 10 tumorus by (normalised) Orphan drug development activity score
top_10_RS <- big_df %>%
  slice_max(order_by = as.numeric(Research_score), n = 10, with_ties = FALSE)
top_10_nRS <- big_df %>%
  slice_max(order_by = as.numeric(Normalised_research_score), n = 10, with_ties = FALSE)
bottom_10_RS <- big_df %>%
  slice_min(order_by = as.numeric(Research_score), n = 10, with_ties = FALSE)
bottom_10_nRS <- big_df %>%
  slice_min(order_by = as.numeric(Normalised_research_score), n = 10, with_ties = FALSE)

# Plot
p1 <- ggplot(top_10_RS, aes(x = reorder(abbreviation, Research_score), y = Research_score)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) +
  labs(
    x = "Tumour",
    y = "Orphan drug development score",
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
    y = "Normalised Orphan drug development score",
    title = "Baseline Analysis"
  ) +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 22),
    plot.title = element_text(size = 26)
  )

p3 <- ggplot(bottom_10_RS, aes(x = reorder(abbreviation, Research_score), y = Research_score)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) +
  labs(
    x = "Tumour",
    y = "Orphan drug development score",
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
    y = "Normalised Orphan drug development score",
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
  select(Tumour, Tumours.ID, Crude.incidence.rate.per.100.000, abbreviation, Research_score, Normalised_research_score) %>%
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

# Retrieve top and bottom 10 tumorus by (normalised) Orphan drug development activity score
top_10_RS1 <- big_df1 %>%
  slice_max(order_by = as.numeric(Research_score), n = 10, with_ties = FALSE)
top_10_nRS1 <- big_df1 %>%
  slice_max(order_by = as.numeric(Normalised_research_score), n = 10, with_ties = FALSE)
bottom_10_RS1 <- big_df1 %>%
  slice_min(order_by = as.numeric(Research_score), n = 10, with_ties = FALSE)
bottom_10_nRS1 <- big_df1 %>%
  slice_min(order_by = as.numeric(Normalised_research_score), n = 10, with_ties = FALSE)

# Plot
p1s1 <- ggplot(top_10_RS1, aes(x = reorder(abbreviation, Research_score), y = Research_score)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) +
  labs(
    x = "Tumour",
    y = "Orphan drug development score",
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
    y = "Normalised Orphan drug development score",
    title = "MA More Important"
  ) +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 22),
    plot.title = element_text(size = 26)
  )

p3s1 <- ggplot(bottom_10_RS1, aes(x = reorder(abbreviation, Research_score), y = Research_score)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) +
  labs(
    x = "Tumour",
    y = "Orphan drug development score",
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
    y = "Normalised Orphan drug development score",
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
  select(Tumour, Tumours.ID, Crude.incidence.rate.per.100.000, abbreviation, Research_score, Normalised_research_score) %>%
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

# Retrieve top and bottom 10 tumorus by (normalised) Orphan drug development activity score
top_10_RS2 <- big_df2 %>%
  slice_max(order_by = as.numeric(Research_score), n = 10, with_ties = FALSE)
top_10_nRS2 <- big_df2 %>%
  slice_max(order_by = as.numeric(Normalised_research_score), n = 10, with_ties = FALSE)
bottom_10_RS2 <- big_df2 %>%
  slice_min(order_by = as.numeric(Research_score), n = 10, with_ties = FALSE)
bottom_10_nRS2 <- big_df2 %>%
  slice_min(order_by = as.numeric(Normalised_research_score), n = 10, with_ties = FALSE)

# Plot
p1s2 <- ggplot(top_10_RS2, aes(x = reorder(abbreviation, Research_score), y = Research_score)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) +
  labs(
    x = "Tumour",
    y = "Orphan drug development score",
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
    y = "Normalised Orphan drug development score",
    title = "MA Removed"
  ) +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 22),
    plot.title = element_text(size = 26)
  )

p3s2 <- ggplot(bottom_10_RS2, aes(x = reorder(abbreviation, Research_score), y = Research_score)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) +
  labs(
    x = "Tumour",
    y = "Orphan drug development score",
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
    y = "Normalised Orphan drug development score",
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
  select(Tumour, Tumours.ID, Crude.incidence.rate.per.100.000, abbreviation, Research_score, Normalised_research_score) %>%
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

# Retrieve top and bottom 10 tumorus by (normalised) Orphan drug development activity score
top_10_RS3 <- big_df3 %>%
  slice_max(order_by = as.numeric(Research_score), n = 10, with_ties = FALSE)
top_10_nRS3 <- big_df3 %>%
  slice_max(order_by = as.numeric(Normalised_research_score), n = 10, with_ties = FALSE)
bottom_10_RS3 <- big_df3 %>%
  slice_min(order_by = as.numeric(Research_score), n = 10, with_ties = FALSE)
bottom_10_nRS3 <- big_df3 %>%
  slice_min(order_by = as.numeric(Normalised_research_score), n = 10, with_ties = FALSE)

# Plot
p1s3 <- ggplot(top_10_RS3, aes(x = reorder(abbreviation, Research_score), y = Research_score)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) +
  labs(
    x = "Tumour",
    y = "Orphan drug development score",
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
    y = "Normalised Orphan drug development activity score",
    title = "MA Less Important"
  ) +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 22),
    plot.title = element_text(size = 26)
  )

p3s3 <- ggplot(bottom_10_RS3, aes(x = reorder(abbreviation, Research_score), y = Research_score)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) +
  labs(
    x = "Tumour",
    y = "Orphan drug development score",
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
    y = "Normalised Orphan drug development activity score",
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
  title = "Sensitivity Analysis: Top 10 Rare Tumours by Orphan Drug Development Score",
  tag_levels = "A",
  theme = theme(
    plot.title = element_text(size = 32, face = "bold", hjust = 0.5, margin = margin(b = 20))
  )
)

plot2 <- (p2/p2s1) | (p2s2/p2s3) 
plot2 +
  plot_annotation(
    title = "Sensitivity Analysis: Top 10 Rare Tumours by Normalised Orphan drug development activity score",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 32, face = "bold", hjust = 0.5, margin = margin(b = 20))
    )
  )

plot3 <- (p3/p3s1) | (p3s2/p3s3) 
plot3 +
  plot_annotation(
    title = "Sensitivity Analysis: Bottom 10 Rare Tumours by Orphan Drug Development Score",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 32, face = "bold", hjust = 0.5, margin = margin(b = 20))
    )
  )

plot4 <- (p4/p4s1) | (p4s2/p4s3) 
plot4 +
  plot_annotation(
    title = "Sensitivity Analysis: Bottom 10 Rare Tumours by Normalised Orphan drug development activity score",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 32, face = "bold", hjust = 0.5, margin = margin(b = 20))
    )
  )

# Function that enables multiple (normalised) Orphan drug development activity scores and group to be added to one dataframe
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

