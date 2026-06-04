library(dplyr)
library(stringr)
library(tidyr)
library(lubridate)
library(purrr)
library(ggplot2)
library(patchwork)
library(readxl)
library(ggrepel)

# Read data and only include trials with products
tumours_raw <- read_excel("E:/data/RARECAREnet_list_of_rare_cancers.xlsx")
tumours_final <- read.csv("E:/data/Splitted2/tumours_final.csv")
products_df <- read.csv("E:/data/Matched/products.csv")
research_scores_df <- read.csv("E:/data/Matched/Research_scores.csv") %>%
  mutate(Tumours.ID = as.character(Tumours.ID))
research_scores_s1_df <- read.csv("E:/data/Matched/Research_scores_s1.csv") %>%
  mutate(Tumours.ID = as.character(Tumours.ID))
research_scores_s2_df <- read.csv("E:/data/Matched/Research_scores_s2.csv") %>%
  mutate(Tumours.ID = as.character(Tumours.ID))
research_scores_s3_df <- read.csv("E:/data/Matched/Research_scores_s3.csv") %>%
  mutate(Tumours.ID = as.character(Tumours.ID))
indications_df <- read.csv("E:/data/Matched/indications.csv") %>%
  select(NCT.ID, CT.EDU.ID, tumours_id) %>%
  semi_join(products_df, by = c("CT.EDU.ID", "NCT.ID"))
orphan_indications_df <- read.csv("E:/data/Matched/orphan_designation_better.csv") %>%
  distinct()
PubMed_df <- read.csv("E:/data/Matched/PubMed.csv") %>%
  select(Tumours.ID, PubMed.Count.General) %>%
  mutate(Tumours.ID = as.character(Tumours.ID)) %>%
  mutate(
    PubMed.Count.General = if_else(Tumours.ID == 197,
                       14365,
                       PubMed.Count.General)
  )
sponsor_df <- read.csv("E:/data/Splitted2/Sponsors.csv") %>%
  semi_join(products_df, by = c("CT.EDU.ID", "NCT.ID")) %>%
  semi_join(indications_df, by = c("CT.EDU.ID", "NCT.ID")) %>%
  mutate(
    Sponsor.Group = case_when(
      Sponsor.Type %in% c("industry", "pharmaceutical company") ~ "industry",
      TRUE ~ "other"
    )
  )
trial_information_df <- read.csv("E:/data/Splitted2/Trial_information.csv") %>%
  semi_join(products_df, by = c("CT.EDU.ID", "NCT.ID")) %>%
  semi_join(indications_df, by = c("CT.EDU.ID", "NCT.ID")) %>%
  separate_rows(Member.State, sep = "\\s*,\\s*") %>%
  # Change member state name
  mutate(
    Member.State = ifelse(
      Member.State == "bulgarian drug agency",
      "bulgaria",
      Member.State
    ))
base_groups2 <- read.csv("E:/data/Matched/Results3.csv", sep = ",") 

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

# Filter other datasets on rare tumour trials
products_df_filtered <- products_df %>%
  semi_join(indications_df, by = c("CT.EDU.ID", "NCT.ID"))

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
trial_information_df_filtered <- trial_information_df %>%
  left_join(trial_information_update, by = "CT.EDU.ID", suffix = c("", ".new")) %>%
  mutate(
    End.of.Trial = ymd(End.of.Trial),
    Trial.Status = coalesce(Trial.Status.new, Trial.Status),
    End.of.Trial = coalesce(End.of.Trial.new, End.of.Trial)
  ) %>%
  select(-Trial.Status.new, -End.of.Trial.new)

# Separate tumours from clinical trial data
indications_long <- indications_df %>%
  separate_rows(tumours_id, sep = ",") %>%
  mutate(Tumours.ID = str_trim(tumours_id),
         Tumours.ID = str_squish(Tumours.ID),
         Trial.ID = coalesce(NCT.ID, CT.EDU.ID)) %>%
  left_join(trial_information_df_filtered, by = c("NCT.ID", "CT.EDU.ID"), relationship = "many-to-many") %>%
  group_by(Tumours.ID, Trial.ID) %>%
  # Get last start date
  slice_max(Start.of.Trial, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(Trial.ID, NCT.ID, CT.EDU.ID, Tumours.ID, Start.of.Trial)

# Separate tumours from orphan designation data
orphan_long <- orphan_indications_df %>%
  separate_rows(Tumours.ID, sep = ",") %>%
  mutate(Tumours.ID = str_squish(Tumours.ID),
         Orphan.ID = coalesce(as.character(FDA.ID), EMA.ID))
  
# Research question 1:
# Determine the number of trials per RARECARE tumour
number_of_trials <- indications_long %>%
  filter(!is.na(Trial.ID)) %>%
  arrange(Tumours.ID, Trial.ID, Start.of.Trial) %>%                             # First trial first (if multiple trials under the same ID)
  #filter(Start.of.Trial >= as.Date("2015-01-01")) %>%                          # Only recent trials
  distinct(Tumours.ID, Trial.ID) %>%                                            # Remove duplicate trials (i.e., those that are done in multiple EU countries)
  group_by(Tumours.ID) %>%
  summarise(n_trials = n(), .groups = "drop")

# Determine the number of orphan designations per RARECARE tumour
number_of_orphan_designations <- orphan_long %>%
  filter(!is.na(Orphan.ID)) %>%
  #filter(Designation.Date >= as.Date("2015-01-01")) %>%                        # Only recent designations
  distinct(Tumours.ID, Orphan.ID) %>%                                           # Remove duplicate trials (i.e., those that are done in multiple EU countries)
  group_by(Tumours.ID) %>%
  summarise(n_orphan_designations = n(), .groups = "drop")

# Determine the number of orphan approvals per RARECARE tumour
number_of_orphan_approvals <- orphan_long %>%
  filter(!is.na(Orphan.ID), Orphan.Designation.Status == "approved") %>%
  distinct(Tumours.ID, Orphan.ID) %>%                                           # Remove duplicate trials (i.e., those that are done in multiple EU countries)
  group_by(Tumours.ID) %>%
  summarise(n_orphan_approvals = n(), .groups = "drop")

# Get total number of PubMed publications per tumour (looking at tumour - orphan drugs)
pubmed <- orphan_long %>%
  filter(!is.na(Drugbank.ID)) %>%
  group_by(Tumours.ID) %>%
  summarise(
    # Get total number of publications and calculate the average number of publications per drug
    total_publications = sum(as.numeric(PubMed.Count), na.rm = TRUE),
    n_drugs = n_distinct(Drugbank.ID),
    avg_publications_per_drug = round(total_publications / n_drugs, 1),
    .groups = "drop"
  )

# Combine tumour names with the number of orphan designations/trials
RQ_1 <- tumour_names %>%
  left_join(number_of_trials, by = "Tumours.ID") %>%
  left_join(number_of_orphan_approvals, by = "Tumours.ID") %>%
  left_join(number_of_orphan_designations, by = "Tumours.ID") %>%
  left_join(pubmed, by = c("Tumours.ID")) %>%
  replace_na(list(n_trials = 0, n_orphan_designations = 0, n_orphan_approvals = 0)) %>%
  mutate(total_activity = n_trials + n_orphan_designations) %>%                 # Add together
  distinct()

# Add PubMed publications
pubmed <- pubmed %>%
  left_join(RQ_1, by = c("Tumours.ID", "total_publications", "n_drugs", "avg_publications_per_drug")) %>%
  select(Tumours.ID, Tumour, Crude.incidence.rate.per.100.000, total_publications, n_drugs, avg_publications_per_drug)
  
# Get top 10 tumours based on number of publications
pubmed_publications <- PubMed_df %>%
  left_join(tumour_names, by = "Tumours.ID") %>%
  slice_max(order_by = as.numeric(PubMed.Count.General), n = 10, with_ties = FALSE)

# Retrieve top 10 tumours by activity  
top_10_MA <- RQ_1 %>%
  slice_max(order_by = as.numeric(n_orphan_approvals), n = 10, with_ties = FALSE)
top_10_OD <- RQ_1 %>%
  slice_max(order_by = as.numeric(n_orphan_designations), n = 10, with_ties = FALSE)
top_10_CT <- RQ_1 %>%
  slice_max(order_by = as.numeric(n_trials), n = 10, with_ties = FALSE)

# Max width of y-label text
wrapper <- function(x, ...) {
  paste(strwrap(x, width = 30), collapse = "\n")
}

# Plot the top 10 tumours per type of activity
p1 <- ggplot(top_10_MA, aes(x = reorder(abbreviation, n_orphan_approvals), y = n_orphan_approvals)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) +
  labs(
    title = "Marketing Authorisations",
    x = "Tumour",
    y = "Number of marketing authorisations"
  ) +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 18),
    plot.title = element_text(size = 26)
  )

p2 <- ggplot(top_10_OD, aes(x = reorder(abbreviation, n_orphan_designations), y = n_orphan_designations)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) + 
  labs(
    title = "Orphan Designations",
    x = "Tumour",
    y = "Number of orphan designations"
  ) +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 18),
    plot.title = element_text(size = 26)
  ) 

p3 <- ggplot(top_10_CT, aes(x = reorder(abbreviation, n_trials), y = n_trials)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) +
  labs(
    title = "Clinical Trials",
    x = "Tumour",
    y = "Number of clinical trials"
  ) +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 18),
    plot.title = element_text(size = 26)
  )

p4 <- ggplot(pubmed_publications, aes(x = reorder(abbreviation, PubMed.Count.General), y = PubMed.Count.General)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) +
  labs(
    title = "Published Literature Density",
    x = "Tumour",
    y = "Published literature density"
  ) +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 18),
    plot.title = element_text(size = 26)
  )

# Combine plots
plot_RQ1 <- (p2 / p1) | (p3 / p4)
plot_RQ1 +
  plot_annotation(
    title = "Top 10 Rare Tumours Across Multiple Orphan Drug Development Activity Metrics",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 28, face = "bold", hjust = 0.5, margin = margin(b = 20))
  )
)

# Only labels above a certain threshold for the scatter plot
RQ_1_label <- RQ_1 %>%
  mutate(
    label_status = case_when(
      as.numeric(n_orphan_designations) > 25 | as.numeric(n_trials) > 450 ~ Tumour,
      TRUE ~ NA_character_
    )
  )

# Make a scatter plot based on number of orphan designations and clinical trials
scatter_plot <- ggplot(RQ_1_label, aes(x = as.numeric(n_orphan_designations), 
                                 y = as.numeric(n_trials),
                                 label = Tumour)) +
  geom_point() +
  geom_jitter(width = 0.2, height = 0.2, alpha = 0.6) + 
  geom_text(aes(label = label_status), vjust = -0.5, hjust = 0.8, size = 4, check_overlap = TRUE, na.rm = TRUE) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") + 
  labs(
    title = "Relationship between orphan designations and clinical trials",
    x = "Number of orphan designations",
    y = "Number of clinical trials"
  ) +
  theme_minimal() 

# Test correlation between orphan designations and clinical trials
cor.test(
  as.numeric(RQ_1_label$n_orphan_designations),
  as.numeric(RQ_1_label$n_trials),
  method = "pearson"
)

# Research question 2:
# Divide RQ_1 results by incidence
RQ_2 <- RQ_1 %>%
  mutate(
    n_trials_100k = round(as.numeric(n_trials / sqrt(Crude.incidence.rate.per.100.000)), 0),
    n_orphan_100k = round(as.numeric(n_orphan_designations / sqrt(Crude.incidence.rate.per.100.000)), 0),
    total_activity_100k = round(as.numeric(total_activity / sqrt(Crude.incidence.rate.per.100.000)), 0),
    n_publications_100k = round(as.numeric(total_publications / sqrt(Crude.incidence.rate.per.100.000)), 0),
    n_approvals_100k = round(as.numeric(n_orphan_approvals / sqrt(Crude.incidence.rate.per.100.000)), 0)
  )
  
# Retrieve top 10 tumours by activity  
top_10_MA <- RQ_2 %>%
  slice_max(order_by = as.numeric(n_approvals_100k), n = 10, with_ties = FALSE)
top_10_OD <- RQ_2 %>%
  slice_max(order_by = as.numeric(n_orphan_100k), n = 10, with_ties = FALSE)
top_10_CT <- RQ_2 %>%
  slice_max(order_by = as.numeric(n_trials_100k), n = 10, with_ties = FALSE)
top_10_PM <- RQ_2 %>%
  slice_max(order_by = as.numeric(n_publications_100k), n = 10, with_ties = FALSE)

# Plot the top 10 tumours per type of activity
p1 <- ggplot(top_10_MA, aes(x = reorder(Tumour, n_approvals_100k), y = n_approvals_100k)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) +
  labs(
    title = "Marketing Authorisations",
    x = "Tumour",
    y = "Number of marketing authorisations"
  ) +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 18),
    plot.title = element_text(size = 26)
  )

p2 <- ggplot(top_10_OD, aes(x = reorder(Tumour, n_orphan_100k), y = n_orphan_100k)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) + 
  labs(
    title = "Orphan Designations",
    x = "Tumour",
    y = "Number of orphan designations"
  ) +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 18),
    plot.title = element_text(size = 26)
  )

p3 <- ggplot(top_10_CT, aes(x = reorder(Tumour, n_trials_100k), y = n_trials_100k)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) +
  labs(
    title = "Clinical Trials",
    x = "Tumour",
    y = "Number of clinical trials"
  ) +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 18),
    plot.title = element_text(size = 26)
  )

p4 <- ggplot(top_10_PM, aes(x = reorder(Tumour, n_publications_100k), y = n_publications_100k)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) +
  labs(
    title = "PubMed Publications",
    x = "Tumour",
    y = "Number of PubMed publications"
  ) +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 18),
    plot.title = element_text(size = 26)
  )

# Combine plots
plot_RQ2 <- (p2 / p1) | (p3 / p4)
plot_RQ2 +
  plot_annotation(
    title = "Top 10 Rare Tumours Across Multiple Incidence-Adjusted Research Activity Metrics",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 28, face = "bold", hjust = 0.5, margin = margin(b = 20))
    )
  )

# Only labels above a certain threshold for the scatter plot
RQ_2_label <- RQ_2 %>%
  mutate(
    label_status = case_when(
      as.numeric(total_activity) > 40 ~ Tumour,
      TRUE ~ NA_character_
    )
  )

# Make a scatter plot based on tumour incidence and total activity
scatter_plot <- ggplot(RQ_2_label, aes(x = as.numeric(Crude.incidence.rate.per.100.000), 
                                 y = as.numeric(total_activity),
                                 label = Tumour)) +
  geom_point() +
  geom_text(aes(label = label_status), vjust = -0.5, hjust = 0.8, size = 4, check_overlap = TRUE, na.rm = TRUE) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") + 
  geom_jitter(width = 0.2, height = 0.2, alpha = 0.6) + 
  labs(
    title = "Relationship between tumour incidence and total research activity",
    x = "Crude incidence rate per 100.000",
    y = "Total research activity"
  ) +
  theme_minimal() 

# Test correlation between orphan designations and clinical trials
cor.test(
  as.numeric(RQ_2_label$Crude.incidence.rate.per.100.000),
  as.numeric(RQ_2_label$total_activity),
  method = "pearson"
)

# Research question 3:
# Retrieve approval and withdrawal rates and combine with RQ_1 
RQ_3 <- orphan_long %>%
  filter(!is.na(Orphan.ID)) %>%
  distinct(Tumours.ID, Orphan.ID, Orphan.Designation.Status) %>%                # Remove duplicate trials (i.e., those that are done in multiple EU countries)
  group_by(Tumours.ID) %>%
  summarise(n_approved = sum(Orphan.Designation.Status == "approved"),
            n_withdrawed = sum(Orphan.Designation.Status == "withdrawn"),
            n_orphan_designations = n(), .groups = "drop") %>%
  mutate(approval_rate = round(n_approved / n_orphan_designations, 2),
         withdrawal_rate = round(n_withdrawed / n_orphan_designations, 2)) %>%
  left_join(RQ_1, by = "Tumours.ID") %>%
  left_join(RQ_2, by = c("Tumours.ID", "Tumour", "abbreviation")) %>%
  select(Tumour, Tumours.ID, abbreviation, Crude.incidence.rate.per.100.000.x, n_orphan_designations.x, n_withdrawed, n_approved, approval_rate, withdrawal_rate, total_activity.x)

# Get top 10 tumours based on approval rate
top_10_approval <- RQ_3 %>%
  slice_max(order_by = approval_rate, n = 10, with_ties = FALSE)
# Adjusted for tumour incidence
top_10_approval_100 <- RQ_3 %>%
  mutate(approval_rate_100k = round(as.numeric(approval_rate / sqrt(Crude.incidence.rate.per.100.000.x)), 0)
  ) %>%
  slice_max(order_by = approval_rate_100k, n = 10, with_ties = FALSE)

# Plot the top 10 tumours based on approval rate
p1 <- ggplot(top_10_approval, aes(x = reorder(Tumour, approval_rate), y = approval_rate)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) +
  labs(
    title = "Top 10 Rare Tumours by Orphan Designation Approval Rate",
    x = "Tumour",
    y = "Approval rate"
  ) +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    plot.title = element_text(size = 24, face = "bold")
  )

# Plot the top 10 tumours based on incidence-adjusted approval rate
p2 <- ggplot(top_10_approval_100, aes(x = reorder(Tumour, approval_rate_100k), y = approval_rate_100k)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) +
  labs(
    title = "Top 10 Rare Tumours by Incidence-Adjusted Orphan Designation Approval Rate",
    x = "Tumour",
    y = "Approval rate"
  ) +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    plot.title = element_text(size = 24, face = "bold")
  )

(p1|p2)

# Only show names above a certain threshold
RQ_3_label <- RQ_3 %>%
  filter(approval_rate > 0.5 |
           n_orphan_designations.x > 100)

# Plot number of designations to approval rate
ggplot(RQ_3, aes(
  x = n_orphan_designations.x,
  y = approval_rate)) +
  geom_point(alpha = 0.7) +
  scale_y_continuous(
    limits = c(0, max(RQ_3$approval_rate * 1.05, na.rm = TRUE)),
    breaks = scales::pretty_breaks(n = 5),
  ) +
  geom_text_repel(
    data = RQ_3_label,
    aes(label = abbreviation),
    size = 6,
    max.overlaps = Inf,
    box.padding = 0.6,
    point.padding = 0.6,
    force = 6,
    segment.alpha = 1.7,
    nudge_y = 0.03
  ) +
  scale_size_continuous(range = c(2, 10)) +
  labs(
    title = "Orphan Designations and Approval Rate across Rare Tumour Types",
    x = "Number of orphan designations",
    y = "Approval rate"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 22),
    plot.title = element_text(size = 32, face = "bold"),
    strip.text = element_text(size = 22),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 18),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

# Make a scatter plot based on approval rate and total activity
scatter_plot <- ggplot(RQ_3, aes(x = as.numeric(total_activity.x), 
                                 y = as.numeric(approval_rate),
                                 label = Tumour)) +
  geom_point() +
  geom_jitter(width = 0.2, height = 0.2, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    title = "Relationship between total research activity and orphan drug approval rate",
    x = "Total research activity",
    y = "Approval rate"
  ) +
  theme_minimal() 

# Test correlation
cor.test(
  as.numeric(RQ_3$total_activity.x),
  as.numeric(RQ_3$approval_rate),
  method = "pearson"
)

# Test correlation
cor.test(
  as.numeric(RQ_3$n_orphan_designations.x),
  as.numeric(RQ_3$approval_rate),
  method = "pearson"
)

# Research question 4:
# Retrieve start dates and trial statuses of clinical trials
trial_df <- trial_information_df_filtered %>%
  mutate(
    Start.Date = ymd(Start.of.Trial),
    End.Date = ymd(End.of.Trial),
    
    is_EU = !is.na(CT.EDU.ID),
    is_US = !is.na(NCT.ID) & is.na(CT.EDU.ID),
    
    Region = case_when(
      is_EU ~ "EU",
      is_US ~ "US",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Start.Date),
         year(Start.Date) <= 2025,
         year(End.Date) <= 2025) %>%
  distinct(CT.EDU.ID, Trial.Status, Start.Date, End.Date, .keep_all = TRUE)

# Determine how often trials are completed/ongoing/etc
#RQ_4_trial_status <- trial_df %>%
#  group_by(Trial.Status, NCT.ID, CT.EDU.ID) %>%
#  summarise(n_status = n(), .groups = "drop")

# Retrieve completed trials and year of completion
RQ_4_completed <- trial_df %>%
  filter(
    Trial.Status == "completed",
    !is.na(End.Date)
  ) %>%
  group_by(Year = year(End.Date), Region) %>%
  summarise(n_completed = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = Region,
    values_from = n_completed,
    values_fill = 0
  ) %>%
  mutate(
    n_completed_EU = EU,
    n_completed_US = US
  ) %>%
  select(Year, n_completed_EU, n_completed_US)

# Determine number of trials per year
RQ_4_trials <- trial_df %>%
  group_by(Year = year(Start.Date), Region) %>%
  summarise(n_trials = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = Region,
    values_from = n_trials,
    values_fill = 0
  ) %>%
  mutate(
    n_trials_EU = EU,
    n_trials_US = US
  ) %>%
  full_join(RQ_4_completed, by = "Year") %>%
  mutate(
    across(everything(), ~replace_na(.x, 0))) %>%
  select(Year, n_trials_EU, n_trials_US, n_completed_EU, n_completed_US)

# Retrieve orphan designation, authorisation, implementation, and withdrawal dates
orphan_dates <- orphan_indications_df %>%
  filter(!is.na(Designation.Date)) %>%
  mutate(
    Year.Designation = year(Designation.Date),
    Year.Withdrawal = year(Date.Designation.Withdrawn.or.Revoked),
    # Combine authorisation dates (US) and implementation dates (EU)
    Year.Authorisation = year(Marketing.Authorisation.Date),
    is_FDA = !is.na(FDA.ID),
    is_EMA = !is.na(EMA.ID)
  )

# Retrieve orphan designation, authorisation, implementation, and withdrawl dates
orphan_dates2 <- orphan_indications_df %>%
  separate_rows(Tumours.ID, sep = ",\\s*") %>%
  filter(!is.na(Designation.Date)) %>%
  mutate(
    Year.Designation = year(Designation.Date),
    Year.Authorisation = year(Marketing.Authorisation.Date),
    
    keep_designation = Year.Designation < 2026,
    keep_authorisation = !is.na(Year.Authorisation) & Year.Authorisation < 2026
  ) %>%
  full_join(
    base_groups2 %>%
      select(Tumours.ID, tumor_family) %>%
      mutate(Tumours.ID = as.character(Tumours.ID)),
    by = "Tumours.ID"
  )

# Designations per year
designations_time <- orphan_dates2 %>%
  filter(keep_designation) %>%
  group_by(Year.Designation, tumor_family) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(type = "Orphan Designations")

# Authorisations per year
authorisations_time <- orphan_dates2 %>%
  filter(!is.na(Year.Authorisation) & keep_authorisation) %>%
  group_by(Year.Authorisation, tumor_family) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(type = "Marketing Authorisations")

# Combine
plot_time_tumour <- bind_rows(
  designations_time %>% rename(Year = Year.Designation),
  authorisations_time %>% rename(Year = Year.Authorisation)
)

# Add important events
events <- data.frame(
  Year = c(1997, 2000, 2012, 2016),
  label = c(
    "FDAMA",
    "EMA",
    "FDASIA",
    "PRIME"
  ))

# Plot results
ggplot(plot_time_tumour,
       aes(x = Year, y = count, color = tumor_family, group = tumor_family)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  facet_wrap(~type, ncol = 1, scales = "free_y") +
  geom_vline(
    data = events,
    aes(xintercept = Year),
    linetype = "longdash",
    colour = "grey40",
    linewidth = 0.8
  ) +
  geom_text(
    data = events,
    aes(
      x = Year + 0.1,
      y = Inf,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1.5,
    size = 5
  ) +
  scale_x_continuous(
    breaks = function(x) {
      seq(
        from = floor(min(x, na.rm = TRUE) / 5) * 5,
        to   = ceiling(max(x, na.rm = TRUE) / 5) * 5,
        by   = 5
      )}) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Orphan Drug Development Activity Over Time Across Tumour Groups",
    x = "Year",
    y = "Number of events",
    color = "Tumour group"
  ) +
  theme_gray() +
  theme(
    plot.margin = margin(
      t = 40,
      r = 20,
      b = 20,
      l = 20
    ),
    text = element_text(size = 14),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    strip.text = element_text(size = 20),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 16),
    plot.title = element_text(size = 26, face = "bold"),
    legend.position = "bottom"
  )

# Determine number of designations/authorisations/implementations/withdrawals per year
RQ_4_orphan <- orphan_dates %>%
  summarise(
    n_designations = n(),
    n_designations_FDA = sum(is_FDA),
    n_designations_EMA = sum(is_EMA),
    
    n_withdrawals = sum(!is.na(Year.Withdrawal)),
    n_withdrawals_FDA = sum(is_FDA & !is.na(Year.Withdrawal)),
    n_withdrawals_EMA = sum(is_EMA & !is.na(Year.Withdrawal)),
    
    n_authorisations = sum(!is.na(Year.Authorisation)),
    n_authorisations_FDA = sum(is_FDA & !is.na(Year.Authorisation)),
    n_authorisations_EMA = sum(is_EMA & !is.na(Year.Authorisation)),
    
    .by = Year.Designation
  ) %>%
  arrange(Year.Designation) %>%
  filter(Year.Designation <= 2025)

# Designations per year
designations_time <- orphan_dates2 %>%
  filter(keep_designation) %>%
  group_by(Year.Designation, tumor_family) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(type = "Orphan Designations")

# Authorisations per year
authorisations_time <- orphan_dates2 %>%
  filter(!is.na(Year.Authorisation) & keep_authorisation) %>%
  group_by(Year.Authorisation, tumor_family) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(type = "Marketing Authorisations")

# Change to long format
RQ_4_orphan_long <- RQ_4_orphan %>%
  select(
    Year.Designation,
    n_designations_FDA,
    n_designations_EMA,
    n_authorisations_FDA,
    n_authorisations_EMA
  ) %>%
  pivot_longer(
    cols = -Year.Designation,
    names_to = "metric",
    values_to = "count"
  ) %>%
  mutate(
    Event = case_when(
      grepl("designations", metric) ~ "Designations",
      grepl("authorisations", metric) ~ "Authorisations"
    ),
    Region = case_when(
      grepl("FDA", metric) ~ "FDA",
      grepl("EMA", metric) ~ "EMA"
    )
  )

# Long format of trial data
trials_long <- RQ_4_trials %>%
  select(Year, n_trials_EU, n_trials_US, n_completed_EU, n_completed_US) %>%
  pivot_longer(
    cols = -Year,
    names_to = "metric",
    values_to = "count"
  ) %>%
  mutate(
    Event = case_when(
      grepl("n_trials", metric) ~ "Clinical Trials Started",
      grepl("n_completed", metric) ~ "Clinical Trials Completed"
    ),
    Region = case_when(
      grepl("EU", metric) ~ "EU",
      grepl("US", metric) ~ "US"
    )
  )

# Events
events <- data.frame(
   Year = c(2000, 2011, 2022),
   label = c(
     "ClinicalTrials.gov",
     "EUCTR",
     "CTIS"
   )
 )

# Combine
plot_time_tumour <- bind_rows(
  designations_time %>% rename(Year = Year.Designation),
  authorisations_time %>% rename(Year = Year.Authorisation)
)

# Combine trial information with tumour information
trial_ind <- trial_df %>%
  mutate(
    EDU.CT.ID = as.character(CT.EDU.ID),
    NCT.ID = as.character(NCT.ID)
  ) %>%
  inner_join(
    indications_df %>%
      mutate(
        EDU.CT.ID = as.character(CT.EDU.ID),
        NCT.ID = as.character(NCT.ID)
      ),
    by = c("EDU.CT.ID", "NCT.ID")
  ) %>%
  rename(Tumours.ID = tumours_id) %>%
  separate_rows(Tumours.ID, sep = ",\\s*") %>%
  full_join(
    base_groups2 %>%
      select(Tumours.ID, tumor_family) %>%
      mutate(Tumours.ID = as.character(Tumours.ID)),
    by = "Tumours.ID"
  )

# Combine with tumour names
trial_ind_tumour <- trial_ind %>%
  mutate(Tumours.ID = as.character(Tumours.ID)) %>%
  left_join(
    tumour_names,
    by = "Tumours.ID"
  ) %>%
  # Get start and completion years
  mutate(
    Start.Year = year(Start.Date),
    End.Year = ifelse(Trial.Status == "completed", year(End.Date), NA)
  )

# Count number of events per year, per tumour group
start_events <- trial_ind_tumour %>%
  filter(!is.na(Start.Year)) %>%
  group_by(Start.Year, tumor_family) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(type = "Trial Initiations")

end_events <- trial_ind_tumour %>%
  filter(!is.na(End.Year)) %>%
  group_by(End.Year, tumor_family) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(type = "Trial Completions")

# Combine
plot_trial_time <- bind_rows(
  start_events %>% rename(Year = Start.Year),
  end_events %>% rename(Year = End.Year)
)

# Plot results
ggplot(plot_trial_time,
       aes(x = Year, y = count, color = tumor_family, group = tumor_family)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  facet_wrap(~type, ncol = 1, scales = "free_y") +
  geom_vline(
    data = events,
    aes(xintercept = Year),
    linetype = "longdash",
    colour = "grey40",
    linewidth = 0.8
  ) +
  geom_text(
    data = events,
    aes(
      x = Year + 0.1,
      y = Inf,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1.5,
    size = 5
  ) +
  scale_x_continuous(
    breaks = function(x) {
      seq(
        from = floor(min(x, na.rm = TRUE) / 5) * 5,
        to   = ceiling(max(x, na.rm = TRUE) / 5) * 5,
        by   = 5
      )}) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Clinical Trial Activity Over Time Across Tumour Groups",
    x = "Year",
    y = "Number of events",
    color = "Tumour group"
  ) +
  theme_gray() +
  theme(
    plot.margin = margin(
      t = 40,
      r = 20,
      b = 20,
      l = 20
    ),
    text = element_text(size = 14),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    strip.text = element_text(size = 20),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 16),
    plot.title = element_text(size = 26, face = "bold"),
    legend.position = "bottom"
  )

# Research question 5:
# Determine time difference between orphan designation and authorisation
Time_difference <- orphan_indications_df %>%
  # Convert to right date format
  mutate(Designation.Date = ymd(Designation.Date),
         Marketing.Authorisation.Date = ymd(Marketing.Authorisation.Date),
         # Add region data field
         region = case_when(
           !is.na(EMA.ID) ~ "EU",
           !is.na(FDA.ID) ~ "US",
           TRUE ~ "Unknown"
           ),
         # Calculate time difference in years
         time_to_aut_years = as.numeric(difftime(Marketing.Authorisation.Date, Designation.Date, units = "days")) / 365.25,
         # Calculate time difference in months
         time_to_aut_months = as.numeric(difftime(Marketing.Authorisation.Date, Designation.Date, units = "days")) / 30.44,
         ) %>%
  # Only include approved/implemented orphan drugs
  filter(!is.na(time_to_aut_years))

# Determine the median/mean/sd
RQ_5 <- Time_difference %>%
  group_by(region) %>%
  summarise(
    n = n(),
    mean_time_years = round(mean(time_to_aut_years, na.rm = TRUE), 1),
    median_time_years = round(median(time_to_aut_years, na.rm = TRUE), 1),
    sd_time_years = round(sd(time_to_aut_years, na.rm = TRUE), 1)
  )

# Add tumour family groups
# base_groups <- Time_difference %>%
#   mutate(tumor_family = case_when(
#     str_detect(Tumour, "squamous") ~ "Squamous cell carcinomas",
#     str_detect(Tumour, "adeno") ~ "Adenocarcinomas",
#     str_detect(Tumour, "undifferentiated") ~ "Undifferentiated carcinomas",
#     str_detect(Tumour, "leukemia|lymphoma|Myeloma|Hodgkin|Myeloproliferative|Myelodysplastic|mast cell|histiocytic") ~ "Haematological malignancies",
#     str_detect(Tumour, "sarcoma|Ewing|Osteogenic|Chondrogenic|rhabdomyo|vascular") ~ "Sarcomas",
#     str_detect(Tumour, "germ cell|seminoma|teratoma") ~ "Germ cell tumours",
#     TRUE ~ "Other epithelial tumours"
#   )) %>%
#   filter(tumor_family != "Other epithelial tumours")

# Boxplot of region differences
p1 <- ggplot(Time_difference, aes(x = region, y = time_to_aut_years, fill = region)) +
  geom_boxplot(alpha = 0.7) +
  labs(
    title = "Time from Designation to Marketing Authorisation: EU verus US",
    x = "Region",
    y = "Time (years)"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("EU" = "coral", "US" = "steelblue")) +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    strip.text = element_text(size = 18),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 18),
    plot.title = element_text(size = 24, face = "bold"),
    legend.position = "none"
  )

# Density plot of region differences
p2 <- ggplot(Time_difference, aes(x = time_to_aut_years, fill = region, color = region)) +
  geom_density(alpha = 0.4) +
  labs(
    title = "Distribution of time to authorisation",
    subtitle = "EU vs US comparison (density plots)",
    x = "Time (years)",
    y = "Density"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    strip.text = element_text(size = 18),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 18),
    plot.title = element_text(size = 24, face = "bold")
  ) +
  scale_fill_manual(values = c("EU" = "coral", "US" = "steelblue")) +
  scale_color_manual(values = c("EU" = "coral", "US" = "steelblue"))

# Identify drugs appearing in both the EU and US
drugs_both_regions <- orphan_indications_df %>%
  filter(!is.na(Drugbank.ID)) %>%                                               # Only drugs with a Drugbank ID
  group_by(Drugbank.ID, Tumours.ID) %>%
  summarise(
    EU = any(!is.na(EMA.ID)),
    US = any(!is.na(FDA.ID)),
    # Retrieve dates
    EU_approval_date = if (any(!is.na(EMA.ID) & !is.na(Marketing.Authorisation.Date))) {
      min(ymd(Marketing.Authorisation.Date[!is.na(EMA.ID) & !is.na(Marketing.Authorisation.Date)]))
    } else {
      as.Date(NA)
    },
    US_approval_date = if (any(!is.na(FDA.ID) & !is.na(Marketing.Authorisation.Date))) {
      min(ymd(Marketing.Authorisation.Date[!is.na(FDA.ID) & !is.na(Marketing.Authorisation.Date)]))
    } else {
      as.Date(NA)
    },
    .groups = "drop"
  ) %>%
  filter(EU == TRUE & US == TRUE) %>%                                           # Only drugs from both regions
  mutate(
    # Calculate the difference in years
    time_diff_years = as.numeric(difftime(US_approval_date, EU_approval_date, units = "days")) / 365.25,
    first_region = case_when(
      EU_approval_date < US_approval_date ~ "EU first",
      US_approval_date < EU_approval_date ~ "US first",
      TRUE ~ "Same time"
    )
  ) %>%
  # Filter on valid time differences
  filter(is.finite(time_diff_years)) %>%
  distinct(Drugbank.ID, time_diff_years, first_region, EU_approval_date, US_approval_date)

# Plot histogram of results
p1 <- ggplot(drugs_both_regions, aes(x = time_diff_years, fill = first_region)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = "Time Difference Between US and EU Orphan Drug Marketing Authorisation",
    subtitle = "Drugs Authorised for the Same Tumour in Both Territories",
    x = "Time difference (years)",
    y = "Count",
    fill = "First territory"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 22),
    strip.text = element_text(size = 22),
    legend.title = element_text(size = 26),
    legend.text = element_text(size = 18),
    plot.title = element_text(size = 28, face = "bold"),
    plot.subtitle = element_text(size = 24)
  ) +
  scale_fill_manual(values = c("EU first" = "coral", "US first" = "steelblue", "Same time" = "gray"))

# Combine approval and implementation date in one column
drugs_long <- drugs_both_regions %>%
  pivot_longer(cols = c(EU_approval_date, US_approval_date),
               names_to = "region",
               values_to = "approval_date") %>%
  mutate(region = ifelse(region == "EU_approval_date", "EU", "US"))

# Plot boxplot of results
p2 <- ggplot(drugs_long, aes(x = region, y = approval_date, fill = region)) +
  geom_boxplot(alpha = 0.7) +
  labs(
    title = "Marketing Authorisation Years: EU versus US",
    x = "Region",
    y = "Year"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    strip.text = element_text(size = 18),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 18),
    plot.title = element_text(size = 24, face = "bold"),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("EU" = "coral", "US" = "steelblue"))

# Research question 6:
# Get all drug-tumour combinations
orphan_drug_tumor <- orphan_long %>%
  filter(!is.na(Drugbank.ID)) %>%
  group_by(Drugbank.ID, Tumours.ID) %>%
  summarise(
    designation_date = min(ymd(Designation.Date), na.rm = TRUE),
    .groups = "drop"
  )

# Combine clinical trial product data with tumour information data and trial information data
trial_indication_drug <- trial_information_df_filtered %>%
  left_join(products_df_filtered, by = c("CT.EDU.ID", "NCT.ID"), relationship = "many-to-many") %>%
  left_join(indications_long, by = c("CT.EDU.ID", "NCT.ID"), relationship = "many-to-many") %>%
  filter(!is.na(Drugbank.ID)) %>%
  select(CT.EDU.ID, NCT.ID, Drugbank.ID, Start.of.Trial.x, Tumours.ID) %>%
  rename(Start.of.Trial = Start.of.Trial.x) %>%
  mutate(
    trial_start_date = ymd(Start.of.Trial)
  ) %>%
  filter(!is.na(trial_start_date))                                              # Remove fields with no trial start date

# Combine with the orphan drug information
RQ_6 <- orphan_drug_tumor %>%
  inner_join(trial_indication_drug, by = c("Drugbank.ID", "Tumours.ID")) %>%
  mutate(
    # Calculate difference between clinical trial initiation and orphan drug designation date
    different_days = as.numeric(difftime(trial_start_date, designation_date, units = "days")),
    different_years = different_days / 365.25,
    # Determine whether clinical trials preceded orphan designations or vice versa
    designation_timing = case_when(
      different_years < 0 ~ "Before orphan designation",
      different_years > 0 ~ "After orphan designation",
      TRUE ~ "Same time"
    )
  ) %>%
  select(Drugbank.ID, Tumours.ID, CT.EDU.ID, NCT.ID, designation_date, trial_start_date, different_years, designation_timing
  ) %>%
  distinct()

# Histogram of results
ggplot(RQ_6, aes(x = different_years, fill = designation_timing)) +
  geom_histogram(bins = 60, alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Timing of Clinical Trial Initiation Relative to Orphan Designation",
    x = "Time (years)",
    y = "Count",
    fill = "Timing of clinical trial initiation"
  ) +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    strip.text = element_text(size = 18),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 18),
    plot.title = element_text(size = 26, face = "bold")
  ) +
  scale_fill_manual(values = c("After orphan designation" = "coral", "Before orphan designation" = "steelblue"))

# Determine the number of orphan approvals per RARECARE tumour
number_of_orphan_approvals <- orphan_long %>%
  filter(!is.na(Orphan.ID), Orphan.Designation.Status == "approved") %>%
  distinct(Tumours.ID, Orphan.ID) %>%                                           # Remove duplicate trials (i.e., those that are done in multiple EU countries)
  group_by(Tumours.ID) %>%
  summarise(n_orphan_approvals = n(), .groups = "drop")

# Get all orphan-tumour combinations
orphan_drugs <- orphan_long %>%
  distinct(Drugbank.ID, Tumours.ID) %>%
  mutate(is_orphan = TRUE)

# Combine drugs with tumour information
trial_drugs <- products_df_filtered %>%
  distinct(CT.EDU.ID, NCT.ID, Drugbank.ID) %>%
  inner_join(
    indications_df %>%
      distinct(CT.EDU.ID, NCT.ID, tumours_id),
    by = c("CT.EDU.ID", "NCT.ID")
  ) %>%
  separate_rows(tumours_id, sep = ",") %>%
  mutate(tumours_id = str_squish(tumours_id)) %>%
  filter(!is.na(Drugbank.ID),
         !is.na(tumours_id))

# Combine with orphan-tumour mapping
trial_drugs_flagged <- trial_drugs %>%
  left_join(orphan_drugs,
            by = c("Drugbank.ID", "tumours_id" ="Tumours.ID")) %>%
  rename("Tumours.ID" = "tumours_id") %>%
  mutate(is_orphan = coalesce(is_orphan, FALSE)) %>%
  distinct(Tumours.ID, CT.EDU.ID, NCT.ID, is_orphan, Drugbank.ID)

# Retrieve number of orphan drugs used in trials per tumour
trial_orphan <- trial_drugs_flagged %>%
  group_by(Tumours.ID) %>%
  summarise(
    n_orphan_true = sum(is_orphan, na.rm = TRUE),
    .groups = "drop"
  )

# Combine all important data into one data frame
big_df <- tumour_names %>%
  # Full joins with number of x and scores
  full_join(number_of_orphan_approvals, by = "Tumours.ID") %>%
  full_join(number_of_trials, by = "Tumours.ID") %>%
  full_join(number_of_orphan_designations, by = "Tumours.ID") %>%
  full_join(trial_orphan, by = "Tumours.ID") %>%
  full_join(research_scores_df, by = "Tumours.ID") %>%
  full_join(RQ_3 %>% select(Tumours.ID, approval_rate, withdrawal_rate, n_withdrawed), by = "Tumours.ID") %>%
  full_join(pubmed %>% select(Tumours.ID, total_publications), by = "Tumours.ID") %>%
  rename(total_publications_od = total_publications) %>%
  full_join(PubMed_df, by = "Tumours.ID") %>% 
  rename(total_publications = "PubMed.Count.General") %>%
  distinct() %>%
  mutate(across(c(n_orphan_designations, n_orphan_approvals, n_withdrawed, n_trials, approval_rate, withdrawal_rate, total_publications, total_publications_od, n_orphan_true),
    ~ coalesce(.x, 0))) %>%
  mutate(Tumour = as.character(Tumour)) %>%
  select(Tumour, abbreviation, Tumours.ID, Crude.incidence.rate.per.100.000, n_orphan_designations, n_orphan_approvals, approval_rate, n_withdrawed, withdrawal_rate, total_publications, total_publications_od, n_trials, n_orphan_true, Research_score, Normalised_research_score) %>%
  arrange(desc(Research_score)) %>%
  # Split scores into tertitles and put tumours in groups
  mutate(
    no_activity = n_trials == 0 & Research_score == 0,
    research_group = case_when(
      no_activity ~ "Non-researched",
      ntile(if_else(!no_activity, Research_score, NA_real_), 3) == 1 ~ "Under-researched",
      ntile(if_else(!no_activity, Research_score, NA_real_), 3) == 2 ~ "Moderately-researched",
      ntile(if_else(!no_activity, Research_score, NA_real_), 3) == 3 ~ "Well-researched"
    ),
    research_group_normalised = case_when(
      no_activity ~ "Non-researched",
      ntile(if_else(!no_activity, Normalised_research_score, NA_real_), 3) == 1 ~ "Under-researched",
      ntile(if_else(!no_activity, Normalised_research_score, NA_real_), 3) == 2 ~ "Moderately-researched",
      ntile(if_else(!no_activity, Normalised_research_score, NA_real_), 3) == 3 ~ "Well-researched"
    ),
    # Rank the tumours
    rank = rank(-Research_score, ties.method = "first"),
    rank_normalised = rank(-Normalised_research_score, ties.method = "first"),
    top_bottom_group = case_when(
      rank <= 10 ~ "Top 10",
      rank > n() - 10 ~ "Bottom 10",
      TRUE ~ NA_character_
    )
  )

# Save data
write.csv2(big_df, "Results.csv", row.names = FALSE, quote = TRUE)

# Calculate time from x to y
time <- trial_drugs_flagged %>%
  left_join(orphan_long, by = c("Tumours.ID", "Drugbank.ID")) %>%
  left_join(trial_information_df_filtered, 
            by = c("NCT.ID", "CT.EDU.ID")) %>%
  left_join(big_df, by = "Tumours.ID") %>%
  distinct() %>%
  mutate(
    CT_event = ifelse(!is.na(Start.of.Trial), 1, 0),
    MA_event = ifelse(!is.na(Marketing.Authorisation.Date), 1, 0),
    time_ODD_CT = as.numeric(as.Date(Start.of.Trial) - as.Date(Designation.Date)) / 365.25,
    time_CT_MA = as.numeric(as.Date(Marketing.Authorisation.Date) - as.Date(Start.of.Trial)) / 365.25,
    time_ODD_MA = as.numeric(as.Date(Marketing.Authorisation.Date) - as.Date(Designation.Date)) / 365.25,
    research_burden = n_trials + total_publications + n_orphan_designations
  ) 

# Calculate the number of trials before obtaining marketing authorisation
trials_to_MA <- time %>%
  filter(!is.na(Marketing.Authorisation.Date)) %>%              
  filter(!is.na(Start.of.Trial)) %>%                            
  filter(as.Date(Start.of.Trial) <= as.Date(Marketing.Authorisation.Date)) %>%  
  filter(as.Date(Start.of.Trial) >= as.Date(Designation.Date)) %>%
  group_by(Tumours.ID, Drugbank.ID) %>%
  summarise(
    trials_to_MA = n_distinct(NCT.ID, CT.EDU.ID),
    .groups = "drop"
  )

# Calculate the number of authorisations achieved after 0 trials
MA_no_prior_trial <- time %>%
  filter(!is.na(Marketing.Authorisation.Date)) %>%
  group_by(Tumours.ID, Drugbank.ID, Marketing.Authorisation.Date, Designation.Date) %>%
  summarise(
    has_prior_trial = any(
      !is.na(Start.of.Trial) &
        as.Date(Start.of.Trial) <= as.Date(Marketing.Authorisation.Date) &
        as.Date(Start.of.Trial) >= as.Date(Designation.Date)
    ),
    .groups = "drop"
  ) %>%
  filter(!has_prior_trial)

# Calculate the number of trials before obtaining orphan designation
trials_to_ODD <- time %>%
  filter(!is.na(Designation.Date)) %>%
  filter(!is.na(Start.of.Trial)) %>%
  filter(as.Date(Start.of.Trial) <= as.Date(Designation.Date)) %>%
  group_by(Tumours.ID, Drugbank.ID) %>%
  summarise(
    trials_to_ODD = n_distinct(NCT.ID, CT.EDU.ID),
    .groups = "drop"
  )

# Calculate the number of designations achieved after 0 trials
OD_no_prior_trial <- time %>%
  filter(!is.na(Designation.Date)) %>%
  group_by(Tumours.ID, Drugbank.ID, Designation.Date) %>%
  summarise(
    has_prior_trial = any(
      !is.na(Start.of.Trial) &
        as.Date(Start.of.Trial) <= as.Date(Designation.Date)
    ),
    .groups = "drop"
  ) %>%
  filter(!has_prior_trial)

# Add zeros
OD_no <- OD_no_prior_trial %>%
  mutate(type = "OD", trials = 0)

MA_no <- MA_no_prior_trial %>%
  mutate(type = "MA", trials = 0)

# Add all information together
df_cum_plot <- trials_to_MA %>%
  mutate(type = "MA", trials = trials_to_MA) %>%
  select(type, trials) %>%
  bind_rows(
    trials_to_ODD %>%
      mutate(type = "OD", trials = trials_to_ODD) %>%
      select(type, trials)
  ) %>%
  bind_rows(
    OD_no,
    MA_no
  ) %>%
  group_by(type) %>%
  arrange(trials, .by_group = TRUE) %>%
  # Calculate fractions
  mutate(
    cum_fraction = row_number() / n()
  ) %>%
  ungroup() 

# Plot the dataframe
ggplot(df_cum_plot, aes(x = trials, y = cum_fraction, color = type)) +
  geom_step(linewidth = 1) +
  coord_cartesian(xlim = c(0, 50)) +
  scale_color_discrete(
    labels = c(
      "Marketing authorisation",
      "Orphan designation"
    )
  ) +
  labs(
    title = "Proportion of Tumours Achieving Orphan Designation or Marketing Authorisation by Number of Clinical Trials",
    x = "Number of clinical trials",
    y = "Proportion reached",
    color = "Regulatory milestone"
  ) +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    strip.text = element_text(size = 18),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 18),
    plot.title = element_text(size = 24, face = "bold")
  ) 

# Change data type
big_df <- big_df %>%
  mutate(
    Crude.incidence.rate.per.100.000 = as.numeric(Crude.incidence.rate.per.100.000),
    Research_score = as.numeric(Research_score)
  )

# Only add tumour names above certain thresholds
big_df_label <- big_df %>%
  select(-top_bottom_group) %>%
  filter(Crude.incidence.rate.per.100.000 > 5 |
           n_orphan_designations > 200)

# Plot crude incidence against number of orphan designations and clinical trials
ggplot(big_df, aes(
  x = Crude.incidence.rate.per.100.000,
  y = n_orphan_designations
)) +
  scale_y_continuous(
    limits = c(0, max(big_df$n_orphan_designations * 1.05, na.rm = TRUE)),
    breaks = scales::pretty_breaks(n = 5),
  ) +
  geom_point(aes(
    size = n_trials,
  ), alpha = 0.7) +
  geom_text_repel(
    data = big_df_label,
    aes(
      x = Crude.incidence.rate.per.100.000,
      y = n_orphan_designations,
      label = abbreviation
    ),
    inherit.aes = FALSE,
    size = 6,
    max.overlaps = Inf,
    box.padding = 0.6,
    point.padding = 0.6,
    force = 6,
    segment.alpha = 0.7,
    nudge_y = 0.03
  ) +
  labs(
    x = "Tumour incidence (per 100,000)",
    y = "Number of orphan designations",
    size = "Number of clinical trials",
    title = "Tumour Incidence, Orphan Designations, and Clinical Trials by Tumour Type"
  )  +
  theme_classic() +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 22),
    plot.title = element_text(size = 32, face = "bold"),
    strip.text = element_text(size = 22),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 18),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

# Only add tumour names above certain thresholds
big_df_label <- big_df %>%
  select(-top_bottom_group) %>%
  filter(Crude.incidence.rate.per.100.000 > 5 |
           n_orphan_approvals > 10)

# Plot crude incidence against number of marketing authorisations and publications
ggplot(big_df, aes(
  x = Crude.incidence.rate.per.100.000,
  y = n_orphan_approvals
)) +
  scale_y_continuous(
    limits = c(0, max(big_df$n_orphan_approvals * 1.05, na.rm = TRUE)),
    breaks = scales::pretty_breaks(n = 5),
  ) +
  geom_point(aes(
    size = total_publications,
  ), alpha = 0.7) +
  geom_text_repel(
    data = big_df_label,
    aes(
      x = Crude.incidence.rate.per.100.000,
      y = n_orphan_approvals,
      label = abbreviation
    ),
    inherit.aes = FALSE,
    size = 6,
    max.overlaps = Inf,
    box.padding = 0.6,
    point.padding = 0.6,
    force = 6,
    segment.alpha = 0.7,
    nudge_y = 0.03
  ) +
  labs(
    x = "Tumour incidence (per 100,000)",
    y = "Number of marketing authorisations",
    size = "Published literature density",
    title = "Tumour Incidence, Marketing Authorisations, and Published Literature Density by Tumour Type"
  )  +
  theme_classic() +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 22),
    plot.title = element_text(size = 32, face = "bold"),
    strip.text = element_text(size = 22),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 18),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

big_df_label <- big_df %>%
  select(-top_bottom_group) %>%
  filter(Crude.incidence.rate.per.100.000 > 5 |
           Research_score > 0.4)

# Plot crude incidence against research score
ggplot(big_df, aes(
  x = Crude.incidence.rate.per.100.000,
  y = Research_score
)) +
  geom_point(alpha = 0.7) +
  scale_y_continuous(
    limits = c(0, max(big_df$Research_score * 1.05, na.rm = TRUE)),
    breaks = scales::pretty_breaks(n = 5),
  ) +
  geom_text_repel(
    data = big_df_label,
    aes(
      x = Crude.incidence.rate.per.100.000,
      y = Research_score,
      label = abbreviation
    ),
    inherit.aes = FALSE,
    size = 6,
    max.overlaps = Inf,
    box.padding = 0.6,
    point.padding = 0.6,
    force = 6,
    segment.alpha = 0.7,
    nudge_y = 0.03
  ) +
  labs(
    x = "Tumour incidence (per 100,000)",
    y = "Orphan drug development score",
    title = "Tumour Incidence and Orphan Drug Development Score by Tumour Type"
  )  +
  theme_classic() +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 22),
    plot.title = element_text(size = 32, face = "bold"),
    strip.text = element_text(size = 22),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 18),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

# Top 10 tumours based on (normalised) research score
top_10_RS <- big_df %>%
  slice_max(order_by = as.numeric(Research_score), n = 10, with_ties = FALSE)
top_10_nRS <- big_df %>%
  slice_max(order_by = as.numeric(Crude.incidence.rate.per.100.000), n = 10, with_ties = FALSE)
bottom_10_RS <- big_df %>%
  slice_min(order_by = as.numeric(Research_score), n = 10, with_ties = FALSE)
bottom_10_nRS <- big_df %>%
  slice_min(order_by = as.numeric(Crude.incidence.rate.per.100.000), n = 10, with_ties = FALSE)

# Find maximum values for later plotting
max_RS <- max(c(top_10_RS$Research_score,
                bottom_10_RS$Research_score))
max_nRS <- max(c(top_10_nRS$Crude.incidence.rate.per.100.000,
                 bottom_10_nRS$Crude.incidence.rate.per.100.000))

# Plot top 10 tumours based on research score
p1 <- ggplot(top_10_RS, aes(x = reorder(abbreviation, Research_score), y = Research_score)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) +
  scale_y_continuous(limits = c(0, max_RS)) +
  labs(
    title = "Orphan Drug Development Score Top 10",
    x = "Tumour",
    y = "Orphan drug development score"
  ) +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    plot.title = element_text(size = 22)
  )

# Plot top 10 tumours based on incidence
p2 <- ggplot(top_10_nRS, aes(x = reorder(abbreviation, Crude.incidence.rate.per.100.000), y = Crude.incidence.rate.per.100.000)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) +
  scale_y_continuous() +
  labs(
    title = "Incidence Top 10",
    x = "Tumour",
    y = "Tumour incidence (per 100,000)"
  ) +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    plot.title = element_text(size = 22)
  )

# Plot bottom 10 tumours based on research score
p3 <- ggplot(bottom_10_RS, aes(x = reorder(abbreviation, Research_score), y = Research_score)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) +
  scale_y_continuous(limits = c(0, max_RS)) +
  labs(
    title = "Orphan Drug Development Score Bottom 10",
    x = "Tumour",
    y = "Orphan drug development score"
  ) +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    plot.title = element_text(size = 22)
  )

# Plot bottom 10 tumours based on incidence
p4 <- ggplot(bottom_10_nRS, aes(x = reorder(abbreviation, Crude.incidence.rate.per.100.000), y = Crude.incidence.rate.per.100.000)) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = function(x) sapply(x, wrapper)) +
  scale_y_continuous() +
  labs(
    title = "Incidence Bottom 10",
    x = "Tumour",
    y = "Tumour incidence (per 100,000)"
  ) +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    plot.title = element_text(size = 22)
  )

# Add plots together
plot_RQ_main <- (p1|p3)/(p2|p4)
plot_RQ_main +
  plot_annotation(
    title = "Top and Bottom 10 Rare Tumours by Orphan Drug Development Score and Tumour Incidence",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 26, face = "bold", hjust = 0.5, margin = margin(b = 20))
    )
  )

# Plot trials against OD Trials
ggplot(big_df, aes(x = n_trials, y = n_orphan_true)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm") +
  labs(
    title = "Orphan Drug Presence in Clinical Trial Activity",
    x = "Total number of clinical trials",
    y = "Number of orphan drug trials"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    strip.text = element_text(size = 18),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 18),
    plot.title = element_text(size = 24, face = "bold")
  ) 

# Determine how many clinical trials it takes for one orphan drug marketing authorisation
# Combine trial data
trials_full <- trial_information_df_filtered %>%
  select(NCT.ID, CT.EDU.ID, Start.of.Trial) %>%
  left_join(products_df_filtered %>% select(NCT.ID, CT.EDU.ID, Drugbank.ID), by = c("NCT.ID", "CT.EDU.ID"),  relationship = "many-to-many") %>%
  left_join(indications_long %>% select(NCT.ID, CT.EDU.ID, Tumours.ID), by = c("NCT.ID", "CT.EDU.ID"),  relationship = "many-to-many") %>%
  # Add regions for later comparison
  mutate(
    region = case_when(
      !is.na(CT.EDU.ID) ~ "EU",
      TRUE ~ "US"
    )) %>%
  distinct()

# Get all US trials that were started before the drug got a marketing authorisation in the US
fda_trials_before <- orphan_long %>%
  filter(!is.na(FDA.ID)) %>%
  select(FDA.ID, Tumours.ID, Drugbank.ID, Marketing.Authorisation.Date) %>%
  distinct() %>%
  left_join(trials_full, by = c("Tumours.ID", "Drugbank.ID"), relationship = "many-to-many") %>%
  distinct() %>%
  filter(region == "US",
         Start.of.Trial < Marketing.Authorisation.Date) %>%
  group_by(Tumours.ID, Drugbank.ID) %>%
  summarise(
    n_trials_before_FDA = n_distinct(NCT.ID),
    .groups = "drop"
  )
  
# Get all EU trials that were started before the drug got a marketing authorisation in the EU
ema_trials_before <- orphan_long %>%
  filter(!is.na(EMA.ID)) %>%
  select(EMA.ID, Tumours.ID, Drugbank.ID, Marketing.Authorisation.Date) %>%
  distinct() %>%
  left_join(trials_full, by = c("Tumours.ID", "Drugbank.ID"), relationship = "many-to-many") %>%
  distinct() %>%
  filter(region == "EU",
         Start.of.Trial < Marketing.Authorisation.Date) %>%
  group_by(Tumours.ID, Drugbank.ID) %>%
  summarise(
    n_trials_before_EMA = n_distinct(CT.EDU.ID),
    .groups = "drop"
  )

# Determine how many trials were needed before a marketing authorisation was given
comparison_df <- ema_trials_before %>%
  full_join(fda_trials_before, by = c("Tumours.ID", "Drugbank.ID"), relationship = "many-to-many") %>%
  distinct() 

# Get all US drugs that were given a market authorisation before clinical trials were started
fda_trials_after <- orphan_long %>%
  filter(!is.na(FDA.ID)) %>%
  select(FDA.ID, Tumours.ID, Drugbank.ID, Marketing.Authorisation.Date) %>%
  distinct() %>%
  left_join(trials_full, by = c("Tumours.ID", "Drugbank.ID"), relationship = "many-to-many") %>%
  distinct() %>%
  filter(region == "US") %>%
  mutate(is_before = Start.of.Trial < Marketing.Authorisation.Date) %>%
  group_by(Tumours.ID, Drugbank.ID) %>%
  summarise(
    n_trials_before_FDA = sum(is_before, na.rm = TRUE),
    n_trials_after_FDA = sum(!is_before & !is.na(Start.of.Trial)),
    .groups = "drop"
  ) %>%
  filter(n_trials_before_FDA == 0) %>%
  select(Tumours.ID, Drugbank.ID, n_trials_after_FDA) %>%
  filter(!is.na(n_trials_after_FDA),
         n_trials_after_FDA != 0)

# Get all EU drugs that were given a market authorisation before clinical trials were started
ema_trials_after <- orphan_long %>%
  filter(!is.na(EMA.ID)) %>%
  select(EMA.ID, Tumours.ID, Drugbank.ID, Marketing.Authorisation.Date) %>%
  distinct() %>%
  left_join(trials_full, by = c("Tumours.ID", "Drugbank.ID"), relationship = "many-to-many") %>%
  distinct() %>%
  filter(region == "EU") %>%
  mutate(is_before = Start.of.Trial < Marketing.Authorisation.Date) %>%
  group_by(Tumours.ID, Drugbank.ID) %>%
  summarise(
    n_trials_before_EMA = sum(is_before, na.rm = TRUE),
    n_trials_after_EMA = sum(!is_before & !is.na(Start.of.Trial)),
    .groups = "drop"
  ) %>%
  filter(n_trials_before_EMA == 0) %>%
  select(Tumours.ID, Drugbank.ID, n_trials_after_EMA) %>%
  filter(!is.na(n_trials_after_EMA),
         n_trials_after_EMA != 0)

# Determine how many trials were started after a marketing authorisation was given
comparison_df_1 <- ema_trials_after %>%
  full_join(fda_trials_after, by = c("Tumours.ID", "Drugbank.ID"), relationship = "many-to-many") %>%
  full_join(comparison_df, by = c("Tumours.ID", "Drugbank.ID"), relationship = "many-to-many") %>%
  distinct() 

# Determine if trials are done multicentre or unicentre
trial_centre_classification <- trial_information_df_filtered %>%
  group_by(CT.EDU.ID) %>%
  filter(!is.na(Member.State)) %>%
  summarise(
    n_countries = n_distinct(Member.State),
    n_status = n_distinct(Trial.Status),
    n_phase = n_distinct(Trial.Phase),
    .groups = "drop"
  )

# Add trial type to EU trials
trial_centres <- trial_centre_classification %>%
  mutate(
    trial_type = case_when(
      n_countries > 1 & n_status == 1 & n_phase == 1 ~ "multicentre",
      n_countries == 1 ~ "unicentre",
      TRUE ~ "multicountry"
    )
  )

# Compare trial phases between multicentre and unicentre
trial_phase_centre <- trial_information_df_filtered %>%
  left_join(trial_centres %>% select(CT.EDU.ID, trial_type),
            by = "CT.EDU.ID") %>%
  filter(!is.na(trial_type),
         trial_type != "multicountry") %>%
  mutate(
    Trial.Phase = ifelse(is.na(Trial.Phase), "NA", Trial.Phase)
  ) %>%
  group_by(trial_type, Trial.Phase) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = Trial.Phase,
    values_from = n,
    values_fill = 0) %>%
  # Normalise the counts
  mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(across(where(is.numeric), ~ round(.x / total, 2), .names = "{.col}_prop")) %>%
  select(-total)

# Change to long format and recode names
trial_phase_long <- trial_phase_centre %>%
  pivot_longer(
    cols = c(I, `II`, `III`, `IV`, `I and II`, `I and III`, `II and III`, `III and IV`),
    names_to = "trial_phase",
    values_to = "value"
  ) %>%
  mutate(trial_phase = recode(trial_phase,
                              "I_prop" = "I",
                              "II_prop" = "II",
                              "III_prop" = "III",
                              "IV_prop" = "IV",
                              "I and II_prop" = "I and II",
                              "I and III_prop" = "I and III",
                              "II and III_prop" = "II and III",
                              "III and IV_prop" = "III and IV"
                              ),
         trial_type = recode(trial_type,
                             "unicentre" = "Unicentre",
                             "multicentre" = "Multicentre"
                             )
         )

# Plot the results
p1 <- ggplot(trial_phase_long, aes(x = trial_type, y = value, fill = trial_phase)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    x = "Centre type",
    y = "Count",
    fill = "Trial phase",
    title = "Trial Phases"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 22),
    strip.text = element_text(size = 22),
    legend.title = element_text(size = 26),
    legend.text = element_text(size = 22),
    plot.title = element_text(size = 26)
  ) 

# Compare trial statuses between multicentre and unicentre
trial_status_centre <- trial_information_df_filtered %>%
  left_join(trial_centres %>% select(CT.EDU.ID, trial_type),
            by = "CT.EDU.ID") %>%
  filter(!is.na(trial_type),
         trial_type != "multicountry") %>%
  mutate(
    Trial.Status = ifelse(is.na(Trial.Status), "NA", Trial.Status)
  ) %>%
  group_by(trial_type, Trial.Status) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = Trial.Status,
    values_from = n,
    values_fill = 0) %>%
  # Normalise the counts
  mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(across(where(is.numeric), ~ round(.x / total, 2), .names = "{.col}_prop")) %>%
  select(-total)

# Change to long format and recode names
trial_status_long <- trial_status_centre %>%
  pivot_longer(
    cols = c(authorised, completed, ongoing, `prematurely ended`, `temporarily halted`, `trial now transitioned`, `eea`),
    names_to = "trial_status",
    values_to = "value"
  ) %>%
  mutate(trial_status = recode(trial_status,
                              "authorised" = "Authorised",
                              "completed" = "Completed",
                              "ongoing" = "Ongoing",
                              "prematurely ended" = "Prematurely ended",
                              "temporarily halted" = "Temporarily halted",
                              "trial now transitioned" = "Transitioned",
                              "eea" = "EEA"
  ),
  trial_type = recode(trial_type,
                      "unicentre" = "Unicentre",
                      "multicentre" = "Multicentre"
  )
  )

# Plot the results
p2 <- ggplot(trial_status_long, aes(x = trial_type, y = value, fill = trial_status)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    x = "Centre type",
    y = "Count",
    fill = "Trial status",
    title = "Trial Statuses"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 22),
    strip.text = element_text(size = 22),
    legend.title = element_text(size = 26),
    legend.text = element_text(size = 22),
    plot.title = element_text(size = 26)
  ) 

# Plot together
plot_centres <- (p1|p2)
plot_centres +
  plot_annotation(
    title = "Counts of Trial Characteristics by Centre Type",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 28, face = "bold", hjust = 0.5, margin = margin(b = 20))
    )
  )

# Determine which EU country has the most trials
trials_per_country <- trial_information_df_filtered %>%
  distinct(CT.EDU.ID, Member.State) %>%
  count(Member.State, sort = TRUE) %>%
  filter(!is.na(Member.State)) %>%
  mutate(
    Member.State = case_when(
      Member.State == "uk" ~ "United Kingdom",
      TRUE ~ str_to_title(Member.State))
  )

# Plot the results
ggplot(trials_per_country, aes(x = reorder(Member.State, n), y = n)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Member state",
    y = "Number of trials",
    title = "Number of Clinical Trials per European Country"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 22),
    strip.text = element_text(size = 22),
    legend.title = element_text(size = 26),
    legend.text = element_text(size = 22),
    plot.title = element_text(size = 28, face = "bold")
  ) 

# Determine trial status per country
status_per_country <- trial_information_df_filtered %>%
  distinct(CT.EDU.ID, Member.State, Trial.Status) %>%
  filter(!is.na(Member.State)) %>%
  mutate(
    Member.State = case_when(
      Member.State == "uk" ~ "United Kingdom",
      TRUE ~ str_to_title(Member.State)),
    Trial.Status = case_when(
      Trial.Status == "eea" ~ "EEA",
      Trial.Status == "trial now transitioned" ~ "Transitioned",
      Trial.Status == "recruitment" ~ "Recruiting",
      TRUE ~ str_to_sentence(Trial.Status))
  ) %>%
  count(Member.State, Trial.Status) %>%
  group_by(Member.State) %>%
  mutate(prop = n / sum(n)) %>%
  group_by(Member.State) %>%
  mutate(total = sum(n)) %>%
  ungroup() %>%
  # Sort results from highest to lowest
  mutate(Member.State = factor(Member.State,
                               levels = names(sort(tapply(n, Member.State, sum),
                                                   decreasing = FALSE))))

# Plot the results
p1 <- ggplot(status_per_country, aes(x = Member.State, y = n, fill = Trial.Status)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Member state",
    y = "Count",
    fill = "Trial status",
    title = "Trial Statuses"
  )  +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    strip.text = element_text(size = 18),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 18),
    plot.title = element_text(size = 22)
  ) 

# Determine trial phase per country
phase_per_country <- trial_information_df_filtered %>%
  distinct(CT.EDU.ID, Member.State, Trial.Phase) %>%
  filter(!is.na(Member.State),
         !is.na(Trial.Phase)) %>%
  mutate(
    Member.State = case_when(
      Member.State == "uk" ~ "United Kingdom",
      TRUE ~ str_to_title(Member.State)
    )
  ) %>%
  count(Member.State, Trial.Phase) %>%
  group_by(Member.State) %>%
  mutate(prop = n / sum(n)) %>%
  group_by(Member.State) %>%
  mutate(total = sum(n)) %>%
  ungroup() %>%
  # Sort results from highest to lowest
  mutate(Member.State = factor(Member.State,
                               levels = names(sort(tapply(n, Member.State, sum),
                                                   decreasing = FALSE))))

# Plot the results
p2 <- ggplot(phase_per_country, aes(x = Member.State, y = n, fill = Trial.Phase)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Member state",
    y = "Count",
    fill = "Trial phase",
    title = "Trial Phases"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    strip.text = element_text(size = 18),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 18),
    plot.title = element_text(size = 22)
  ) 

# Plot together
(p2 | p1) +
  plot_annotation(
    title = "Counts of Trial Characteristics across European Countries",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 30, face = "bold", hjust = 0.5, margin = margin(b = 20))
    )
  )

# Merge sponsor information with trial information
sponsor_trial <- sponsor_df %>%
  inner_join(
    trial_information_df_filtered,
    by = c("NCT.ID", "CT.EDU.ID"),
    relationship = "many-to-many"
  )

# Determine number of trials per sponsor group
trials_per_sponsor <- sponsor_trial %>%
  distinct(NCT.ID, CT.EDU.ID, Sponsor.Group) %>%
  count(Sponsor.Group, sort = TRUE)

# Determine trial phase per sponsor group
phase_per_sponsor <- sponsor_trial %>%
  distinct(NCT.ID, CT.EDU.ID, Sponsor.Group, Trial.Phase) %>%
  filter(!is.na(Trial.Phase)) %>%
  count(Sponsor.Group, Trial.Phase) %>%
  mutate(Sponsor.Group = str_to_sentence(Sponsor.Group)) %>%
  group_by(Sponsor.Group) %>%
  mutate(prop = n / sum(n))

# Plot the results
p1 <- ggplot(phase_per_sponsor, aes(x = Sponsor.Group, y = n, fill = Trial.Phase)) +
  geom_col() +
  labs(
    x = "Sponsor type",
    y = "Count",
    fill = "Trial phase",
    title = "Trial Phases"
  )  +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 22),
    strip.text = element_text(size = 22),
    legend.title = element_text(size = 26),
    legend.text = element_text(size = 22),
    plot.title = element_text(size = 26)
  ) 

# Determine trial status per sponsor group
status_per_sponsor <- sponsor_trial %>%
  distinct(NCT.ID, CT.EDU.ID, Sponsor.Group, Trial.Status) %>%
  mutate(
    Trial.Status = case_when(
    Trial.Status == "eea" ~ "EEA",
    TRUE ~ str_to_sentence(Trial.Status)),
    Sponsor.Group = str_to_sentence(Sponsor.Group)
    ) %>%
  count(Sponsor.Group, Trial.Status) %>%
  group_by(Sponsor.Group) %>%
  mutate(prop = n / sum(n))

# Plot the results
p2 <- ggplot(status_per_sponsor, aes(x = Sponsor.Group, y = n, fill = Trial.Status)) +
  geom_col() +
  labs(
    x = "Sponsor type",
    y = "Count",
    fill = "Trial status",
    title = "Trial Statuses"
  )  +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 22),
    strip.text = element_text(size = 22),
    legend.title = element_text(size = 26),
    legend.text = element_text(size = 22),
    plot.title = element_text(size = 26)
  ) 

# Plot together
(p1|p2) +
  plot_annotation(
    title = "Counts of Trial Characteristics across Sponsor Type",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 30, face = "bold", hjust = 0.5, margin = margin(b = 20))
    )
  )

# Calculate number of trials per trial phase and research group
trial_phase_plot_df <- trial_information_df_filtered %>%
  left_join(
    indications_long %>%
      select(CT.EDU.ID, NCT.ID, Tumours.ID) %>%
      distinct(),
    by = c("CT.EDU.ID", "NCT.ID")
    ) %>%
  left_join(
      big_df %>% select(Tumours.ID, research_group),
      by = "Tumours.ID"
      ) %>%
  distinct() %>%
  filter(
    !is.na(research_group),
    research_group != "Non-researched"
    ) %>%
  mutate(
    Trial.Phase = ifelse(is.na(Trial.Phase), "Unknown", Trial.Phase),
    research_group = recode(research_group,
                            "Moderately-researched" = "Moderately",
                            "Well-researched" = "Well",
                            "Under-researched" = "Under")
    ) %>%
  count(research_group, Trial.Phase) %>%
  group_by(research_group) %>%
  mutate(prop = n / sum(n))

# Plot the results
p1 <- ggplot(trial_phase_plot_df,
       aes(x = research_group, y = n, fill = Trial.Phase)) +
  geom_col(width = 0.7) +
  labs(
    x = "Research intensity group",
    y = "Count",
    fill = "Trial phase",
    title = "Trial Phases"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    strip.text = element_text(size = 18),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 18),
    plot.title = element_text(size = 22)
  ) 

# Calculate number of trials per trial status and research group
trial_status_plot_df <- trial_information_df_filtered %>%
  left_join(
    indications_long %>%
      select(CT.EDU.ID, NCT.ID, Tumours.ID) %>%
      distinct(),
    by = c("CT.EDU.ID", "NCT.ID")
  ) %>%
  left_join(
    big_df %>% select(Tumours.ID, research_group),
    by = "Tumours.ID"
  ) %>%
  distinct() %>%
  filter(
    !is.na(research_group),
    research_group != "Non-researched"
  ) %>%
  mutate(
    Trial.Status = ifelse(is.na(Trial.Status), "unknown", Trial.Status),
    Trial.Status = recode(Trial.Status,
                                 "authorised" = "Authorised",
                                 "completed" = "Completed",
                                 "ongoing" = "Ongoing",
                                 "prematurely ended" = "Prematurely ended",
                                 "recruitment" = "Recruiting",
                                 "prohibited" = "Prohibited",
                                 "restarted" = "Restarted",
                                 "temporarily halted" = "Temporarily halted",
                                 "trial now transitioned" = "Transitioned",
                                 "eea" = "EEA",
                                 "unknown" = "Unknown"),
    research_group = recode(research_group,
                          "Moderately-researched" = "Moderately",
                          "Well-researched" = "Well",
                          "Under-researched" = "Under")
  ) %>%
  count(research_group, Trial.Status) %>%
  group_by(research_group) %>%
  mutate(prop = n / sum(n))

# Plot the results
p2 <- ggplot(trial_status_plot_df,
       aes(x = research_group, y = n, fill = Trial.Status)) +
  geom_col(width = 0.7) +
  labs(
    x = "Research intensity group",
    y = "Count",
    fill = "Trial status",
    title = "Trial Statuses"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    strip.text = element_text(size = 18),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 18),
    plot.title = element_text(size = 22)
  ) 

# Plot together
(p1|p2) +
  plot_annotation(
    title = "Counts of Trial Characteristics by Research Intensity Group",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 28, face = "bold", hjust = 0.5, margin = margin(b = 20))
    )
  )
