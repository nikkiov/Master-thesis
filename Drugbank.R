library(xml2)
library(dplyr)
library(purrr)
library(XML)

# Read data
FDA_df <- read.csv("E:/data/First Filtering2/FDA.csv")
EMA_df <- read.csv("E:/data/First Filtering2/EMA.csv")

xsd <- read_xml("C:/Users/Nikki/Downloads/drugbank(1).xsd")
xml_structure <- xml_children(xsd)
xml_name(xml_structure)
xml_find_all(xsd, ".//xs:element[contains(translate(@name,'ABCDEFGHIJKLMNOPQRSTUVWXYZ','abcdefghijklmnopqrstuvwxyz'),'fda')]", xml_ns(xsd)) |>
  xml_attr("name")
doc <- read_xml("E:/Master thesis Data/drugbank_all_full_database.xml/full database.xml")
ns <- xml_ns(doc)

# Select specific nodes
drugs <- xml_find_all(doc, ".//d1:drug", ns)

# Helper function
collapse_nodes <- function(node, xpath, ns){
  xml_find_all(node, xpath, ns) |>
    xml_text() |>
    paste(collapse = ";")
}

# Extract data frame
drug_df <- data.frame(
  drugbank_id = sapply(drugs, function(x)
    xml_text(xml_find_first(x, ".//d1:drugbank-id[@primary='true']", ns))
  ),
  name = sapply(drugs, function(x)
    xml_text(xml_find_first(x, ".//d1:name", ns))
  ),
  cas_number = sapply(drugs, function(x)
    xml_text(xml_find_first(x, ".//d1:cas-number", ns))
  ),
  indication = sapply(drugs, function(x)
    xml_text(xml_find_first(x, ".//d1:indication", ns))
  ),
  classification = sapply(drugs, function(x)
    xml_text(xml_find_first(x, ".//d1:classification/d1:description", ns))
  ),
  synonyms = sapply(drugs, function(x)
    collapse_nodes(x, ".//d1:synonym", ns)
  ),
  atc_codes = sapply(drugs, function(x)
    collapse_nodes(x, ".//d1:atc-code/d1:code", ns)
  ),
  pubmed_ids = sapply(drugs, function(x)
    collapse_nodes(x, ".//d1:pubmed-id", ns)
  ),
  external_ids = sapply(drugs, function(x)
    collapse_nodes(x, ".//d1:external-identifier/d1:identifier", ns)
  ),
  products = sapply(drugs, function(x)
    collapse_nodes(x, ".//d1:product/d1:name", ns)
  ),
  stringsAsFactors = FALSE
)

# Add NAs
drug_df <- drug_df %>%
  mutate(across(everything(), ~na_if(.x, "")))
