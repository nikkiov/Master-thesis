library(xml2)
library(dplyr)
library(purrr)
library(XML)
library(stringr)
library(stringdist)

# Read data
doc <- read_xml("Master Thesis/full database.xml")
xml_structure <- xml_children(doc)

ns <- xml_ns(doc)
xsd <- read_xml("Master Thesis/drugbank.xsd")
xml_find_all(xsd, ".//xs:element", xml_ns(xsd)) |>
  xml_attr("name")

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
  mech_action = sapply(drugs, function(x)
    xml_text(xml_find_first(x, ".//d1:action", ns))
  ),
  target = sapply(drugs, function(x)
    xml_text(xml_find_first(x, ".//d1:targets", ns))
  ),
  indication = sapply(drugs, function(x)
    xml_text(xml_find_first(x, ".//d1:indication", ns))
  ),
  classification = sapply(drugs, function(x)
    xml_text(xml_find_first(x, ".//d1:classification/d1:class", ns))
  ),
  type = sapply(drugs, function(x)
    xml_attr(x, "type")
  ),
  categories = sapply(drugs, function(x)
    xml_text(xml_find_first(x, ".//d1:category", ns))
  ),
  groups = sapply(drugs, function(x)
    xml_text(xml_find_first(x, ".//d1:group", ns))
  ),
  classification2 = sapply(drugs, function(x)
    xml_text(xml_find_first(x, ".//d1:groups/d1:group", ns))
  ),
  synonyms = sapply(drugs, function(x)
    collapse_nodes(x, ".//d1:synonym", ns)
  ),
  products = sapply(drugs, function(x)
    collapse_nodes(x, ".//d1:product/d1:name", ns)
  ),
  international_brands = sapply(drugs, function(x)
    collapse_nodes(x, ".//d1:international-brand", ns)
  ),
  fda_label = sapply(drugs, function(x)
    collapse_nodes(x, ".//d1:fda-label", ns)
  ),
  pubmed_id = sapply(drugs, function(x)
    collapse_nodes(x, ".//d1:pubmed-id", ns)
  ),
  ema_product_code = sapply(drugs, function(x)
    collapse_nodes(x, ".//d1:ema-product-code", ns)
  ),
  ema_ma_number = sapply(drugs, function(x)
    collapse_nodes(x, ".//d1:ema-ma-number", ns)
  ),
  fda_application_number = sapply(drugs, function(x)
    collapse_nodes(x, ".//d1:fda-application-number", ns)
  ),
  stringsAsFactors = FALSE
)

# Add NAs
drug_df <- drug_df %>%
  mutate(across(everything(), ~na_if(.x, "")))

# Save data
write.csv(drug_df, "Drugbank.csv", row.names = FALSE)
