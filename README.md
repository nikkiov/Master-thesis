# Master-thesis

# Cleaning.R: A script that cleans database fields and harmonises them between databases. The output is a filtered version of the databases.
# Splitting.R: A script that splits the complex clinical trial databases into tumour information, drug information and trial information dataframes. The output is these three different databases.
# DrugBank.R: A script that reads the DrugBank database data and retrieves the important fields. The output is part of the DrugBank database.
# Drug_Matching.R: A script that loops through all drugs used in the orphan drug designation and clinical trial databases and matches them to DrugBank IDs, using string matching techniques.
# Tumour_Matching.R: A script that loops through all tumours from in the orphan drug designation and clinical trial databases and matches them to the RARECARE-Orphanet nomenclature, using string matching techniques.
# Pubmed.R: A script used to retrieve PubMed results using a custom search query and API. 
# Matches.R: A script used to manually edit inaccurate scores from the matching techniques and determine the threshold value of a good-bad match. The output is the orphan drug designation and clinical trial databases with standardised tumour/drug names and identifiable IDs. 

# For the analysis, the following scripts were used:
# Main_Question.R: A script used to calculate the orphan drug development score.
# Analysis.R: A script used to do all analyses and create visualisations of the results. 
# Sensitivity.R: A script used to calculate the orphan drug development scores with differing component weights.
# Groups.R: A script used to determine which tumours are part of which tumour family groups. The output is Analysis_File.csv.
# Analysis_File.csv: The final file consisting of all tumour scores and tumour family groups.
