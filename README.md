# Blood resistome

Blood_dataset.tsv: file of blood isolates information derived from PATRIC Database. This table was used to download proteins from blood isolates.

table_filter_blood_review.py: Python script used to filter all PATRIC Database resulting in Blood_dataset.tsv

  Usage: python3 table_filter_blood_review input_table_from_patric_database.tsv

resistome.py: Python script used to analyse results from hmmsearch from blood isolates against Resfams HMM profile 
