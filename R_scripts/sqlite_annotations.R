# An R script to take sqlite databases from Trinotate and filter by bit scores and E values for annotation reports.
library(DBI)
library(tidyverse)

# A function to convert sqlite databases and xls annotation reports from Trinotate to tables, and filter by both bit scores and E values
sqlite_to_table <- function(sqlite_db, xls_annotations, bit_score, E_value){
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = sqlite_db)
  dbListTables(con)
  
  con_db <- tbl(con, "BlastDbase")
  con_db <- as_tibble(con_db)
  filtered_con <- con_db %>% 
    filter(BitScore > bit_score & Evalue < E_value)

  sqlite_filtered <- filtered_con %>% 
    select(TrinityID, Evalue, BitScore) %>% 
    collect() 
  
  dbDisconnect(con)
  
  xls_trinotate <- read_tsv(xls_annotations)
  
  print("Number of unique genes")
  print(length(unique(xls_trinotate$`#gene_id`))) # Number of unique genes
  print("Number of unique transcripts")
  print(length(unique(xls_trinotate$transcript_id))) # Number of unique transcripts
  
  sqlite_filtered <- left_join(sqlite_filtered, xls_trinotate, by = c("TrinityID" = "transcript_id"))
  
  
  # Separate th sprot top blastx hit column into various other columns ID1, ID2, Q, perctID, and so on. Separated by the ^
  # While some transcripts have many, many annotations associated with them, they are sorted with the lowest E-values reported first. Therefore splitting the blastx string and keeping only the first one will allow us to keep the annotation we are most confident in.
  blastx <- sqlite_filtered %>% 
    separate(sprot_Top_BLASTX_hit, c("sp.BX_UniProt.ID.1", "sp.BX_UniProt.ID.2", "sp.BX_Q.H", "sp.BX_perctID", 
                                     "sp.BX_E", "sp.BX_RecName", "sp.BX_Tax.Lineage"), sep = "\\^") %>% 
    select(-sp.BX_UniProt.ID.2)  %>% # remove column since same as ...ID.1
    dplyr::rename(sp.BX_UniProt.ID = sp.BX_UniProt.ID.1) %>%  # rename UniProt_ID column
    separate(sp.BX_E, c("E", "sp.BX_E"), sep = "\\:") %>% # Separating the E's apart from E values to we can sort by quality
    select(-E)
  blastx$sp.BX_E <- as.numeric(blastx$sp.BX_E)
  
  # Cleaning up the transcript names and percent ID columns by separating the useful information from the blastx headers
  blastx <- blastx %>% 
    separate(sp.BX_RecName, c("recname", "sp.BX_transcript_name"), sep = "=") %>% 
    select(-recname)
  # Removing the semicolons from transcript descriptions just because they don't look good
  blastx$sp.BX_transcript_name <- str_replace(blastx$sp.BX_transcript_name, ";", "")
  
  # Return the blastx table as our final output
  return(blastx)
}

# Use the filtering and conversion function for the 13 tissues.
# Brain
brain_blastx <- sqlite_to_table(sqlite_db = "Trinotate_brain.sqlite", xls_annotations = "trinotate_brain.xls", bit_score = 50, E_value = 1e-6)
# Anterior intestine
ant_intestine_blastx <- sqlite_to_table(sqlite_db = "Trinotate.AW19_ST26_anterior_intestine.sqlite", xls_annotations = "trinotate_anterior_intestine.xls", bit_score = 50, E_value = 1e-6)
# Pyloric ceca
pyloric_ceca_blastx <- sqlite_to_table(sqlite_db = "Trinotate.ST_pyloric_ceca.sqlite", xls_annotations = "trinotate_pyloric_ceca.xls", bit_score = 50, E_value = 1e-6)
# Rectum
rectum_blastx <- sqlite_to_table(sqlite_db = "Trinotate.ST_rectum.sqlite", xls_annotations = "trinotate_rectum.xls", bit_score = 50, E_value = 1e-6)
# Spiral valve
spiral_valve_blastx <- sqlite_to_table(sqlite_db = "Trinotate.ST_spiral_valve.sqlite", xls_annotations = "trinotate_spiral_valve.xls", bit_score = 50, E_value = 1e-6)
# Distal stomach
stomach_distal_blastx <- sqlite_to_table(sqlite_db = "Trinotate.ST_stomach_distal.sqlite", xls_annotations = "trinotate_stomach_distal.xls", bit_score = 50, E_value = 1e-6)
# Proximal stomach
stomach_proximal_blastx <- sqlite_to_table(sqlite_db = "Trinotate.ST_stomach_proximal.sqlite", xls_annotations = "trinotate_stomach_proximal.xls", bit_score = 50, E_value = 1e-6)
# Gill
gill_blastx <- sqlite_to_table(sqlite_db = "Trinotate_gill.sqlite", xls_annotations = "trinotate_gill.xls", bit_score = 50, E_value = 1e-6)
# Head kidney
head_kidney_blastx <- sqlite_to_table(sqlite_db = "Trinotate_head_kidney.sqlite", xls_annotations = "trinotate_head_kidney.xls", bit_score = 50, E_value = 1e-6)
# Heart
heart_blastx <- sqlite_to_table(sqlite_db = "Trinotate_heart.sqlite", xls_annotations = "trinotate_heart.xls", bit_score = 50, E_value = 1e-6)
# Liver
liver_blastx <- sqlite_to_table(sqlite_db = "Trinotate_liver.sqlite", xls_annotations = "trinotate_liver.xls", bit_score = 50, E_value = 1e-6)
# White muscle
white_muscle_blastx <- sqlite_to_table(sqlite_db = "Trinotate_white_muscle.sqlite", xls_annotations = "trinotate_white_muscle.xls", bit_score = 50, E_value = 1e-6)
# Esophagus
esophagus_blastx <- sqlite_to_table(sqlite_db = "Trinotate.ST_esophagus.sqlite", xls_annotations = "trinotate_esophagus.xls", bit_score = 50, E_value = 1e-6)

# Head kidney 2, the re-done head kidney annotation
head_kidney2_blastx <- sqlite_to_table(sqlite_db = "Trinotate_head_kidney2.sqlite", xls_annotations = "head_kidney2_report_filtered.xls", bit_score = 50, E_value = 1e-6)

# The overall transcriptome with data from all 13 tissues
overall_blastx <- sqlite_to_table(sqlite_db = "Trinotate_overall.sqlite", xls_annotations = "lkst_overall_report.xls", bit_score = 50, E_value = 1e-6)


# Writing the uniprot gene IDs to the uniprot retrieve ID/mapping tool for conversion to gene names, which are much more likely to be found in the EnrichR database down the line
# https://www.uniprot.org/uploadlists/
write_delim(as.data.frame(unique(brain_blastx$sp.BX_UniProt.ID)), "brain_uniprot_IDs.txt", col_names = FALSE, delim = "\t")
# Read in the uniprot gene names with gene IDs from above, for attachment to the annotation spreadsheet
brain_gene_names <- read_delim("brain_uniprot2genename.txt", delim = "\t", col_names = TRUE) %>% 
  rename(sp.BX_UniProt.ID = From, gene_name = To)
# Attach the Uniprot gene names to the annotation spreadsheet
brain_blastx <- left_join(brain_blastx, brain_gene_names) %>% 
  relocate(gene_name, .before = sp.BX_transcript_name) # Move the gene name column before its description
# Remove the gene name conversion table to keep the global environment clean(er)
rm(brain_gene_names)

# Repeat the process
write_delim(as.data.frame(unique(ant_intestine_blastx$sp.BX_UniProt.ID)), "anterior_intestine_uniprot_IDs.txt", col_names = FALSE, delim = "\t")
# Read in the uniprot gene names with gene IDs from above, for attachment to the annotation spreadsheet
ant_intestine_gene_names <- read_delim("anterior_intestine_uniprot2genename.txt", delim = "\t", col_names = TRUE) %>% 
  rename(sp.BX_UniProt.ID = From, gene_name = To)
# Attach the Uniprot gene names to the annotation spreadsheet
ant_intestine_blastx <- left_join(ant_intestine_blastx, ant_intestine_gene_names) %>% 
  relocate(gene_name, .before = sp.BX_transcript_name) # Move the gene name column before its description
# Remove the gene name conversion table to keep the global environment clean(er)
rm(ant_intestine_gene_names)

# Repeat the process
write_delim(as.data.frame(unique(pyloric_ceca_blastx$sp.BX_UniProt.ID)), "pyloric_ceca_uniprot_IDs.txt", col_names = FALSE, delim = "\t")
# Read in the uniprot gene names with gene IDs from above, for attachment to the annotation spreadsheet
pyloric_ceca_gene_names <- read_delim("pyloric_ceca_uniprot2genename.txt", delim = "\t", col_names = TRUE) %>% 
  rename(sp.BX_UniProt.ID = From, gene_name = To)
# Attach the Uniprot gene names to the annotation spreadsheet
pyloric_ceca_blastx <- left_join(pyloric_ceca_blastx, pyloric_ceca_gene_names) %>% 
  relocate(gene_name, .before = sp.BX_transcript_name) # Move the gene name column before its description
# Remove the gene name conversion table to keep the global environment clean(er)
rm(pyloric_ceca_gene_names)

# Repeat the process
write_delim(as.data.frame(unique(rectum_blastx$sp.BX_UniProt.ID)), "rectum_uniprot_IDs.txt", col_names = FALSE, delim = "\t")
# Read in the uniprot gene names with gene IDs from above, for attachment to the annotation spreadsheet
rectum_gene_names <- read_delim("rectum_uniprot2genename.txt", delim = "\t", col_names = TRUE) %>% 
  rename(sp.BX_UniProt.ID = From, gene_name = To)
# Attach the Uniprot gene names to the annotation spreadsheet
rectum_blastx <- left_join(rectum_blastx, rectum_gene_names) %>% 
  relocate(gene_name, .before = sp.BX_transcript_name) # Move the gene name column before its description
# Remove the gene name conversion table to keep the global environment clean(er)
rm(rectum_gene_names)

# Repeat the process
write_delim(as.data.frame(unique(spiral_valve_blastx$sp.BX_UniProt.ID)), "spiral_valve_uniprot_IDs.txt", col_names = FALSE, delim = "\t")
# Read in the uniprot gene names with gene IDs from above, for attachment to the annotation spreadsheet
spiral_valve_gene_names <- read_delim("spiral_valve_uniprot2genename.txt", delim = "\t", col_names = TRUE) %>% 
  rename(sp.BX_UniProt.ID = From, gene_name = To)
# Attach the Uniprot gene names to the annotation spreadsheet
spiral_valve_blastx <- left_join(spiral_valve_blastx, spiral_valve_gene_names) %>% 
  relocate(gene_name, .before = sp.BX_transcript_name) # Move the gene name column before its description
# Remove the gene name conversion table to keep the global environment clean(er)
rm(spiral_valve_gene_names)

# Repeat the process
write_delim(as.data.frame(unique(stomach_distal_blastx$sp.BX_UniProt.ID)), "stomach_distal_uniprot_IDs.txt", col_names = FALSE, delim = "\t")
# Read in the uniprot gene names with gene IDs from above, for attachment to the annotation spreadsheet
stomach_distal_gene_names <- read_delim("stomach_distal_uniprot2genename.txt", delim = "\t", col_names = TRUE) %>% 
  rename(sp.BX_UniProt.ID = From, gene_name = To)
# Attach the Uniprot gene names to the annotation spreadsheet
stomach_distal_blastx <- left_join(stomach_distal_blastx, stomach_distal_gene_names) %>% 
  relocate(gene_name, .before = sp.BX_transcript_name) # Move the gene name column before its description
# Remove the gene name conversion table to keep the global environment clean(er)
rm(stomach_distal_gene_names)

# Repeat the process
write_delim(as.data.frame(unique(stomach_proximal_blastx$sp.BX_UniProt.ID)), "stomach_proximal_uniprot_IDs.txt", col_names = FALSE, delim = "\t")
# Read in the uniprot gene names with gene IDs from above, for attachment to the annotation spreadsheet
stomach_proximal_gene_names <- read_delim("stomach_proximal_uniprot2genename.txt", delim = "\t", col_names = TRUE) %>% 
  rename(sp.BX_UniProt.ID = From, gene_name = To)
# Attach the Uniprot gene names to the annotation spreadsheet
stomach_proximal_blastx <- left_join(stomach_proximal_blastx, stomach_proximal_gene_names) %>% 
  relocate(gene_name, .before = sp.BX_transcript_name) # Move the gene name column before its description
# Remove the gene name conversion table to keep the global environment clean(er)
rm(stomach_proximal_gene_names)

# Repeat the process
write_delim(as.data.frame(unique(gill_blastx$sp.BX_UniProt.ID)), "gill_uniprot_IDs.txt", col_names = FALSE, delim = "\t")
# Read in the uniprot gene names with gene IDs from above, for attachment to the annotation spreadsheet
gill_gene_names <- read_delim("gill_uniprot2genename.txt", delim = "\t", col_names = TRUE) %>% 
  rename(sp.BX_UniProt.ID = From, gene_name = To)
# Attach the Uniprot gene names to the annotation spreadsheet
gill_blastx <- left_join(gill_blastx, gill_gene_names) %>% 
  relocate(gene_name, .before = sp.BX_transcript_name) # Move the gene name column before its description
# Remove the gene name conversion table to keep the global environment clean(er)
rm(gill_gene_names)

# Repeat the process
write_delim(as.data.frame(unique(head_kidney_blastx$sp.BX_UniProt.ID)), "head_kidney_uniprot_IDs.txt", col_names = FALSE, delim = "\t")
# Read in the uniprot gene names with gene IDs from above, for attachment to the annotation spreadsheet
head_kidney_gene_names <- read_delim("head_kidney_uniprot2genename.txt", delim = "\t", col_names = TRUE) %>% 
  rename(sp.BX_UniProt.ID = From, gene_name = To)
# Attach the Uniprot gene names to the annotation spreadsheet
head_kidney_blastx <- left_join(head_kidney_blastx, head_kidney_gene_names) %>% 
  relocate(gene_name, .before = sp.BX_transcript_name) # Move the gene name column before its description
# Remove the gene name conversion table to keep the global environment clean(er)
rm(head_kidney_gene_names)

# Repeat the process
write_delim(as.data.frame(unique(heart_blastx$sp.BX_UniProt.ID)), "heart_uniprot_IDs.txt", col_names = FALSE, delim = "\t")
# Read in the uniprot gene names with gene IDs from above, for attachment to the annotation spreadsheet
heart_gene_names <- read_delim("heart_uniprot2genename.txt", delim = "\t", col_names = TRUE) %>% 
  rename(sp.BX_UniProt.ID = From, gene_name = To)
# Attach the Uniprot gene names to the annotation spreadsheet
heart_blastx <- left_join(heart_blastx, heart_gene_names) %>% 
  relocate(gene_name, .before = sp.BX_transcript_name) # Move the gene name column before its description
# Remove the gene name conversion table to keep the global environment clean(er)
rm(heart_gene_names)

# Repeat the process
write_delim(as.data.frame(unique(liver_blastx$sp.BX_UniProt.ID)), "liver_uniprot_IDs.txt", col_names = FALSE, delim = "\t")
# Read in the uniprot gene names with gene IDs from above, for attachment to the annotation spreadsheet
liver_gene_names <- read_delim("liver_uniprot2genename.txt", delim = "\t", col_names = TRUE) %>% 
  rename(sp.BX_UniProt.ID = From, gene_name = To)
# Attach the Uniprot gene names to the annotation spreadsheet
liver_blastx <- left_join(liver_blastx, liver_gene_names) %>% 
  relocate(gene_name, .before = sp.BX_transcript_name) # Move the gene name column before its description
# Remove the gene name conversion table to keep the global environment clean(er)
rm(liver_gene_names)

# Repeat the process
write_delim(as.data.frame(unique(white_muscle_blastx$sp.BX_UniProt.ID)), "white_muscle_uniprot_IDs.txt", col_names = FALSE, delim = "\t")
# Read in the uniprot gene names with gene IDs from above, for attachment to the annotation spreadsheet
white_muscle_gene_names <- read_delim("white_muscle_uniprot2genename.txt", delim = "\t", col_names = TRUE) %>% 
  rename(sp.BX_UniProt.ID = From, gene_name = To)
# Attach the Uniprot gene names to the annotation spreadsheet
white_muscle_blastx <- left_join(white_muscle_blastx, white_muscle_gene_names) %>% 
  relocate(gene_name, .before = sp.BX_transcript_name) # Move the gene name column before its description
# Remove the gene name conversion table to keep the global environment clean(er)
rm(white_muscle_gene_names)

# Repeat the process
write_delim(as.data.frame(unique(esophagus_blastx$sp.BX_UniProt.ID)), "esophagus_uniprot_IDs.txt", col_names = FALSE, delim = "\t")
# Read in the uniprot gene names with gene IDs from above, for attachment to the annotation spreadsheet
esophagus_gene_names <- read_delim("esophagus_uniprot2genename.txt", delim = "\t", col_names = TRUE) %>% 
  rename(sp.BX_UniProt.ID = From, gene_name = To)
# Attach the Uniprot gene names to the annotation spreadsheet
esophagus_blastx <- left_join(esophagus_blastx, esophagus_gene_names) %>% 
  relocate(gene_name, .before = sp.BX_transcript_name) # Move the gene name column before its description
# Remove the gene name conversion table to keep the global environment clean(er)
rm(esophagus_gene_names)

# Repeat the process
write_delim(as.data.frame(unique(head_kidney2_blastx$sp.BX_UniProt.ID)), "head_kidney2_uniprot_IDs.txt", col_names = FALSE, delim = "\t")
# Read in the uniprot gene names with gene IDs from above, for attachment to the annotation spreadsheet
head_kidney2_gene_names <- read_delim("head_kidney2_uniprot2genename.txt", delim = "\t", col_names = TRUE) %>% 
  rename(sp.BX_UniProt.ID = From, gene_name = To)
# Attach the Uniprot gene names to the annotation spreadsheet
head_kidney2_blastx <- left_join(head_kidney2_blastx, head_kidney2_gene_names) %>% 
  relocate(gene_name, .before = sp.BX_transcript_name) # Move the gene name column before its description
# Remove the gene name conversion table to keep the global environment clean(er)
rm(head_kidney2_gene_names)

# Repeat the process
write_delim(as.data.frame(unique(overall_blastx$sp.BX_UniProt.ID)), "overall_uniprot_IDs.txt", col_names = FALSE, delim = "\t")
# Read in the uniprot gene names with gene IDs from above, for attachment to the annotation spreadsheet
overall_gene_names <- read_delim("overall_uniprot-download_true_fields_accession_2Creviewed_2Cid_2Cgene_names_-2022.08.29-13.58.55.44.tsv", delim = "\t", col_names = TRUE) %>% 
  rename(sp.BX_UniProt.ID = From, gene_name = `Gene Names (primary)`) %>% 
  mutate(gene_name = toupper(gene_name)) %>% 
  select(sp.BX_UniProt.ID, gene_name)
# Attach the Uniprot gene names to the annotation spreadsheet
overall_blastx <- left_join(overall_blastx, overall_gene_names) %>% 
  relocate(gene_name, .before = sp.BX_transcript_name) # Move the gene name column before its description
# Remove the gene name conversion table to keep the global environment clean(er)
rm(overall_gene_names)

# Write out the different annotation tables for downstream analyses
write_tsv(x = ant_intestine_blastx, file = "anterior_intestine_annotation_bit50_E1e-6.txt")
write_tsv(x = brain_blastx, file = "brain_annotation_bit50_E1e-6.txt")
write_tsv(x = gill_blastx, file = "gill_annotation_bit50_E1e-6.txt")
write_tsv(x = head_kidney_blastx, file = "head_kidney_annotation_bit50_E1e-6.txt")
write_tsv(x = heart_blastx, file = "heart_annotation_bit50_E1e-6.txt")
write_tsv(x = liver_blastx, file = "liver_annotation_bit50_E1e-6.txt")
write_tsv(x = pyloric_ceca_blastx, file = "pyloric_ceca_annotation_bit50_E1e-6.txt")
write_tsv(x = rectum_blastx, file = "rectum_annotation_bit50_E1e-6.txt")
write_tsv(x = spiral_valve_blastx, file = "spiral_valve_annotation_bit50_E1e-6.txt")
write_tsv(x = stomach_distal_blastx, file = "stomach_distal_annotation_bit50_E1e-6.txt")
write_tsv(x = stomach_proximal_blastx, file = "stomach_proximal_annotation_bit50_E1e-6.txt")
write_tsv(x = white_muscle_blastx, file = "white_muscle_annotation_bit50_E1e-6.txt")
write_tsv(x = esophagus_blastx, file = "esophagus_annotation_bit50_E1e-6.txt")

write_tsv(x = head_kidney2_blastx, file = "head_kidney2_annotation_bit50_E1e-6.txt")

write_tsv(x = overall_blastx, file = "overall_annotation_bit50_E1e-6.txt")
