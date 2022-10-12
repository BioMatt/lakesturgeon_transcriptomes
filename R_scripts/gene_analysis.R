# A script to analyze genes from the 13 different tissue-specific transcriptomes for lake sturgeon.

library(tidyverse)
library(enrichR)
library(UpSetR)
library(SuperExactTest)
library(ggrepel)
library(patchwork)
library(fishualize) # Colour palettes based on fish species

# A function to take an EnrichrR table and split the GO term and definition for input into Revigo
split_go <- function(x) {
  # Taking apart the GO descriptions and GO ID terms
  x <- x %>%
    select(-starts_with("Old")) %>%
    separate(Term, c("GO_term", "GO_ID"), sep = "GO")
  
  # Changing all the floating :#### GO terms to be GO:####
  x <- mutate_if(x, is.character, str_replace_all, pattern = ":", replacement = "GO:")
  
  # Removing the last parentheses in the GO ID 
  x$GO_ID <- x$GO_ID %>% 
    str_replace("\\)", "")
  
  # Removing the last parentheses in the GO term to have clean looking cells 
  x$GO_term <- x$GO_term %>% 
    str_replace("\\($", "")
  return(x)
}

# A function to read in a transcriptome report and create a list of the Biological Process 2021, Molecular Function 2021, and Cellular Component 2021 results filtered for q < 0.05
enrichR_tissue <- function(annotation_report){
  # Read in the gene annotations for the annotation report. Capitalize gene names to keep them consistent
  report <- read_tsv(annotation_report) %>% 
    mutate(gene_name = toupper(gene_name))
  # Run enrichR
  enrichr_results <- enrichr(unique(report$gene_name), databases = c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021"))
  
  # Pull the results from one of the databases, biological process first
  bio_process <- as_tibble(enrichr_results[["GO_Biological_Process_2021"]])
  bio_process <- dplyr::filter(bio_process, Adjusted.P.value < 0.05)
  # Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
  bio_process <- split_go(bio_process)
  
  # Pull the results from one of the databases, now molecular function
  mol_function <- as_tibble(enrichr_results[["GO_Molecular_Function_2021"]])
  mol_function <- dplyr::filter(mol_function, Adjusted.P.value < 0.05)
  # Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
  mol_function <- split_go(mol_function)
  
  # Pull the results from one of the databases, last cellular component
  cell_comp <- as_tibble(enrichr_results[["GO_Cellular_Component_2021"]])
  cell_comp <- dplyr::filter(cell_comp, Adjusted.P.value < 0.05)
  # Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
  cell_comp <- split_go(cell_comp)
  
  tissue_list <- list(bio = bio_process, mol = mol_function, cell = cell_comp)
  return(tissue_list)
}

# A simpler function to read the annotation report and capitalize gene names.
# Filter out NA values for gene names, only keep distinct gene names and only return that column to save on memory.
read_annotation <- function(annotation_report){
  # Read in the gene annotations for the annotation report. Capitalize gene names to keep them consistent
  report <- read_tsv(annotation_report) %>% 
    mutate(gene_name = toupper(gene_name)) %>% 
    filter(!is.na(gene_name)) %>% 
    distinct(gene_name)
  return(report)
}

# Run the custom function to find GO terms for all 13 tissues
ant_intestine_GO <- enrichR_tissue("anterior_intestine_annotation_bit50_E1e-6.txt")
brain_GO <- enrichR_tissue("brain_annotation_bit50_E1e-6.txt")
esophagus_GO <- enrichR_tissue("esophagus_annotation_bit50_E1e-6.txt")
gill_GO <- enrichR_tissue("gill_annotation_bit50_E1e-6.txt")
#head_kidney_GO <- enrichR_tissue("head_kidney_annotation_bit50_E1e-6.txt")
heart_GO <- enrichR_tissue("heart_annotation_bit50_E1e-6.txt")
liver_GO <- enrichR_tissue("liver_annotation_bit50_E1e-6.txt")
pyloric_ceca_GO <- enrichR_tissue("pyloric_ceca_annotation_bit50_E1e-6.txt")
rectum_GO <- enrichR_tissue("rectum_annotation_bit50_E1e-6.txt")
spiral_valve_GO <- enrichR_tissue("spiral_valve_annotation_bit50_E1e-6.txt")
stomach_distal_GO <- enrichR_tissue("stomach_distal_annotation_bit50_E1e-6.txt")
stomach_proximal_GO <- enrichR_tissue("stomach_proximal_annotation_bit50_E1e-6.txt")
white_muscle_GO <- enrichR_tissue("white_muscle_annotation_bit50_E1e-6.txt")

head_kidney2_GO <- enrichR_tissue("head_kidney2_annotation_bit50_E1e-6.txt")

# For each of the annotations, count the number of GO terms from each category
print(nrow(ant_intestine_GO$bio)) # Biological process
print(nrow(ant_intestine_GO$mol)) # Molecular function
print(nrow(ant_intestine_GO$cell)) # Cellular component

print(nrow(brain_GO$bio)) # Biological process
print(nrow(brain_GO$mol)) # Molecular function
print(nrow(brain_GO$cell)) # Cellular component

print(nrow(esophagus_GO$bio)) # Biological process
print(nrow(esophagus_GO$mol)) # Molecular function
print(nrow(esophagus_GO$cell)) # Cellular component

print(nrow(gill_GO$bio)) # Biological process
print(nrow(gill_GO$mol)) # Molecular function
print(nrow(gill_GO$cell)) # Cellular component

print(nrow(head_kidney_GO$bio)) # Biological process
print(nrow(head_kidney_GO$mol)) # Molecular function
print(nrow(head_kidney_GO$cell)) # Cellular component

print(nrow(heart_GO$bio)) # Biological process
print(nrow(heart_GO$mol)) # Molecular function
print(nrow(heart_GO$cell)) # Cellular component

print(nrow(liver_GO$bio)) # Biological process
print(nrow(liver_GO$mol)) # Molecular function
print(nrow(liver_GO$cell)) # Cellular component

print(nrow(pyloric_ceca_GO$bio)) # Biological process
print(nrow(pyloric_ceca_GO$mol)) # Molecular function
print(nrow(pyloric_ceca_GO$cell)) # Cellular component

print(nrow(rectum_GO$bio)) # Biological process
print(nrow(rectum_GO$mol)) # Molecular function
print(nrow(rectum_GO$cell)) # Cellular component

print(nrow(spiral_valve_GO$bio)) # Biological process
print(nrow(spiral_valve_GO$mol)) # Molecular function
print(nrow(spiral_valve_GO$cell)) # Cellular component

print(nrow(stomach_distal_GO$bio)) # Biological process
print(nrow(stomach_distal_GO$mol)) # Molecular function
print(nrow(stomach_distal_GO$cell)) # Cellular component

print(nrow(stomach_proximal_GO$bio)) # Biological process
print(nrow(stomach_proximal_GO$mol)) # Molecular function
print(nrow(stomach_proximal_GO$cell)) # Cellular component

print(nrow(white_muscle_GO$bio)) # Biological process
print(nrow(white_muscle_GO$mol)) # Molecular function
print(nrow(white_muscle_GO$cell)) # Cellular component

print(nrow(head_kidney2_GO$bio)) # Biological process
print(nrow(head_kidney2_GO$mol)) # Molecular function
print(nrow(head_kidney2_GO$cell)) # Cellular component

# Trying an upsetR plot to look at overlap in genes among the different data sets
# First, the biological process GO terms
bio_GO_list <- list(ant_intestine_GO$bio$GO_ID, brain_GO$bio$GO_ID, esophagus_GO$bio$GO_ID, gill_GO$bio$GO_ID, head_kidney2_GO$bio$GO_ID, heart_GO$bio$GO_ID, liver_GO$bio$GO_ID, pyloric_ceca_GO$bio$GO_ID, rectum_GO$bio$GO_ID, spiral_valve_GO$bio$GO_ID, stomach_distal_GO$bio$GO_ID, stomach_proximal_GO$bio$GO_ID, white_muscle_GO$bio$GO_ID)
bio_upset_data <- reshape2::melt(bio_GO_list) %>%
  magrittr::set_colnames(c("dataset", "GO_term")) %>%
  dplyr::group_by(dataset, GO_term) %>%
  dplyr::mutate(value = 1)  %>%
  tidyr::spread(GO_term, value, fill = 0) %>% 
  rename(GO_term = 1, `Anterior Intestine` = 2, `Brain` = 3, `Esophagus` = 4, `Gill` = 5, `Head Kidney` = 6, `Heart` = 7, `Liver` = 8, `Pyloric Caecum` = 9, `Rectum` = 10, `Spiral Valve` = 11, `Muscular Stomach` = 12, `Glandular Stomach` = 13, `White Muscle` = 14) 
bio_upset_data <- as.data.frame(bio_upset_data)
for(i in 2:ncol(bio_upset_data)){ bio_upset_data[ , i] <- as.integer(bio_upset_data[ , i]) }
# View the upset plot
upset(bio_upset_data, nsets = ncol(bio_upset_data) - 1, order.by = "freq", sets.x.label = "Biological Process 2021 Set Size")
# Which bio process terms are shared among all tissues?
View(filter(bio_upset_data, across(is.numeric, ~ .x == 1)))


# The molecular function GO terms
mol_GO_list <- list(ant_intestine_GO$mol$GO_ID, brain_GO$mol$GO_ID, esophagus_GO$mol$GO_ID, gill_GO$mol$GO_ID, head_kidney2_GO$mol$GO_ID, heart_GO$mol$GO_ID, liver_GO$mol$GO_ID, pyloric_ceca_GO$mol$GO_ID, rectum_GO$mol$GO_ID, spiral_valve_GO$mol$GO_ID, stomach_distal_GO$mol$GO_ID, stomach_proximal_GO$mol$GO_ID, white_muscle_GO$mol$GO_ID)
mol_upset_data <- reshape2::melt(mol_GO_list) %>%
  magrittr::set_colnames(c("dataset", "GO_term")) %>%
  dplyr::group_by(dataset, GO_term) %>%
  dplyr::mutate(value = 1)  %>%
  tidyr::spread(GO_term, value, fill = 0) %>% 
  rename(GO_term = 1, `Anterior Intestine` = 2, `Brain` = 3, `Esophagus` = 4, `Gill` = 5, `Head Kidney` = 6, `Heart` = 7, `Liver` = 8, `Pyloric Caecum` = 9, `Rectum` = 10, `Spiral Valve` = 11, `Muscular Stomach` = 12, `Glandular Stomach` = 13, `White Muscle` = 14)
mol_upset_data <- as.data.frame(mol_upset_data)
for(i in 2:ncol(mol_upset_data)){ mol_upset_data[ , i] <- as.integer(mol_upset_data[ , i]) }
# View the upset plot
upset(mol_upset_data, nsets = ncol(mol_upset_data) - 1, order.by = "freq", sets.x.label = "Molecular Function 2021 Set Size")

# The cellular component GO terms
cell_GO_list <- list(ant_intestine_GO$cell$GO_ID, brain_GO$cell$GO_ID, esophagus_GO$cell$GO_ID, gill_GO$cell$GO_ID, head_kidney2_GO$cell$GO_ID, heart_GO$cell$GO_ID, liver_GO$cell$GO_ID, pyloric_ceca_GO$cell$GO_ID, rectum_GO$cell$GO_ID, spiral_valve_GO$cell$GO_ID, stomach_distal_GO$cell$GO_ID, stomach_proximal_GO$cell$GO_ID, white_muscle_GO$cell$GO_ID)
cell_upset_data <- reshape2::melt(cell_GO_list) %>%
  magrittr::set_colnames(c("dataset", "GO_term")) %>%
  dplyr::group_by(dataset, GO_term) %>%
  dplyr::mutate(value = 1)  %>%
  tidyr::spread(GO_term, value, fill = 0) %>% 
  rename(GO_term = 1, `Anterior Intestine` = 2, `Brain` = 3, `Esophagus` = 4, `Gill` = 5, `Head Kidney` = 6, `Heart` = 7, `Liver` = 8, `Pyloric Caecum` = 9, `Rectum` = 10, `Spiral Valve` = 11, `Muscular Stomach` = 12, `Glandular Stomach` = 13, `White Muscle` = 14)
cell_upset_data <- as.data.frame(cell_upset_data)
for(i in 2:ncol(cell_upset_data)){ cell_upset_data[ , i] <- as.integer(cell_upset_data[ , i]) }
# View the upset plot
upset(cell_upset_data, nsets = ncol(cell_upset_data) - 1, order.by = "freq", sets.x.label = "Cellular Component 2021 Set Size")


# Combine all three upset dataframes into a list
upset_list <- list(bio = bio_upset_data, mol = mol_upset_data, cell = cell_upset_data)

# A function to attach tissue-specificity information to each GO term, from enrichR results
tissue_GO <- function(tissue_data, tissue, UpSet_data){
  # Start with biological process results
  bio_GO <- tissue_data$bio # Pull the GO database results to be worked on
  bio_GO$unique <- NA # Add an empty column for uniqueness for the particular tissue
  bio_GO$database <- rep("Biological_Process_2021", nrow(bio_GO)) # Add information for the database used
  bio_GO <- left_join(bio_GO, UpSet_data$bio, by = c("GO_ID" = "GO_term")) # Attach UpSetR results to the GO results
  tissue_index <- match(tissue, names(bio_GO)) # Take the index from the UpSetR results where the tissue data is
  # Identify GO terms that are *only* present in the tissue of focus and fill in 1 in the 'unique' column if it is unique. If not unique, fill in 0
  bio_GO <- bio_GO %>% 
    mutate(unique = dplyr::if_else(bio_GO[,tissue_index] == 1 & rowSums(bio_GO[,-c(1:10, tissue_index)]) == 0, 1, 0))
  
  # Repeat the process for molecular function GO terms
  mol_GO <- tissue_data$mol # Pull the GO database results to be worked on
  mol_GO$unique <- NA # Add an empty column for uniqueness for the particular tissue
  mol_GO$database <- rep("Molecular_Function_2021", nrow(mol_GO)) # Add information for the database used
  mol_GO <- left_join(mol_GO, UpSet_data$mol, by = c("GO_ID" = "GO_term")) # Attach UpSetR results to the GO results
  tissue_index <- match(tissue, names(mol_GO)) # Take the index from the UpSetR results where the tissue data is
  # Identify GO terms that are *only* present in the tissue of focus and fill in 1 in the 'unique' column if it is unique. If not unique, fill in 0
  mol_GO <- mol_GO %>% 
    mutate(unique = dplyr::if_else(mol_GO[,tissue_index] == 1 & rowSums(mol_GO[,-c(1:10, tissue_index)]) == 0, 1, 0))
  
  # Repeat the process for cellular component GO terms
  cell_GO <- tissue_data$cell # Pull the GO database results to be worked on
  cell_GO$unique <- NA # Add an empty column for uniqueness for the particular tissue
  cell_GO$database <- rep("Cellular_Component_2021", nrow(cell_GO)) # Add information for the database used
  cell_GO <- left_join(cell_GO, UpSet_data$cell, by = c("GO_ID" = "GO_term")) # Attach UpSetR results to the GO results
  tissue_index <- match(tissue, names(cell_GO)) # Take the index from the UpSetR results where the tissue data is
  # Identify GO terms that are *only* present in the tissue of focus and fill in 1 in the 'unique' column if it is unique. If not unique, fill in 0
  cell_GO <- cell_GO %>% 
    mutate(unique = dplyr::if_else(cell_GO[,tissue_index] == 1 & rowSums(cell_GO[,-c(1:10, tissue_index)]) == 0, 1, 0))
  
  # Combine all three new dataframes with uniqueness into a list, and export
  return(rbind(bio_GO, mol_GO, cell_GO))
}

# Run the custom function for GO term uniqueness for each tissue
brain_unique <- tissue_GO(brain_GO, "brain", upset_list)
esophagus_unique <- tissue_GO(esophagus_GO, "esophagus", upset_list)
gill_unique <- tissue_GO(gill_GO, "gill", upset_list)
head_kidney_unique <- tissue_GO(head_kidney2_GO, "head_kidney", upset_list)
heart_unique <- tissue_GO(heart_GO, "heart", upset_list)
liver_unique <- tissue_GO(liver_GO, "liver", upset_list)
ant_intestine_unique <- tissue_GO(ant_intestine_GO, "ant_intestine", upset_list)
pyloric_ceca_unique <- tissue_GO(pyloric_ceca_GO, "pyloric_ceca", upset_list)
rectum_unique <- tissue_GO(rectum_GO, "rectum", upset_list)
spiral_valve_unique <- tissue_GO(spiral_valve_GO, "spiral_valve", upset_list)
stomach_distal_unique <- tissue_GO(stomach_distal_GO, "stomach_distal", upset_list)
stomach_proximal_unique <- tissue_GO(stomach_proximal_GO, "stomach_proximal", upset_list)
white_muscle_unique <- tissue_GO(white_muscle_GO, "white_muscle", upset_list)



# Write out the different GO term uniqueness results
write_csv(x = ant_intestine_unique, file = "ls_transcriptome_GOunique/ant_intestine_GO_uniquness.csv")

write_csv(x = brain_unique, file = "ls_transcriptome_GOunique/brain_GO_uniquness.csv")

write_csv(x = esophagus_unique, file = "ls_transcriptome_GOunique/esophagus_GO_uniquness.csv")

write_csv(x = gill_unique, file = "ls_transcriptome_GOunique/gill_GO_uniquness.csv")

write_csv(x = head_kidney_unique, file = "ls_transcriptome_GOunique/head_kidney_GO_uniquness.csv")

write_csv(x = heart_unique, file = "ls_transcriptome_GOunique/heart_GO_uniquness.csv")

write_csv(x = liver_unique, file = "ls_transcriptome_GOunique/liver_GO_uniquness.csv")

write_csv(x = pyloric_ceca_unique, file = "ls_transcriptome_GOunique/pyloric_ceca_GO_uniquness.csv")

write_csv(x = rectum_unique, file = "ls_transcriptome_GOunique/rectum_GO_uniquness.csv")

write_csv(x = spiral_valve_unique, file = "ls_transcriptome_GOunique/spiral_valve_GO_uniquness.csv")

write_csv(x = stomach_distal_unique, file = "ls_transcriptome_GOunique/stomach_distal_GO_uniquness.csv")

write_csv(x = stomach_proximal_unique, file = "ls_transcriptome_GOunique/stomach_proximal_GO_uniquness.csv")

write_csv(x = white_muscle_unique, file = "ls_transcriptome_GOunique/white_muscle_GO_uniquness.csv")

########################################################################################################
# Plot the unique GO terms for each tissue
# This is a custom function to take only unique GO terms and plot them as a sideways barplot
plot_Unique_GO <- function(input_GO, tissue){
  input_GO$database <- factor(input_GO$database, levels = c("Biological_Process_2021", "Molecular_Function_2021", "Cellular_Component_2021"))
  filtered_GO <- dplyr::filter(input_GO, unique == 1) %>% 
    arrange(desc(database), Combined.Score)
  filtered_GO$GO_ID <- factor(filtered_GO$GO_ID, levels = filtered_GO$GO_ID)
  filtered_GO$GO_term <- factor(filtered_GO$GO_term, levels = filtered_GO$GO_term)
  if (nrow(filtered_GO) == 0){
    stop("No unique GO terms for tissue!")
  }
  plot <- ggplot(data = filtered_GO, aes(x = GO_ID, y = log(Combined.Score), colour = database, fill = database, group = database)) +
    geom_col() +
    scale_x_discrete(labels = filtered_GO$GO_term) +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    theme_classic() +
    xlab(element_blank()) +
    ylab(bquote('log'[10]~'Combined Score')) +
    scale_colour_manual(name = "Enrichment Database", labels = c("Biological_Process_2021" = "Biological Process 2021", "Molecular_Function_2021" = "Molecular Function 2021", "Cellular_Component_2021" = "Cellular Component 2021"), values = c("#4575b4", "#fdae61", "#d73027"), aesthetics = c("colour", "fill"), breaks = c("Biological_Process_2021", "Molecular_Function_2021", "Cellular_Component_2021")) +
    theme(text=element_text(size=16), axis.ticks = element_blank()) +
    coord_flip() +
    ggtitle(tissue)
  return(plot)
}

ggsave(filename = "ls_transcriptome_GOunique/heart_uniqueGO.png", plot = plot_Unique_GO(heart_unique, "Heart"), dpi = 900)
ggsave(filename = "ls_transcriptome_GOunique/ant_intestine_uniqueGO.png", plot = plot_Unique_GO(ant_intestine_unique, "Anterior Intestine"), dpi = 900)
ggsave(filename = "ls_transcriptome_GOunique/brain_uniqueGO.png", plot = plot_Unique_GO(brain_unique, "Brain"), dpi = 900, height = 25, width = 25)
ggsave(filename = "ls_transcriptome_GOunique/esophagus_uniqueGO.png", plot = plot_Unique_GO(esophagus_unique, "Esophagus"), dpi = 900)
ggsave(filename = "ls_transcriptome_GOunique/gill_uniqueGO.png", plot = plot_Unique_GO(gill_unique, "Gill"), dpi = 900, height = 10, width = 15)
ggsave(filename = "ls_transcriptome_GOunique/head_kidney_uniqueGO.png", plot = plot_Unique_GO(head_kidney_unique, "Head Kidney"), dpi = 900)
ggsave(filename = "ls_transcriptome_GOunique/heart_uniqueGO.png", plot = plot_Unique_GO(heart_unique, "Heart"), dpi = 900, height = 25, width = 25)
ggsave(filename = "ls_transcriptome_GOunique/liver_uniqueGO.png", plot = plot_Unique_GO(liver_unique, "Liver"), dpi = 900, height = 25, width = 25)
ggsave(filename = "ls_transcriptome_GOunique/pyloric_ceca_uniqueGO.png", plot = plot_Unique_GO(pyloric_ceca_unique, "Pyloric Ceca"), dpi = 900)
ggsave(filename = "ls_transcriptome_GOunique/rectum_uniqueGO.png", plot = plot_Unique_GO(rectum_unique, "Rectum"), dpi = 900)
ggsave(filename = "ls_transcriptome_GOunique/spiral_valve_uniqueGO.png", plot = plot_Unique_GO(spiral_valve_unique, "Spiral Valve"), dpi = 900)
ggsave(filename = "ls_transcriptome_GOunique/stomach_distal_uniqueGO.png", plot = plot_Unique_GO(stomach_distal_unique, "Distal Stomach"), dpi = 900)
ggsave(filename = "ls_transcriptome_GOunique/stomach_proximal_uniqueGO.png", plot = plot_Unique_GO(stomach_proximal_unique, "Proximal Stomach"), dpi = 900)
ggsave(filename = "ls_transcriptome_GOunique/white_muscle_uniqueGO.png", plot = plot_Unique_GO(white_muscle_unique, "White Muscle"), dpi = 900, height = 25, width = 25)


###########################################################################################################
# Read in the gene annotations for all 13 transcriptomes. Capitalize gene names to keep them consistent
anterior_intestine <- read_tsv("anterior_intestine_annotation_bit50_E1e-6.txt") %>% 
  mutate(gene_name = toupper(gene_name)) 
nrow(filter(anterior_intestine, !is.na(gene_name))) # Number of transcripts with annotations
length(unique(anterior_intestine$gene_name)) # Number of unique genes annotated in the transcriptome

brain <- read_tsv("brain_annotation_bit50_E1e-6.txt") %>% 
  mutate(gene_name = toupper(gene_name))
nrow(filter(brain, !is.na(gene_name))) # Number of transcripts with annotations
length(unique(brain$gene_name)) # Number of unique genes annotated in the transcriptome

esophagus <- read_tsv("esophagus_annotation_bit50_E1e-6.txt") %>% 
  mutate(gene_name = toupper(gene_name))
nrow(filter(esophagus, !is.na(gene_name))) # Number of transcripts with annotations
length(unique(esophagus$gene_name)) # Number of unique genes annotated in the transcriptome

gill <- read_tsv("gill_annotation_bit50_E1e-6.txt") %>% 
  mutate(gene_name = toupper(gene_name))
nrow(filter(gill, !is.na(gene_name))) # Number of transcripts with annotations
length(unique(gill$gene_name)) # Number of unique genes annotated in the transcriptome

head_kidney <- read_tsv("head_kidney_annotation_bit50_E1e-6.txt") %>% 
  mutate(gene_name = toupper(gene_name))
nrow(filter(head_kidney, !is.na(gene_name))) # Number of transcripts with annotations
length(unique(head_kidney$gene_name)) # Number of unique genes annotated in the transcriptome

heart <- read_tsv("heart_annotation_bit50_E1e-6.txt") %>% 
  mutate(gene_name = toupper(gene_name))
nrow(filter(heart, !is.na(gene_name))) # Number of transcripts with annotations
length(unique(heart$gene_name)) # Number of unique genes annotated in the transcriptome

liver <- read_tsv("liver_annotation_bit50_E1e-6.txt") %>% 
  mutate(gene_name = toupper(gene_name))
nrow(filter(liver, !is.na(gene_name))) # Number of transcripts with annotations
length(unique(liver$gene_name)) # Number of unique genes annotated in the transcriptome

pyloric_ceca <- read_tsv("pyloric_ceca_annotation_bit50_E1e-6.txt") %>% 
  mutate(gene_name = toupper(gene_name))
nrow(filter(pyloric_ceca, !is.na(gene_name))) # Number of transcripts with annotations
length(unique(pyloric_ceca$gene_name)) # Number of unique genes annotated in the transcriptome

rectum <- read_tsv("rectum_annotation_bit50_E1e-6.txt") %>% 
  mutate(gene_name = toupper(gene_name))
nrow(filter(rectum, !is.na(gene_name))) # Number of transcripts with annotations
length(unique(rectum$gene_name)) # Number of unique genes annotated in the transcriptome

spiral_valve <- read_tsv("spiral_valve_annotation_bit50_E1e-6.txt") %>% 
  mutate(gene_name = toupper(gene_name))
nrow(filter(spiral_valve, !is.na(gene_name))) # Number of transcripts with annotations
length(unique(spiral_valve$gene_name)) # Number of unique genes annotated in the transcriptome

stomach_distal <- read_tsv("stomach_distal_annotation_bit50_E1e-6.txt") %>% 
  mutate(gene_name = toupper(gene_name))
nrow(filter(stomach_distal, !is.na(gene_name))) # Number of transcripts with annotations
length(unique(stomach_distal$gene_name)) # Number of unique genes annotated in the transcriptome

stomach_proximal <- read_tsv("stomach_proximal_annotation_bit50_E1e-6.txt") %>% 
  mutate(gene_name = toupper(gene_name))
nrow(filter(stomach_proximal, !is.na(gene_name))) # Number of transcripts with annotations
length(unique(stomach_proximal$gene_name)) # Number of unique genes annotated in the transcriptome

white_muscle <- read_tsv("white_muscle_annotation_bit50_E1e-6.txt") %>% 
  mutate(gene_name = toupper(gene_name))
nrow(filter(white_muscle, !is.na(gene_name))) # Number of transcripts with annotations
length(unique(white_muscle$gene_name)) # Number of unique genes annotated in the transcriptome

head_kidney2 <- read_tsv("head_kidney2_annotation_bit50_E1e-6.txt") %>% 
  mutate(gene_name = toupper(gene_name))
nrow(filter(head_kidney2, !is.na(gene_name))) # Number of transcripts with annotations
length(unique(head_kidney2$gene_name)) # Number of unique genes annotated in the transcriptome

# Take all of the annotations and collect a table of the number of total transcripts annotated and unique genes annotated to those transcripts

#############################################################################################################
# Gut-specific analyses. 
# First, an upsetR plot to look at overlap in genes among the different data sets in just the gut
# First, the biological process GO terms
gut_bio_GO_list <- list(ant_intestine_GO$bio$GO_ID, esophagus_GO$bio$GO_ID, pyloric_ceca_GO$bio$GO_ID, rectum_GO$bio$GO_ID, spiral_valve_GO$bio$GO_ID, stomach_distal_GO$bio$GO_ID, stomach_proximal_GO$bio$GO_ID)
gut_bio_upset_data <- reshape2::melt(gut_bio_GO_list) %>%
  magrittr::set_colnames(c("dataset", "GO_term")) %>%
  dplyr::group_by(dataset, GO_term) %>%
  dplyr::mutate(value = 1)  %>%
  tidyr::spread(GO_term, value, fill = 0) %>% 
  rename(GO_term = 1, `Anterior Intestine` = 2, `Esophagus` = 3, `Pyloric Caecum` = 4, `Rectum` = 5, `Spiral Valve` = 6, `Muscular Stomach` = 7, `Glandular Stomach` = 8) 
gut_bio_upset_data <- as.data.frame(gut_bio_upset_data)
for(i in 2:ncol(gut_bio_upset_data)){ gut_bio_upset_data[ , i] <- as.integer(gut_bio_upset_data[ , i]) }
# View the upset plot
upset(gut_bio_upset_data, nsets = ncol(gut_bio_upset_data) - 1, order.by = "freq", sets.x.label = "Gut Biological Process 2021 Set Size")

# The molecular function GO terms
gut_mol_GO_list <- list(ant_intestine_GO$mol$GO_ID, esophagus_GO$mol$GO_ID, pyloric_ceca_GO$mol$GO_ID, rectum_GO$mol$GO_ID, spiral_valve_GO$mol$GO_ID, stomach_distal_GO$mol$GO_ID, stomach_proximal_GO$mol$GO_ID)
gut_mol_upset_data <- reshape2::melt(gut_mol_GO_list) %>%
  magrittr::set_colnames(c("dataset", "GO_term")) %>%
  dplyr::group_by(dataset, GO_term) %>%
  dplyr::mutate(value = 1)  %>%
  tidyr::spread(GO_term, value, fill = 0) %>% 
  rename(GO_term = 1, `Anterior Intestine` = 2, `Esophagus` = 3, `Pyloric Caecum` = 4, `Rectum` = 5, `Spiral Valve` = 6, `Muscular Stomach` = 7, `Glandular Stomach` = 8) 
gut_mol_upset_data <- as.data.frame(gut_mol_upset_data)
for(i in 2:ncol(gut_mol_upset_data)){ gut_mol_upset_data[ , i] <- as.integer(gut_mol_upset_data[ , i]) }
# View the upset plot
upset(gut_mol_upset_data, nsets = ncol(gut_mol_upset_data) - 1, order.by = "freq", sets.x.label = "Gut Molecular Function 2021 Set Size")

# The cellular component GO terms
gut_cell_GO_list <- list(ant_intestine_GO$cell$GO_ID, esophagus_GO$cell$GO_ID, pyloric_ceca_GO$cell$GO_ID, rectum_GO$cell$GO_ID, spiral_valve_GO$cell$GO_ID, stomach_distal_GO$cell$GO_ID, stomach_proximal_GO$cell$GO_ID)
gut_cell_upset_data <- reshape2::melt(gut_cell_GO_list) %>%
  magrittr::set_colnames(c("dataset", "GO_term")) %>%
  dplyr::group_by(dataset, GO_term) %>%
  dplyr::mutate(value = 1)  %>%
  tidyr::spread(GO_term, value, fill = 0) %>% 
  rename(GO_term = 1, `Anterior Intestine` = 2, `Esophagus` = 3, `Pyloric Caecum` = 4, `Rectum` = 5, `Spiral Valve` = 6, `Muscular Stomach` = 7, `Glandular Stomach` = 8) 
gut_cell_upset_data <- as.data.frame(gut_cell_upset_data)
for(i in 2:ncol(gut_cell_upset_data)){ gut_cell_upset_data[ , i] <- as.integer(gut_cell_upset_data[ , i]) }
# View the upset plot
upset(gut_cell_upset_data, nsets = ncol(gut_cell_upset_data) - 1, order.by = "freq", sets.x.label = "Gut Cellular Component 2021 Set Size")

# Combine all three upset gut dataframes into a list
gut_upset_list <- list(bio = gut_bio_upset_data, mol = gut_mol_upset_data, cell = gut_cell_upset_data)


# A function to attach tissue-specificity information to each GO term, from enrichR results. This one is modified to only identify if a tissue within the gut is unique to other tissues in the gut, or shared. 
gut_GO <- function(tissue_data, tissue, UpSet_data){
  # Start with biological process results
  bio_GO <- tissue_data$bio # Pull the GO database results to be worked on
  bio_GO$unique <- NA # Add an empty column for uniqueness for the particular tissue
  bio_GO$n_tissues <- NA # An empty column to represent the number of tissues a GO term is present in. Useful if it is not unique
  bio_GO$database <- rep("Biological_Process_2021", nrow(bio_GO)) # Add information for the database used
  bio_GO <- left_join(bio_GO, UpSet_data$bio, by = c("GO_ID" = "GO_term")) # Attach UpSetR results to the GO results
  bio_tissue_index <- as.integer(match(tissue, names(bio_GO))) # Take the index from the UpSetR results where the tissue data is
  # Identify GO terms that are *only* present in the tissue of focus and fill in 1 in the 'unique' column if it is unique. If not unique, fill in 0
  bio_GO <- bio_GO %>% 
    mutate(unique = dplyr::if_else(bio_GO[,bio_tissue_index] == 1 & rowSums(bio_GO[,-c(1:11, bio_tissue_index)]) == 0, 1, 0)) %>% 
    mutate(n_tissues = rowSums(bio_GO[,-c(1:11)]))
  
  # Repeat the process for molecular function GO terms
  mol_GO <- tissue_data$mol # Pull the GO database results to be worked on
  mol_GO$unique <- NA # Add an empty column for uniqueness for the particular tissue
  mol_GO$n_tissues <- NA # An empty column to represent the number of tissues a GO term is present in. Useful if it is not unique
  mol_GO$database <- rep("Molecular_Function_2021", nrow(mol_GO)) # Add information for the database used
  mol_GO <- left_join(mol_GO, UpSet_data$mol, by = c("GO_ID" = "GO_term")) # Attach UpSetR results to the GO results
  mol_tissue_index <- as.integer(match(tissue, names(mol_GO))) # Take the index from the UpSetR results where the tissue data is
  # Identify GO terms that are *only* present in the tissue of focus and fill in 1 in the 'unique' column if it is unique. If not unique, fill in 0
  mol_GO <- mol_GO %>% 
    mutate(unique = dplyr::if_else(mol_GO[,mol_tissue_index] == 1 & rowSums(mol_GO[,-c(1:11, mol_tissue_index)]) == 0, 1, 0)) %>% 
    mutate(n_tissues = rowSums(mol_GO[,-c(1:11)]))

  
  # Repeat the process for cellular component GO terms
  cell_GO <- tissue_data$cell # Pull the GO database results to be worked on
  cell_GO$unique <- NA # Add an empty column for uniqueness for the particular tissue
  cell_GO$n_tissues <- NA # An empty column to represent the number of tissues a GO term is present in. Useful if it is not unique
  cell_GO$database <- rep("Cellular_Component_2021", nrow(cell_GO)) # Add information for the database used
  cell_GO <- left_join(cell_GO, UpSet_data$cell, by = c("GO_ID" = "GO_term")) # Attach UpSetR results to the GO results
  cell_tissue_index <- as.integer(match(tissue, names(cell_GO))) # Take the index from the UpSetR results where the tissue data is
  # Identify GO terms that are *only* present in the tissue of focus and fill in 1 in the 'unique' column if it is unique. If not unique, fill in 0
  cell_GO <- cell_GO %>% 
    mutate(unique = dplyr::if_else(cell_GO[,cell_tissue_index] == 1 & rowSums(cell_GO[,-c(1:11, cell_tissue_index)]) == 0, 1, 0)) %>% 
    mutate(n_tissues = rowSums(cell_GO[,-c(1:11)]))
  
  # Combine all three new dataframes with uniqueness into a list, and export
  return(rbind(bio_GO, mol_GO, cell_GO))
}

# Use the new gut-specific function on all the gut tissues
ant_intestine_gut_unique <- gut_GO(ant_intestine_GO, "ant_intestine", gut_upset_list)
esophagus_gut_unique <- gut_GO(brain_GO, "esophagus", gut_upset_list)
pyloric_ceca_gut_unique <- gut_GO(pyloric_ceca_GO, "pyloric_ceca", gut_upset_list)
rectum_gut_unique <- gut_GO(rectum_GO, "rectum", gut_upset_list)
spiral_valve_gut_unique <- gut_GO(spiral_valve_GO, "spiral_valve", gut_upset_list)
stomach_distal_gut_unique <- gut_GO(stomach_distal_GO, "stomach_distal", gut_upset_list)
stomach_proximal_gut_unique <- gut_GO(stomach_proximal_GO, "stomach_proximal", gut_upset_list)

# Write out the uniqueness results
write_csv(x = ant_intestine_gut_unique, file = "digestion_GOunique/ant_intestine_GO_gut_uniquness.csv")
write_csv(x = esophagus_gut_unique, file = "digestion_GOunique/esophagus_GO_gut_uniquness.csv")
write_csv(x = pyloric_ceca_gut_unique, file = "digestion_GOunique/pyloric_ceca_GO_gut_uniquness.csv")
write_csv(x = rectum_gut_unique, file = "digestion_GOunique/rectum_GO_gut_uniquness.csv")
write_csv(x = spiral_valve_gut_unique, file = "digestion_GOunique/spiral_valve_GO_gut_uniquness.csv")
write_csv(x = stomach_distal_gut_unique, file = "digestion_GOunique/stomach_distal_GO_gut_uniquness.csv")
write_csv(x = stomach_proximal_gut_unique, file = "digestion_GOunique/stomach_proximal_GO_gut_uniquness.csv")

# Use the plotting function to get gut-specific plots
ggsave(filename = "digestion_GOunique/ant_intestine_uniqueGO.png", plot = plot_Unique_GO(ant_intestine_gut_unique, "Anterior Intestine"), dpi = 900, height = 25, width = 25)
ggsave(filename = "digestion_GOunique/esophagus_uniqueGO.png", plot = plot_Unique_GO(esophagus_gut_unique, "Esophagus"), dpi = 900)
ggsave(filename = "digestion_GOunique/pyloric_ceca_uniqueGO.png", plot = plot_Unique_GO(pyloric_ceca_gut_unique, "Pyloric Ceca"), dpi = 900, height = 20, width = 20)
ggsave(filename = "digestion_GOunique/rectum_uniqueGO.png", plot = plot_Unique_GO(rectum_gut_unique, "Rectum"), dpi = 900)
ggsave(filename = "digestion_GOunique/spiral_valve_uniqueGO.png", plot = plot_Unique_GO(spiral_valve_gut_unique, "Spiral Valve"), dpi = 900)
ggsave(filename = "digestion_GOunique/stomach_distal_uniqueGO.png", plot = plot_Unique_GO(stomach_distal_gut_unique, "Distal Stomach"), dpi = 900)
ggsave(filename = "digestion_GOunique/stomach_proximal_uniqueGO.png", plot = plot_Unique_GO(stomach_proximal_gut_unique, "Proximal Stomach"), dpi = 900)

#############################################################################################################################
# Code to investigate microbiota in the different tissues

# A function to read in the gene annotations for the annotation report. Capitalize gene names to keep them consistent
# Filter for only genes with annotations, and filter out eukaryotes and viruses to leave bacteria and archaea
# Keep only the Trinity ID with the highest bit score to avoid redundant annotations, and in the case of ties, only keep one
microbiome_annotation <- function(ann){
  annotation <- read_tsv(ann) %>% 
    mutate(gene_name = toupper(gene_name)) %>% 
    filter(!is.na(gene_name), !is.na(sp.BX_transcript_name)) %>% 
    filter(!grepl("Eukaryota|Viruses", sp.BX_Tax.Lineage)) %>% 
    group_by(TrinityID) %>% 
    slice_max(BitScore, n = 1, with_ties = FALSE) %>% 
    separate(sp.BX_Tax.Lineage, into = c("Lineage", "Genus"), sep = ";(?!(.|\n)*;)") %>% 
    mutate(Genus = str_squish(Genus)) %>% # Remove extra whitespace from genera names
    separate(Genus, into = c("Genus", "errata"), sep = "`") # Some genera have additional information attached, which must be separated
  return(annotation)
}

# Step through all the tissues with this microbiome-specific function
pyloric_microbiota <- microbiome_annotation("pyloric_ceca_annotation_bit50_E1e-6.txt")

ant_intestine_microbiota <- microbiome_annotation("anterior_intestine_annotation_bit50_E1e-6.txt")
brain_microbiota <- microbiome_annotation("brain_annotation_bit50_E1e-6.txt")
esophagus_microbiota <- microbiome_annotation("esophagus_annotation_bit50_E1e-6.txt")
gill_microbiota <- microbiome_annotation("gill_annotation_bit50_E1e-6.txt")
head_kidney_microbiota <- microbiome_annotation("head_kidney_annotation_bit50_E1e-6.txt") # Somehow the head kidney report is formatted differently than the others, so no taxonomic lineages are available.
heart_microbiota <- microbiome_annotation("heart_annotation_bit50_E1e-6.txt")
liver_microbiota <- microbiome_annotation("liver_annotation_bit50_E1e-6.txt")
rectum_microbiota <- microbiome_annotation("rectum_annotation_bit50_E1e-6.txt")
spiral_valve_microbiota <- microbiome_annotation("spiral_valve_annotation_bit50_E1e-6.txt")
stomach_distal_microbiota <- microbiome_annotation("stomach_distal_annotation_bit50_E1e-6.txt")
stomach_proximal_microbiota <- microbiome_annotation("stomach_proximal_annotation_bit50_E1e-6.txt")
white_muscle_microbiota <- microbiome_annotation("white_muscle_annotation_bit50_E1e-6.txt")

head_kidney2_microbiota <- microbiome_annotation("head_kidney2_annotation_bit50_E1e-6.txt")


# Take the various microbiome annotations and make a stacked barplot
# Make a tally of all the tissues
pyloric_tally <- pyloric_microbiota %>% 
  group_by(Genus) %>% 
  tally() %>% 
  mutate(tissue = "Pyloric Caeca")

ant_intestine_tally <- ant_intestine_microbiota %>% 
  group_by(Genus) %>% 
  tally() %>% 
  mutate(tissue = "Anterior Intestine")

brain_tally <- brain_microbiota %>% 
  group_by(Genus) %>% 
  tally() %>% 
  mutate(tissue = "Brain")

esophagus_tally <- esophagus_microbiota %>% 
  group_by(Genus) %>% 
  tally() %>% 
  mutate(tissue = "Esophagus")

gill_tally <- gill_microbiota %>% 
  group_by(Genus) %>% 
  tally() %>% 
  mutate(tissue = "Gill")

heart_tally <- heart_microbiota %>% 
  group_by(Genus) %>% 
  tally() %>% 
  mutate(tissue = "Heart")

liver_tally <- liver_microbiota %>% 
  group_by(Genus) %>% 
  tally() %>% 
  mutate(tissue = "Liver")

rectum_tally <- rectum_microbiota %>% 
  group_by(Genus) %>% 
  tally() %>% 
  mutate(tissue = "Rectum")

spiral_valve_tally <- spiral_valve_microbiota %>% 
  group_by(Genus) %>% 
  tally() %>% 
  mutate(tissue = "Spiral Valve")

stomach_proximal_tally <- stomach_proximal_microbiota %>% 
  group_by(Genus) %>% 
  tally() %>% 
  mutate(tissue = "Proximal Stomach")

stomach_distal_tally <- stomach_distal_microbiota %>% 
  group_by(Genus) %>% 
  tally() %>% 
  mutate(tissue = "Distal Stomach")

white_muscle_tally <- white_muscle_microbiota %>% 
  group_by(Genus) %>% 
  tally() %>% 
  mutate(tissue = "White Muscle")

head_kidney2_tally <- head_kidney2_microbiota %>% 
  group_by(Genus) %>% 
  tally() %>% 
  mutate(tissue = "Head Kidney")
###################################################################################################
# Read in lists of anaerobic and aerobic bacteria from bacdive.dsmz.de
# First the anaerobes, re-formatted from the website by Will Bugg
anaerobe_genera <- read_csv("Anerobe List BAC Drive_wsb.csv", skip = 2) %>% 
  select(`Final anerobic genus list no duplicates`) %>% 
  rename(Genus = `Final anerobic genus list no duplicates`) %>% 
  mutate(O2_tolerance = "anaerobe")

# Read in the aerobes
aerobic_genera <- read_csv("aerobe_export_bacdive_adv_search_table.csv", skip = 2) %>% 
  separate(species, into = c("Genus", "species"), sep = " ") %>% 
  select(Genus) %>% 
  distinct() %>% 
  mutate(O2_tolerance = "aerobe")

# Combine the two tables, call this genera metadata
#genera_metadata <- rbind(anaerobe_genera, aerobic_genera)
genera_metadata <- anti_join(anaerobe_genera, aerobic_genera, by = "Genus") %>% 
  filter(!is.na(Genus))

# Search a table of bacterial genera against the metadata, and create a new column of anaerobic or aerobic bacteria
#ant_intestine_microbiota <- left_join(ant_intestine_microbiota, genera_metadata)

##############################################################################
# Plot the bacterial patterns
# Combine all the tallies
# Head kidney was skipped because its bacterial transcripts were missing from the annotation
# Note, head kidney 2 represents a re-annotated head kidney transcriptome that had microbial nformation
all_tissue_microbe_tally <- rbind(esophagus_tally, stomach_proximal_tally, stomach_distal_tally, ant_intestine_tally, spiral_valve_tally, pyloric_tally, rectum_tally, brain_tally, gill_tally, heart_tally, liver_tally, white_muscle_tally, head_kidney2_tally)

# Lock in the order of these tissues by specifying tissues as factors
all_tissue_microbe_tally$tissue <- factor(all_tissue_microbe_tally$tissue, levels = c("Esophagus", "Proximal Stomach", "Distal Stomach", "Anterior Intestine", "Spiral Valve", "Pyloric Caeca", "Rectum", "Brain", "Gill", "Heart", "Liver", "White Muscle", "Head Kidney"))

# Look at a stacked bar plot
tissue_micro_barplot <- ggplot(all_tissue_microbe_tally, aes(fill=Genus, y=n, x=tissue)) + 
  geom_bar(position="stack", stat="identity") +
  xlab("Tissue") +
  ylab("Number of Genera Weighted by Unique Transcripts Present") +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1), text = element_text(size = 14))
#tissue_micro_barplot
ggsave(filename = "genera_barplot.png", plot = tissue_micro_barplot, dpi = 900, width = 22)

# Count up the number of unique genera per tissue
all_tissue_unique_genera_tally <- all_tissue_microbe_tally %>% 
  group_by(tissue, Genus) %>% 
  tally()

# Add in the metadata of aerobic or anaerobic bacteria to the tissue tally
all_tissue_unique_genera_tally <- left_join(all_tissue_unique_genera_tally, genera_metadata)

# Make a similar stacked barplot, this one *not* weighted by number of genes
tissue_micro_uniquegenera_barplot <- ggplot(all_tissue_unique_genera_tally, aes(y=n, x=tissue)) + 
  geom_bar(stat="identity") +
  xlab("Tissue") +
  ylab("Number of Unique Bacterial Genera") +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1), text = element_text(size = 20))
tissue_micro_uniquegenera_barplot
ggsave(filename = "unique_genera_barplot.png", plot = tissue_micro_uniquegenera_barplot, dpi = 900)

# A stacked barplot, coloured by anaerobe vs aerobe
o2_barplot <- ggplot(all_tissue_unique_genera_tally, aes(fill=O2_tolerance, y=n, x=tissue)) + 
  geom_bar(position="stack", stat="identity") +
  xlab("Tissue") +
  ylab("Number of Unique Bacterial Genera") +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust = 1), text = element_text(size = 20))
o2_barplot


# An upsetR plot to look at overlap in bacterial genes among the different data sets in just the gut
microbiome_list <- list(unique(esophagus_microbiota$gene_name), unique(stomach_proximal_microbiota$gene_name), unique(stomach_distal_microbiota$gene_name), unique(ant_intestine_microbiota$gene_name), unique(spiral_valve_microbiota$gene_name), unique(pyloric_microbiota$gene_name), unique(rectum_microbiota$gene_name), unique(brain_microbiota$gene_name), unique(gill_microbiota$gene_name), unique(heart_microbiota$gene_name), unique(liver_microbiota$gene_name), unique(white_muscle_microbiota$gene_name), unique(head_kidney2_microbiota$gene_name))
microbiome_gene_upset_data <- reshape2::melt(microbiome_list) %>%
  magrittr::set_colnames(c("dataset", "gene_name")) %>%
  dplyr::group_by(dataset, gene_name) %>%
  dplyr::mutate(value = 1)  %>%
  tidyr::spread(gene_name, value, fill = 0) %>% 
  rename(gene_name = 1, Esophagus = 2, `Glandular Stomach` = 3, `Muscular Stomach` = 4, `Anterior Intestine` = 5, `Spiral Valve` = 6, `Pyloric Caecum` = 7, Rectum = 8, Brain = 9, Gill = 10, Heart = 11, Liver = 12, `White Muscle` = 13, `Head Kidney` = 14) 

microbiome_gene_upset_data <- as.data.frame(microbiome_gene_upset_data)
for(i in 2:ncol(microbiome_gene_upset_data)){ microbiome_gene_upset_data[ , i] <- as.integer(microbiome_gene_upset_data[ , i]) }
# View the upset plot
upset(microbiome_gene_upset_data, nsets = ncol(microbiome_gene_upset_data) - 1, order.by = "freq", sets.x.label = "Number of Microbial Genes")

View(filter(microbiome_gene_upset_data, across(is.numeric, ~ .x == 1))) # The microbial-annotated genes present among all tissues
# Look at the number of microbial genes present in each tissue
microbiome_gene_upset_data %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))


# An upset plot of bacterial genera
microbiome_genera_list <- list(unique(esophagus_microbiota$Genus), unique(stomach_proximal_microbiota$Genus), unique(stomach_distal_microbiota$Genus), unique(ant_intestine_microbiota$Genus), unique(spiral_valve_microbiota$Genus), unique(pyloric_microbiota$Genus), unique(rectum_microbiota$Genus), unique(brain_microbiota$Genus), unique(gill_microbiota$Genus), unique(heart_microbiota$Genus), unique(liver_microbiota$Genus), unique(white_muscle_microbiota$Genus), unique(head_kidney2_microbiota$Genus))
microbiome_genera_upset_data <- reshape2::melt(microbiome_genera_list) %>%
  magrittr::set_colnames(c("dataset", "Genus")) %>%
  dplyr::group_by(dataset, Genus) %>%
  dplyr::mutate(value = 1)  %>%
  tidyr::spread(Genus, value, fill = 0) %>% 
  rename(Genus = 1, Esophagus = 2, `Glandular Stomach` = 3, `Muscular Stomach` = 4, `Anterior Intestine` = 5, `Spiral Valve` = 6, `Pyloric Caecum` = 7, Rectum = 8, Brain = 9, Gill = 10, Heart = 11, Liver = 12, `White Muscle` = 13, `Head Kidney` = 14) 

microbiome_genera_upset_data <- as.data.frame(microbiome_genera_upset_data)
for(i in 2:ncol(microbiome_genera_upset_data)){ microbiome_genera_upset_data[ , i] <- as.integer(microbiome_genera_upset_data[ , i]) }
# View the upset plot
upset(microbiome_genera_upset_data, nsets = ncol(microbiome_genera_upset_data) - 1, order.by = "freq", sets.x.label = "Number of Microbial Genera")

filter(microbiome_genera_upset_data, across(is.numeric, ~ .x == 1))

# Look at the number of microbial genera present in each tissue
microbiome_genera_upset_data %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

upset_genera_o2 <- left_join(microbiome_genera_upset_data, genera_metadata)
View(filter(upset_genera_o2, across(is.numeric, ~ .x == 1)))

#######################################################

# Similar to the previous enrichR function, but this one does not read in the annotation report and only searches 3 databases and collates results
enrichr_search <- function(report){
# Run enrichR
enrichr_results <- enrichr(unique(report$gene_name), databases = c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021"))

# Pull the results from one of the databases, biological process first
bio_process <- as_tibble(enrichr_results[["GO_Biological_Process_2021"]])
bio_process <- dplyr::filter(bio_process, Adjusted.P.value < 0.05)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
bio_process <- split_go(bio_process)

# Pull the results from one of the databases, now molecular function
mol_function <- as_tibble(enrichr_results[["GO_Molecular_Function_2021"]])
mol_function <- dplyr::filter(mol_function, Adjusted.P.value < 0.05)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
mol_function <- split_go(mol_function)

# Pull the results from one of the databases, last cellular component
cell_comp <- as_tibble(enrichr_results[["GO_Cellular_Component_2021"]])
cell_comp <- dplyr::filter(cell_comp, Adjusted.P.value < 0.05)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
cell_comp <- split_go(cell_comp)

tissue_list <- list(bio = bio_process, mol = mol_function, cell = cell_comp)
return(tissue_list)
}
pyloric_micro_enrichr <- enrichr_search(pyloric_microbiota)
nrow(pyloric_micro_enrichr$bio)

ant_intestine_micro_enrichr <- enrichr_search(ant_intestine_microbiota)
nrow(ant_intestine_micro_enrichr$bio)

brain_micro_enrichr <- enrichr_search(brain_microbiota)
nrow(brain_micro_enrichr$bio)

esophagus_micro_enrichr <- enrichr_search(esophagus_microbiota)
nrow(esophagus_micro_enrichr$bio)

gill_micro_enrichr <- enrichr_search(gill_microbiota)
nrow(gill_micro_enrichr$bio)

heart_micro_enrichr <- enrichr_search(heart_microbiota)
nrow(heart_micro_enrichr$bio)

liver_micro_enrichr <- enrichr_search(liver_microbiota)
nrow(liver_micro_enrichr$bio)

rectum_micro_enrichr <- enrichr_search(rectum_microbiota)
nrow(rectum_micro_enrichr$bio)

spiral_valve_micro_enrichr <- enrichr_search(spiral_valve_microbiota)
nrow(spiral_valve_micro_enrichr$bio)

stomach_distal_micro_enrichr <- enrichr_search(stomach_distal_microbiota)
nrow(stomach_distal_micro_enrichr$bio)

# Nothing came up as significant

##########################################################################################################################
# Trying a superexacttest to look at possible differences among tissues
# data("eqtls")
# (length.gene.sets=sapply(cis.eqtls,length))
# total=18196
# (num.expcted.overlap=total*do.call(prod,as.list(length.gene.sets/total)))
# (p=sapply(0:101,function(i) dpsets(i, length.gene.sets, n=total)))
# fit=MSET(cis.eqtls, n=total, lower.tail=FALSE)
# fit$FE
# fit$p.value

# The super exact test applied to just the genes annotated to microbes
exact_fit <- MSET(microbiome_list, n = nrow(microbiome_gene_upset_data), lower.tail = FALSE)
exact_fit$FE
exact_fit$p.value

#(length.gene.sets=sapply(microbiome_list,length))

#(num.expcted.overlap=nrow(microbiome_gene_upset_data)*do.call(prod,as.list(length.gene.sets/nrow(microbiome_gene_upset_data))))

#(p=sapply(0:101,function(i) dpsets(i, length.gene.sets, n=nrow(microbiome_gene_upset_data))))


#res=supertest(microbiome_list, n=nrow(microbiome_gene_upset_data))
#microbial_gene_significance <- summary(res)

#microbial_gene_significance <- microbial_gene_significance$Table
#microbial_gene_significance
#plot(res, sort.by="size", margin=c(2,2,2,2), color.scale.pos=c(0.85,1), legend.pos=c(0.9,0.15))
#plot(res, Layout="landscape", degree=2:4, sort.by="size", margin=c(0.5,5,1,2))

###################################################################################
# Use PCAs to visualize how different the transcriptomes are

# Run the PCA, transposed so the transcriptomes and not the GO terms are analyzed
bio_pca <- prcomp(t(bio_upset_data[,-1]))
# Get variance explained from the PCA
summary(bio_pca)
# Look at an initial plot
plot(bio_pca$x)
# Pull the PCA loading information
bio_pca_loadings <- as_tibble(bio_pca$x)
# Check the row names for the PCA, manually set transcriptome names for the plot using those row names
rownames(bio_pca$x)
bio_pca_loadings$transcriptome <- c("Anterior Intestine", "Brain", "Esophagus", "Gill", "Head Kidney", "Heart", "Liver", "Pyloric Caecum", "Rectum", "Spiral Valve", "Muscular Stomach", "Glandular Stomach", "White Muscle")
  
# Create a manuscript-ready figure
bio_pca_plot <- ggplot(data = bio_pca_loadings, aes(x = PC1, y = PC2, label = transcriptome)) +
  geom_point() +
  geom_text_repel() +
  xlab("PC1 (22.6% variance)") +
  ylab("PC2 (12.5% variance)") +
  theme_bw() +
  theme(text = element_text(size = 18), axis.text.x = element_blank(), axis.text.y = element_blank())

# Repeat the process for the molecular function GO terms
mol_pca <- prcomp(t(mol_upset_data[,-1]))
# Get variance explained from the PCA
summary(mol_pca)
# Look at an initial plot
plot(mol_pca$x)
# Pull the PCA loading information
mol_pca_loadings <- as_tibble(mol_pca$x)
# Check the row names for the PCA, manually set transcriptome names for the plot using those row names
rownames(mol_pca$x)
mol_pca_loadings$transcriptome <- c("Anterior Intestine", "Brain", "Esophagus", "Gill", "Head Kidney", "Heart", "Liver", "Pyloric Caecum", "Rectum", "Spiral Valve", "Muscular Stomach", "Glandular Stomach", "White Muscle")

# Create a manuscript-ready figure
mol_pca_plot <- ggplot(data = mol_pca_loadings, aes(x = PC1, y = PC2, label = transcriptome)) +
  geom_point() +
  geom_text_repel() +
  xlab("PC1 (22.4% variance)") +
  ylab("PC2 (14.6% variance)") +
  theme_bw() +
  theme(text = element_text(size = 18), axis.text.x = element_blank(), axis.text.y = element_blank())

# Repeat the process for the cellular component GO terms
cell_pca <- prcomp(t(cell_upset_data[,-1]))
# Get variance explained from the PCA
summary(cell_pca)
# Look at an initial plot
plot(cell_pca$x)
# Pull the PCA loading information
cell_pca_loadings <- as_tibble(cell_pca$x)
# Check the row names for the PCA, manually set transcriptome names for the plot using those row names
rownames(cell_pca$x)
cell_pca_loadings$transcriptome <- c("Anterior Intestine", "Brain", "Esophagus", "Gill", "Head Kidney", "Heart", "Liver", "Pyloric Caecum", "Rectum", "Spiral Valve", "Muscular Stomach", "Glandular Stomach", "White Muscle")

# Create a manuscript-ready figure
cell_pca_plot <- ggplot(data = cell_pca_loadings, aes(x = PC1, y = PC2, label = transcriptome)) +
  geom_point() +
  geom_text_repel() +
  xlab("PC1 (29.5% variance)") +
  ylab("PC2 (12.5% variance)") +
  theme_bw() +
  theme(text = element_text(size = 18), axis.text.x = element_blank(), axis.text.y = element_blank())

# Combine the 3 different GO term PCA plots with patchwork
combined_pca_plot <- bio_pca_plot / mol_pca_plot / cell_pca_plot + plot_annotation(tag_levels = 'A')
# Write out the resulting figure
ggsave(plot = combined_pca_plot, filename = "combined_pca.pdf", height = 15, width = 10)

########################################################################
# Look at a PCA with genes among all 13 tissues
# Run the custom function to find GO terms for all 13 tissues
ant_intestine_ann <- read_annotation("anterior_intestine_annotation_bit50_E1e-6.txt")
brain_ann <- read_annotation("brain_annotation_bit50_E1e-6.txt")
esophagus_ann <- read_annotation("esophagus_annotation_bit50_E1e-6.txt")
gill_ann <- read_annotation("gill_annotation_bit50_E1e-6.txt")
#head_kidney_ann <- read_annotation("head_kidney_annotation_bit50_E1e-6.txt")
heart_ann <- read_annotation("heart_annotation_bit50_E1e-6.txt")
liver_ann <- read_annotation("liver_annotation_bit50_E1e-6.txt")
pyloric_ceca_ann <- read_annotation("pyloric_ceca_annotation_bit50_E1e-6.txt")
rectum_ann <- read_annotation("rectum_annotation_bit50_E1e-6.txt")
spiral_valve_ann <- read_annotation("spiral_valve_annotation_bit50_E1e-6.txt")
stomach_distal_ann <- read_annotation("stomach_distal_annotation_bit50_E1e-6.txt")
stomach_proximal_ann <- read_annotation("stomach_proximal_annotation_bit50_E1e-6.txt")
white_muscle_ann <- read_annotation("white_muscle_annotation_bit50_E1e-6.txt")

head_kidney2_ann <- read_annotation("head_kidney2_annotation_bit50_E1e-6.txt")

# Take genes from the annotation reports, and format them for an UpSet plot
gene_annotation_list <- list(ant_intestine_ann$gene_name, brain_ann$gene_name, esophagus_ann$gene_name, gill_ann$gene_name, head_kidney2_ann$gene_name, heart_ann$gene_name, liver_ann$gene_name, pyloric_ceca_ann$gene_name, rectum_ann$gene_name, spiral_valve_ann$gene_name, stomach_distal_ann$gene_name, stomach_proximal_ann$gene_name, white_muscle_ann$gene_name)

gene_upset_data <- reshape2::melt(gene_annotation_list) %>%
  magrittr::set_colnames(c("dataset", "gene")) %>%
  dplyr::group_by(dataset, gene) %>%
  dplyr::mutate(value = 1)  %>%
  tidyr::spread(gene, value, fill = 0) %>% 
  rename(gene = 1, ant_intestine = 2, brain = 3, esophagus = 4, gill = 5, head_kidney = 6, heart = 7, liver = 8, pyloric_ceca = 9, rectum = 10, spiral_valve = 11, stomach_distal = 12, stomach_proximal = 13, white_muscle = 14) 
gene_upset_data <- as.data.frame(gene_upset_data)

for(i in 2:ncol(gene_upset_data)){ gene_upset_data[ , i] <- as.integer(gene_upset_data[ , i]) }
# View the upset plot
upset(gene_upset_data, nsets = ncol(gene_upset_data) - 1, order.by = "freq", sets.x.label = "Gene Set Size")


# Take the upset data and run a PCA with all the transcriptomes
gene_pca <- prcomp(t(gene_upset_data[,-1]))
# Get variance explained from the PCA
summary(gene_pca)
# Look at an initial plot
plot(gene_pca$x)
# Pull the PCA loading information
gene_pca_loadings <- as_tibble(gene_pca$x)
# Check the row names for the PCA, manually set transcriptome names for the plot using those row names
rownames(gene_pca$x)
gene_pca_loadings$transcriptome <- c("Anterior Intestine", "Brain", "Esophagus", "Gill", "Head Kidney", "Heart", "Liver", "Pyloric Caecum", "Rectum", "Spiral Valve", "Muscular Stomach", "Glandular Stomach", "White Muscle")

gene_pca_loadings$Category <- c("Gut", "Peripheral", "Gut", "Peripheral", "Peripheral", "Peripheral", "Peripheral", "Gut", "Gut", "Gut", "Gut", "Gut", "Peripheral")

# Create a manuscript-ready figure from the PCA
gene_pca_plot <- ggplot(data = gene_pca_loadings, aes(x = PC1, y = PC2, label = transcriptome, col = Category)) +
  geom_point(size = 2) +
  geom_text_repel(size = 8, show.legend = FALSE) +
  xlab("PC1 (23.7% variance)") +
  ylab("PC2 (12.8% variance)") +
  scale_colour_fish_d(option = "Epinephelus_lanceolatus", direction = -1) +
  theme_bw() +
  theme(text = element_text(size = 26), axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = c(.88, .9), legend.background = element_blank())
gene_pca_plot

ggsave(filename = "gene_pca.pdf", plot = gene_pca_plot, dpi = 4000, scale = 1.5, height = 5, width = 5)
ggsave(filename = "gene_pca2.pdf", plot = gene_pca_plot, dpi = 4000, scale = 1.5, height = 6.24, width = 8.55)

########################################################################
# Read in the overall transcriptome, created much later than the others and repeat the process of looking at GO terms
overall_GO <- enrichR_tissue("overall_annotation_bit50_E1e-6.txt")

overall_enrichr_report <- read_tsv("overall_annotation_bit50_E1e-6.txt") %>% 
  mutate(gene_name = toupper(gene_name))
# Run enrichR
overall_enrichr_results <- enrichr(unique(overall_enrichr_report$gene_name), databases = c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021"))

# Pull the results from one of the databases, biological process first
overall_bio_process <- as_tibble(overall_enrichr_results[["GO_Biological_Process_2021"]])
overall_bio_process <- dplyr::filter(overall_bio_process, Adjusted.P.value < 0.05)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
overall_bio_process <- split_go(overall_bio_process)

overall_mol_function <- as_tibble(overall_enrichr_results[["GO_Molecular_Function_2021"]])
overall_mol_function <- dplyr::filter(overall_mol_function, Adjusted.P.value < 0.05)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
overall_mol_function <- split_go(overall_mol_function)

overall_cell_comp <- as_tibble(overall_enrichr_results[["GO_Cellular_Component_2021"]])
overall_cell_comp <- dplyr::filter(overall_cell_comp, Adjusted.P.value < 0.05)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
overall_cell_comp <- split_go(overall_cell_comp)

# For each of the annotations, count the number of GO terms from each category
print(nrow(overall_enrichr_results$GO_Biological_Process_2021)) # Biological process
print(nrow(overall_enrichr_results$GO_Molecular_Function_2021)) # Molecular function
print(nrow(overall_enrichr_results$GO_Cellular_Component_2021)) # Cellular component

overall_ann <- read_annotation("overall_annotation_bit50_E1e-6.txt")

# Make a PCA that includes the overall transcriptome
gene_annotation_list_overall <- list(ant_intestine_ann$gene_name, brain_ann$gene_name, esophagus_ann$gene_name, gill_ann$gene_name, head_kidney2_ann$gene_name, heart_ann$gene_name, liver_ann$gene_name, pyloric_ceca_ann$gene_name, rectum_ann$gene_name, spiral_valve_ann$gene_name, stomach_distal_ann$gene_name, stomach_proximal_ann$gene_name, white_muscle_ann$gene_name, overall_ann$gene_name)

gene_upset_data_overall <- reshape2::melt(gene_annotation_list_overall) %>%
  magrittr::set_colnames(c("dataset", "gene")) %>%
  dplyr::group_by(dataset, gene) %>%
  dplyr::mutate(value = 1)  %>%
  tidyr::spread(gene, value, fill = 0) %>% 
  rename(gene = 1, ant_intestine = 2, brain = 3, esophagus = 4, gill = 5, head_kidney = 6, heart = 7, liver = 8, pyloric_ceca = 9, rectum = 10, spiral_valve = 11, stomach_distal = 12, stomach_proximal = 13, white_muscle = 14, overall = 15) 
gene_upset_data_overall <- as.data.frame(gene_upset_data_overall)

for(i in 2:ncol(gene_upset_data_overall)){ gene_upset_data_overall[ , i] <- as.integer(gene_upset_data_overall[ , i]) }
# View the upset plot. This is a supplementary figure in the manuscript
upset(gene_upset_data_overall, nsets = ncol(gene_upset_data_overall) - 1, order.by = "freq", sets.x.label = "Gene Set Size")


# Pull the group of genes unique to the overall transcriptome. First, sum presence across the 13 tissue-specific transcriptomes. Filter for genes that are not present in the tissue-specific transcriptomes, but are present in the overall transcriptome.
overall_uniqueness <- gene_upset_data_overall %>% 
  mutate(tissue_sum = rowSums(.[2:14])) %>% 
  filter(tissue_sum == 0 & overall == 1)

# Run enrichR on the 3,065 genes unique to the overall transcriptome
overall_unique_enrichr <- enrichr(unique(overall_uniqueness$gene), databases = c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021"))

# Pull the results from one of the databases, biological process first
overall_bio_process <- as_tibble(overall_unique_enrichr[["GO_Biological_Process_2021"]])
overall_bio_process <- dplyr::filter(overall_bio_process, Adjusted.P.value < 0.05)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
overall_bio_process <- split_go(overall_bio_process)
nrow(overall_bio_process)

overall_mol_function <- as_tibble(overall_unique_enrichr[["GO_Molecular_Function_2021"]])
overall_mol_function <- dplyr::filter(overall_mol_function, Adjusted.P.value < 0.05)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
overall_mol_function <- split_go(overall_mol_function)
nrow(overall_mol_function)

overall_cell_comp <- as_tibble(overall_unique_enrichr[["GO_Cellular_Component_2021"]])
overall_cell_comp <- dplyr::filter(overall_cell_comp, Adjusted.P.value < 0.05)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
overall_cell_comp <- split_go(overall_cell_comp)
nrow(overall_cell_comp)

# Take the upset data and run a PCA with all the transcriptomes
gene_pca_overall <- prcomp(t(gene_upset_data_overall[,-1]))
# Get variance explained from the PCA
summary(gene_pca_overall)
# Look at an initial plot
plot(gene_pca_overall$x)
# Pull the PCA loading information
gene_pca_loadings_overall <- as_tibble(gene_pca_overall$x)
# Check the row names for the PCA, manually set transcriptome names for the plot using those row names
rownames(gene_pca_overall$x)
gene_pca_loadings_overall$transcriptome <- c("Anterior Intestine", "Brain", "Esophagus", "Gill", "Head Kidney", "Heart", "Liver", "Pyloric Caecum", "Rectum", "Spiral Valve", "Muscular Stomach", "Glandular Stomach", "White Muscle", "Overall")

gene_pca_loadings_overall$Category <- c("Gut", "Peripheral", "Gut", "Peripheral", "Peripheral", "Peripheral", "Peripheral", "Gut", "Gut", "Gut", "Gut", "Gut", "Peripheral", "Overall")

# Create a manuscript-ready figure from the PCA
gene_pca_plot_overall <- ggplot(data = gene_pca_loadings_overall, aes(x = PC1, y = PC2, label = transcriptome, col = Category)) +
  geom_point(size = 2) +
  geom_text_repel(size = 8, show.legend = FALSE) +
  xlab("PC1 (23.7% variance)") +
  ylab("PC2 (12.8% variance)") +
  scale_colour_fish_d(option = "Epinephelus_lanceolatus", direction = -1) +
  theme_bw() +
  theme(text = element_text(size = 26), axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = c(.88, .9), legend.background = element_blank())
gene_pca_plot_overall
