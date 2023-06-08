library(tidyverse)

gff <- read_tsv("ncbi-genomes-2023-06-07/GCF_010645085.2_ASM1064508v2_genomic.gff", skip = 8, col_names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")) |>
  separate(col = attribute, into = c("sseqid", "attribute"), sep = ";", extra = "merge") |>
  mutate(sseqid = str_replace_all(sseqid, "ID=gene-|ID=rna-|ID=exon-|ID=cds-|ID=id-", "")) |>
  mutate(sseqid = str_replace_all(sseqid, "ID=", ""))



# Read in the different transcriptomes, add column names based on the diamond/blast default settings for blastn, then fitler for e values less than 1e-6 and bit scores greater than 50 consistent with the rest of the manuscript
overall_res <- read_tsv("sterlet_blast_res/overall_out.tsv", col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")) |> 
  filter(evalue < 1e-6, bitscore > 50) |>
  full_join(gff, by = "sseqid") |>
  group_by(qseqid) |>
  filter(evalue == min(evalue), bitscore == max(bitscore), pident == max(pident)) |>
  slice_head(n = 1) |>
  rename(LKST_transcript_ID = qseqid, Sterlet_ID = qseqid, perc_identity = pident, LKST_start = qstart, LKST_end = qend, sterlet_start = sstart, sterlet_end = send, sterlet_chrom = seqname) # Rename column names from the defaults to things more informative for supplementary tables

write_tsv(x = overall_res, file = "LKST_overall_sterlet.txt")

# Repeat the process for the rest of the transcriptomes
antintestine_res <- read_tsv("sterlet_blast_res/Trinity_anterior_intestine.fasta_out.tsv", col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")) |> 
  filter(evalue < 1e-6, bitscore > 50) |>
  full_join(gff, by = "sseqid") |>
  group_by(qseqid) |>
  filter(evalue == min(evalue), bitscore == max(bitscore), pident == max(pident)) |>
  slice_head(n = 1) |>
  rename(LKST_transcript_ID = qseqid, Sterlet_ID = qseqid, perc_identity = pident, LKST_start = qstart, LKST_end = qend, sterlet_start = sstart, sterlet_end = send, sterlet_chrom = seqname) # Rename column names from the defaults to things more informative for supplementary tables

write_tsv(x = antintestine_res, file = "LKST_ant_intestine_sterlet.txt")

brain_res <- read_tsv("sterlet_blast_res/Trinity_brain.fasta_out.tsv", col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")) |> 
  filter(evalue < 1e-6, bitscore > 50) |>
  full_join(gff, by = "sseqid") |>
  group_by(qseqid) |>
  filter(evalue == min(evalue), bitscore == max(bitscore), pident == max(pident)) |>
  slice_head(n = 1) |>
  rename(LKST_transcript_ID = qseqid, Sterlet_ID = qseqid, perc_identity = pident, LKST_start = qstart, LKST_end = qend, sterlet_start = sstart, sterlet_end = send, sterlet_chrom = seqname) # Rename column names from the defaults to things more informative for supplementary tables

write_tsv(x = brain_res, file = "LKST_brain_sterlet.txt")


esophagus_res <- read_tsv("sterlet_blast_res/Trinity_esophagus.fasta_out.tsv", col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")) |> 
  filter(evalue < 1e-6, bitscore > 50) |>
  full_join(gff, by = "sseqid") |>
  group_by(qseqid) |>
  filter(evalue == min(evalue), bitscore == max(bitscore), pident == max(pident)) |>
  slice_head(n = 1) |>
  rename(LKST_transcript_ID = qseqid, Sterlet_ID = qseqid, perc_identity = pident, LKST_start = qstart, LKST_end = qend, sterlet_start = sstart, sterlet_end = send, sterlet_chrom = seqname) # Rename column names from the defaults to things more informative for supplementary tables

write_tsv(x = esophagus_res, file = "LKST_esophagus_sterlet.txt")


gill_res <- read_tsv("sterlet_blast_res/Trinity_gill.fasta_out.tsv", col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")) |> 
  filter(evalue < 1e-6, bitscore > 50) |>
  full_join(gff, by = "sseqid") |>
  group_by(qseqid) |>
  filter(evalue == min(evalue), bitscore == max(bitscore), pident == max(pident)) |>
  slice_head(n = 1) |>
  rename(LKST_transcript_ID = qseqid, Sterlet_ID = qseqid, perc_identity = pident, LKST_start = qstart, LKST_end = qend, sterlet_start = sstart, sterlet_end = send, sterlet_chrom = seqname) # Rename column names from the defaults to things more informative for supplementary tables

write_tsv(x = gill_res, file = "LKST_gill_sterlet.txt")


head_kidney_res <- read_tsv("sterlet_blast_res/Trinity_head_kidney.fasta_out.tsv", col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")) |> 
  filter(evalue < 1e-6, bitscore > 50) |>
  full_join(gff, by = "sseqid") |>
  group_by(qseqid) |>
  filter(evalue == min(evalue), bitscore == max(bitscore), pident == max(pident)) |>
  slice_head(n = 1) |>
  rename(LKST_transcript_ID = qseqid, Sterlet_ID = qseqid, perc_identity = pident, LKST_start = qstart, LKST_end = qend, sterlet_start = sstart, sterlet_end = send, sterlet_chrom = seqname) # Rename column names from the defaults to things more informative for supplementary tables

write_tsv(x = head_kidney_res, file = "LKST_head_kidney_sterlet.txt")


heart_res <- read_tsv("sterlet_blast_res/Trinity_heart.fasta_out.tsv", col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")) |> 
  filter(evalue < 1e-6, bitscore > 50) |>
  full_join(gff, by = "sseqid") |>
  group_by(qseqid) |>
  filter(evalue == min(evalue), bitscore == max(bitscore), pident == max(pident)) |>
  slice_head(n = 1) |>
  rename(LKST_transcript_ID = qseqid, Sterlet_ID = qseqid, perc_identity = pident, LKST_start = qstart, LKST_end = qend, sterlet_start = sstart, sterlet_end = send, sterlet_chrom = seqname) # Rename column names from the defaults to things more informative for supplementary tables

write_tsv(x = heart_res, file = "LKST_heart_sterlet.txt")


liver_res <- read_tsv("sterlet_blast_res/Trinity_liver.fasta_out.tsv", col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")) |> 
  filter(evalue < 1e-6, bitscore > 50) |>
  full_join(gff, by = "sseqid") |>
  group_by(qseqid) |>
  filter(evalue == min(evalue), bitscore == max(bitscore), pident == max(pident)) |>
  slice_head(n = 1) |>
  rename(LKST_transcript_ID = qseqid, Sterlet_ID = qseqid, perc_identity = pident, LKST_start = qstart, LKST_end = qend, sterlet_start = sstart, sterlet_end = send, sterlet_chrom = seqname) # Rename column names from the defaults to things more informative for supplementary tables

write_tsv(x = liver_res, file = "LKST_liver_sterlet.txt")


pyloric_ceca_res <- read_tsv("sterlet_blast_res/Trinity_pyloric_ceca.fasta_out.tsv", col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")) |> 
  filter(evalue < 1e-6, bitscore > 50) |>
  full_join(gff, by = "sseqid") |>
  group_by(qseqid) |>
  filter(evalue == min(evalue), bitscore == max(bitscore), pident == max(pident)) |>
  slice_head(n = 1) |>
  rename(LKST_transcript_ID = qseqid, Sterlet_ID = qseqid, perc_identity = pident, LKST_start = qstart, LKST_end = qend, sterlet_start = sstart, sterlet_end = send, sterlet_chrom = seqname) # Rename column names from the defaults to things more informative for supplementary tables

write_tsv(x = pyloric_ceca_res, file = "LKST_pyloric_ceca_sterlet.txt")


rectum_res <- read_tsv("sterlet_blast_res/Trinity_rectum.fasta_out.tsv", col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")) |> 
  filter(evalue < 1e-6, bitscore > 50) |>
  full_join(gff, by = "sseqid") |>
  group_by(qseqid) |>
  filter(evalue == min(evalue), bitscore == max(bitscore), pident == max(pident)) |>
  slice_head(n = 1) |>
  rename(LKST_transcript_ID = qseqid, Sterlet_ID = qseqid, perc_identity = pident, LKST_start = qstart, LKST_end = qend, sterlet_start = sstart, sterlet_end = send, sterlet_chrom = seqname) # Rename column names from the defaults to things more informative for supplementary tables

write_tsv(x = rectum_res, file = "LKST_rectum_sterlet.txt")


spiral_valve_res <- read_tsv("sterlet_blast_res/Trinity_spiral_valve.fasta_out.tsv", col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")) |> 
  filter(evalue < 1e-6, bitscore > 50) |>
  full_join(gff, by = "sseqid") |>
  group_by(qseqid) |>
  filter(evalue == min(evalue), bitscore == max(bitscore), pident == max(pident)) |>
  slice_head(n = 1) |>
  rename(LKST_transcript_ID = qseqid, Sterlet_ID = qseqid, perc_identity = pident, LKST_start = qstart, LKST_end = qend, sterlet_start = sstart, sterlet_end = send, sterlet_chrom = seqname) # Rename column names from the defaults to things more informative for supplementary tables

write_tsv(x = spiral_valve_res, file = "LKST_spiral_valve_sterlet.txt")


muscular_stomach_res <- read_tsv("sterlet_blast_res/Trinity_stomach_distal.fasta_out.tsv", col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")) |> 
  filter(evalue < 1e-6, bitscore > 50) |>
  full_join(gff, by = "sseqid") |>
  group_by(qseqid) |>
  filter(evalue == min(evalue), bitscore == max(bitscore), pident == max(pident)) |>
  slice_head(n = 1) |>
  rename(LKST_transcript_ID = qseqid, Sterlet_ID = qseqid, perc_identity = pident, LKST_start = qstart, LKST_end = qend, sterlet_start = sstart, sterlet_end = send, sterlet_chrom = seqname) # Rename column names from the defaults to things more informative for supplementary tables

write_tsv(x = muscular_stomach_res, file = "LKST_muscular_stomach_sterlet.txt")


glandular_stomach_res <- read_tsv("sterlet_blast_res/Trinity_stomach_proximal.fasta_out.tsv", col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")) |> 
  filter(evalue < 1e-6, bitscore > 50) |>
  full_join(gff, by = "sseqid") |>
  group_by(qseqid) |>
  filter(evalue == min(evalue), bitscore == max(bitscore), pident == max(pident)) |>
  slice_head(n = 1) |>
  rename(LKST_transcript_ID = qseqid, Sterlet_ID = qseqid, perc_identity = pident, LKST_start = qstart, LKST_end = qend, sterlet_start = sstart, sterlet_end = send, sterlet_chrom = seqname) # Rename column names from the defaults to things more informative for supplementary tables

write_tsv(x = glandular_stomach_res, file = "LKST_glandular_stomach_sterlet.txt")


white_muscle_res <- read_tsv("sterlet_blast_res/Trinity_white_muscle.fasta_out.tsv", col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")) |> 
  filter(evalue < 1e-6, bitscore > 50) |>
  full_join(gff, by = "sseqid") |>
  group_by(qseqid) |>
  filter(evalue == min(evalue), bitscore == max(bitscore), pident == max(pident)) |>
  slice_head(n = 1) |>
  rename(LKST_transcript_ID = qseqid, Sterlet_ID = qseqid, perc_identity = pident, LKST_start = qstart, LKST_end = qend, sterlet_start = sstart, sterlet_end = send, sterlet_chrom = seqname) # Rename column names from the defaults to things more informative for supplementary tables

write_tsv(x = white_muscle_res, file = "LKST_white_muslce_sterlet.txt")
