library(viralORF)
library(dplyr)
library(magrittr)
library(stringr)

# --- Input --- 
all_ORFs <- viralORF::find_startstop_ORFs("your_sequence_file.fasta") # find all ORFs across reading frames of any sequence you want 
kozak_reference <- read.csv("kozak_strength_table.csv") # kozak_reference should include: dinuc, kozak_strength 

# PRE-PROCESSING 
all_ORFs <- all_ORFs %>%
  mutate(frame = case_when(
    frame == "plus_0" ~ "1", 
    frame == "plus_1" ~ "2", 
    frame == "plus_2" ~ "3", 
    TRUE ~ frame
  )) 

# Isolate the +1 and +2 ORFs (candidate AcORFs)
plus_1_2_orfs <- all_ORFs %>% 
  filter(frame %in% c("2", "3")) %>%
  group_by(seq, end) %>%
  arrange(desc(length)) %>%
  slice(1) %>%
  ungroup() %>%
  select(seq, ORF_ID, start, end, length, nuc_sequence, kozak_dinuc, frame, sense, kozak_strength)

# Isolate the +0 ORFs (reference gene) 
plus_0_orfs <- all_ORFs %>%
  filter(frame == "1") %>%
  group_by(seq) %>%
  arrange(desc(length)) %>%
  slice(1) %>%
  ungroup() %>%
  select(seq, gene_start = start, gene_length = length, gene_kozak_strength = kozak_strength, gene_seq = nuc_sequence, frame)

# Combine reference gene and acORFs 
combined_data <- plus_0_orfs %>%
  left_join(plus_1_2_orfs, by = "seq", suffix = c("_gene", "_orf")) %>%
  filter(!is.na(start)) %>%
  ungroup() 

# Get intervening ATG Kozak strengths 
get_iATG_kozaks <- function(seq_region, kozak_reference, gene_strength) {
  kozak_strengths <- c(if (is.na(gene_strength)) "unknown" else gene_strength) 
  if (is.na(seq_region) || nchar(seq_region) < 6) return(kozak_strengths) 
  atg_positions <- unlist(gregexpr("ATG", seq_region, fixed = TRUE)) 
  if (length(atg_positions) == 0 || atg_positions[1] == -1) return(kozak_strengths) 
  intervening_positions <- atg_positions[atg_positions > 1] 
  if (length(intervening_positions) > 0) { 
    for (pos in intervening_positions) { 
      pre <- if (pos >= 3) substr(seq_region, pos - 3, pos - 3) else "_" 
      post <- if ((pos + 3) <= nchar(seq_region)) substr(seq_region, pos + 3, pos + 3) else "_" 
      dinuc <- paste0(pre, post) 
      kozak_strength <- kozak_reference$kozak_strength[match(dinuc, kozak_reference$dinuc)] 
      if (length(kozak_strength) == 0 || is.na(kozak_strength)) kozak_strength <- "unknown" 
      kozak_strengths <- c(kozak_strengths, kozak_strength) 
    } 
  } 
  return(kozak_strengths) 
} 

# Convert Kozak qualitative strength to skipping probability 
convert_kozak_to_prob <- function(k) { 
  if (is.na(k)) return(1) 
  k <- tolower(k) 
  if (k == "strong") return(0.15) 
  if (k == "moderate") return(0.453333333) 
  if (k == "weak") return(0.85) 
  return(1) 
} 

# Distance scoring 
get_distance_score <- function(distance) { 
  if (is.na(distance)) return(0) 
  if (distance < 100) return(1) 
  if (distance > 300) return(0) 
  return(0.5) 
} 

# Length scoring 
get_length_score <- function(length) { 
  if (is.na(length)) return(0) 
  if (length >= 120) return(1) 
  return(0) 
} 

# Main computation 
combined_data_with_iATGs <- combined_data %>% 
  rowwise() %>% 
  mutate( 
    unaligned_5prime = if_else(
      !is.na(gene_start) & !is.na(start) & (as.numeric(start) > as.numeric(gene_start)), 
      substr(gene_seq, 1, (as.numeric(start) - as.numeric(gene_start))), 
      "" 
    ), 
    iATG_Vs = list(get_iATG_kozaks(unaligned_5prime, kozak_reference, gene_kozak_strength)), 
    unaligned_atg_count = length(iATG_Vs) 
  ) %>% 
  ungroup() 

# Expand iATG columns 
max_iATGs <- max(lengths(combined_data_with_iATGs$iATG_Vs), na.rm = TRUE) 
iATG_cols <- paste0("iATG_V", seq_len(max_iATGs)) 
for (col in iATG_cols) { 
  combined_data_with_iATGs[[col]] <- NA_character_ 
} 
for (i in seq_len(nrow(combined_data_with_iATGs))) { 
  vset <- combined_data_with_iATGs$iATG_Vs[[i]] 
  if (length(vset) > 0) { 
    for (j in seq_along(vset)) { 
      combined_data_with_iATGs[[iATG_cols[j]]][i] <- vset[j] 
    } 
  } 
} 
combined_data_with_iATGs$iATG_Vs <- NULL 

# Compute final numeric scores and initiation potential 
combined_data_with_scores <- combined_data_with_iATGs %>% 
  rowwise() %>% 
  mutate( 
    distance = if_else(!is.na(start) & !is.na(gene_start), as.numeric(start) - as.numeric(gene_start), NA_real_), 
    kozak_score = convert_kozak_to_prob(kozak_strength), 
    initiation_potential = 1 - kozak_score, 
    distance_score = get_distance_score(distance), 
    length_score = get_length_score(length), 
    iATG_prob_product = { 
      raw_data <- c_across(all_of(iATG_cols)) 
      numeric_vals <- vapply(na.omit(raw_data), function(k) { 
        as.numeric(convert_kozak_to_prob(k)) 
      }, FUN.VALUE = numeric(1)) 
      if (length(numeric_vals) == 0) 1 else prod(numeric_vals) 
    }, 
    final_score = initiation_potential * distance_score * length_score * iATG_prob_product 
  ) %>% 
  ungroup() 

# Save output 
output_scored_name <- "gene_scored_ORFs.csv" 
write.csv(combined_data_with_scores, file = output_scored_name, row.names = FALSE) 
cat("Scoring complete.\n") 
cat("File saved as:", output_scored_name, "\n")