
library(tidyverse)
library(GetFeatistics)

restuls_BEHPOS_imported <- read_tsv("metab_stat_w_features_known_BEH_POS.txt")
restuls_C18NEG_imported <- read_tsv("metab_stat_w_features_known_C18_NEG.txt")
restuls_C18POS_imported <- read_tsv("metab_stat_w_features_known_C18_POS.txt")


sample_names <- c("x017_MIC_LAC_AMYL_rbiol2_rtech1", "x022_MIC_LAC_KEFI_rbiol2_rtech1", "x026_MIC_BIF_BREV_rbiol2_rtech1", "x028_MIC_LAC_KEFI_rbiol1_rtech1", "x032_MIC_LAC_SALI_rbiol1_rtech1", "x033_MIC_LAC_PLAN_rbiol1_rtech1",
                  "x038_MIC_LAC_AMYL_rbiol1_rtech1", "x041_MIC_LAC_SALI_rbiol2_rtech1", "x042_MIC_BIF_BREV_rbiol1_rtech1", "x043_MIC_LAC_PLAN_rbiol2_rtech1", "x047_MIC_BIF_BREV_rbiol1_rtech2", "x048_MIC_LAC_AMYL_rbiol2_rtech2",
                  "x049_MIC_BIF_BREV_rbiol2_rtech2", "x053_MIC_LAC_PLAN_rbiol2_rtech2", "x057_MIC_LAC_KEFI_rbiol1_rtech2", "x058_MIC_LAC_PLAN_rbiol1_rtech2", "x059_MIC_LAC_KEFI_rbiol2_rtech2", "x067_MIC_LAC_SALI_rbiol1_rtech2",
                  "x068_MIC_LAC_AMYL_rbiol1_rtech2", "x073_MIC_LAC_SALI_rbiol2_rtech2", "x078_MIC_LAC_SALI_rbiol1_rtech3", "x079_MIC_BIF_BREV_rbiol1_rtech3", "x080_MIC_LAC_KEFI_rbiol2_rtech3", "x082_MIC_LAC_KEFI_rbiol1_rtech3",
                  "x087_MIC_LAC_SALI_rbiol2_rtech3", "x090_MIC_LAC_AMYL_rbiol2_rtech3", "x092_MIC_LAC_PLAN_rbiol2_rtech3", "x097_MIC_LAC_AMYL_rbiol1_rtech3", "x100_MIC_LAC_PLAN_rbiol1_rtech3", "x101_MIC_BIF_BREV_rbiol2_rtech3")


improve_imported_table <- function(imported_table, sample_name_suffix) {
  
  improved_table <- imported_table %>%
    add_column(Name = as.character(NA), .before = "Name_bis.x") %>%
    add_column(mean_intensity = as.numeric(NA), median_intensity = as.numeric(NA), .before = "MSDialName")
  
  for (i in 1:length(pull(improved_table, 1))) {
    improved_table[i, "Name"] <- improved_table[i, "Name_bis.y"]
    
    this_row_intensities <- as.vector(improved_table[i, paste0(sample_names, sample_name_suffix)], mode = "numeric")
    if (length(this_row_intensities)!= 30) {stop("The intensities are not exactly 30")}
    if (any(is.na(this_row_intensities))) {warning("There are still some missing values in the intensities")}
    
    
    improved_table[i, "mean_intensity"] <- mean(this_row_intensities, na.rm = TRUE)
    improved_table[i, "median_intensity"] <- median(this_row_intensities, na.rm = TRUE)
  }
  
  return(improved_table)
}


restuls_BEHPOS <- improve_imported_table(restuls_BEHPOS_imported, sample_name_suffix = "_HILIC_pos")
restuls_C18NEG <- improve_imported_table(restuls_C18NEG_imported, sample_name_suffix = "_RP_neg")
restuls_C18POS <- improve_imported_table(restuls_C18POS_imported, sample_name_suffix = "_RP_pos")


results_list <- list(BEH_POS = restuls_BEHPOS,
                     C18_NEG = restuls_C18NEG,
                     C18_POS = restuls_C18POS)



results_combined <- merge_results(results_list = results_list,
                                  strings_in_colnames_to_remove = c("_HILIC_pos", "_RP_neg", "_RP_pos"),
                                  by_column = "Name",
                                  if_duplicated_consider_columns = c("mean_intensity", "median_intensity"),
                                  decreasing = c(TRUE, TRUE),
                                  name_new_column = "chromatographic_run")


export_the_table(results_combined, "results_combined_known_even_non_sign", "txt")
export_the_table(results_combined, "results_combined_known_even_non_sign", "csv")
export_the_table(results_combined, "results_combined_known_even_non_sign", "xlsx")




