
library(tidyverse)
library(GetFeatistics)

restuls_C8NEG_imported <- read_tsv("metab_stat_w_features_w_diff_more_C8NEG.txt")
restuls_C8POS_imported <- read_tsv("metab_stat_w_features_w_diff_more_C8POS.txt")


results_list <- list(C8_NEG = restuls_C8NEG_imported,
                     C8_POS = restuls_C8POS_imported)



results_combined <- merge_results(results_list = results_list,
                                  strings_in_colnames_to_remove = "",
                                  by_column = "MSDialName",
                                  if_duplicated_consider_columns = NULL,
                                  decreasing = FALSE,
                                  name_new_column = "chromatographic_run")


export_the_table(results_combined, "results_combined_known_even_non_sign", "txt")
export_the_table(results_combined, "results_combined_known_even_non_sign", "csv")
export_the_table(results_combined, "results_combined_known_even_non_sign", "xlsx")




