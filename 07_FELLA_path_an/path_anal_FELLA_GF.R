
#

library(tidyverse)
library(GetFeatistics)
library(metaboliteIDmapping) # for converting from CID to KEGG or others
library(FELLA) # for enrichment analyses
library(igraph)


names_samples <- c("x017_MIC_LAC_AMYL_rbiol2_rtech1", "x022_MIC_LAC_KEFI_rbiol2_rtech1", "x026_MIC_BIF_BREV_rbiol2_rtech1", "x028_MIC_LAC_KEFI_rbiol1_rtech1", "x032_MIC_LAC_SALI_rbiol1_rtech1", "x033_MIC_LAC_PLAN_rbiol1_rtech1",
                  "x038_MIC_LAC_AMYL_rbiol1_rtech1", "x041_MIC_LAC_SALI_rbiol2_rtech1", "x042_MIC_BIF_BREV_rbiol1_rtech1", "x043_MIC_LAC_PLAN_rbiol2_rtech1", "x047_MIC_BIF_BREV_rbiol1_rtech2", "x048_MIC_LAC_AMYL_rbiol2_rtech2",
                  "x049_MIC_BIF_BREV_rbiol2_rtech2", "x053_MIC_LAC_PLAN_rbiol2_rtech2", "x057_MIC_LAC_KEFI_rbiol1_rtech2", "x058_MIC_LAC_PLAN_rbiol1_rtech2", "x059_MIC_LAC_KEFI_rbiol2_rtech2", "x067_MIC_LAC_SALI_rbiol1_rtech2",
                  "x068_MIC_LAC_AMYL_rbiol1_rtech2", "x073_MIC_LAC_SALI_rbiol2_rtech2", "x078_MIC_LAC_SALI_rbiol1_rtech3", "x079_MIC_BIF_BREV_rbiol1_rtech3", "x080_MIC_LAC_KEFI_rbiol2_rtech3", "x082_MIC_LAC_KEFI_rbiol1_rtech3",
                  "x087_MIC_LAC_SALI_rbiol2_rtech3", "x090_MIC_LAC_AMYL_rbiol2_rtech3", "x092_MIC_LAC_PLAN_rbiol2_rtech3", "x097_MIC_LAC_AMYL_rbiol1_rtech3", "x100_MIC_LAC_PLAN_rbiol1_rtech3", "x101_MIC_BIF_BREV_rbiol2_rtech3") 


##

results_combined <- read_csv("results_combined.csv")

results_combined_known_even_non_sign <- read_csv("results_combined_known_even_non_sign.csv")




# adding KEGG code: 

results_combined_withcodes <- results_combined %>% 
  add_column(PubChemCID0 = results_combined$PubChemCID,
             HMDBcode = as.character(rep(NA, length(pull(results_combined, 1)))),
             KEGGcode = as.character(rep(NA, length(pull(results_combined, 1)))),
             .after = "PubChemCID") %>%
  separate(PubChemCID0,
           into = paste0("PubChemCID", 1:max(map_int(strsplit(results_combined$PubChemCID, ";"), length))),
           sep = ";",
           fill = "right")





for (i in 1:length(pull(results_combined_withcodes, 1))) {
  
  if (!is.na(pull(results_combined_withcodes, "PubChemCID1")[i]) & pull(results_combined_withcodes, "PubChemCID1")[i] %in% metabolitesMapping$CID) {
    
    results_combined_withcodes[i, "HMDBcode"] <- pull(filter(metabolitesMapping, CID == pull(results_combined_withcodes, "PubChemCID1")[i]), "HMDB")[1]
    results_combined_withcodes[i, "KEGGcode"] <- pull(filter(metabolitesMapping, CID == pull(results_combined_withcodes, "PubChemCID1")[i]), "KEGG")[1]
    
  } else if (!is.na(pull(results_combined_withcodes, "PubChemCID2")[i]) & pull(results_combined_withcodes, "PubChemCID2")[i] %in% metabolitesMapping$CID) {
    
    results_combined_withcodes[i, "HMDBcode"] <- pull(filter(metabolitesMapping, CID == pull(results_combined_withcodes, "PubChemCID2")[i]), "HMDB")[1]
    results_combined_withcodes[i, "KEGGcode"] <- pull(filter(metabolitesMapping, CID == pull(results_combined_withcodes, "PubChemCID2")[i]), "KEGG")[1]
    
  } else if (!is.na(pull(results_combined_withcodes, "PubChemCID3")[i]) & pull(results_combined_withcodes, "PubChemCID3")[i] %in% metabolitesMapping$CID) {
    
    results_combined_withcodes[i, "HMDBcode"] <- pull(filter(metabolitesMapping, CID == pull(results_combined_withcodes, "PubChemCID3")[i]), "HMDB")[1]
    results_combined_withcodes[i, "KEGGcode"] <- pull(filter(metabolitesMapping, CID == pull(results_combined_withcodes, "PubChemCID3")[i]), "KEGG")[1]
    
  } else if (!is.na(pull(results_combined_withcodes, "PubChemCID4")[i]) & pull(results_combined_withcodes, "PubChemCID4")[i] %in% metabolitesMapping$CID) {
    
    results_combined_withcodes[i, "HMDBcode"] <- pull(filter(metabolitesMapping, CID == pull(results_combined_withcodes, "PubChemCID4")[i]), "HMDB")[1]
    results_combined_withcodes[i, "KEGGcode"] <- pull(filter(metabolitesMapping, CID == pull(results_combined_withcodes, "PubChemCID4")[i]), "KEGG")[1]
    
  } else {
    
    print(paste0("CID ", pull(results_combined_withcodes, "PubChemCID")[i], ", ", pull(results_combined_withcodes, "Name")[i]  ,", is not present in the metabolitesMapping library"))
    
  }
}





write_csv(results_combined_withcodes, "results_combined_withcodes.csv")
write_tsv(results_combined_withcodes, "results_combined_withcodes.txt")


####

# FOUND NEW KEGG CODES THAT WERE ABSENT THROUGH METABOANALYST.
# THEREFORE, LOADED THIS NEW MANUALLY UPDATED TABLE:


results_combined_withcodes_metAn <- read_csv("results_combined_withcodes_MetAn.csv")


results_combined_withcodes_fin <- mutate(results_combined_withcodes_metAn,
                                         KEGGcode = ifelse(is.na(KEGGcode_metMapp), KEGGcode_metAn, KEGGcode_metMapp))

write_csv(results_combined_withcodes_fin, "results_combined_withcodes_fin.csv")
write_tsv(results_combined_withcodes_fin, "results_combined_withcodes_fin.txt")

#


Significative_molecules_KEGG_codes <- results_combined_withcodes_fin$KEGGcode[!is.na(results_combined_withcodes_fin$KEGGcode)]



### pathway analysis with the FELLA package con human come database:

## Human:

# hsa	KGB	Homo sapiens (human)

## with the new function:

do_FELLA_enrichment_analysis(organism_code = "hsa", KEGG_codes = Significative_molecules_KEGG_codes,
                             path_databases = "C:/databases/FELLA/", output_prefix = "EA00_", output_suffix = "human_all_GF")





### getting only codes for molecules that are significative in a certain group:

get_sign_met_in_strain <- function(strain) {
  
  if (!strain %in% c("MIC_LAC_AMYL", "MIC_LAC_KEFI", "MIC_LAC_SALI", "MIC_LAC_PLAN", "MIC_BIF_BREV")) {warning("Not exactly a name of a strain!")}
  
  
  rcwff <- filter(results_combined_withcodes_fin,
                  grepl(strain, MIC_LAC_AMYL_vs_MIC_BIF_BREV) | grepl(strain, MIC_LAC_KEFI_vs_MIC_BIF_BREV) | grepl(strain, MIC_LAC_PLAN_vs_MIC_BIF_BREV) | grepl(strain, MIC_LAC_SALI_vs_MIC_BIF_BREV) | grepl(strain, MIC_LAC_KEFI_vs_MIC_LAC_AMYL) | 
                    grepl(strain, MIC_LAC_PLAN_vs_MIC_LAC_AMYL) | grepl(strain, MIC_LAC_SALI_vs_MIC_LAC_AMYL) | grepl(strain, MIC_LAC_PLAN_vs_MIC_LAC_KEFI) | grepl(strain, MIC_LAC_SALI_vs_MIC_LAC_KEFI) | grepl(strain, MIC_LAC_SALI_vs_MIC_LAC_PLAN))
  
  return(rcwff)
  
}





### pathway analysis with the FELLA package using different strains as database:


## MIC_LAC_AMYL = L. Amylovorus SGL14:

# lam	KGB	Lactobacillus amylovorus GRL 1112
# lai	KGB	Lactobacillus amylovorus 30SC
# lay	KGB	Lactobacillus amylovorus GRL1118


sign_MIC_LAC_AMYL <- get_sign_met_in_strain("MIC_LAC_AMYL")
write_tsv(sign_MIC_LAC_AMYL, "sign_MIC_LAC_AMYL.txt")
write_csv(sign_MIC_LAC_AMYL, "sign_MIC_LAC_AMYL.csv")

sign_MIC_LAC_AMYL_KEGGcodes <- sign_MIC_LAC_AMYL$KEGGcode[!is.na(sign_MIC_LAC_AMYL$KEGGcode)]

do_FELLA_enrichment_analysis(organism_code = "lam", KEGG_codes = sign_MIC_LAC_AMYL_KEGGcodes,
                             path_databases = "C:/databases/FELLA/", output_prefix = "EA01", output_suffix = "_MIC_LAC_AMYL_GF")

do_FELLA_enrichment_analysis(organism_code = "lai", KEGG_codes = sign_MIC_LAC_AMYL_KEGGcodes,
                             path_databases = "C:/databases/FELLA/", output_prefix = "EA02", output_suffix = "_MIC_LAC_AMYL_GF")

do_FELLA_enrichment_analysis(organism_code = "lay", KEGG_codes = sign_MIC_LAC_AMYL_KEGGcodes,
                             path_databases = "C:/databases/FELLA/", output_prefix = "EA03", output_suffix = "_MIC_LAC_AMYL_GF")






## MIC_LAC_KEFI = L. Kefiri SGL13:

# lkf	KGB	Lentilactobacillus kefiri

sign_MIC_LAC_KEFI <- get_sign_met_in_strain("MIC_LAC_KEFI")
write_tsv(sign_MIC_LAC_KEFI, "sign_MIC_LAC_KEFI.txt")
write_csv(sign_MIC_LAC_KEFI, "sign_MIC_LAC_KEFI.csv")

sign_MIC_LAC_KEFI_KEGGcodes <- sign_MIC_LAC_KEFI$KEGGcode[!is.na(sign_MIC_LAC_KEFI$KEGGcode)]


do_FELLA_enrichment_analysis(organism_code = "lkf", KEGG_codes = sign_MIC_LAC_KEFI_KEGGcodes,
                             path_databases = "C:/databases/FELLA/", output_prefix = "EA04", output_suffix = "_MIC_LAC_KEFI_GF")





## MIC_LAC_SALI = L. Salivarus SGL03:

# lsl	KGB	Ligilactobacillus salivarius UCC118
# lsi	KGB	Ligilactobacillus salivarius CECT 5713
# lsj	KGB	Ligilactobacillus salivarius JCM1046


sign_MIC_LAC_SALI <- get_sign_met_in_strain("MIC_LAC_SALI")
write_tsv(sign_MIC_LAC_SALI, "sign_MIC_LAC_SALI.txt")
write_csv(sign_MIC_LAC_SALI, "sign_MIC_LAC_SALI.csv")

sign_MIC_LAC_SALI_KEGGcodes <- sign_MIC_LAC_SALI$KEGGcode[!is.na(sign_MIC_LAC_SALI$KEGGcode)]

do_FELLA_enrichment_analysis(organism_code = "lsl", KEGG_codes = sign_MIC_LAC_SALI_KEGGcodes,
                             path_databases = "C:/databases/FELLA/", output_prefix = "EA05", output_suffix = "_MIC_LAC_SALI_GF")

do_FELLA_enrichment_analysis(organism_code = "lsi", KEGG_codes = sign_MIC_LAC_SALI_KEGGcodes,
                             path_databases = "C:/databases/FELLA/", output_prefix = "EA06", output_suffix = "_MIC_LAC_SALI_GF")

do_FELLA_enrichment_analysis(organism_code = "lsj", KEGG_codes = sign_MIC_LAC_SALI_KEGGcodes,
                             path_databases = "C:/databases/FELLA/", output_prefix = "EA07", output_suffix = "_MIC_LAC_SALI_GF")





## MIC_LAC_PLAN = L. Plantarum SGL07

# lpl	KGB	Lactiplantibacillus plantarum WCFS1
# lpj	KGB	Lactiplantibacillus plantarum JDM1
# lpt	KGB	Lactiplantibacillus plantarum ZJ316
# lps	KGB	Lactiplantibacillus plantarum ST-III
# lpr	KGB	Lactiplantibacillus plantarum subsp. plantarum P-8
# lpz	KGB	Lactiplantibacillus plantarum 16	2013
# lpb	KGB	Lactiplantibacillus plantarum B21	2015
# lpx	KGB	Lactiplantibacillus paraplantarum	2016

sign_MIC_LAC_PLAN <- get_sign_met_in_strain("MIC_LAC_PLAN")
write_tsv(sign_MIC_LAC_PLAN, "sign_MIC_LAC_PLAN.txt")
write_csv(sign_MIC_LAC_PLAN, "sign_MIC_LAC_PLAN.csv")

sign_MIC_LAC_PLAN_KEGGcodes <- sign_MIC_LAC_PLAN$KEGGcode[!is.na(sign_MIC_LAC_PLAN$KEGGcode)]

do_FELLA_enrichment_analysis(organism_code = "lpl", KEGG_codes = sign_MIC_LAC_PLAN_KEGGcodes,
                             path_databases = "C:/databases/FELLA/", output_prefix = "EA08", output_suffix = "_MIC_LAC_PLAN_GF")

do_FELLA_enrichment_analysis(organism_code = "lpj", KEGG_codes = sign_MIC_LAC_PLAN_KEGGcodes,
                             path_databases = "C:/databases/FELLA/", output_prefix = "EA09", output_suffix = "_MIC_LAC_PLAN_GF")

do_FELLA_enrichment_analysis(organism_code = "lpt", KEGG_codes = sign_MIC_LAC_PLAN_KEGGcodes,
                             path_databases = "C:/databases/FELLA/", output_prefix = "EA10", output_suffix = "_MIC_LAC_PLAN_GF")

do_FELLA_enrichment_analysis(organism_code = "lps", KEGG_codes = sign_MIC_LAC_PLAN_KEGGcodes,
                             path_databases = "C:/databases/FELLA/", output_prefix = "EA11", output_suffix = "_MIC_LAC_PLAN_GF")

do_FELLA_enrichment_analysis(organism_code = "lpr", KEGG_codes = sign_MIC_LAC_PLAN_KEGGcodes,
                             path_databases = "C:/databases/FELLA/", output_prefix = "EA12", output_suffix = "_MIC_LAC_PLAN_GF")

do_FELLA_enrichment_analysis(organism_code = "lpz", KEGG_codes = sign_MIC_LAC_PLAN_KEGGcodes,
                             path_databases = "C:/databases/FELLA/", output_prefix = "EA13", output_suffix = "_MIC_LAC_PLAN_GF")

do_FELLA_enrichment_analysis(organism_code = "lpb", KEGG_codes = sign_MIC_LAC_PLAN_KEGGcodes,
                             path_databases = "C:/databases/FELLA/", output_prefix = "EA14", output_suffix = "_MIC_LAC_PLAN_GF")

do_FELLA_enrichment_analysis(organism_code = "lpx", KEGG_codes = sign_MIC_LAC_PLAN_KEGGcodes,
                             path_databases = "C:/databases/FELLA/", output_prefix = "EA15", output_suffix = "_MIC_LAC_PLAN_GF")







## MIC_BIF_BREV = B. Breve SGB01

# bbv 	KGB	Bifidobacterium breve ACS-071-V-Sch8b
# bbru	KGB	Bifidobacterium breve UCC2003
# bbre	KGB	Bifidobacterium breve 12L
# bbrv	KGB	Bifidobacterium breve 689b
# bbrj	KGB	Bifidobacterium breve JCM 7017
# bbrc	KGB	Bifidobacterium breve JCM 7019
# bbrn	KGB	Bifidobacterium breve NCFB 2258
# bbrs	KGB	Bifidobacterium breve S27
# bbrd	KGB	Bifidobacterium breve DSM 20213 = JCM 1192


sign_MIC_BIF_BREV <- get_sign_met_in_strain("MIC_BIF_BREV")
write_tsv(sign_MIC_BIF_BREV, "sign_MIC_BIF_BREV.txt")
write_csv(sign_MIC_BIF_BREV, "sign_MIC_BIF_BREV.csv")

sign_MIC_BIF_BREV_KEGGcodes <- sign_MIC_BIF_BREV$KEGGcode[!is.na(sign_MIC_BIF_BREV$KEGGcode)]

do_FELLA_enrichment_analysis(organism_code = "bbv", KEGG_codes = sign_MIC_BIF_BREV_KEGGcodes,
                             path_databases = "C:/databases/FELLA/", output_prefix = "EA16", output_suffix = "_MIC_BIF_BREV_GF")

do_FELLA_enrichment_analysis(organism_code = "bbru", KEGG_codes = sign_MIC_BIF_BREV_KEGGcodes,
                             path_databases = "C:/databases/FELLA/", output_prefix = "EA17", output_suffix = "_MIC_BIF_BREV_GF")

do_FELLA_enrichment_analysis(organism_code = "bbre", KEGG_codes = sign_MIC_BIF_BREV_KEGGcodes,
                             path_databases = "C:/databases/FELLA/", output_prefix = "EA18", output_suffix = "_MIC_BIF_BREV_GF")

do_FELLA_enrichment_analysis(organism_code = "bbrv", KEGG_codes = sign_MIC_BIF_BREV_KEGGcodes,
                             path_databases = "C:/databases/FELLA/", output_prefix = "EA19", output_suffix = "_MIC_BIF_BREV_GF")

do_FELLA_enrichment_analysis(organism_code = "bbrj", KEGG_codes = sign_MIC_BIF_BREV_KEGGcodes,
                             path_databases = "C:/databases/FELLA/", output_prefix = "EA20", output_suffix = "_MIC_BIF_BREV_GF")

do_FELLA_enrichment_analysis(organism_code = "bbrc", KEGG_codes = sign_MIC_BIF_BREV_KEGGcodes,
                             path_databases = "C:/databases/FELLA/", output_prefix = "EA21", output_suffix = "_MIC_BIF_BREV_GF")

do_FELLA_enrichment_analysis(organism_code = "bbrn", KEGG_codes = sign_MIC_BIF_BREV_KEGGcodes,
                             path_databases = "C:/databases/FELLA/", output_prefix = "EA22", output_suffix = "_MIC_BIF_BREV_GF")

do_FELLA_enrichment_analysis(organism_code = "bbrs", KEGG_codes = sign_MIC_BIF_BREV_KEGGcodes,
                             path_databases = "C:/databases/FELLA/", output_prefix = "EA23", output_suffix = "_MIC_BIF_BREV_GF")

do_FELLA_enrichment_analysis(organism_code = "bbrd", KEGG_codes = sign_MIC_BIF_BREV_KEGGcodes,
                             path_databases = "C:/databases/FELLA/", output_prefix = "EA24", output_suffix = "_MIC_BIF_BREV_GF")






## MIC_BIF_BIF (ma NON lo ABBIAMO TRA I CAMPIONI) = B. Bifidum SGB02

# bbp	KGB	Bifidobacterium bifidum PRL2010
# bbi	KGB	Bifidobacterium bifidum S17
# bbf	KGB	Bifidobacterium bifidum BGN4



##### put together pathway analyses for the same strain

combine_enrich_analyses <- function(enrich_tables_file_names, output_prefix, output_suffix) {
  the_object_name <- paste0(output_prefix, "enrich_analysis_run_resultsGraph_", output_suffix)
  the_output_name <- paste0(the_object_name, ".csv")
  
  the_tibble_combined <- read_csv(enrich_tables_file_names[1])[0,]
  
  for (a in enrich_tables_file_names) {
    the_tibble_combined <- bind_rows(the_tibble_combined,
                                     read_csv(a))
  }
  
  the_right_order <- unique(the_tibble_combined$Entry.type)
  
  the_tibble_combined <- the_tibble_combined %>%
    mutate(Entry.type = factor(Entry.type, levels = the_right_order)) %>%
    arrange(Entry.type) %>%
    filter(!duplicated(KEGG.id))
  
  
  assign(the_object_name, the_tibble_combined, envir = .GlobalEnv)
  
  write_csv(the_tibble_combined, the_output_name)
}



### put together MIC_LAC_AMYL:

combine_enrich_analyses(enrich_tables_file_names = c("EA01enrich_analysis_run_resultsGraph_lam_MIC_LAC_AMYL_GF.csv",
                                                     "EA02enrich_analysis_run_resultsGraph_lai_MIC_LAC_AMYL_GF.csv",
                                                     "EA03enrich_analysis_run_resultsGraph_lay_MIC_LAC_AMYL_GF.csv"),
                        output_prefix = "EAC01",
                        output_suffix = "MIC_LAC_AMYL_GF")

  


### put together MIC_LAC_KEFI:

combine_enrich_analyses(enrich_tables_file_names = c("EA04enrich_analysis_run_resultsGraph_lkf_MIC_LAC_KEFI_GF.csv"),
                        output_prefix = "EAC02",
                        output_suffix = "MIC_LAC_KEFI_GF")




### put together MIC_LAC_SALI:

combine_enrich_analyses(enrich_tables_file_names = c("EA05enrich_analysis_run_resultsGraph_lsl_MIC_LAC_SALI_GF.csv",
                                                     "EA06enrich_analysis_run_resultsGraph_lsi_MIC_LAC_SALI_GF.csv",
                                                     "EA07enrich_analysis_run_resultsGraph_lsj_MIC_LAC_SALI_GF.csv"),
                        output_prefix = "EAC03",
                        output_suffix = "MIC_LAC_SALI_GF")






### put together MIC_LAC_PLAN:

combine_enrich_analyses(enrich_tables_file_names = c("EA08enrich_analysis_run_resultsGraph_lpl_MIC_LAC_PLAN_GF.csv",
                                                     "EA09enrich_analysis_run_resultsGraph_lpj_MIC_LAC_PLAN_GF.csv",
                                                     "EA10enrich_analysis_run_resultsGraph_lpt_MIC_LAC_PLAN_GF.csv",
                                                     "EA11enrich_analysis_run_resultsGraph_lps_MIC_LAC_PLAN_GF.csv",
                                                     "EA12enrich_analysis_run_resultsGraph_lpr_MIC_LAC_PLAN_GF.csv",
                                                     "EA13enrich_analysis_run_resultsGraph_lpz_MIC_LAC_PLAN_GF.csv",
                                                     "EA14enrich_analysis_run_resultsGraph_lpb_MIC_LAC_PLAN_GF.csv",
                                                     "EA15enrich_analysis_run_resultsGraph_lpx_MIC_LAC_PLAN_GF.csv"),
                        output_prefix = "EAC04",
                        output_suffix = "MIC_LAC_PLAN_GF")




### put together MIC_BIF_BREV:

combine_enrich_analyses(enrich_tables_file_names = c("EA16enrich_analysis_run_resultsGraph_bbv_MIC_BIF_BREV_GF.csv",
                                                     "EA17enrich_analysis_run_resultsGraph_bbru_MIC_BIF_BREV_GF.csv",
                                                     "EA18enrich_analysis_run_resultsGraph_bbre_MIC_BIF_BREV_GF.csv",
                                                     "EA19enrich_analysis_run_resultsGraph_bbrv_MIC_BIF_BREV_GF.csv",
                                                     "EA20enrich_analysis_run_resultsGraph_bbrj_MIC_BIF_BREV_GF.csv",
                                                     "EA21enrich_analysis_run_resultsGraph_bbrc_MIC_BIF_BREV_GF.csv",
                                                     "EA22enrich_analysis_run_resultsGraph_bbrn_MIC_BIF_BREV_GF.csv",
                                                     "EA23enrich_analysis_run_resultsGraph_bbrs_MIC_BIF_BREV_GF.csv",
                                                     "EA24enrich_analysis_run_resultsGraph_bbrd_MIC_BIF_BREV_GF.csv"),
                        output_prefix = "EAC05",
                        output_suffix = "MIC_BIF_BREV_GF")



## and indicate the significant metabolites that produced them.
