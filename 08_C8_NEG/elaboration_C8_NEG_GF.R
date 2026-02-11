

#

library(LSD)
library(pals)
library(pcaMethods)
library(ComplexHeatmap)
library(pals)
library(ComplexHeatmap)
library(notame)
library(Hmisc)
library(clusterProfiler)
library(tidyverse)
library(margheRita)


###
## 3 Input data
###
mL <- read_input_file(feature_file = "Normalized_0_202112281320_NoBiFBIFI_OnlyKnown_Reduced_2.txt",sample_file = "meta.txt")



## blank subtraction performed by Gianfranco Frigerio:

if (!all(colnames(mL[["data"]]) == mL[["sample_ann"]][["id"]])) {stop("col names and sample_ann id are not identical!")}

for (a in colnames(mL[["data"]])[grepl("MIC", colnames(mL[["data"]]))]) {
  
  if (grepl("MIC_BIF_BREV", a)) {
    
    blank_samples <- colnames(mL[["data"]])[grepl("BLK_MDE_BREV", colnames(mL[["data"]]))]
    
  } else if (grepl("MIC_LAC_AMYL", a)) {
    
    blank_samples <- colnames(mL[["data"]])[grepl("BLK_TRIL_AMYL", colnames(mL[["data"]]))]
    
  } else if (grepl("MIC_LAC_KEFI", a)) {
    
    blank_samples <- colnames(mL[["data"]])[grepl("BLK_KEF_KEFI", colnames(mL[["data"]]))]
    
  } else if (grepl("MIC_LAC_PLAN", a)) {
    
    blank_samples <- colnames(mL[["data"]])[grepl("BLK_MDE_PLAN", colnames(mL[["data"]]))]
    
  } else if (grepl("MIC_LAC_SALI", a)) {
    
    blank_samples <- colnames(mL[["data"]])[grepl("BLK_MDE_SALI", colnames(mL[["data"]]))]
    
  } else {stop("Missing match with any blanks!")}
  
  
  for (i in 1:length(pull(mL[["data"]], a))) {
    
    mL[["data"]][i,a] <- as.numeric(mL[["data"]][i,a])-mean(as.numeric(mL[["data"]][i, blank_samples]))
    
  }
}

# removing outliers

mL <- remove_samples(mRList = mL, ids = colnames(mL[["data"]])[grepl("BLK", colnames(mL[["data"]])) | grepl("solvent", colnames(mL[["data"]]))])
mL <- remove_samples(mRList = mL, ids = c("x011_QC1_C8_pos","x012_QC2_C8_pos", "100_MIC_LAC_PLAN_rbiol1_rtech3_C8_neg"), column = "id")



write_csv(bind_cols(tibble(feature = row.names(mL[["data"]])), as_tibble(mL[["data"]])), "database_blank_subtracted.csv")


####

###
## 5 Filtering, imputation and normalization
###



heatscatter_chromatography(mL)



mL_norm <- normalize_profiles(mL, method = "pqn")

write_csv(bind_cols(tibble(feature = row.names(mL_norm[["data"]])), as_tibble(mL_norm[["data"]])), "database_blank_subtracted_PQN_normalised.csv")





###
## 6 Principal Component Analysis
###


mL_norm <- mR_pca(mRList = mL_norm, nPcs=5, scaling="pareto", include_QC=TRUE, dirout = "PCA")

Plot2DPCA(mRList = mL_norm, pcx=1, pcy=2, col_by="class", include_QC=TRUE, dirout = "PCA", name_samples = TRUE)



### PCA di nuovo dopo specificando include_QC come FALSE

mL_norm <- mR_pca(mRList = mL_norm, nPcs=5, scaling="pareto", include_QC=FALSE, dirout = "PCA_bis")

Plot2DPCA(mRList = mL_norm, pcx=1, pcy=2, col_by="class", col_pal = pals::trubetskoy(5), include_QC=FALSE, dirout = "PCA_bis", name_samples = TRUE)


### PCA di nuovo dopo specificando include_QC come FALSE NOMI FALSE

mL_norm <- mR_pca(mRList = mL_norm, nPcs=5, scaling="pareto", include_QC=FALSE, dirout = "PCA_tris")

Plot2DPCA(mRList = mL_norm, pcx=1, pcy=2, col_by="class", col_pal = pals::trubetskoy(5), include_QC=FALSE, dirout = "PCA_tris", name_samples = FALSE)





###
##8 Statistical analysis"
###

mean_median_stdev_samples(mL_norm)

mL_norm <- univariate(mL_norm, test_method="anova", exp.levels = c("MIC_BIF_BREV","MIC_LAC_AMYL","MIC_LAC_KEFI","MIC_LAC_PLAN","MIC_LAC_SALI"), exp.factor = "class")



###
## 9 Metabolite identification
###

mR_library <- select_library(source = "MS-DIAL", column = "LipC8", mode = "NEG", accept_RI=10)

mL_norm <- metabolite_identification(mL_norm, library_list = mR_library)

# just to check better the tables:
write.table(mL_norm$metabolite_identification$associations,
            file = "metabolite_identification_associations.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

write.table(mL_norm$metabolite_identification$associations_summary,
            file = "metabolite_identification_associations_summary.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

write.table(mL_norm$metab_ann,
            file = "metabolite_identification_metab_ann.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

write.table(mL_norm$data,
            file = "metabolite_identification_data.txt",
            sep = "\t",
            row.names = TRUE,
            col.names = NA)

write.table(mL_norm$data_ann,
            file = "metabolite_identification_data_ann.txt",
            sep = "\t",
            row.names = TRUE,
            col.names = NA)



#### heatmap with metabolite names

significant_features <- select_sign_features(mRList=mL_norm, test_method="anova", test_value = "q", cutoff_value = 0.05, feature_id="Name", values=FALSE)

png(file = "heat_map.png", width = 200, height = 320, res = 300, units = "mm")
h_map(mL_norm, scale_features=TRUE, column_ann="class", data.use = "data_ann" , col_ann= pals::trubetskoy(5), col=NULL, features = significant_features, top = Inf, show_column_names=FALSE)
dev.off()



###
## 10 Retriving statistics for identified metabolites
###

metab_stat <- annotate_univariate_results(mRList = mL_norm, feature_stats="anova", add.metabolite.levels=TRUE, dirout=NULL)

#

merge_stats_with_features <- function(mRList = NULL, feature_stats=NULL){
  
  if(is.null(feature_stats)){
    stop("feature_stats argument can not be NULL.\n")
  }
  if(!is.data.frame(feature_stats)){
    feature_stats <- mRList[[feature_stats]]
  }
  
  ans <- merge(mRList$metab_ann, feature_stats, by.x = "Feature_ID", by.y=0)
  ans <- ans[, ! colnames(ans) %in% c("rt", "mz", "MS1.isotopic.spectrum", "MS_MS_spectrum", "quality", "reference")]
  
  return(ans)
  
}
metab_stat <- merge_stats_with_features(mRList = mL_norm, feature_stats = "anova")
metab_stat <- left_join(x = metab_stat, y = mL_norm[["metabolite_identification"]][["associations"]][order(as.numeric(mL_norm[["metabolite_identification"]][["associations"]]$Level)),], by = "Feature_ID", suffix = c("", "_bis"), multiple = "first")
write.table(metab_stat,
            file = "metab_stat.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)



## merge metab_stat with features data


data_intensities <- cbind(data.frame(Feature_ID = row.names(mL_norm[["data"]])),
                          as.data.frame(mL_norm[["data"]]))


metab_stat_w_features <- left_join(x = metab_stat, y = data_intensities, by = "Feature_ID")
write.table(metab_stat_w_features,
            file = "metab_stat_w_features.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)



## same but only known:


data_intensities_known <- cbind(data.frame(Name = row.names(mL_norm[["data_ann"]])),
                                as.data.frame(mL_norm[["data_ann"]]))

metab_stat$Level <- as.numeric(metab_stat$Level)
metab_stat <- arrange(metab_stat, Level, q)
colnames(data_intensities_known)[duplicated(colnames(data_intensities_known)) & colnames(data_intensities_known)== "Name"] <- "Name_bis"
metab_stat_w_features_known <- left_join(x = data_intensities_known, y = metab_stat, by = "Name", multiple = "first")
write.table(metab_stat_w_features_known,
            file = "metab_stat_w_features_known.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)




## same but only known and significative:



data_intensities_known_only_sign <- data_intensities_known[which(data_intensities_known$Name %in% significant_features),]



metab_stat_w_features_known_only_sign <- left_join(x = data_intensities_known_only_sign , y = metab_stat, by = "Name", multiple = "first")
write.table(metab_stat_w_features_known_only_sign,
            file = "metab_stat_w_features_known_only_sign.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)




## 2024-08-30 adding differences from ANOVA GF

library(GetFeatistics)

data_per_ANOVA_GF <- transpose_feat_table(metab_stat_w_features[,c("Feature_ID",
                                                                   "x017_MIC_LAC_AMYL_rbiol2_rtech1_C8_neg", "x022_MIC_LAC_KEFI_rbiol2_rtech1_C8_neg", 
                                                                   "x026_MIC_BIF_BREV_rbiol2_rtech1_C8_neg", "x028_MIC_LAC_KEFI_rbiol1_rtech1_C8_neg",
                                                                   "x032_MIC_LAC_SALI_rbiol1_rtech1_C8_neg", "x033_MIC_LAC_PLAN_rbiol1_rtech1_C8_neg",
                                                                   "x038_MIC_LAC_AMYL_rbiol1_rtech1_C8_neg", "x041_MIC_LAC_SALI_rbiol2_rtech1_C8_neg",
                                                                   "x042_MIC_BIF_BREV_rbiol1_rtech1_C8_neg", "x043_MIC_LAC_PLAN_rbiol2_rtech1_C8_neg",
                                                                   "x047_MIC_BIF_BREV_rbiol1_rtech2_C8_neg", "x048_MIC_LAC_AMYL_rbiol2_rtech2_C8_neg",
                                                                   "x049_MIC_BIF_BREV_rbiol2_rtech2_C8_neg", "x053_MIC_LAC_PLAN_rbiol2_rtech2_C8_neg",
                                                                   "x057_MIC_LAC_KEFI_rbiol1_rtech2_C8_neg", "x058_MIC_LAC_PLAN_rbiol1_rtech2_C8_neg",
                                                                   "x059_MIC_LAC_KEFI_rbiol2_rtech2_C8_neg", "x067_MIC_LAC_SALI_rbiol1_rtech2_C8_neg",
                                                                   "x068_MIC_LAC_AMYL_rbiol1_rtech2_C8_neg", "x073_MIC_LAC_SALI_rbiol2_rtech2_C8_neg",
                                                                   "x078_MIC_LAC_SALI_rbiol1_rtech3_C8_neg", "x079_MIC_BIF_BREV_rbiol1_rtech3_C8_neg",
                                                                   "x080_MIC_LAC_KEFI_rbiol2_rtech3_C8_neg", "x082_MIC_LAC_KEFI_rbiol1_rtech3_C8_neg",
                                                                   "x087_MIC_LAC_SALI_rbiol2_rtech3_C8_neg", "x090_MIC_LAC_AMYL_rbiol2_rtech3_C8_neg",
                                                                   "x092_MIC_LAC_PLAN_rbiol2_rtech3_C8_neg", "x097_MIC_LAC_AMYL_rbiol1_rtech3_C8_neg",
                                                                   "x101_MIC_BIF_BREV_rbiol2_rtech3_C8_neg", "x100_MIC_LAC_PLAN_rbiol1_rtech3_C8_neg" )]) %>%
  add_column(group = as.factor((mL_norm[["sample_ann"]][["class"]])), .after = "samples")

colnames(data_per_ANOVA_GF) <- fix_names(colnames(data_per_ANOVA_GF))


ANOVA_GF <- gentab_P.1wayANOVA_posthocTukeyHSD(DF= data_per_ANOVA_GF,
                                               v= colnames(data_per_ANOVA_GF)[3:length(colnames(data_per_ANOVA_GF))],
                                               f = "group",
                                               FDR = TRUE,
                                               groupdiff = TRUE,
                                               pcutoff = 0.05)



metab_stat_w_features_w_diff <- cbind(metab_stat_w_features, ANOVA_GF)

export_the_table(metab_stat_w_features_w_diff, "metab_stat_w_features_w_diff", "txt")

## 
## preparing the table each group vs the others

metab_stat_w_features_w_diff_more <- as_tibble(metab_stat_w_features_w_diff) %>%
  mutate(MIC_LAC_AMYL_vs_others = character(length(pull(metab_stat_w_features_w_diff, 1))),
         MIC_LAC_KEFI_vs_others = character(length(pull(metab_stat_w_features_w_diff, 1))),
         MIC_LAC_PLAN_vs_others = character(length(pull(metab_stat_w_features_w_diff, 1))),
         MIC_LAC_SALI_vs_others = character(length(pull(metab_stat_w_features_w_diff, 1))),
         MIC_BIF_BREV_vs_others = character(length(pull(metab_stat_w_features_w_diff, 1))))

# 

for (i in 1:length(pull(metab_stat_w_features_w_diff_more, 1))) {
  
  # MIC_LAC_AMYL
  
  df_MIC_LAC_AMYL <- metab_stat_w_features_w_diff_more %>%
    select(MIC_LAC_AMYL_vs_MIC_BIF_BREV, MIC_LAC_KEFI_vs_MIC_LAC_AMYL, MIC_LAC_PLAN_vs_MIC_LAC_AMYL, MIC_LAC_SALI_vs_MIC_LAC_AMYL)
  
  df_MIC_LAC_AMYL_fil <- df_MIC_LAC_AMYL[i,]
  
  df_MIC_LAC_AMYL_fil_vect <- as.character(df_MIC_LAC_AMYL_fil)
  
  df_MIC_LAC_AMYL_fil_vect_noNA <- df_MIC_LAC_AMYL_fil_vect[which(!is.na(df_MIC_LAC_AMYL_fil_vect))]
  
  if (length(df_MIC_LAC_AMYL_fil_vect_noNA)>0) {
    for (u in 1:length(df_MIC_LAC_AMYL_fil_vect_noNA)) {
      if (df_MIC_LAC_AMYL_fil_vect_noNA[u] == "MIC_LAC_AMYL > MIC_LAC_KEFI") {df_MIC_LAC_AMYL_fil_vect_noNA[u] <- "> MIC_LAC_KEFI"}
      if (df_MIC_LAC_AMYL_fil_vect_noNA[u] == "MIC_LAC_KEFI > MIC_LAC_AMYL") {df_MIC_LAC_AMYL_fil_vect_noNA[u] <- "< MIC_LAC_KEFI"}
      if (df_MIC_LAC_AMYL_fil_vect_noNA[u] == "MIC_LAC_AMYL > MIC_LAC_PLAN") {df_MIC_LAC_AMYL_fil_vect_noNA[u] <- "> MIC_LAC_PLAN"}
      if (df_MIC_LAC_AMYL_fil_vect_noNA[u] == "MIC_LAC_PLAN > MIC_LAC_AMYL") {df_MIC_LAC_AMYL_fil_vect_noNA[u] <- "< MIC_LAC_PLAN"}
      if (df_MIC_LAC_AMYL_fil_vect_noNA[u] == "MIC_LAC_AMYL > MIC_LAC_SALI") {df_MIC_LAC_AMYL_fil_vect_noNA[u] <- "> MIC_LAC_SALI"}
      if (df_MIC_LAC_AMYL_fil_vect_noNA[u] == "MIC_LAC_SALI > MIC_LAC_AMYL") {df_MIC_LAC_AMYL_fil_vect_noNA[u] <- "< MIC_LAC_SALI"}
      if (df_MIC_LAC_AMYL_fil_vect_noNA[u] == "MIC_LAC_AMYL > MIC_BIF_BREV") {df_MIC_LAC_AMYL_fil_vect_noNA[u] <- "> MIC_BIF_BREV"}
      if (df_MIC_LAC_AMYL_fil_vect_noNA[u] == "MIC_BIF_BREV > MIC_LAC_AMYL") {df_MIC_LAC_AMYL_fil_vect_noNA[u] <- "< MIC_BIF_BREV"}
    }
  }
  
  metab_stat_w_features_w_diff_more[i, "MIC_LAC_AMYL_vs_others"] <- paste0(df_MIC_LAC_AMYL_fil_vect_noNA, collapse = "; ")
  
  
  # MIC_LAC_KEFI
  
  df_MIC_LAC_KEFI <- metab_stat_w_features_w_diff_more %>%
    select(MIC_LAC_KEFI_vs_MIC_BIF_BREV, MIC_LAC_KEFI_vs_MIC_LAC_AMYL, MIC_LAC_PLAN_vs_MIC_LAC_KEFI, MIC_LAC_SALI_vs_MIC_LAC_KEFI)
  
  df_MIC_LAC_KEFI_fil <- df_MIC_LAC_KEFI[i,]
  
  df_MIC_LAC_KEFI_fil_vect <- as.character(df_MIC_LAC_KEFI_fil)
  
  df_MIC_LAC_KEFI_fil_vect_noNA <- df_MIC_LAC_KEFI_fil_vect[which(!is.na(df_MIC_LAC_KEFI_fil_vect))]
  
  if (length(df_MIC_LAC_KEFI_fil_vect_noNA)>0) {
    for (u in 1:length(df_MIC_LAC_KEFI_fil_vect_noNA)) {
      if (df_MIC_LAC_KEFI_fil_vect_noNA[u] == "MIC_LAC_KEFI > MIC_LAC_AMYL") {df_MIC_LAC_KEFI_fil_vect_noNA[u] <- "> MIC_LAC_AMYL"}
      if (df_MIC_LAC_KEFI_fil_vect_noNA[u] == "MIC_LAC_AMYL > MIC_LAC_KEFI") {df_MIC_LAC_KEFI_fil_vect_noNA[u] <- "< MIC_LAC_AMYL"}
      if (df_MIC_LAC_KEFI_fil_vect_noNA[u] == "MIC_LAC_KEFI > MIC_LAC_PLAN") {df_MIC_LAC_KEFI_fil_vect_noNA[u] <- "> MIC_LAC_PLAN"}
      if (df_MIC_LAC_KEFI_fil_vect_noNA[u] == "MIC_LAC_PLAN > MIC_LAC_KEFI") {df_MIC_LAC_KEFI_fil_vect_noNA[u] <- "< MIC_LAC_PLAN"}
      if (df_MIC_LAC_KEFI_fil_vect_noNA[u] == "MIC_LAC_KEFI > MIC_LAC_SALI") {df_MIC_LAC_KEFI_fil_vect_noNA[u] <- "> MIC_LAC_SALI"}
      if (df_MIC_LAC_KEFI_fil_vect_noNA[u] == "MIC_LAC_SALI > MIC_LAC_KEFI") {df_MIC_LAC_KEFI_fil_vect_noNA[u] <- "< MIC_LAC_SALI"}
      if (df_MIC_LAC_KEFI_fil_vect_noNA[u] == "MIC_LAC_KEFI > MIC_BIF_BREV") {df_MIC_LAC_KEFI_fil_vect_noNA[u] <- "> MIC_BIF_BREV"}
      if (df_MIC_LAC_KEFI_fil_vect_noNA[u] == "MIC_BIF_BREV > MIC_LAC_KEFI") {df_MIC_LAC_KEFI_fil_vect_noNA[u] <- "< MIC_BIF_BREV"}
    }
  }
  
  metab_stat_w_features_w_diff_more[i, "MIC_LAC_KEFI_vs_others"] <- paste0(df_MIC_LAC_KEFI_fil_vect_noNA, collapse = "; ")
  
  
  
  # MIC_LAC_PLAN
  
  df_MIC_LAC_PLAN <- metab_stat_w_features_w_diff_more %>%
    select(MIC_LAC_PLAN_vs_MIC_BIF_BREV, MIC_LAC_PLAN_vs_MIC_LAC_AMYL, MIC_LAC_PLAN_vs_MIC_LAC_KEFI, MIC_LAC_SALI_vs_MIC_LAC_PLAN)
  
  df_MIC_LAC_PLAN_fil <- df_MIC_LAC_PLAN[i,]
  
  df_MIC_LAC_PLAN_fil_vect <- as.character(df_MIC_LAC_PLAN_fil)
  
  df_MIC_LAC_PLAN_fil_vect_noNA <- df_MIC_LAC_PLAN_fil_vect[which(!is.na(df_MIC_LAC_PLAN_fil_vect))]
  
  if (length(df_MIC_LAC_PLAN_fil_vect_noNA)>0) {
    for (u in 1:length(df_MIC_LAC_PLAN_fil_vect_noNA)) {
      if (df_MIC_LAC_PLAN_fil_vect_noNA[u] == "MIC_LAC_PLAN > MIC_LAC_AMYL") {df_MIC_LAC_PLAN_fil_vect_noNA[u] <- "> MIC_LAC_AMYL"}
      if (df_MIC_LAC_PLAN_fil_vect_noNA[u] == "MIC_LAC_AMYL > MIC_LAC_PLAN") {df_MIC_LAC_PLAN_fil_vect_noNA[u] <- "< MIC_LAC_AMYL"}
      if (df_MIC_LAC_PLAN_fil_vect_noNA[u] == "MIC_LAC_PLAN > MIC_LAC_KEFI") {df_MIC_LAC_PLAN_fil_vect_noNA[u] <- "> MIC_LAC_KEFI"}
      if (df_MIC_LAC_PLAN_fil_vect_noNA[u] == "MIC_LAC_KEFI > MIC_LAC_PLAN") {df_MIC_LAC_PLAN_fil_vect_noNA[u] <- "< MIC_LAC_KEFI"}
      if (df_MIC_LAC_PLAN_fil_vect_noNA[u] == "MIC_LAC_PLAN > MIC_LAC_SALI") {df_MIC_LAC_PLAN_fil_vect_noNA[u] <- "> MIC_LAC_SALI"}
      if (df_MIC_LAC_PLAN_fil_vect_noNA[u] == "MIC_LAC_SALI > MIC_LAC_PLAN") {df_MIC_LAC_PLAN_fil_vect_noNA[u] <- "< MIC_LAC_SALI"}
      if (df_MIC_LAC_PLAN_fil_vect_noNA[u] == "MIC_LAC_PLAN > MIC_BIF_BREV") {df_MIC_LAC_PLAN_fil_vect_noNA[u] <- "> MIC_BIF_BREV"}
      if (df_MIC_LAC_PLAN_fil_vect_noNA[u] == "MIC_BIF_BREV > MIC_LAC_PLAN") {df_MIC_LAC_PLAN_fil_vect_noNA[u] <- "< MIC_BIF_BREV"}
    }
  }
  
  metab_stat_w_features_w_diff_more[i, "MIC_LAC_PLAN_vs_others"] <- paste0(df_MIC_LAC_PLAN_fil_vect_noNA, collapse = "; ")
  
  
  
  # MIC_LAC_SALI
  
  df_MIC_LAC_SALI <- metab_stat_w_features_w_diff_more %>%
    select(MIC_LAC_SALI_vs_MIC_BIF_BREV, MIC_LAC_SALI_vs_MIC_LAC_AMYL, MIC_LAC_SALI_vs_MIC_LAC_KEFI, MIC_LAC_SALI_vs_MIC_LAC_PLAN)
  
  df_MIC_LAC_SALI_fil <- df_MIC_LAC_SALI[i,]
  
  df_MIC_LAC_SALI_fil_vect <- as.character(df_MIC_LAC_SALI_fil)
  
  df_MIC_LAC_SALI_fil_vect_noNA <- df_MIC_LAC_SALI_fil_vect[which(!is.na(df_MIC_LAC_SALI_fil_vect))]
  
  if (length(df_MIC_LAC_SALI_fil_vect_noNA)>0) {
    for (u in 1:length(df_MIC_LAC_SALI_fil_vect_noNA)) {
      if (df_MIC_LAC_SALI_fil_vect_noNA[u] == "MIC_LAC_SALI > MIC_LAC_AMYL") {df_MIC_LAC_SALI_fil_vect_noNA[u] <- "> MIC_LAC_AMYL"}
      if (df_MIC_LAC_SALI_fil_vect_noNA[u] == "MIC_LAC_AMYL > MIC_LAC_SALI") {df_MIC_LAC_SALI_fil_vect_noNA[u] <- "< MIC_LAC_AMYL"}
      if (df_MIC_LAC_SALI_fil_vect_noNA[u] == "MIC_LAC_SALI > MIC_LAC_KEFI") {df_MIC_LAC_SALI_fil_vect_noNA[u] <- "> MIC_LAC_KEFI"}
      if (df_MIC_LAC_SALI_fil_vect_noNA[u] == "MIC_LAC_KEFI > MIC_LAC_SALI") {df_MIC_LAC_SALI_fil_vect_noNA[u] <- "< MIC_LAC_KEFI"}
      if (df_MIC_LAC_SALI_fil_vect_noNA[u] == "MIC_LAC_SALI > MIC_LAC_PLAN") {df_MIC_LAC_SALI_fil_vect_noNA[u] <- "> MIC_LAC_PLAN"}
      if (df_MIC_LAC_SALI_fil_vect_noNA[u] == "MIC_LAC_PLAN > MIC_LAC_SALI") {df_MIC_LAC_SALI_fil_vect_noNA[u] <- "< MIC_LAC_PLAN"}
      if (df_MIC_LAC_SALI_fil_vect_noNA[u] == "MIC_LAC_SALI > MIC_BIF_BREV") {df_MIC_LAC_SALI_fil_vect_noNA[u] <- "> MIC_BIF_BREV"}
      if (df_MIC_LAC_SALI_fil_vect_noNA[u] == "MIC_BIF_BREV > MIC_LAC_SALI") {df_MIC_LAC_SALI_fil_vect_noNA[u] <- "< MIC_BIF_BREV"}
    }
  }
  
  metab_stat_w_features_w_diff_more[i, "MIC_LAC_SALI_vs_others"] <- paste0(df_MIC_LAC_SALI_fil_vect_noNA, collapse = "; ")
  
  
  
  # MIC_BIF_BREV
  
  df_MIC_BIF_BREV <- metab_stat_w_features_w_diff_more %>%
    select(MIC_LAC_AMYL_vs_MIC_BIF_BREV, MIC_LAC_KEFI_vs_MIC_BIF_BREV, MIC_LAC_PLAN_vs_MIC_BIF_BREV, MIC_LAC_SALI_vs_MIC_BIF_BREV)
  
  df_MIC_BIF_BREV_fil <- df_MIC_BIF_BREV[i,]
  
  df_MIC_BIF_BREV_fil_vect <- as.character(df_MIC_BIF_BREV_fil)
  
  df_MIC_BIF_BREV_fil_vect_noNA <- df_MIC_BIF_BREV_fil_vect[which(!is.na(df_MIC_BIF_BREV_fil_vect))]
  
  if (length(df_MIC_BIF_BREV_fil_vect_noNA)>0) {
    for (u in 1:length(df_MIC_BIF_BREV_fil_vect_noNA)) {
      if (df_MIC_BIF_BREV_fil_vect_noNA[u] == "MIC_BIF_BREV > MIC_LAC_AMYL") {df_MIC_BIF_BREV_fil_vect_noNA[u] <- "> MIC_LAC_AMYL"}
      if (df_MIC_BIF_BREV_fil_vect_noNA[u] == "MIC_LAC_AMYL > MIC_BIF_BREV") {df_MIC_BIF_BREV_fil_vect_noNA[u] <- "< MIC_LAC_AMYL"}
      if (df_MIC_BIF_BREV_fil_vect_noNA[u] == "MIC_BIF_BREV > MIC_LAC_KEFI") {df_MIC_BIF_BREV_fil_vect_noNA[u] <- "> MIC_LAC_KEFI"}
      if (df_MIC_BIF_BREV_fil_vect_noNA[u] == "MIC_LAC_KEFI > MIC_BIF_BREV") {df_MIC_BIF_BREV_fil_vect_noNA[u] <- "< MIC_LAC_KEFI"}
      if (df_MIC_BIF_BREV_fil_vect_noNA[u] == "MIC_BIF_BREV > MIC_LAC_PLAN") {df_MIC_BIF_BREV_fil_vect_noNA[u] <- "> MIC_LAC_PLAN"}
      if (df_MIC_BIF_BREV_fil_vect_noNA[u] == "MIC_LAC_PLAN > MIC_BIF_BREV") {df_MIC_BIF_BREV_fil_vect_noNA[u] <- "< MIC_LAC_PLAN"}
      if (df_MIC_BIF_BREV_fil_vect_noNA[u] == "MIC_BIF_BREV > MIC_LAC_SALI") {df_MIC_BIF_BREV_fil_vect_noNA[u] <- "> MIC_LAC_SALI"}
      if (df_MIC_BIF_BREV_fil_vect_noNA[u] == "MIC_LAC_SALI > MIC_BIF_BREV") {df_MIC_BIF_BREV_fil_vect_noNA[u] <- "< MIC_LAC_SALI"}
    }
  }
  
  metab_stat_w_features_w_diff_more[i, "MIC_BIF_BREV_vs_others"] <- paste0(df_MIC_BIF_BREV_fil_vect_noNA, collapse = "; ")
  
  
}

write_tsv(metab_stat_w_features_w_diff_more, "metab_stat_w_features_w_diff_more.txt")






