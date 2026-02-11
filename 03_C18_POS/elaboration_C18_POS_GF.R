

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
library(openxlsx)
library(tidyverse)
library(margheRita)



#################################################



h_map <- function(mRList=NULL, column_ann="class", data.use = c("data", "data_ann"), col_ann=NULL, col=NULL, scale_features=TRUE, features=NULL, feature_id = "Feature_ID", samples=NULL, top=200, use_top_annoteted_metab = FALSE, test_method="anova", test_value = "q", ...){

  data.use <- match.arg(data.use, c("data", "data_ann"))
  
  data <- mRList[[data.use]]
  stopifnot(length(data)>0)
  if(data.use == "data_ann"){
    
    if (feature_id == "Feature_ID") {
      rownames(data) <- data$Feature_ID
    } else if (feature_id == "Name") {
      rownames(data) <- data$Name
    } else {
      stop('feature_id must be "Feature_ID" or "Name"')
    }
    
    data <- data[, ! colnames(data) %in% c("Feature_ID", "Name")]
  }
  
  if(!is.null(samples)){
    data <- data[, colnames(data) %in% samples]
  }
  stopifnot(length(data)>0)
  
  if(scale_features){
    data<-t(scale(t(data)))
    hm_name <- "z"
  }else{
    hm_name <- "y"
  }
  
  if(is.null(features)){
    if(length(top_var)<top) {top <- length(top_var)}
    top_var <- apply(data, 1, var)
    top_var <- order(top_var, decreasing = T)[1:top]
    data <- data[top_var, ]
  }else{
    if(length(features)>top) {features <- features[1:top]}
    data <- data[rownames(data) %in% features, ]
  }
  
  
  column_ann <- as.factor(mRList$sample_ann[match(colnames(data), rownames(mRList$sample_ann)), column_ann])
  
  if(is.null(col_ann)){
    col_ann <- setNames(polychrome(length(levels(column_ann))), levels(column_ann))
  } else {
    if (length(col_ann) != length(levels(column_ann))) {stop("the length of col_ann must be identical to the length of groups")}
    names(col_ann) <- levels(column_ann)
  }
  
  column_ha <- HeatmapAnnotation(class=column_ann, col = list(class=col_ann))
  
  
  
  if(is.null(col)){
    if(scale_features){
      col <- rev(brewer.rdylbu(7))
    }else{
      col <- brewer.purples(7)
    }
  }
  
  #GF, selecting the best annotation metabolite as a name, also considering the statistics:
  
  if (use_top_annoteted_metab) {
    if (data.use == "data_ann") {
      
      new_row_names_data <- row.names(data)
      
      statistical_test_table <- mRList[[test_method]][order(as.vector(mRList[[test_method]][,test_value])),]
      statistical_test_table[,"Feature_ID"] <- row.names(statistical_test_table)
      statistical_test_table <- statistical_test_table[which(!duplicated(statistical_test_table$Feature_ID)),]
      
      annotation_table <- mL_norm[["metabolite_identification"]][["associations"]][order(as.numeric(as.vector(mL_norm[["metabolite_identification"]][["associations"]][,"Level"]))),]
      annotation_table <- annotation_table[which(!duplicated(annotation_table$Feature_ID)),]
      
      combined_table <- merge(x = statistical_test_table, y = annotation_table, by = "Feature_ID")
      
      
      
      for (i in 1:length(new_row_names_data)) {
        this_feat_id <- as.vector(mRList[[data.use]][,"Feature_ID"])[which(as.vector(mRList[[data.use]][,"Name"]) ==  row.names(data)[i])]
        
        if (length(this_feat_id) ==1) {
          if (this_feat_id %in% combined_table$Feature_ID) {
            new_row_names_data[i] <- combined_table$Name[which(combined_table$Feature_ID==this_feat_id)]
          }
        } else if (length(this_feat_id) >1) {
          this_feat_id <- this_feat_id[which(this_feat_id %in% combined_table$Feature_ID)]
          if (length(this_feat_id)>0) {
            new_row_names_data[i] <- combined_table$Name[which(combined_table$Feature_ID==this_feat_id)][1]
          }
        }
      }
      
      row.names(data) <- new_row_names_data
    }
  }
  
  
  row.names(data) <- sub(";.*", "", row.names(data))
  
  #png(filename = file.path(dirout, "Heatmap.png"), width = 200, height = 200, units = "mm", res=300)
  
  print(Heatmap(as.matrix(data), top_annotation = column_ha, col=col, name = hm_name, ...))
  
  #dev.off()

}


select_sign_features <- function(mRList=NULL, test_method="anova", test_value = "q", cutoff_value = 0.05, feature_id="Feature_ID", values=FALSE, split_by_semicolon = FALSE) {
  
  #data.use <- match.arg(data.use, c("data", "data_ann"))
  #data <- mRList[[data.use]]
  
  stopifnot(feature_id %in% colnames(mRList$data_ann))  
  
  if(values){
    
    ans <- setNames(mRList[[test_method]][, test_value], rownames(mRList[[test_method]]))
    ans <- sort(ans[ans < cutoff_value])
    
    if(feature_id == "Name"){
      
      ans <- merge(data.frame(Feature_ID=names(ans), as.numeric(ans), stringsAsFactors = F), mRList$data_ann[, c("Feature_ID", feature_id)], by="Feature_ID")
      colnames(ans)[2] <- test_value
      
      if (split_by_semicolon) {
        if(any(grepl(";", ans[, feature_id]))){
          unl <- strsplit(ans[, feature_id], ";")
          unl <- data.frame(old=rep(ans[, feature_id], times=lengths(unl)), unlist(unl), stringsAsFactors = F)
          ans <- merge(ans, unl, by.x=3, by.y=1)
          ans <- ans[!is.na(ans[, 4]), ]
          ans <- ans[order(ans[, colnames(ans) == test_value]), ]
          idx_keep <- !duplicated(ans[, 4])
          ans <- setNames(ans[idx_keep, test_value], ans[idx_keep, 4])
        }
      }
    } else if (feature_id != "Name" & feature_id != "Feature_ID") {
      stop('feature_id must be "Feature_ID" or "Name"')
    }
    
    
  }else{
    
    mRList[[test_method]] <- mRList[[test_method]][order(mRList[[test_method]][ , test_value]),]
    ans <- row.names(mRList[[test_method]])[which(mRList[[test_method]][ , test_value] < cutoff_value)]
    
    
    if(feature_id == "Name"){
      
      ans <- as.character(mRList$data_ann[mRList$data_ann$Feature_ID %in% ans, feature_id])
      
      # 1 -> n mapping results in multiple identifiers for each feature
      
      if (split_by_semicolon) {
        if(any(grepl(";", ans))){
          ans <- unlist(strsplit(ans, ";"))
        }
      }
      
      ans <- unique(ans[!is.na(ans)])
      
    } else if (feature_id != "Name" & feature_id != "Feature_ID") {
      stop('feature_id must be "Feature_ID" or "Name"')
    }
    
  }
  
  return(ans)
  
}


merge_stats_w_features <- function(mRList = NULL, test_method = "anova", test_value = "q", cutoff_value = 0.05){
  
  if(is.null(test_method)){
    stop("feature_stats argument can not be NULL.\n")
  }
  if(!is.data.frame(test_method)){
    tab_test_method <- mRList[[test_method]]
  }
  
  ans <- merge(mRList$metab_ann, tab_test_method, by.x = "Feature_ID", by.y=0)
  ans <- ans[, ! colnames(ans) %in% c("rt", "mz", "MS1.isotopic.spectrum", "MS_MS_spectrum", "quality", "reference")]
  
  data_annotation <- as_tibble(mRList[["metabolite_identification"]][["associations"]])
  # order data annotation not only per Level, but also for the other characteristic, to better prioritize:
  data_annotation$Level <- as.numeric(data_annotation$Level)
  data_annotation$mass_status <- factor(data_annotation$mass_status, levels = c("suffer", "acceptable", "super"))
  data_annotation$peaks_found_ppm_RI <- as.numeric(data_annotation$peaks_found_ppm_RI)
  data_annotation$matched_peaks_ratio <- as.numeric(data_annotation$matched_peaks_ratio)
  data_annotation$RT_class <- factor(data_annotation$RT_class, levels = c("unacceptable", "acceptable", "super"))
  data_annotation$ppm_error <- as.numeric(data_annotation$ppm_error)
  data_annotation$RT_err <- as.numeric(data_annotation$RT_err)
  
  data_annotation_ordered <- arrange(data_annotation, Level, desc(mass_status), desc(peaks_found_ppm_RI), desc(matched_peaks_ratio), desc(RT_class), ppm_error, RT_err)
  
  
  
  metab_stat <- left_join(x = ans, y = data_annotation_ordered, by = "Feature_ID", suffix = c("", "_bis"), multiple = "first")
  write.table(metab_stat,
              file = "metab_stat.txt",
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE)
  assign("metab_stat", metab_stat, envir = .GlobalEnv)
  
  
  ## merge metab_stat with features data
  
  
  data_intensities <- cbind(data.frame(Feature_ID = row.names(mRList[["data"]])),
                            as.data.frame(mRList[["data"]]))
  
  
  metab_stat_w_features <- left_join(x = metab_stat, y = data_intensities, by = "Feature_ID")
  write.table(metab_stat_w_features,
              file = "metab_stat_w_features.txt",
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE)
  assign("metab_stat_w_features", metab_stat_w_features, envir = .GlobalEnv)
  
  
  
  
  ## same but only known:
  
  
  data_intensities_known <- cbind(data.frame(Name = row.names(mRList[["data_ann"]])),
                                  as.data.frame(mRList[["data_ann"]]))
  colnames(data_intensities_known)[duplicated(colnames(data_intensities_known)) & colnames(data_intensities_known)== "Name"] <- "Name_bis"
  
  
  metab_stat$Level <- as.numeric(metab_stat$Level)
  metab_stat[,test_value] <- as.numeric(pull(metab_stat, test_value))
  
  metab_stat <- arrange(metab_stat, Level, !!sym(test_value))
  metab_stat_w_features_known <- left_join(x = data_intensities_known, y = metab_stat, by = "Feature_ID", multiple = "first")
  write.table(metab_stat_w_features_known,
              file = "metab_stat_w_features_known.txt",
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE)
  assign("metab_stat_w_features_known", metab_stat_w_features_known, envir = .GlobalEnv)
  
  
  ## same but only significative:
  
  sign_feat_featID <- select_sign_features(mRList = mRList, test_method = test_method, test_value = test_value, cutoff_value = cutoff_value, feature_id = "Feature_ID")
  
  data_intensities_only_sign <- data_intensities[which(data_intensities$Feature_ID %in% sign_feat_featID),]
  
  
  metab_stat_w_features_only_sign <- left_join(x = data_intensities_only_sign , y = metab_stat, by = "Feature_ID", multiple = "first")
  write.table(metab_stat_w_features_only_sign,
              file = "metab_stat_w_features_only_sign.txt",
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE)
  assign("metab_stat_w_features_only_sign", metab_stat_w_features_only_sign, envir = .GlobalEnv)
  
  
  
  ## same but only known and significative:
  
  
  data_intensities_known_only_sign <- data_intensities_known[which(data_intensities_known$Feature_ID %in% data_intensities_only_sign$Feature_ID),]
  
  
  metab_stat_w_features_known_only_sign <- left_join(x = data_intensities_known_only_sign , y = metab_stat, by = "Feature_ID", multiple = "first")
  write.table(metab_stat_w_features_known_only_sign,
              file = "metab_stat_w_features_known_only_sign.txt",
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE)
  assign("metab_stat_w_features_known_only_sign", metab_stat_w_features_known_only_sign, envir = .GlobalEnv)
  
  
}


filter_metabolite_associations <- function(mRList=NULL){
  
  out_levels <- mRList$metabolite_identification$associations
  
  cat("Initial associations:", nrow(out_levels), "\n")
  
  out_levels$ID_Feature_ID <- apply(out_levels[, c("ID", "Feature_ID")], 1, paste0, collapse="_")
  out_levels$ID_Feature_ID <- gsub(" +", "", out_levels$ID_Feature_ID)
  
  
  ### 1) To resolve the one-to-many associations between (f_i,m_j) pairs and MS/MS spectra s_jk, for each pair (f_i, m_j), select the s_jk that has:
  ### - highest number of peaks;
  ### - highest peaks ratio;
  cat("1) Filtering redundant MS/MS spectra for the same metabolite-feature...\n")
  cat("\tby highest number of peaks\n")
  cat("\tby highest peaks ratio...")
  
  out_levels_NA <- out_levels[is.na(out_levels$peaks_found_ppm_RI), ]
  out_levels_tofilter <- out_levels[!is.na(out_levels$peaks_found_ppm_RI), ]
  
  out_levels_tofilter <- out_levels_tofilter[order(-out_levels_tofilter$peaks_found_ppm_RI, -out_levels_tofilter$matched_peaks_ratio), ]
  
  idx_rm <- which(duplicated(out_levels_tofilter$ID_Feature_ID))
  if(length(idx_rm) > 0){
    #cat("Removing", length(idx_rm), "\n")
    out_levels_tofilter <- out_levels_tofilter[-idx_rm, ]
  }
  cat(nrow(out_levels_NA) + nrow(out_levels_tofilter), "\n")
  
  #2) To resolve the one-to-many associations between m_j and f_i, for each m_j, select the f_i that has:
  #a.	the highest ratio between #{feature MS/MS peaks that match metabolite MS/MS peaks} and #{metabolite MS/MS peaks};
  #b.	the best annotation level;
  #c.	the lowest ppm error
  cat("2) 1:n associations between metabolites and features...\n")
  cat("\tby highest peaks ratio\n")
  cat("\tby best annotaiton level\n")
  cat("\tby lowest ppm error...")
  
  out_levels_tofilter <- out_levels_tofilter[order(-out_levels_tofilter$matched_peaks_ratio, out_levels_tofilter$Level, out_levels_tofilter$ppm_error), ]
  
  idx_rm <- which(duplicated(out_levels_tofilter$ID))
  if(length(idx_rm) > 0){
    #cat("Removing", length(idx_rm), "\n")
    out_levels_tofilter <- out_levels_tofilter[-idx_rm, ]
  }
  out_levels <- rbind(out_levels_tofilter, out_levels_NA)
  out_levels$ID_Feature_ID <- NULL
  
  cat(nrow(out_levels), "\n")
  
  #	3 To resolve the one-to-many associations between m_j and f_i not supported by MS/MS spectra, for each m_j, select the f_i that has: 
  # the best annotation level;
  # the best mass match;
  # the best retention time match
  # Lowest rt error
  # Lowest ppm error
  
  cat("3) 1:n associations between metabolites and features, not supported by MS/MS...\n")
  # cat("\tby best annotaiton level\n")
  # cat("\tby best ppm match\n")
  # cat("\tby best RT match\n")
  # cat("\tby lowest RT error\n")
  # cat("\tby lowest ppm error\n")
  
  ####### by LEVEL 	######
  cat("\tby level... ")
  for(lev in c(1, 2)){
    lev_m <- out_levels$ID[out_levels$Level==lev]
    if(length(lev_m) > 0){
      idx_rm <- which(out_levels$ID %in% lev_m & out_levels$Level > lev)
      if(length(idx_rm) > 0){
        #cat(unique(out_levels$Level[idx_rm]), "\n")
        out_levels <- out_levels[-idx_rm, ]
      }
    }
  }
  
  lev_m <- out_levels$ID[out_levels$Level==3 & out_levels$Level_note ==  "mz, rt"]
  if(length(lev_m) > 0){
    idx_rm <- which(out_levels$ID %in% lev_m & out_levels$Level_note == "mz")
    if(length(idx_rm) > 0){
      #cat(unique(out_levels$Level[idx_rm]), "\n")
      out_levels <- out_levels[-idx_rm, ]
    }
  }
  cat(nrow(out_levels), "\n")
  
  ####### by mass status 	######
  if(any(duplicated(out_levels$ID))){
    cat("\tby mass (category)... ")
    mass_status_class <- as.numeric(factor(out_levels$mass_status, levels = c("super", "acceptable", "suffer")))
    ms_num_set <- sort(unique(mass_status_class))
    if(length(ms_num_set)>1){
      for(ms in ms_num_set[-length(ms_num_set)]){ #the highest value does have higher values
        feat_selected <- out_levels$ID[mass_status_class==ms] #IDs that have level ms
        if(length(feat_selected) > 0){
          idx_rm <- which(out_levels$ID %in% feat_selected & mass_status_class > ms) #remove mass_status higher than ms for the same IDs
          if(length(idx_rm) > 0){
            #cat(unique(out_levels$Level[idx_rm]), "\n")
            out_levels <- out_levels[-idx_rm, ]
            mass_status_class <- mass_status_class[-idx_rm]
          }
        }
      }
    }
  }
  cat(nrow(out_levels), "\n")
  
  ##### RT_class #####
  if(any(duplicated(out_levels$ID))){
    
    cat("\tby RT (category)... ")
    
    rt_class <- as.numeric(factor(out_levels$RT_class, levels = c("super", "acceptable")))
    rt_class_set <- sort(unique(rt_class))
    
    if(length(rt_class_set)>1){
      for(rt_i in rt_class_set[-length(rt_class_set)]){ #the highest value does have higher values
        feat_selected <- out_levels$ID[rt_class==rt_i] #IDs that have level ms
        if(length(feat_selected) > 0){
          idx_rm <- which(out_levels$ID %in% feat_selected & rt_class > rt_i) #remove mass_status higher than ms for the same IDs
          if(length(idx_rm) > 0){
            #cat(unique(out_levels$Level[idx_rm]), "\n")
            out_levels <- out_levels[-idx_rm, ]
            rt_class <- rt_class[-idx_rm]
          }
        }
      }
    }
  }
  
  #duplicated IDs
  m_id <- unique(out_levels$ID[duplicated(out_levels$ID)])
  
  out_levels$selected <- TRUE
  
  for(i in 1:length(m_id)){
    
    #by lowest RT error
    if(sum(out_levels$selected[out_levels$ID == m_id[i]])>1){
      
      rt_err_min <- min(out_levels$RT_err[out_levels$ID == m_id[i] & out_levels$selected])
      
      if(any(out_levels$RT_class[out_levels$ID ==  m_id[i] & out_levels$selected] != "unacceptable")){
        
        if(any(abs(out_levels$RT_err[out_levels$ID ==  m_id[i] & out_levels$selected] - rt_err_min) > 1e-4)){
          
          cat("\tby RT... ")
          idx_rm <- which(out_levels$ID == m_id[i] & out_levels$selected & abs(out_levels$RT_err - rt_err_min) > 1e-4)
          out_levels$selected[idx_rm] <- FALSE
          #out_levels <- out_levels[-idx_rm, ]
          #print(out_levels[out_levels$ID == m_id[i], ])
          cat(sum(out_levels$selected), "\n")
          
        }
        
      }
    }
    
    #by lowest ppm error
    if(sum(out_levels$selected[out_levels$ID == m_id[i]])>1){
      
      ppm_err_min <- min(out_levels$ppm_error[out_levels$ID == m_id[i] & out_levels$selected])
      
      if(any(abs(out_levels$ppm_error[out_levels$ID ==  m_id[i] & out_levels$selected] - ppm_err_min) > 1e-4)){
        
        cat("\tby ppm... ")
        idx_rm <- which(out_levels$ID == m_id[i] & out_levels$selected & abs(out_levels$ppm_error - ppm_err_min) > 1e-4)
        out_levels$selected[idx_rm] <- FALSE
        #out_levels <- out_levels[-idx_rm, ]
        #print(out_levels[out_levels$ID == m_id[i], ])
        cat(sum(out_levels$selected), "\n")
        
      }
    }
    
    
  }
  
  out_levels <- out_levels[out_levels$selected, ]
  out_levels$selected <- NULL
  
  #cat(nrow(out_levels), "\n")
  
  # To resolve the one-to-many associations between f_i and m_j, for each f_i keep the m_j that has:
  #   the best level;
  # the best mass match;
  # the highest number of matching MS/MS peaks;
  # the best retention time.
  # the highest ratio between #{feature MS/MS peaks that match metabolite MS/MS peaks} and #{metabolite MS/MS peaks}
  # Lowest rt error
  # Lowest ppm error
  
  cat("4) 1:n associations between features and metabolites...\n")
  # cat("\tby best annotaiton level\n")
  # cat("\tby best ppm match\n")
  # cat("\tby highest number of matching MS/MS peaks\n")
  # cat("\tby best RT match\n")
  # cat("\tby highest ratio of MS/MS peaks\n")
  # cat("\tby lowest RT error\n")
  # cat("\tby lowest ppm error\n")
  
  ####### by LEVEL 	######
  cat("\tby level... ")
  for(lev in c(1, 2)){
    lev_m <- out_levels$Feature_ID[out_levels$Level==lev]
    if(length(lev_m) > 0){
      idx_rm <- which(out_levels$Feature_ID %in% lev_m & out_levels$Level > lev)
      if(length(idx_rm) > 0){
        #cat(unique(out_levels$Level[idx_rm]), "\n")
        out_levels <- out_levels[-idx_rm, ]
      }
    }
  }
  
  #level 3 mz,rt annotated also as level 3 mz
  lev_m <- out_levels$Feature_ID[out_levels$Level==3 & out_levels$Level_note ==  "mz, rt"]
  if(length(lev_m) > 0){
    idx_rm <- which(out_levels$Feature_ID %in% lev_m & out_levels$Level_note == "mz")
    if(length(idx_rm) > 0){
      #cat(unique(out_levels$Level[idx_rm]), "\n")
      out_levels <- out_levels[-idx_rm, ]
    }
  }
  cat(nrow(out_levels), "\n")
  
  
  ####### by mass status 	######
  if(any(duplicated(out_levels$Feature_ID))){
    cat("\tby mass (category)... ")
    mass_status_class <- as.numeric(factor(out_levels$mass_status, levels = c("super", "acceptable", "suffer")))
    ms_num_set <- sort(unique(mass_status_class))
    if(length(ms_num_set)>1){
      for(ms in ms_num_set[-length(ms_num_set)]){ #the highest value does have higher values
        feat_selected <- out_levels$Feature_ID[mass_status_class==ms] #IDs that have level ms
        if(length(feat_selected) > 0){
          idx_rm <- which(out_levels$Feature_ID %in% feat_selected & mass_status_class > ms) #remove mass_status higher than ms for the same IDs
          if(length(idx_rm) > 0){
            #cat(unique(out_levels$Level[idx_rm]), "\n")
            out_levels <- out_levels[-idx_rm, ]
            mass_status_class <- mass_status_class[-idx_rm]
          }
        }
      }
    }
  }
  cat(nrow(out_levels), "\n")
  
  ####### by Peaks_found_ppm_RI ####
  if(any(duplicated(out_levels$Feature_ID))){
    cat("\tby number of matching MS/MS peaks... ")
    num_peaks <- out_levels$peaks_found_ppm_RI
    num_peaks[is.na(num_peaks)] <- 0
    num_peaks_set <- unique(sort(num_peaks, decreasing = T))
    if(any(num_peaks>0)){
      for(ps in num_peaks_set[-length(num_peaks_set)]){ #the highest value does have less peaks
        feat_selected <- out_levels$Feature_ID[num_peaks==ps] #IDs with ps peaks
        if(length(feat_selected) > 0){
          idx_rm <- which(out_levels$Feature_ID %in% feat_selected & num_peaks < ps) #remove mass_status higher than ms for the same IDs >
          if(length(idx_rm) > 0){
            #cat(unique(out_levels$Level[idx_rm]), "\n")
            out_levels <- out_levels[-idx_rm, ]
            num_peaks <- num_peaks[-idx_rm]
          }
        }
      }
    }
  }
  cat(nrow(out_levels), "\n")
  
  if(any(duplicated(out_levels$Feature_ID))){
    
    ##### RT_class #####
    cat("\tby RT (category)... ")
    rt_class <- as.numeric(factor(out_levels$RT_class, levels = c("super", "acceptable")))
    rt_class_set <- sort(unique(rt_class))
    
    if(length(rt_class_set)>1){
      for(rt_i in rt_class_set[-length(rt_class_set)]){ #the highest value does have higher values
        feat_selected <- out_levels$Feature_ID[rt_class==rt_i] #IDs that have level ms
        if(length(feat_selected) > 0){
          idx_rm <- which(out_levels$Feature_ID %in% feat_selected & rt_class > rt_i) #remove mass_status higher than ms for the same IDs
          if(length(idx_rm) > 0){
            #cat(unique(out_levels$Level[idx_rm]), "\n")
            out_levels <- out_levels[-idx_rm, ]
            rt_class <- rt_class[-idx_rm]
          }
        }
      }
    }
  }
  cat(nrow(out_levels), "\n")
  
  f_id <- unique(out_levels$Feature_ID[duplicated(out_levels$Feature_ID)])
  
  out_levels$selected <- TRUE
  
  for(i in 1:length(f_id)){
    
    
    if(sum(out_levels$selected[out_levels$Feature_ID == f_id[i]])>1){
      
      ### by matched_peaks_ratio
      
      if(any(!is.na(out_levels$matched_peaks_ratio[out_levels$Feature_ID == f_id[i]]))){
        pr_max <- max(out_levels$matched_peaks_ratio[out_levels$Feature_ID == f_id[i]], na.rm = T)
        
        if(!is.na(pr_max)){
          
          if(any(abs(out_levels$matched_peaks_ratio[out_levels$Feature_ID ==  f_id[i]] - pr_max) > 1e-4)){
            
            cat("\tby matched peak ratio... ")
            idx_rm <- which(out_levels$Feature_ID == f_id[i] & abs(out_levels$matched_peaks_ratio - pr_max) > 1e-4)
            out_levels$selected[idx_rm] <- FALSE
            cat(sum(out_levels$selected), "\n")
            
          }
          
        }
      }
      
    }
    
    if(sum(out_levels$selected[out_levels$Feature_ID == f_id[i]])>1){
      
      rt_err_min <- min(out_levels$RT_err[out_levels$Feature_ID == f_id[i] & out_levels$selected])
      
      if(any(out_levels$RT_class[out_levels$Feature_ID ==  f_id[i] & out_levels$selected] != "unacceptable")){
        
        if(any(abs(out_levels$RT_err[out_levels$Feature_ID ==  f_id[i] & out_levels$selected] - rt_err_min) > 1e-4)){
          
          cat("\tby RT... ")
          idx_rm <- which(out_levels$Feature_ID == f_id[i] & out_levels$selected & abs(out_levels$RT_err - rt_err_min) > 1e-4)
          out_levels$selected[idx_rm] <- FALSE
          #out_levels <- out_levels[-idx_rm, ]
          #print(out_levels[out_levels$Feature_ID == f_id[i], ])
          cat(sum(out_levels$selected), "\n")
          
        }
        
      }
    }
    
    
    if(sum(out_levels$selected[out_levels$Feature_ID == f_id[i]])>1){
      
      
      ppm_err_min <- min(out_levels$ppm_error[out_levels$Feature_ID == f_id[i] & out_levels$selected])
      
      if(any(abs(out_levels$ppm_error[out_levels$Feature_ID ==  f_id[i] & out_levels$selected] - ppm_err_min) > 1e-4)){
        
        cat("\tby ppm... ")
        idx_rm <- which(out_levels$Feature_ID == f_id[i] & out_levels$selected & abs(out_levels$ppm_error - ppm_err_min) > 1e-4)
        out_levels$selected[idx_rm] <- FALSE
        #out_levels <- out_levels[-idx_rm, ]
        #print(out_levels[out_levels$Feature_ID == f_id[i], ])
        cat(sum(out_levels$selected), "\n")
        
      }
    }
    
    
  }
  
  out_levels <- out_levels[out_levels$selected, ]
  out_levels$selected <- NULL
  
  mRList$metabolite_identification$associations <- out_levels
  
  mRList$metabolite_identification$associations_summary <- unique(out_levels[, c("Feature_ID", "ID", "Name", "Level", "Level_note")])
  
  return(mRList)
  
}


metabolite_identification <- function(mRList = NULL, features = NULL, library_list = NULL, rt_err=1, rt_best_thr=0.5, unaccept_flag=20, accept_flag=5, suffer_flag=10, min_RI = 10, ppm_err = 20, RI_err=20, RI_err_type="rel", filter_ann=FALSE, lib_ann_fields=c("ID", "Name", "PubChemCID"), dirout=NULL){
  
  if(is.null(features)){
    idx <- 1:nrow(mRList$metab_ann) 
  }else{
    idx <- which(mRList$metab_ann$Feature_ID %in% features) 
  }
  
  if (!is.null(dirout)) {
    dir.create(dirout, showWarnings = F)
  }else{
    dirout <- getwd()
  }
  
  
  cat("# Features considered: ", length(idx), ".\n")
  
  #extract feature information from metabolite metadata
  out_levels <- unique(mRList$metab_ann[idx, c("Feature_ID", "rt", "mz", "MS_MS_spectrum")])
  out_levels$MS_MS_spectrum[out_levels$MS_MS_spectrum == ""] <- NA
  
  ### check retention time and mass
  if(length(library_list$lib_precursor$rt)>0){
    cat("Checking RT... ")
    RT <-  check_RT(feature_data = out_levels, reference = library_list$lib_precursor, rt_err_thr=rt_err, rt_best_thr = rt_best_thr)
    cat("done\n")
  }
  
  cat("Checking precursor mz... ")
  mass <- check_mass(feature_data = out_levels, reference = library_list$lib_precursor, unaccept_flag=unaccept_flag, accept_flag=accept_flag, suffer_flag=suffer_flag)
  cat("done\n")
  
  if(length(library_list$lib_precursor$rt)>0){
    cat("Integrating RT and mz results... ")
    RT_mass <- check_RT_mass(RT = RT, mass = mass)
    cat("done\n")
  }
  
  cat("Assembling output for precursors...")
  ### features with appropriate mass
  mass_ok <- stack(lapply(mass, function(x) x$Feature_ID))
  mass_ok <- merge(library_list$lib_precursor, mass_ok, by.x="ID", by.y="ind", all.y=T)
  colnames(mass_ok) <- replace(colnames(mass_ok), list = colnames(mass_ok)=="values", "Feauture_ID")
  
  out_levels <- merge(out_levels, data.frame(ID=rep(names(mass), unlist(lapply(mass, nrow))), do.call(rbind, lapply(mass, function(x) x[, colnames(x) != "mz"])), stringsAsFactors = F), by="Feature_ID", all=T)
  out_levels <- merge(library_list$lib_precursor, out_levels, by="ID", all = T, suffixes = c("_lib", ""))
  
  ### features with appropriate mass and RT
  if(length(library_list$lib_precursor$rt)>0){
    
    RT_mass_ok <- stack(lapply(RT_mass, function(x) x$Feature_ID))
    RT_mass_ok <- merge(library_list$lib_precursor, RT_mass_ok, by.x="ID", by.y="ind", all.y=T)
    
    out_levels <- merge(out_levels, data.frame(ID=rep(names(RT_mass), unlist(lapply(RT_mass, nrow))), do.call(rbind, lapply(RT_mass, function(x) x[, c("Feature_ID", "rt", "RT_err", "RT_class", "RT_flag")])), stringsAsFactors = F), by=c("ID", "Feature_ID", "rt"), all=T)
    
  }else{
    out_levels$RT_err <- NA
    out_levels$RT_class <- NA
    out_levels$RT_flag <- NA
    out_levels$rt_lib <- NA
  }
  
  cat("done\n")
  
  
  ###MS/MS matching - features that match at mass level and have MSMSspectra
  cat("MS/MS analysis...\n")
  feature_spectra_list <- unique(out_levels$Feature_ID[which(!is.na(out_levels$MS_MS_spectrum) & out_levels$mass_flag)]) ## features already associated at least with mass and MS/MS available
  
  feature_spectra_list <- get_spectra_list_from_vector(spectra = setNames(mRList$metab_ann$MS_MS_spectrum[mRList$metab_ann$Feature_ID %in% feature_spectra_list], mRList$metab_ann$Feature_ID[mRList$metab_ann$Feature_ID %in% feature_spectra_list])) #list of features and their spectra
  
  RI_sample <- calc_RI(feature_spectra=feature_spectra_list, accept_RI = min_RI) #define RI
  
  #remove from mass those features that do not appear in RI_sample, because do not have MSMS spectra
  mass <- lapply(mass, function(x) x[x$Feature_ID %in% names(RI_sample), ])
  mass <- mass[unlist(lapply(mass, nrow)) > 0]
  
  intense_peak <- peak_matching(library_list = library_list, library_matched_features = mass, RI_sample = RI_sample, ppm_err = ppm_err, intensity = RI_err, RI_diff_type = RI_err_type)
  cat("done\n")
  
  
  cat("Assembling final output... ")
  
  if(is.null(intense_peak$matched_peaks)) {
    intense_peak$matched_peaks <- data.frame(ID = character(),
                                             Feature_ID = character(),
                                             mz = numeric(),
                                             ppm_error = numeric(),
                                             mass_flag = logical(),
                                             mass_status = character(),
                                             ID_peaks = character(),
                                             peaks_found_ppm_RI = numeric(),
                                             matched_peaks_ratio = numeric(),
                                             precursor_in_MSMS = logical(),
                                             CAS = character(),
                                             Name = character(),
                                             rt = numeric(),
                                             mz_lib = numeric(),
                                             PubChemCID = character(),
                                             SMILES = character())
  }
  
  intense_peak_unq <- unique(intense_peak$matched_peaks[, c("ID", "Feature_ID")])
  intense_peak_unq <- merge(library_list$lib_precursor, intense_peak_unq, by="ID", all.y=T)
  
  
  out_levels <- merge(out_levels, intense_peak$matched_peaks[, c("ID", "Feature_ID", "ID_peaks", "peaks_found_ppm_RI", "matched_peaks_ratio", "precursor_in_MSMS")], by=c("ID", "Feature_ID"), all=T)
  
  out_levels$Level <- ""
  out_levels$Level_note <- ""
  
  idx <- which(out_levels$mass_flag)
  out_levels$Level[idx] <- 3
  out_levels$Level_note[idx] <- "mz"
  
  idx <- which(out_levels$mass_flag & out_levels$RT_flag)
  out_levels$Level[idx] <- 3
  out_levels$Level_note[idx] <- "mz, rt"
  
  idx <- which(!is.na(out_levels$peaks_found_ppm_RI))
  out_levels$Level[idx] <- 2
  out_levels$Level_note[idx] <- ""
  
  idx <- which(!is.na(out_levels$peaks_found_ppm_RI) & out_levels$mass_flag & out_levels$RT_flag)
  out_levels$Level[idx] <- 1
  out_levels$Level_note[idx] <- ""
  
  
  ### overall Reliability
  #out_levels$Reliability <- ""
  #out_levels$Reliability
  
  ### summary
  out_levels <- out_levels[out_levels$Level != "", ]
  
  #out_levels <- out_levels[, c("Feature_ID", "rt", "mz", "MS_MS_spectrum", "RT_err", "RT_flag", "ppm_error", "mass_flag", "mass_status", "ID_peaks", "peaks_found_ppm_RI", "precursor_in_MSMS", "ID", "Name", "rt_lib", "mz_lib", "Level", "Level_note", "CAS", "PubChemCID")]
  out_levels <- out_levels[, c("Feature_ID", "rt", "mz", "RT_err", "RT_class", "RT_flag", "ppm_error", "mass_flag", "mass_status", "ID_peaks", "peaks_found_ppm_RI", "matched_peaks_ratio", "precursor_in_MSMS", "ID", "Name", "rt_lib", "mz_lib", "SMILES", "Level", "Level_note")]
  
  #add the results to mRList
  mRList$metabolite_identification <- list(associations=out_levels, RI_sample=RI_sample, MS_MS_info=intense_peak$matrices)
  
  ## FILTERING
  #metabolite-feature association of level 1 imply that these do not appear in other pairs, and so on for other levels
  if(filter_ann){
    mRList <- filter_metabolite_associations(mRList)
  }else{
    ### summary associations and Feature ID assignment
    #summary associations
    mRList$metabolite_identification$associations_summary <- unique(out_levels[, c("Feature_ID", "ID", "Name", "Level", "Level_note")])
  }
  
  
  
  cat("Adding annotations to metab_ann...\n")
  
  if(any(lib_ann_fields %in% colnames(mRList$metab_ann))){
    
    cat("Any of ", lib_ann_fields, "is already in metab_ann table; these columns will be replaced.\n")
    mRList$metab_ann[, which(colnames(mRList$metab_ann) %in% lib_ann_fields)] <- NULL
    
  }
  #add name to metab_ann
  feat_name <- split(mRList$metabolite_identification$associations_summary$ID, mRList$metabolite_identification$associations_summary$Feature_ID)
  
  for(i_col in lib_ann_fields){
    mapped_feat <- as.data.frame(unlist(lapply(feat_name, function(x) paste0(sort(library_list$lib_precursor[match(x, library_list$lib_precursor$ID), i_col]), collapse = ";"))))
    colnames(mapped_feat) <- i_col
    mRList$metab_ann <- merge(mRList$metab_ann, mapped_feat, by.x="Feature_ID", by.y=0, all.x = T, sort = F)
  }
  
  mRList$metab_ann <- mRList$metab_ann[match(rownames(mRList$data), mRList$metab_ann$Feature_ID), ]
  
  ###annotated data
  mRList$data_ann <- merge(unique(mRList$metab_ann[!is.na(mRList$metab_ann$Name), c("Feature_ID", "Name")]), mRList$data, by.x=1, by.y=0, sort=F)
  
  if(any(duplicated(mRList$data_ann$Name))){
    mRList$data_ann <- keep_strongest_representative(X = mRList$data_ann)
  }
  
  wb <- createWorkbook()
  addWorksheet(wb, "associations")
  addWorksheet(wb, "summary")
  addWorksheet(wb, "metab_ann")
  addWorksheet(wb, "data_ann")
  writeDataTable(wb, "associations", mRList$metabolite_identification$associations)
  writeDataTable(wb, "summary", mRList$metabolite_identification$associations_summary)
  writeDataTable(wb, "metab_ann", mRList$metab_ann)
  writeDataTable(wb, "data_ann", mRList$data_ann)
  saveWorkbook(wb, file = file.path(dirout, "metabolite_identification.xlsx"), overwrite = T)
  
  cat("done\n")
  
  return(mRList)
  
}


barplot_ora <- function(pa_res, p = c("pvalue", "p.adjust"), filename = "pathway_analysis_ora") {
  
  p <- match.arg(p, c("pvalue", "p.adjust"))
  
  pa_res_ora_table <- pa_res[["res"]]@result
  
  the_bar_plot <- ggplot(data = pa_res_ora_table, aes(x = Description, y = Count, fill = !!sym(p))) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = "red", high = "blue", limits = c(0, 1)) +
    coord_flip() + 
    ggtitle("Over Representation Analysis (ORA)") +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
          plot.title = element_text(hjust = 0.5))
  
  ggsave(filename = paste0(filename, ".png"), plot = the_bar_plot, device = "png")
  
}


#################################################







###
## 3 Input data
###
mL <- read_input_file(feature_file = "Normalized_0_202112311138_No_BIF_BIFI.txt", sample_file = "meta.txt")



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


mL <- remove_samples(mRList = mL, ids = colnames(mL[["data"]])[grepl("BLK", colnames(mL[["data"]])) | grepl("solvent", colnames(mL[["data"]]))])


write_csv(bind_cols(tibble(feature = row.names(mL[["data"]])), as_tibble(mL[["data"]])), "database_blank_subtracted.csv")


####

###
## 5 Filtering, imputation and normalization
###

mL <- filtering(mL)




heatscatter_chromatography(mL)



mL_norm <- normalize_profiles(mL, method = "pqn")

write_csv(bind_cols(tibble(feature = row.names(mL_norm[["data"]])), as_tibble(mL_norm[["data"]])), "database_blank_subtracted_PQN_normalised.csv")




mL_norm <- CV_ratio(mRList = mL_norm, ratioCV=2)



mL <- RLA(mRList = mL, las=2)
mL <- RLA(mRList = mL_norm, las=2)



###
## 6 Principal Component Analysis
###


mL_norm <- mR_pca(mRList = mL_norm, nPcs=5, scaling="pareto", include_QC=TRUE, dirout = "PCA")

Plot2DPCA(mRList = mL_norm, pcx=1, pcy=2, col_by="class", col_pal = pals::trubetskoy(5), include_QC=TRUE, dirout = "PCA", name_samples = FALSE)



### PCA di nuovo dopo specificando include_QC come FALSE

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

mR_library <- select_library(source = "margheRita", column = "RPLong", mode = "POS", accept_RI=10)



mL_norm <- metabolite_identification(mL_norm, library_list = mR_library, filter_ann = TRUE)



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

png(file = "heat_map.png", width = 200, height = 400, res = 300, units = "mm")
h_map(mL_norm, scale_features=TRUE, column_ann="class", data.use = "data_ann" , col_ann= pals::trubetskoy(5), col=NULL, features = significant_features, feature_id = "Name", top = Inf, use_top_annoteted_metab = TRUE, test_method="anova", test_value = "q", show_column_names=FALSE)
dev.off()

write_tsv(tibble(Significant_features = significant_features), "significant_features.txt")

# metab_boxplot(mRList = mL_norm, col_by="class", col_pal = pals::trubetskoy(5),  group="class", features = "F1048")



#### Merge stat with features:


merge_stats_w_features(mRList = mL_norm, test_method = "anova", test_value = "q", cutoff_value = 0.05)








## 2024-08-27 adding differences from ANOVA

library(GetFeatistics)

data_per_ANOVA_GF <- transpose_feat_table(metab_stat_w_features_known_only_sign[,c(1,4:33)]) %>%
  add_column(group = as.factor((mL_norm[["sample_ann"]][["class"]])), .after = "samples")

colnames(data_per_ANOVA_GF) <- fix_names(colnames(data_per_ANOVA_GF))


ANOVA_GF <- gentab_P.1wayANOVA_posthocTukeyHSD(DF= data_per_ANOVA_GF,
                                               v= colnames(data_per_ANOVA_GF)[3:length(colnames(data_per_ANOVA_GF))],
                                               f = "group",
                                               FDR = TRUE,
                                               groupdiff = TRUE,
                                               pcutoff = 0.05)



metab_stat_w_features_known_only_sign_w_diff <- cbind(metab_stat_w_features_known_only_sign, ANOVA_GF)

export_the_table(metab_stat_w_features_known_only_sign_w_diff, "metab_stat_w_features_known_only_sign_w_diff", "txt")


## preparing the table each group vs the others

metab_stat_w_features_known_only_sign_w_diff_more <- as_tibble(metab_stat_w_features_known_only_sign_w_diff) %>%
  mutate(MIC_LAC_AMYL_vs_others = character(length(pull(metab_stat_w_features_known_only_sign_w_diff, 1))),
         MIC_LAC_KEFI_vs_others = character(length(pull(metab_stat_w_features_known_only_sign_w_diff, 1))),
         MIC_LAC_PLAN_vs_others = character(length(pull(metab_stat_w_features_known_only_sign_w_diff, 1))),
         MIC_LAC_SALI_vs_others = character(length(pull(metab_stat_w_features_known_only_sign_w_diff, 1))),
         MIC_BIF_BREV_vs_others = character(length(pull(metab_stat_w_features_known_only_sign_w_diff, 1))))

# 

for (i in 1:length(pull(metab_stat_w_features_known_only_sign_w_diff_more, 1))) {
  
  # MIC_LAC_AMYL
  
  df_MIC_LAC_AMYL <- metab_stat_w_features_known_only_sign_w_diff_more %>%
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
  
  metab_stat_w_features_known_only_sign_w_diff_more[i, "MIC_LAC_AMYL_vs_others"] <- paste0(df_MIC_LAC_AMYL_fil_vect_noNA, collapse = "; ")
 
  
  # MIC_LAC_KEFI
  
  df_MIC_LAC_KEFI <- metab_stat_w_features_known_only_sign_w_diff_more %>%
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
  
  metab_stat_w_features_known_only_sign_w_diff_more[i, "MIC_LAC_KEFI_vs_others"] <- paste0(df_MIC_LAC_KEFI_fil_vect_noNA, collapse = "; ")
  
  
  
  # MIC_LAC_PLAN
  
  df_MIC_LAC_PLAN <- metab_stat_w_features_known_only_sign_w_diff_more %>%
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
  
  metab_stat_w_features_known_only_sign_w_diff_more[i, "MIC_LAC_PLAN_vs_others"] <- paste0(df_MIC_LAC_PLAN_fil_vect_noNA, collapse = "; ")
  
  
  
  # MIC_LAC_SALI
  
  df_MIC_LAC_SALI <- metab_stat_w_features_known_only_sign_w_diff_more %>%
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
  
  metab_stat_w_features_known_only_sign_w_diff_more[i, "MIC_LAC_SALI_vs_others"] <- paste0(df_MIC_LAC_SALI_fil_vect_noNA, collapse = "; ")
  
  
  
  # MIC_BIF_BREV
  
  df_MIC_BIF_BREV <- metab_stat_w_features_known_only_sign_w_diff_more %>%
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
  
  metab_stat_w_features_known_only_sign_w_diff_more[i, "MIC_BIF_BREV_vs_others"] <- paste0(df_MIC_BIF_BREV_fil_vect_noNA, collapse = "; ")
  
   
}

write_tsv(metab_stat_w_features_known_only_sign_w_diff_more, "metab_stat_w_features_known_only_sign_w_diff_more.txt")
write_csv(metab_stat_w_features_known_only_sign_w_diff_more, "metab_stat_w_features_known_only_sign_w_diff_more.csv")

