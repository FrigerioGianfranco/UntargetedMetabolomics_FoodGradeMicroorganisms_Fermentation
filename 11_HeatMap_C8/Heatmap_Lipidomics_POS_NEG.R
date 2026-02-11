####Heatmap with family
library(pheatmap)
library(pals)
library(readxl)



#Here you can change only the colors for intensity scale
ht  <- function(data, sample_color, condition_color, output, separatore, condition_col=2, sample_col=2, numero_di_separatori=1, ...) {
  df <- as.data.frame(read_excel(data))
  rownames(df) <- df[,1]
  df[1] <- NULL
  meta_r <- data.frame(row.names = rownames(df), 
                       class = df$Family,
                       stringsAsFactors = T)
  df <- df[order(meta_r$class),]
  df[1] <- NULL 
  metadata <- data.frame(condition = stringr::str_split_fixed(names(df), separatore, numero_di_separatori)[,condition_col],
                         sample = stringr::str_split_fixed(names(df), separatore, numero_di_separatori)[,sample_col],
                         row.names = names(df),
                         stringsAsFactors = T)
  annotation_column <- metadata[,1:(dim(metadata)[2])]
  mycolors_s <- sample_color; names(mycolors_s) = levels(annotation_column$sample)
  mycolors_c <- condition_color; names(mycolors_c) = levels(annotation_column$condition)
  mycolors_r <- rainbow(length(unique(meta_r$class))); names(mycolors_r) = levels(meta_r$class)
  ann_colors = list(sample = mycolors_s, condition = mycolors_c, class=mycolors_r)
  colors <- colorRampPalette(c('#2A52BE', 'white','#F400A1'))(255) #change color scale here 
  ht <- pheatmap(df, annotation_col = annotation_column,
                 annotation_colors = ann_colors,
                 filename = output,
                 color = colors,
                 annotation_row = meta_r,
                 labels_row = row.names(df),
                 #cellheight = 7,
                 #cellwidth = 20,
                 ...)
}


ht(data="SIGN_metab_stat_w_features_w_diff_HeatMap_C8_NEG&POS-1.xlsx",
   #sample_color= c("violet", "blue", "red", "black"),#Use this code for colored samples
   sample_color= rep("white", 1, 20), #Use this code for white samples and type the number of replicates
   condition_color = trubetskoy(5),
   output="SIGN_metab_stat_w_features_w_diff_HeatMap_C8_NEG&POS-1.png",##Change name here
   separatore="_",
   condition_col=1,# gruppi principali Ctrl e HSWWP1
   sample_col=2, #nomi dei campioni 1-395 e 1-544 ect
   #height = 100,
   numero_di_separatori=3,#quante collone si formano dopo aver separato con _
   cluster_col=TRUE,
   cluster_row=FALSE, 
   scale="row", # per sfumature che indicano la intesita'
   colors=colors,
   border_color="white",
   show_rownames=TRUE,
   show_colnames=TRUE,
   cellwidth=6,
   cellheight=8, 
   fontsize_row=7, 
   fontsize = 7, 
   )


ht(data="SIGN_metab_stat_w_features_w_diff_HeatMap_C8_NEG&POS-2.xlsx",
   #sample_color= c("violet", "blue", "red", "black"),#Use this code for colored samples
   sample_color= rep("white", 1, 20), #Use this code for white samples and type the number of replicates
   condition_color = trubetskoy(5),
   output="SIGN_metab_stat_w_features_w_diff_HeatMap_C8_NEG&POS-2.png",##Change name here
   separatore="_",
   condition_col=1,# gruppi principali Ctrl e HSWWP1
   sample_col=2, #nomi dei campioni 1-395 e 1-544 ect
   #height = 100,
   numero_di_separatori=3,#quante collone si formano dopo aver separato con _
   cluster_col=TRUE,
   cluster_row=FALSE, 
   scale="row", # per sfumature che indicano la intesita'
   colors=colors,
   border_color="white",
   show_rownames=TRUE,
   show_colnames=TRUE,
   cellwidth=6,
   cellheight=8, 
   fontsize_row=7, 
   fontsize = 7, 
)


