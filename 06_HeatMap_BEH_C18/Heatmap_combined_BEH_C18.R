library(pheatmap)
library(pals)
library(readxl)

ht  <- function(data, sample_color, condition_color, output, separatore, condition_col=1, sample_col=2, numero_di_separatori=1, ...) {
  df <- as.data.frame(read_excel(data))
  rownames(df) <- df[,1]
  df[1] <- NULL
  metadata <- data.frame(condition = stringr::str_split_fixed(names(df), separatore, numero_di_separatori)[,condition_col],
                     sample = stringr::str_split_fixed(names(df), separatore, numero_di_separatori)[,sample_col],
                     row.names = names(df),
                     stringsAsFactors = T)
  annotation_column <- metadata[,1:(dim(metadata)[2])]
  mycolors_s <- sample_color; names(mycolors_s) = levels(annotation_column$sample)
  mycolors_c <- condition_color; names(mycolors_c) = levels(annotation_column$condition)
  ann_colors = list(sample = mycolors_s, condition = mycolors_c)
  colors <- colorRampPalette(c('#2A52BE', 'white','#F400A1'))(255) #change color scale here 
  ht <- pheatmap(df, annotation_col = annotation_column,
                         annotation_colors = ann_colors,
                         filename = output,
                         color = colors,
                         labels_row = row.names(df),
                          ...)
}

ht(data="SINTALresults_combined_fin.xlsx",
 sample_color= rep("white",1,20),
 condition_color = c(trubetskoy(5)),
 output="Heatmap_BEH&C18_merged.png",##Change name here
 separatore="_",
 condition_col=1,
 sample_col=2,
 numero_di_separatori=2,
 cluster_col=TRUE,
 cluster_row=TRUE, 
 scale="row",
 colors=colors,
 border_color="white",
 show_rownames=TRUE, 
 cellwidth = 6,
 cellheight = 8,
 show_colnames=T, 
 fontsize_row=7, 
 fontsize = 7
 )
