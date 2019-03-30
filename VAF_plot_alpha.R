# import pkgs
library(ggplot2)
library(maftools)
library(ggridges)
library(ggsci)

############ Using maftools ###############
# Usage: VAF_plot("maf_file_dir", "Sample_Option")

# VAF drawer
VAF_OFA <- function(cluster_all, theme_option){
  VAF_ofa_cha = paste("ggplot(cluster_all, aes(x=VAF, y=Tumor_Sample_Barcode)) +
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(), axis.line=element_line(size=0.25)) + 
    geom_point(aes(x=VAF, y=Tumor_Sample_Barcode, color = cluster), alpha = 0.5) +
    geom_density_ridges(color = \"cadetblue3\", fill = \"whitesmoke\", calc_ecdf = TRUE, alpha = 0.5) + ", 
    "scale_color_", theme_option, "() + scale_fill_", theme_option, "()", sep="")
  eval(parse(text = VAF_ofa_cha))
}

VAF_draw <- function(cluster_mt, theme_option){
  picv <- ggplot(cluster_mt, aes(x = VAF)) + geom_line(size = 1, colour = "cadetblue3", stat = "density")
  VAF_draw_cha = paste("ggplot(cluster_mt, aes(x = VAF)) + 
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(), axis.line=element_line(size=0.25)) + 
    geom_line(size = 1, colour = \"cadetblue3\", stat = \"density\") + geom_rug(aes(y = 0, colour = cluster), sides = \"b\") + ", 
    VAF_vline(cluster_mt, picv),
    "scale_color_", theme_option, "() + scale_fill_", theme_option, "()", sep="")
  eval(parse(text = VAF_draw_cha))
}

# VAF draw vlines
VAF_vline <- function(cluster_mt, pic){
  VAF_vline_cha <- ""
  cluster_ls <- unique(cluster_mt$cluster)
  density_info <- data.frame(layer_data(pic))
  for (cluster_name in cluster_ls){
    x_end <- max(cluster_mt[which(cluster_mt$cluster == cluster_name)]$VAF)
    x_end_alter <- density_info$x[which.min(abs(outer(density_info$x,x_end,FUN="-")))]
    y_end <- density_info$y[which(density_info$x == x_end_alter)]
    VAF_vline_cha <- paste(VAF_vline_cha, "geom_segment(aes(x = ", x_end_alter,", xend = ", x_end_alter, ", y = 0, yend = ", y_end,"), size = 0.3, colour=\"cadetblue3\", linetype=\"dashed\") + ",sep="")
  }
  VAF_vline_cha
}

# maftools Clusters(could not be used because of the specific position of the point)
VAF_plot <-function(maf_file, sample_option = "OFA", theme_option = "aaas", file_format = "png"){
  # read .maf
  maf_input <- read.table(maf_file, header = TRUE, fill = TRUE, sep = '\t', quote = "")
  laml <- read.maf(maf=maf_file)
  samples <- data.frame(maf_input[,ncol(maf_input)])
  cluster_all <- data.frame()
  n <- length(maf_input[,1])
  vaf_input_mt <- data.frame(maf_input[,1], maf_input[,ncol(maf_input)-3], samples)
  colnames(vaf_input_mt) <- c("Hugo_Symbol", "VAF", "Samples")
  tsb_ls <- as.data.frame(as.data.frame(table(samples))["samples"][which(as.data.frame(table(samples))["samples"]$samples != ""),])
  colnames(tsb_ls) <- c("samples")
  # sample_option
  if (sample_option == "All"){
  # multi samples
  for (counter_mt in 1:length(tsb_ls[,1])){
    for (sample_name_mt in tsb_ls){
      sample_mt <- vaf_input_mt[which(vaf_input_mt$Samples %in% as.character(sample_name_mt)[counter_mt]),]
      cluster_mt = inferHeterogeneity(maf = laml, tsb = as.character(sample_mt[1,3]), vafCol = 'VAF', useSyn = TRUE)$"clusterData"
      colnames(cluster_mt)[6] = "VAF"
      pic <- VAF_draw(cluster_mt, theme_option)
      ggsave(pic, filename =  paste(as.character(sample_name_mt)[counter_mt], "_VAF_Cluster", ".", file_format,sep=""), width = 12, height = 9)
    }
  }
  } else if (sample_option == "OFA"){
    # one pic for all sample
    patientID = strsplit(as.character(tsb_ls[1,1]), "-")[[1]][1]
    # collect cluster results
    for (counter_mt in 1:length(tsb_ls[,1])){
      for (sample_name_mt in tsb_ls){
        sample_mt <- vaf_input_mt[which(vaf_input_mt$Samples %in% as.character(sample_name_mt)[counter_mt]),]
        cluster_mt = inferHeterogeneity(maf = laml, tsb = as.character(sample_mt[1,3]), vafCol = 'VAF', useSyn = TRUE)$"clusterData"
        cluster_all <- rbind(cluster_all, cluster_mt)
        }
    }
    colnames(cluster_all)[6] = "VAF"
    pic <- VAF_OFA(cluster_all, theme_option)
    ggsave(pic, filename =  paste(patientID, "_VAF_Cluster",".", file_format, sep=""), width = 12, height = 9)
  } else {
  # specific sample
    sample_mt <- vaf_input_mt[which(vaf_input_mt$Samples %in% sample_option),]
    cluster_mt = inferHeterogeneity(maf = laml, tsb = as.character(sample_mt[1,3]), vafCol = 'VAF', useSyn = TRUE)$"clusterData"
    colnames(cluster_mt)[6] = "VAF"
    pic <- VAF_draw(cluster_mt, theme_option)
    ggsave(pic, filename =  paste(sample_option,"_VAF_Cluster",".", file_format,sep=""), width = 12, height = 9)
  }
}


########## Directory #######
# maf_dir = "/home/ninomoriaty/R_Project/patients_snv_indel.imputed.maf"
maf_file = "/home/ninomoriaty/R_Project/311252_snv_indel.imputed.maf"
# sample_option = "311252-S"


