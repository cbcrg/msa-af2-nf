library("reshape2")
library("ggplot2")
library("viridis")
library("superheat")
library("pheatmap")
library("ggpubr")
library("ComplexHeatmap")
library("RColorBrewer")
library("tidyverse")
library("SDMTools")
library("scales")
library("ggpmisc")
library("dplyr")
library("ggExtra")

args = commandArgs(trailingOnly=TRUE)
#if (length(args)==0) {
#  stop("At least one argument must be supplied (input dataset).n", call.=FALSE)
#}

intra_element_diff = function(myvector){
  collect_diff = c()
  for(i in 1:length(myvector)) {
    if (i != length(myvector)) {
      for(j in (i+1):length(myvector)){
       collect_diff = c(collect_diff,(myvector[i]-myvector[j]))
      }
    }
  }
  return(collect_diff)
}

rows_comb = function(mydf){
  cmb = combn(seq_len(nrow(mydf)), 2)
  comb_df = cbind(mydf[cmb[1,],], mydf[cmb[2,],])
  return(comb_df)
}

pairwise_plots = function(gdt_ts_df,sop_df,tcs_df,plddt_df,ref_aligner,included_struct){
  for (fam in levels(factor(gdt_ts_df$Family))) {
  tmpdf = gdt_ts_df[which(gdt_ts_df$Family == fam),]
  new_tmpdf = rows_comb(tmpdf)
  if (exists('pair_df') && is.data.frame(get('pair_df'))){
      pair_df = rbind(pair_df,new_tmpdf)
    } else {
      pair_df = new_tmpdf
    }
  }
  for (fam in levels(factor(tcs_df$Family))) {
  tmpdf_tcs = tcs_df[which(tcs_df$Family == fam),]
  new_tmpdf_tcs = rows_comb(tmpdf_tcs)
  if (exists('pair_df_tcs') && is.data.frame(get('pair_df_tcs'))){
      pair_df_tcs = rbind(pair_df_tcs,new_tmpdf_tcs)
    } else {
      pair_df_tcs = new_tmpdf_tcs
    }
  }
  for (fam in levels(factor(plddt_df$Family))) {
  tmpdf_plddt = plddt_df[which(plddt_df$Family == fam),]
  new_tmpdf_plddt = rows_comb(tmpdf_plddt)
  if (exists('pair_df_plddt') && is.data.frame(get('pair_df_plddt'))){
      pair_df_plddt = rbind(pair_df_plddt,new_tmpdf_plddt)
    } else {
      pair_df_plddt = new_tmpdf_plddt
    }
  }

  pair_df = pair_df[!duplicated(as.list(pair_df))]
  colnames(pair_df) = c("Sequence_1","GDT_TS_1","Family","Sequence_2","GDT_TS_2")
  pair_df$GDT_TS_1 = pair_df$GDT_TS_1*100
  pair_df$GDT_TS_2 = pair_df$GDT_TS_2*100

  pair_df_tcs = pair_df_tcs[!duplicated(as.list(pair_df_tcs))]
  colnames(pair_df_tcs) = c("Sequence_1","Family","TCS_1","Sequence_2","TCS_2")
  pair_df_tcs$TCS_1 = pair_df_tcs$TCS_1/10
  pair_df_tcs$TCS_2 = pair_df_tcs$TCS_2/10

  pair_df_plddt = pair_df_plddt[!duplicated(as.list(pair_df_plddt))]
  colnames(pair_df_plddt) = c("Sequence_1","pLDDT_1","Family","Sequence_2","pLDDT_2")
  merged_df_left = merge(merge(merge(sop_df,pair_df,by=c("Sequence_1","Sequence_2","Family")),pair_df_tcs,by=c("Sequence_1","Sequence_2","Family")),pair_df_plddt,by=c("Sequence_1","Sequence_2","Family"))
  colnames(pair_df_plddt) = c("Sequence_2","pLDDT_2","Family","Sequence_1","pLDDT_1")
  merged_df_right = merge(merge(merge(sop_df,pair_df,by=c("Sequence_1","Sequence_2","Family")),pair_df_tcs,by=c("Sequence_1","Sequence_2","Family")),pair_df_plddt,by=c("Sequence_1","Sequence_2","Family"))
  merged_df = rbind(merged_df_left,merged_df_right)

  #merged_df$MIN_GDT_TS = apply(merged_df[,c(5,6)],1, function(x) min(x))
  #merged_df$MEAN_GDT_TS = apply(merged_df[,c(5,6)],1, function(x) mean(x))
  merged_df$GM_GDT_TS = apply(merged_df[,c(5,6)],1, function(x) exp(mean(log(x))))
  #merged_df$MIN_TCS = apply(merged_df[,c(7,8)],1, function(x) min(x))
  #merged_df$MEAN_TCS = apply(merged_df[,c(7,8)],1, function(x) mean(x))
  merged_df$GM_TCS = apply(merged_df[,c(7,8)],1, function(x) exp(mean(log(x))))
  merged_df$GM_pLDDT = apply(merged_df[,c(9,10)],1, function(x) exp(mean(log(x))))


  #sop_title = paste0("SoP score on pairs of sequences ",ref_aligner,"_AF2 vs ",ref_aligner,"_NAT")
  sop_title = paste0("SoP score on pairs of sequences MSA-AF2 vs MSA-PDB")
  
  #gdt_title = paste0("Minimum GDT_TS scores per pair of sequences ")
  #gdt_title_avg = paste0("Average GDT_TS scores per pair of sequences ")
  gdt_title_gm = paste0("Geometric mean of GDT-TS scores on pairs of sequences")
  #tcs_title = paste0("Minimum TCS scores per pair of sequences ")
  #tcs_title_avg = paste0("Average TCS scores per pair of sequences ")
  tcs_title_gm = paste0("Geometric mean of TCS scores per pair of sequences ")
  plddt_title_gm = paste0("Geometric mean of pLDDT scores per pair of sequences ")
  #p = ggplot(merged_df,aes(x=SoP,y=MIN_GDT_TS,color=Family)) + geom_point(shape=1) + stat_cor(inherit.aes = FALSE,data=merged_df,aes(x=SoP,y=MIN_GDT_TS),method = "pearson") + geom_smooth(method='lm',inherit.aes = FALSE,data=merged_df,aes(x=SoP,y=MIN_GDT_TS)) + xlim(0,100) + ylim(0,100) + theme_light() + xlab(sop_title) + ylab(gdt_title)
  #ggsave(paste0("sop_vs_gdt_ts_pairwise_",included_struct,"_alphafold.png"),dpi="retina")
  #p = ggplot(merged_df,aes(x=SoP,y=MEAN_GDT_TS,color=Family)) + geom_point(shape=1) + stat_cor(inherit.aes = FALSE,data=merged_df,aes(x=SoP,y=MEAN_GDT_TS),method = "pearson") + geom_smooth(method='lm',inherit.aes = FALSE,data=merged_df,aes(x=SoP,y=MEAN_GDT_TS)) + xlim(0,100) + ylim(0,100) + theme_light() + xlab(sop_title) + ylab(gdt_title_avg)
  #ggsave(paste0("sop_vs_gdt_ts_pairwise_with_mean_",included_struct,"_alphafold.png"),dpi="retina")
  p = ggplot(merged_df,aes(x=SoP,y=GM_GDT_TS,color=GM_TCS)) + geom_point(alpha = 0.5) + stat_cor(inherit.aes = FALSE,data=merged_df,aes(x=SoP,y=GM_GDT_TS),method = "pearson") + geom_smooth(method='lm',inherit.aes = FALSE,data=merged_df,aes(x=SoP,y=GM_GDT_TS)) + xlim(0,100) + ylim(0,100) + scale_color_gradientn(colors = c("#1fa187","#fde725"), limits=c(75,100),na.value = "black") + theme_light() + xlab(sop_title) + ylab(gdt_title_gm)
  ggsave(paste0("sop_vs_gdt_ts_pairwise_with_gm_",included_struct,"_alphafold.png"),dpi=700)

  merged_df_sop_above_95 = merged_df[which(merged_df$SoP > 95),]
  merged_df_sop_above_95_seq_1 = merged_df_sop_above_95[,c("Sequence_1","GDT_TS_1","Family")]
  colnames(merged_df_sop_above_95_seq_1) = c("Sequence","GDT_TS","Family")
  merged_df_sop_above_95_seq_2 = merged_df_sop_above_95[,c("Sequence_2","GDT_TS_2","Family")]
  colnames(merged_df_sop_above_95_seq_2) = c("Sequence","GDT_TS","Family")
  merged_df_sop_above_95_seq_1_2 = rbind(merged_df_sop_above_95_seq_1,merged_df_sop_above_95_seq_2)
  merged_df_sop_above_95_seq_1_2 = merged_df_sop_above_95_seq_1_2[!duplicated(merged_df_sop_above_95_seq_1_2),]

  p = ggplot(merged_df_sop_above_95_seq_1_2, aes(x=GDT_TS,fill=Family)) + geom_histogram(color="black",position="stack") + theme_light() + xlab("AF2 GDT_TS with pairwise SoP > 95")
  ggsave(paste0("gdt_ts_with_pairwise_sop_above_95_",included_struct,"_alphafold.png"),dpi="retina")

  p = ggplot(merged_df_sop_above_95, aes(x=GM_GDT_TS,fill=Family)) + geom_histogram(color="black",position="stack") + theme_light() + xlab("AF2 geometric mean of GDT_TS pairs with pairwise SoP > 95")
  ggsave(paste0("gdt_ts_gm_with_pairwise_sop_above_95_",included_struct,"_alphafold.png"),dpi="retina")


  #
  # Deltas per family
  #

  for (fam in levels(factor(merged_df$Family))) {
    tmp_df = merged_df[which(merged_df$Family == fam),]
    delta_tcs = intra_element_diff(tmp_df$GM_TCS) #apply(combn(tmp_df$TCS_mTMalign_AF2,2), 2, diff)
    delta_sop = intra_element_diff(tmp_df$SoP) #apply(combn(tmp_df$SoP_mTMalign_AF2,2), 2, diff)
    delta_gdt_ts = intra_element_diff(tmp_df$GM_GDT_TS) #apply(combn(tmp_df$GDT_TS,2), 2, diff)
    delta_plddt = intra_element_diff(tmp_df$GM_pLDDT) #apply(combn(tmp_df$pLDDT,2), 2, diff)
    #delta_nirmsd = intra_element_diff(tmp_df$niRMSD)*(-1)
    
    if (exists('delta_df') && is.data.frame(get('delta_df'))){
      delta_df = rbind(delta_df,data.frame(delta_tcs,delta_sop,delta_gdt_ts,delta_plddt,Family=rep(fam,length(delta_gdt_ts))))
    } else {
      delta_df = data.frame(delta_tcs,delta_sop,delta_gdt_ts,delta_plddt,Family=rep(fam,length(delta_gdt_ts)))
    }
  }

  delta_sop_title = paste0("Δ Sum-of-Pairs scores between ", included_struct," pairs within family - ",ref_aligner,"_AF2 vs ",ref_aligner,"_NAT")
  delta_tcs_title = paste0("Δ TCS scores between ",included_struct, " pairs within family - ",ref_aligner,"_AF2")
  delta_gdt_title = paste0("Δ GDT_TS scores between ",included_struct, " pairs within family")
  delta_plddt_title = paste0("Δ pLDDT scores between ", included_struct," pairs within family")
  #nirmsd_title = paste0("Δ niRMSD scores between ", included_struct, " pairs within family - ",ref_aligner,"_AF2 with AF2 structures")

  p = ggplot(delta_df,aes(x=delta_sop,y=delta_tcs,color=Family)) + geom_point(shape=1,alpha=0.6) + stat_quadrant_counts() + xlim(-100,100) + ylim(-100,100) + theme_light() + xlab(delta_sop_title) + ylab(delta_tcs_title) + theme(axis.title = element_text(size = 9)) #+ stat_cor(inherit.aes = FALSE,data=delta_df,aes(x=delta_sop,y=delta_tcs),method = "pearson")
  ggsave(paste0("delta_sop_vs_delta_tcs_",included_struct,"_ref_",ref_aligner,"_alphafold.png"),dpi="retina")
  test = table(delta_df$delta_sop > 0, delta_df$delta_tcs > 0)

  p = ggplot(delta_df,aes(x=delta_sop,y=delta_plddt,color=Family)) + geom_point(shape=1,alpha=0.6) + stat_quadrant_counts() + xlim(-100,100) + ylim(-100,100) + theme_light() + xlab(delta_sop_title) + ylab(delta_plddt_title) + theme(axis.title = element_text(size = 9)) #+ stat_cor(inherit.aes = FALSE,data=delta_df,aes(x=delta_sop,y=delta_plddt),method = "pearson")
  ggsave(paste0("delta_sop_vs_delta_plddt_",included_struct,"_ref_",ref_aligner,"_alphafold.png"),dpi="retina")
  p = ggplot(delta_df,aes(x=delta_sop,y=delta_gdt_ts,color=Family)) + geom_point(shape=1,alpha=0.6) + stat_quadrant_counts() + xlim(-100,100) + ylim(-100,100) + theme_light() + xlab(delta_sop_title) + ylab(delta_gdt_title) + theme(axis.title = element_text(size = 9)) #+ stat_cor(inherit.aes = FALSE,data=delta_df,aes(x=delta_sop,y=delta_gdt_ts),method = "pearson")
  ggsave(paste0("delta_sop_vs_delta_gdt_ts_",included_struct,"_ref_",ref_aligner,"_alphafold.png"),dpi="retina")

  return(merged_df)

}


pairwise_plots_dmpfold = function(gdt_ts_df,sop_df,tcs_df,ref_aligner,included_struct){
  for (fam in levels(factor(gdt_ts_df$Family))) {
  tmpdf = gdt_ts_df[which(gdt_ts_df$Family == fam),]
  new_tmpdf = rows_comb(tmpdf)
  if (exists('pair_df') && is.data.frame(get('pair_df'))){
      pair_df = rbind(pair_df,new_tmpdf)
    } else {
      pair_df = new_tmpdf
    }
  }
  for (fam in levels(factor(tcs_df$Family))) {
  tmpdf_tcs = tcs_df[which(tcs_df$Family == fam),]
  new_tmpdf_tcs = rows_comb(tmpdf_tcs)
  if (exists('pair_df_tcs') && is.data.frame(get('pair_df_tcs'))){
      pair_df_tcs = rbind(pair_df_tcs,new_tmpdf_tcs)
    } else {
      pair_df_tcs = new_tmpdf_tcs
    }
  }

  pair_df = pair_df[!duplicated(as.list(pair_df))]
  colnames(pair_df) = c("Sequence_1","GDT_TS_1","Family","Sequence_2","GDT_TS_2")
  pair_df$GDT_TS_1 = pair_df$GDT_TS_1*100
  pair_df$GDT_TS_2 = pair_df$GDT_TS_2*100

  pair_df_tcs = pair_df_tcs[!duplicated(as.list(pair_df_tcs))]
  colnames(pair_df_tcs) = c("Sequence_1","Family","TCS_1","Sequence_2","TCS_2")
  pair_df_tcs$TCS_1 = pair_df_tcs$TCS_1/10
  pair_df_tcs$TCS_2 = pair_df_tcs$TCS_2/10

  merged_df = merge(merge(sop_df,pair_df,by=c("Sequence_1","Sequence_2","Family")),pair_df_tcs,by=c("Sequence_1","Sequence_2","Family"))
  #merged_df$MIN_GDT_TS = apply(merged_df[,c(5,6)],1, function(x) min(x))
  #merged_df$MEAN_GDT_TS = apply(merged_df[,c(5,6)],1, function(x) mean(x))
  merged_df$GM_GDT_TS = apply(merged_df[,c(5,6)],1, function(x) exp(mean(log(x))))
  #merged_df$MIN_TCS = apply(merged_df[,c(7,8)],1, function(x) min(x))
  #merged_df$MEAN_TCS = apply(merged_df[,c(7,8)],1, function(x) mean(x))
  merged_df$GM_TCS = apply(merged_df[,c(7,8)],1, function(x) exp(mean(log(x))))

  sop_title = paste0("Sum-of-Pairs scores per pair of sequences ",ref_aligner,"_DMPFOLD vs ",ref_aligner,"_NAT")
  #gdt_title = paste0("Minimum GDT_TS scores per pair of sequences ")
  #gdt_title_avg = paste0("Average GDT_TS scores per pair of sequences ")
  gdt_title_gm = paste0("Geometric mean of GDT_TS scores per pair of sequences ")
  #tcs_title = paste0("Minimum TCS scores per pair of sequences ")
  #tcs_title_avg = paste0("Average TCS scores per pair of sequences ")
  tcs_title_gm = paste0("Geometric mean of TCS scores per pair of sequences ")
  #p = ggplot(merged_df,aes(x=SoP,y=MIN_GDT_TS,color=Family)) + geom_point(shape=1) + stat_cor(inherit.aes = FALSE,data=merged_df,aes(x=SoP,y=MIN_GDT_TS),method = "pearson") + geom_smooth(method='lm',inherit.aes = FALSE,data=merged_df,aes(x=SoP,y=MIN_GDT_TS)) + xlim(0,100) + ylim(0,100) + theme_light() + xlab(sop_title) + ylab(gdt_title)
  #ggsave(paste0("sop_vs_gdt_ts_pairwise_",included_struct,"_dmpfold.png"),dpi="retina")
  #p = ggplot(merged_df,aes(x=SoP,y=MEAN_GDT_TS,color=Family)) + geom_point(shape=1) + stat_cor(inherit.aes = FALSE,data=merged_df,aes(x=SoP,y=MEAN_GDT_TS),method = "pearson") + geom_smooth(method='lm',inherit.aes = FALSE,data=merged_df,aes(x=SoP,y=MEAN_GDT_TS)) + xlim(0,100) + ylim(0,100) + theme_light() + xlab(sop_title) + ylab(gdt_title_avg)
  #ggsave(paste0("sop_vs_gdt_ts_pairwise_with_mean_",included_struct,"_dmpfold.png"),dpi="retina")
  p = ggplot(merged_df,aes(x=SoP,y=GM_GDT_TS,color=GM_TCS)) + geom_point(alpha=0.6) + stat_cor(inherit.aes = FALSE,data=merged_df,aes(x=SoP,y=GM_GDT_TS),method = "pearson") + geom_smooth(method='lm',inherit.aes = FALSE,data=merged_df,aes(x=SoP,y=GM_GDT_TS)) + xlim(0,100) + ylim(0,100) + scale_color_gradientn(colors = c("#1fa187","#fde725"), limits=c(75,100),na.value = "black") + theme_light() + xlab(sop_title) + ylab(gdt_title_gm)
  ggsave(paste0("sop_vs_gdt_ts_pairwise_with_gm_",included_struct,"_dmpfold.png"),dpi="retina")


  for (fam in levels(factor(merged_df$Family))) {
    tmp_df = merged_df[which(merged_df$Family == fam),]
    delta_tcs = intra_element_diff(tmp_df$GM_TCS)
    delta_sop = intra_element_diff(tmp_df$SoP)
    delta_gdt_ts = intra_element_diff(tmp_df$GM_GDT_TS)
    #delta_nirmsd = intra_element_diff(tmp_df$niRMSD)*(-1)
    if (exists('delta_df') && is.data.frame(get('delta_df'))){
      delta_df = rbind(delta_df,data.frame(delta_tcs,delta_sop,delta_gdt_ts,Family=rep(fam,length(delta_gdt_ts))))
    } else {
      delta_df = data.frame(delta_tcs,delta_sop,delta_gdt_ts,Family=rep(fam,length(delta_gdt_ts)))
    }
  }

  delta_sop_title = paste0("Δ Sum-of-Pairs scores between ", included_struct," pairs within family - ",ref_aligner,"_DMPFOLD vs ",ref_aligner,"_NAT")
  delta_tcs_title = paste0("Δ TCS scores between ",included_struct, " pairs within family - ",ref_aligner,"_DMPFOLD")
  delta_gdt_title = paste0("Δ GDT_TS scores between ",included_struct, " pairs within family")
  #delta_nirmsd_title = paste0("Δ niRMSD scores between ", included_struct, " pairs within family - ",ref_aligner,"_DMPFOLD with DMPFOLD structures")

  p = ggplot(delta_df,aes(x=delta_sop,y=delta_tcs,color=Family)) + geom_point(shape=1,alpha=0.6) + stat_quadrant_counts() + xlim(-100,100) + ylim(-100,100) + theme_light() + xlab(delta_sop_title) + ylab(delta_tcs_title) + theme(axis.title = element_text(size = 9))
  ggsave(paste0("delta_sop_vs_delta_tcs_",included_struct,"_ref_",ref_aligner,"_dmpfold.png"),dpi="retina")
  p = ggplot(delta_df,aes(x=delta_sop,y=delta_gdt_ts,color=Family)) + geom_point(shape=1,alpha=0.6) + stat_quadrant_counts() + xlim(-100,100) + ylim(-100,100) + theme_light() + xlab(delta_sop_title) + ylab(delta_gdt_title) + theme(axis.title = element_text(size = 9))
  ggsave(paste0("delta_sop_vs_delta_gdt_ts_",included_struct,"_ref_",ref_aligner,"_dmpfold.png"),dpi="retina")

  return(merged_df)


}


circular_plot = function(tcs_cutoff_list_and_sp,first_label,second_label,output_name,sel_colors) {
  empty_bar <- 6
  to_add <- data.frame( matrix(NA, empty_bar*nlevels(as.factor(tcs_cutoff_list_and_sp$Family)), ncol(tcs_cutoff_list_and_sp) ))
  colnames(to_add) <- colnames(tcs_cutoff_list_and_sp)
  to_add$Family <- rep(levels(as.factor(tcs_cutoff_list_and_sp$Family)), each=empty_bar)
  tcs_cutoff_list_and_sp <- rbind(tcs_cutoff_list_and_sp, to_add)
  tcs_cutoff_list_and_sp <- tcs_cutoff_list_and_sp %>% arrange(Family)
  #tcs_cutoff_list_and_sp$id <- seq(1, nrow(tcs_cutoff_list_and_sp))
  tcs_cutoff_list_and_sp$id <- sort(rep(seq(1, nrow(tcs_cutoff_list_and_sp)/2),times=2))

  tcs_cutoff_list_and_sp$SP = tcs_cutoff_list_and_sp$SP
  tcs_cutoff_list_and_sp$percent_of_seq = -1/2 * tcs_cutoff_list_and_sp$percent_of_seq

  label_data <- tcs_cutoff_list_and_sp %>%
  group_by(id) %>%
  summarize(sp=max(SP),tcs=TCS_mode)
  label_data = label_data[duplicated(label_data),]
  number_of_bar <- nrow(label_data)
  angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  label_data$hjust <- ifelse( angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle+180, angle)

  base_data <- tcs_cutoff_list_and_sp %>% 
  group_by(Family) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

  grid_data <- base_data
  grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
  grid_data$start <- grid_data$start - 1
  grid_data <- grid_data[-1,]


  nb.cols = sel_colors
  mycolors = c(sel_colors,"red","blue")

  p = ggplot(tcs_cutoff_list_and_sp, aes(x=as.factor(id), y=SP)) + 
  geom_bar(aes(x=as.factor(id), y=SP, fill=group),position="identity", stat="identity", alpha=0.6, width=0.8) +
  geom_segment(data=grid_data, aes(x = start-2, y = 80, xend = start, yend = 80), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_bar(aes(x=as.factor(id), y=SP, fill=group),position="identity", stat="identity", alpha=0.6, width=0.8) +
  geom_segment(data=grid_data, aes(x = start-2, y = 80, xend = start, yend = 80), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = start-2, y = 60, xend = start, yend = 60), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = start-2, y = 40, xend = start, yend = 40), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = start-2, y = 20, xend = start, yend = 20), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  annotate("text", x = rep(max(tcs_cutoff_list_and_sp$id),4), y = c(20, 40, 60, 80), label = c("20", "40", "60", "80") , color="grey", size=2 , angle=0, fontface="bold", hjust=1) +
  #geom_bar(aes(x=as.factor(id), y=SP, fill=group),position="identity", stat="identity", alpha=0.9, width=0.5) +
  geom_bar(aes(x=as.factor(id), y=percent_of_seq, fill=Family),position="identity", stat="identity", alpha=0.3, width=0.8) + scale_fill_manual(name = "group", values=mycolors,breaks = c(first_label,second_label)) +
  theme_minimal() +
  theme(
    legend.position = c(0.5,0.5),
    legend.title = element_blank(),
    legend.text=element_text(size=5),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm"),
  ) +
  ylim(-120,140) +
  coord_polar() +  
  geom_text(data=label_data, aes(x=id, y=sp+10, label=tcs, hjust=hjust), color="black", fontface="bold", alpha=0.8, size=1.3, angle= label_data$angle, inherit.aes = FALSE ) +
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end+2, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title+1, y = -18, label=Family), colour = "black", alpha=0.8, size=1.5, fontface="bold", inherit.aes = FALSE)
  ggsave(output_name,dpi="retina")
}

circular_plot_dodge = function(tcs_cutoff_list_and_sp,first_label,second_label,output_name) {
  empty_bar <- 3
  to_add <- data.frame( matrix(NA, empty_bar*nlevels(as.factor(tcs_cutoff_list_and_sp$Family)), ncol(tcs_cutoff_list_and_sp) ))
  colnames(to_add) <- colnames(tcs_cutoff_list_and_sp)
  to_add$Family <- rep(levels(as.factor(tcs_cutoff_list_and_sp$Family)), each=empty_bar)
  tcs_cutoff_list_and_sp <- rbind(tcs_cutoff_list_and_sp, to_add)
  tcs_cutoff_list_and_sp <- tcs_cutoff_list_and_sp %>% arrange(Family)
  tcs_cutoff_list_and_sp$id <- seq(1, nrow(tcs_cutoff_list_and_sp))

  label_data <- tcs_cutoff_list_and_sp
  number_of_bar <- nrow(label_data)
  angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  label_data$hjust <- ifelse( angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle+180, angle)

  base_data <- tcs_cutoff_list_and_sp %>% 
  group_by(Family) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

  grid_data <- base_data
  grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
  grid_data$start <- grid_data$start - 1
  grid_data <- grid_data[-1,]

  tcs_cutoff_list_and_sp$percent_of_seq = -0.5 * tcs_cutoff_list_and_sp$percent_of_seq

  nb.cols = nlevels(as.factor(tcs_cutoff_list_and_sp$Family)) + 2
  mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)
#position_dodge(3)
  p = ggplot(tcs_cutoff_list_and_sp, aes(x=as.factor(id), y=SP)) + geom_bar(aes(x=as.factor(id), y=SP, fill=group),position=position_dodge(3), stat="identity", alpha=0.5) + geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_bar(aes(x=as.factor(id), y=SP, fill=group),position=position_dodge(3), stat="identity", alpha=0.5) +
  geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  annotate("text", x = rep(max(tcs_cutoff_list_and_sp$id),4), y = c(20, 40, 60, 80), label = c("20", "40", "60", "80") , color="grey", size=2 , angle=0, fontface="bold", hjust=1) +
  geom_bar(aes(x=as.factor(id), y=SP, fill=group),position=position_dodge(3), stat="identity", alpha=0.5) +
  geom_bar(aes(x=as.factor(id), y=percent_of_seq, fill=Family),position=position_dodge(3), stat="identity", alpha=0.5) + scale_fill_manual(name = "group", values=mycolors,breaks = c(first_label,second_label)) +
  theme_minimal() +
  theme(
    legend.position = c(0.5,0.5),
    legend.title = element_blank(),
    legend.text=element_text(size=5),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  ylim(-120,120) +
  coord_polar() +  
  geom_text(data=label_data, aes(x=id, y=SP+10, label=TCS_mode, hjust=hjust), color="black", fontface="bold", alpha=0.8, size=1.3, angle= label_data$angle, inherit.aes = FALSE ) +
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -18, label=Family), colour = "black", alpha=0.8, size=2, fontface="bold", inherit.aes = FALSE)
  ggsave(output_name,dpi="retina")
}



delta_analysis_AF2 = function (tcs_df,gdt_ts_df,plddts_df,sop_df,nirmsd_df,ref_aligner,included_struct) {
  merged_df = merge(merge(merge(merge(tcs_df,gdt_ts_df,by=c("Sequence","Family")),plddts_df,by=c("Sequence","Family")),sop_df,by=c("Sequence","Family")),nirmsd_df,by=c("Sequence","Family"))
  merged_df = unique(merged_df)
  colnames(merged_df)[4] = "GDT_TS"
  merged_df$GDT_TS = merged_df$GDT_TS*100
  merged_df$TCS = merged_df$TCS/10
  merged_df = merged_df[order(merged_df$Family),]
  
  color_codes <- colorRampPalette(brewer.pal(8, "Set1"))(12)
  if(FALSE){
    gather_val = c()
    for (weight in seq(1, 100, by=1)) {
    val = mean(((weight*merged_df$TCS) + ((100-weight)*merged_df$pLDDT))/100)
    gather_val = c(gather_val,val)
    }
    mix_tcs_and_plddt = data.frame(seq(1, 100, by=1),gather_val)
    colnames(mix_tcs_and_plddt)[1] = "weight"
    p = ggplot(mix_tcs_and_plddt, aes(x=weight,y=gather_val)) + geom_point() + ylim(c(0,100)) + ylab("mean( (weight*TCS + (100-weight)*pLDDT)/100 )")
    ggsave("mix_tcs_and_plddt.png",dpi="retina")
  }
  

  #
  # Deltas per family
  #

  for (fam in levels(factor(merged_df$Family))) {
    tmp_df = merged_df[which(merged_df$Family == fam),]
    delta_tcs = intra_element_diff(tmp_df$TCS) #apply(combn(tmp_df$TCS_mTMalign_AF2,2), 2, diff)
    delta_sop = intra_element_diff(tmp_df$SoP) #apply(combn(tmp_df$SoP_mTMalign_AF2,2), 2, diff)
    delta_gdt_ts = intra_element_diff(tmp_df$GDT_TS) #apply(combn(tmp_df$GDT_TS,2), 2, diff)
    delta_plddt = intra_element_diff(tmp_df$pLDDT) #apply(combn(tmp_df$pLDDT,2), 2, diff)
    delta_nirmsd = intra_element_diff(tmp_df$niRMSD)*(-1)
    #mult_tcs_and_plddt = delta_tcs * delta_plddt
    #subtr_tcs_and_plddt = delta_tcs - delta_plddt
    if (exists('delta_df') && is.data.frame(get('delta_df'))){
      delta_df = rbind(delta_df,data.frame(delta_tcs,delta_sop,delta_gdt_ts,delta_plddt,delta_nirmsd,Family=rep(fam,length(delta_gdt_ts))))
    } else {
      delta_df = data.frame(delta_tcs,delta_sop,delta_gdt_ts,delta_plddt,delta_nirmsd,Family=rep(fam,length(delta_gdt_ts)))
    }
  }

  sop_title = paste0("Δ Sum-of-Pairs scores between ", included_struct," pairs within family - ",ref_aligner,"_AF2 vs ",ref_aligner,"_NAT")
  tcs_title = paste0("Δ TCS scores between ",included_struct, " pairs within family - ",ref_aligner,"_AF2")
  gdt_title = paste0("Δ GDT_TS scores between ",included_struct, " pairs within family")
  plddt_title = paste0("Δ pLDDT scores between ", included_struct," pairs within family")
  nirmsd_title = paste0("Δ niRMSD scores between ", included_struct, " pairs within family - ",ref_aligner,"_AF2 with AF2 structures")

  sop_title_raw = paste0("Sum-of-Pairs scores per sequence ",ref_aligner,"_AF2 vs ",ref_aligner,"_NAT")
  tcs_title_raw = paste0("TCS scores per sequence ",ref_aligner,"_AF2")
  gdt_title_raw = paste0("GDT_TS scores per sequence ")
  plddt_title_raw = paste0("pLDDT scores per sequence ")
  nirmsd_title_raw = paste0("niRMSD scores per sequence ",ref_aligner,"_AF2 with AF2 structures")


  #p = ggplot(delta_df,aes(x=delta_sop,y=delta_tcs,color=Family)) + geom_point(shape=1) + stat_quadrant_counts() + xlim(-100,100) + ylim(-100,100) + theme_light() + xlab(sop_title) + ylab(tcs_title) + theme(axis.title = element_text(size = 9)) #+ stat_cor(inherit.aes = FALSE,data=delta_df,aes(x=delta_sop,y=delta_tcs),method = "pearson")
  #p = ggplot(delta_df,aes(x=delta_sop,y=delta_tcs)) + geom_density_2d_filled(contour_var = "ndensity") + geom_hline(yintercept=0,size=0.1) + geom_vline(xintercept=0,size=0.1) + stat_quadrant_counts() + scale_x_continuous(limits = c(-100,100),expand = c(0, 0)) + scale_y_continuous(limits = c(-100,100),expand = c(0, 0)) + theme_light() + xlab(sop_title) + ylab(tcs_title) + theme(axis.title = element_text(size = 9)) # + stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) + scale_fill_distiller(palette="Oranges", direction=1) + stat_cor(inherit.aes = FALSE,data=delta_df,aes(x=delta_sop,y=delta_tcs),method = "pearson")
  #ggsave(paste0("delta_sop_vs_delta_tcs_",included_struct,"_ref_",ref_aligner,"_alphafold.png"),dpi="retina")
  
  p = ggplot(delta_df,aes(x=delta_sop,y=delta_plddt,color=Family)) + geom_point(shape=1) + stat_quadrant_counts() + xlim(-100,100) + ylim(-100,100) + theme_light() + xlab(sop_title) + ylab(plddt_title) + theme(axis.title = element_text(size = 9)) #+ stat_cor(inherit.aes = FALSE,data=delta_df,aes(x=delta_sop,y=delta_plddt),method = "pearson")
  ggsave(paste0("delta_sop_vs_delta_plddt_",included_struct,"_ref_",ref_aligner,"_alphafold.png"),dpi="retina")
  p = ggplot(delta_df,aes(x=delta_sop,y=delta_gdt_ts,color=Family)) + geom_point(shape=1) + stat_quadrant_counts() + xlim(-100,100) + ylim(-100,100) + theme_light() + xlab(sop_title) + ylab(gdt_title) + theme(axis.title = element_text(size = 9)) #+ stat_cor(inherit.aes = FALSE,data=delta_df,aes(x=delta_sop,y=delta_gdt_ts),method = "pearson")
  ggsave(paste0("delta_sop_vs_delta_gdt_ts_",included_struct,"_ref_",ref_aligner,"_alphafold.png"),dpi="retina")

  p = ggplot(delta_df,aes(x=delta_sop,y=delta_nirmsd,color=Family)) + geom_point(shape=1) + stat_quadrant_counts() + xlim(-100,100) + ylim(-3,3) + theme_light() + xlab(sop_title) + ylab(nirmsd_title) + theme(axis.title = element_text(size = 9)) #+ stat_cor(inherit.aes = FALSE,data=delta_df,aes(x=delta_sop,y=delta_tcs),method = "pearson")
  ggsave(paste0("delta_sop_vs_delta_nirmsd_",included_struct,"_ref_",ref_aligner,"_alphafold.png"),dpi="retina")

  p = ggplot(delta_df,aes(x=delta_nirmsd,y=delta_tcs,color=Family)) + geom_point(shape=1) + stat_quadrant_counts() + xlim(-3,3) + ylim(-100,100) + theme_light() + xlab(nirmsd_title) + ylab(tcs_title) + theme(axis.title = element_text(size = 9)) # #+ xlim(-3,3) + stat_cor(inherit.aes = FALSE,data=delta_df,aes(x=delta_sop,y=delta_tcs),method = "pearson")
  ggsave(paste0("delta_nirmsd_vs_delta_tcs_",included_struct,"_ref_",ref_aligner,"_alphafold.png"),dpi="retina")
  p = ggplot(delta_df,aes(x=delta_nirmsd,y=delta_plddt,color=Family)) + geom_point(shape=1) + stat_quadrant_counts() + xlim(-3,3) + ylim(-100,100) + theme_light() + xlab(nirmsd_title) + ylab(plddt_title) + theme(axis.title = element_text(size = 9)) # + xlim(-80,80) + stat_cor(inherit.aes = FALSE,data=delta_df,aes(x=delta_sop,y=delta_plddt),method = "pearson")
  ggsave(paste0("delta_nirmsd_vs_delta_plddt_",included_struct,"_ref_",ref_aligner,"_alphafold.png"),dpi="retina")
  p = ggplot(delta_df,aes(x=delta_nirmsd,y=delta_gdt_ts,color=Family)) + geom_point(shape=1) + stat_quadrant_counts() + xlim(-3,3) + ylim(-100,100) + theme_light() + xlab(nirmsd_title) + ylab(gdt_title) + theme(axis.title = element_text(size = 9)) #+ xlim(-80,80) + stat_cor(inherit.aes = FALSE,data=delta_df,aes(x=delta_sop,y=delta_gdt_ts),method = "pearson")
  ggsave(paste0("delta_nirmsd_vs_delta_gdt_ts_",included_struct,"_ref_",ref_aligner,"_alphafold.png"),dpi="retina")

  
  p = ggplot(delta_df,aes(x=delta_gdt_ts,y=delta_tcs,color=Family)) + geom_point(shape=1) + stat_quadrant_counts() + xlim(-100,100) + ylim(-100,100) + theme_light() + xlab(gdt_title) + ylab(tcs_title) + theme(axis.title = element_text(size = 9)) #+ stat_cor(inherit.aes = FALSE,data=delta_df,aes(x=delta_gdt_ts,y=delta_tcs),method = "pearson")
  ggsave(paste0("delta_gdt_ts_vs_delta_tcs_",included_struct,"_ref_",ref_aligner,"_alphafold.png"),dpi="retina")
  p = ggplot(delta_df,aes(x=delta_gdt_ts,y=delta_plddt,color=Family)) + geom_point(shape=1) + stat_quadrant_counts() + xlim(-100,100) + ylim(-100,100) + theme_light() + xlab(gdt_title) + ylab(plddt_title) + theme(axis.title = element_text(size = 9)) #+ stat_cor(inherit.aes = FALSE,data=delta_df,aes(x=delta_gdt_ts,y=delta_plddt),method = "pearson")
  ggsave(paste0("delta_gdt_ts_vs_delta_plddt_",included_struct,"_alphafold.png"),dpi="retina")

  p = ggplot(merged_df,aes(x=SoP,y=niRMSD,color=Family)) + geom_point(shape=1) + stat_cor(inherit.aes = FALSE,data=merged_df,aes(x=SoP,y=niRMSD),method = "pearson") + geom_smooth(method='lm',inherit.aes = FALSE,data=merged_df,aes(x=SoP,y=niRMSD)) + xlim(0,100) + theme_light() + xlab(sop_title_raw) + ylab(nirmsd_title_raw)
  ggsave(paste0("sop_vs_nirmsd_",included_struct,"_alphafold.png"),dpi="retina")

  p = ggplot(merged_df,aes(x=SoP,y=TCS,color=Family)) + geom_point(shape=1) + stat_cor(inherit.aes = FALSE,data=merged_df,aes(x=SoP,y=TCS),method = "pearson") + geom_smooth(method='lm',inherit.aes = FALSE,data=merged_df,aes(x=SoP,y=TCS)) + xlim(0,100) + ylim(0,100) + theme_light() + xlab(sop_title_raw) + ylab(tcs_title_raw)
  ggsave(paste0("sop_vs_tcs_",included_struct,"_alphafold.png"),dpi="retina")

  p = ggplot(merged_df,aes(x=SoP,y=pLDDT,color=Family)) + geom_point(shape=1) + stat_cor(inherit.aes = FALSE,data=merged_df,aes(x=SoP,y=pLDDT),method = "pearson") + geom_smooth(method='lm',inherit.aes = FALSE,data=merged_df,aes(x=SoP,y=pLDDT)) + xlim(0,100) + ylim(0,100) + theme_light() + xlab(sop_title_raw) + ylab(plddt_title_raw)
  ggsave(paste0("sop_vs_plddt_",included_struct,"_alphafold.png"),dpi="retina")

  p = ggplot(merged_df,aes(x=SoP,y=GDT_TS,color=Family)) + geom_point(shape=1) + stat_cor(inherit.aes = FALSE,data=merged_df,aes(x=SoP,y=GDT_TS),method = "pearson") + geom_smooth(method='lm',inherit.aes = FALSE,data=merged_df,aes(x=SoP,y=GDT_TS)) + xlim(0,100) + ylim(0,100) + theme_light() + xlab(sop_title_raw) + ylab(gdt_title_raw)
  ggsave(paste0("sop_vs_gdt_ts_",included_struct,"_alphafold.png"),dpi="retina")


  return(delta_df)

  #p = ggscatter(tmp_df,x="SoP",y="niRMSD",add = "reg.line",conf.int = TRUE,color="Family",palette=color_codes, add.params = list(color = "blue",fill = "lightgray")) + stat_cor(method = "pearson") + xlab(sop_title) + ylab(nirmsd_title) + font("xlab", size = 12) + font("ylab", size = 12)
  #ggsave(paste0("sop_vs_nirmsd_",included_struct,"_alphafold.png"),dpi="retina")

  #p = ggplot(delta_df,aes(x=delta_gdt_ts,y=mult_tcs_and_plddt,color=Family)) + geom_point(shape=1) + stat_cor(inherit.aes = FALSE,data=delta_df,aes(x=delta_gdt_ts,y=mult_tcs_and_plddt),method = "pearson") + xlim(-80,80) + theme_light() + xlab("Δ GDT_TS scores between all pairs within family") + ylab(paste0("Δ TCS scores multiplied by delta pLDDT scores between all pairs within family - ",ref_aligner,"_AF2")) # + stat_quadrant_counts() + ylim(-80,80)
  #ggsave(paste0("delta_gdt_ts_vs_delta_tcs_multiplied_by_plddt_ref_",ref_aligner,"_alphafold.png"),dpi="retina")

  #p = ggplot(delta_df,aes(x=delta_gdt_ts,y=subtr_tcs_and_plddt,color=Family)) + geom_point(shape=1) + stat_cor(inherit.aes = FALSE,data=delta_df,aes(x=delta_gdt_ts,y=subtr_tcs_and_plddt),method = "pearson") + xlim(-80,80) + theme_light() + xlab("Δ GDT_TS scores between all pairs within family") + ylab(paste0("Δ TCS scores subtracted by delta pLDDT scores between all pairs within family - ",ref_aligner,"_AF2")) # + stat_quadrant_counts() + ylim(-80,80)
  #ggsave(paste0("delta_gdt_ts_vs_delta_tcs_subtracted_by_plddt_ref_",ref_aligner,"_alphafold.png"),dpi="retina")

  #
  # Deltas all against all
  #

  delta_tcs_all_vs_all = intra_element_diff(merged_df$TCS) #apply(combn(tmp_df$TCS_mTMalign_AF2,2), 2, diff)
  delta_sop_all_vs_all = intra_element_diff(merged_df$SoP) #apply(combn(tmp_df$SoP_mTMalign_AF2,2), 2, diff)
  delta_gdt_ts_all_vs_all = intra_element_diff(merged_df$GDT_TS) #apply(combn(tmp_df$GDT_TS,2), 2, diff)
  delta_plddt_all_vs_all = intra_element_diff(merged_df$pLDDT) #apply(combn(tmp_df$pLDDT,2), 2, diff)
  delta_nirmsd_all_vs_all = intra_element_diff(merged_df$niRMSD)*(-1)
  delta_df_all_vs_all = data.frame(delta_tcs_all_vs_all,delta_sop_all_vs_all,delta_gdt_ts_all_vs_all,delta_plddt_all_vs_all,delta_nirmsd_all_vs_all)

  sop_all_vs_all_title = paste0("Δ Sum-of-Pairs scores between ", included_struct," pairs - ",ref_aligner,"_AF2 vs ",ref_aligner,"_NAT")
  tcs_all_vs_all_title = paste0("Δ TCS scores between ",included_struct, " pairs - ",ref_aligner,"_AF2")
  gdt_all_vs_all_title = paste0("Δ GDT_TS scores between ",included_struct, " pairs")
  plddt_all_vs_all_title = paste0("Δ pLDDT scores between ", included_struct," pairs")
  nirmsd_all_vs_all_title = paste0("Δ niRMSD scores between ", included_struct, " pairs - ",ref_aligner,"_AF2 with AF2 structures")

  p = ggplot(delta_df_all_vs_all,aes(x=delta_sop_all_vs_all,y=delta_tcs_all_vs_all)) + geom_point(color="darkblue",shape=1) + stat_quadrant_counts() + xlim(-100,100) + ylim(-100,100) + theme_light() + xlab(sop_all_vs_all_title) + ylab(tcs_all_vs_all_title) + theme(axis.title = element_text(size = 9)) #+ stat_cor(inherit.aes = FALSE,data=delta_df_all_vs_all,aes(x=delta_sop_all_vs_all,y=delta_tcs_all_vs_all),method = "pearson")
  #ggsave(paste0("delta_sop_vs_delta_tcs_",included_struct,"_ref_",ref_aligner,"_all_vs_all_alphafold.png"),dpi="retina")
  p = ggplot(delta_df_all_vs_all,aes(x=delta_sop_all_vs_all,y=delta_plddt_all_vs_all)) + geom_point(color="darkblue",shape=1) + stat_quadrant_counts() + xlim(-100,100) + ylim(-100,100) + theme_light() + xlab(sop_all_vs_all_title) + ylab(plddt_all_vs_all_title) + theme(axis.title = element_text(size = 9)) #+ stat_cor(inherit.aes = FALSE,data=delta_df_all_vs_all,aes(x=delta_sop_all_vs_all,y=delta_plddt_all_vs_all),method = "pearson")
  #ggsave(paste0("delta_sop_vs_delta_plddt_",included_struct,"_ref_",ref_aligner,"_all_vs_all_alphafold.png"),dpi="retina")
  p = ggplot(delta_df_all_vs_all,aes(x=delta_sop_all_vs_all,y=delta_gdt_ts_all_vs_all)) + geom_point(color="darkblue",shape=1) + stat_quadrant_counts() + xlim(-100,100) + ylim(-100,100) + theme_light() + xlab(sop_all_vs_all_title) + ylab(gdt_all_vs_all_title) + theme(axis.title = element_text(size = 9)) #+ stat_cor(inherit.aes = FALSE,data=delta_df_all_vs_all,aes(x=delta_sop_all_vs_all,y=delta_gdt_ts_all_vs_all),method = "pearson")
  #ggsave(paste0("delta_sop_vs_delta_gdt_ts_",included_struct,"_ref_",ref_aligner,"_all_vs_all_alphafold.png"),dpi="retina")

  p = ggplot(delta_df_all_vs_all,aes(x=delta_nirmsd_all_vs_all,y=delta_tcs_all_vs_all)) + geom_point(color="darkblue",shape=1) + stat_quadrant_counts() + xlim(-3,3) + ylim(-100,100) + theme_light() + xlab(nirmsd_all_vs_all_title) + ylab(tcs_all_vs_all_title) + theme(axis.title = element_text(size = 9)) #+ xlim(-100,100) + stat_cor(inherit.aes = FALSE,data=delta_df_all_vs_all,aes(x=delta_sop_all_vs_all,y=delta_tcs_all_vs_all),method = "pearson")
  #ggsave(paste0("delta_nirmsd_vs_delta_tcs_",included_struct,"_ref_",ref_aligner,"_all_vs_all_alphafold.png"),dpi="retina")
  p = ggplot(delta_df_all_vs_all,aes(x=delta_nirmsd_all_vs_all,y=delta_plddt_all_vs_all)) + geom_point(color="darkblue",shape=1) + stat_quadrant_counts() + xlim(-3,3) + ylim(-100,100) + theme_light() + xlab(nirmsd_all_vs_all_title) + ylab(plddt_all_vs_all_title) + theme(axis.title = element_text(size = 9)) #+ xlim(-100,100) + stat_cor(inherit.aes = FALSE,data=delta_df_all_vs_all,aes(x=delta_sop_all_vs_all,y=delta_plddt_all_vs_all),method = "pearson")
  #ggsave(paste0("delta_nirmsd_vs_delta_plddt_",included_struct,"_ref_",ref_aligner,"_all_vs_all_alphafold.png"),dpi="retina")
  p = ggplot(delta_df_all_vs_all,aes(x=delta_nirmsd_all_vs_all,y=delta_gdt_ts_all_vs_all)) + geom_point(color="darkblue",shape=1) + stat_quadrant_counts() + xlim(-3,3) + ylim(-100,100) + theme_light() + xlab(nirmsd_all_vs_all_title) + ylab(gdt_all_vs_all_title) + theme(axis.title = element_text(size = 9)) #+ xlim(-100,100) + stat_cor(inherit.aes = FALSE,data=delta_df_all_vs_all,aes(x=delta_sop_all_vs_all,y=delta_gdt_ts_all_vs_all),method = "pearson")
  #ggsave(paste0("delta_nirmsd_vs_delta_gdt_ts_",included_struct,"_ref_",ref_aligner,"_all_vs_all_alphafold.png"),dpi="retina")

  p = ggplot(delta_df_all_vs_all,aes(x=delta_gdt_ts_all_vs_all,y=delta_tcs_all_vs_all)) + geom_point(color="darkblue",shape=1) + stat_quadrant_counts() + xlim(-100,100) + ylim(-100,100) + theme_light() + xlab(gdt_all_vs_all_title) + ylab(tcs_all_vs_all_title) + theme(axis.title = element_text(size = 9)) #+ stat_cor(inherit.aes = FALSE,data=delta_df_all_vs_all,aes(x=delta_gdt_ts_all_vs_all,y=delta_tcs_all_vs_all),method = "pearson")
  #ggsave(paste0("delta_gdt_ts_vs_delta_tcs_",included_struct,"_ref_",ref_aligner,"_all_vs_all_alphafold.png"),dpi="retina")

  p = ggplot(delta_df_all_vs_all,aes(x=delta_gdt_ts_all_vs_all,y=delta_plddt_all_vs_all)) + geom_point(color="darkblue",shape=1) + stat_quadrant_counts() + xlim(-100,100) + ylim(-100,100) + theme_light() + xlab(gdt_all_vs_all_title) + ylab(plddt_all_vs_all_title) + theme(axis.title = element_text(size = 9)) #+ stat_cor(inherit.aes = FALSE,data=delta_df_all_vs_all,aes(x=delta_gdt_ts_all_vs_all,y=delta_plddt_all_vs_all),method = "pearson")
  #ggsave(paste0("delta_gdt_ts_vs_delta_plddt_",included_struct,"_all_vs_all_alphafold.png"),dpi="retina")

}


delta_analysis_DMPFOLD = function (tcs_df,gdt_ts_df,sop_df,nirmsd_df,ref_aligner,included_struct) {
  merged_df = merge(merge(merge(tcs_df,gdt_ts_df,by=c("Sequence","Family")),sop_df,by=c("Sequence","Family")),nirmsd_df,by=c("Sequence","Family"))
  merged_df = unique(merged_df)
  colnames(merged_df)[4] = "GDT_TS"
  merged_df$GDT_TS = merged_df$GDT_TS*100
  merged_df$TCS = merged_df$TCS/10
  merged_df = merged_df[order(merged_df$Family),]
  

  #
  # Deltas per family
  #

  for (fam in levels(factor(merged_df$Family))) {
    tmp_df = merged_df[which(merged_df$Family == fam),]
    delta_tcs = intra_element_diff(tmp_df$TCS)
    delta_sop = intra_element_diff(tmp_df$SoP)
    delta_gdt_ts = intra_element_diff(tmp_df$GDT_TS)
    delta_nirmsd = intra_element_diff(tmp_df$niRMSD)*(-1)
    if (exists('delta_df') && is.data.frame(get('delta_df'))){
      delta_df = rbind(delta_df,data.frame(delta_tcs,delta_sop,delta_gdt_ts,delta_nirmsd,Family=rep(fam,length(delta_gdt_ts))))
    } else {
      delta_df = data.frame(delta_tcs,delta_sop,delta_gdt_ts,delta_nirmsd,Family=rep(fam,length(delta_gdt_ts)))
    }
  }


  sop_title = paste0("Δ Sum-of-Pairs scores between ", included_struct," pairs within family - ",ref_aligner,"_DMPFOLD vs ",ref_aligner,"_NAT")
  tcs_title = paste0("Δ TCS scores between ",included_struct, " pairs within family - ",ref_aligner,"_DMPFOLD")
  gdt_title = paste0("Δ GDT_TS scores between ",included_struct, " pairs within family")
  nirmsd_title = paste0("Δ niRMSD scores between ", included_struct, " pairs within family - ",ref_aligner,"_DMPFOLD with DMPFOLD structures")

  p = ggplot(delta_df,aes(x=delta_sop,y=delta_tcs,color=Family)) + geom_point(shape=1) + stat_quadrant_counts() + xlim(-100,100) + ylim(-100,100) + theme_light() + xlab(sop_title) + ylab(tcs_title) + theme(axis.title = element_text(size = 9))
  ggsave(paste0("delta_sop_vs_delta_tcs_",included_struct,"_ref_",ref_aligner,"_dmpfold.png"),dpi="retina")
  p = ggplot(delta_df,aes(x=delta_sop,y=delta_gdt_ts,color=Family)) + geom_point(shape=1) + stat_quadrant_counts() + xlim(-100,100) + ylim(-100,100) + theme_light() + xlab(sop_title) + ylab(gdt_title) + theme(axis.title = element_text(size = 9))
  ggsave(paste0("delta_sop_vs_delta_gdt_ts_",included_struct,"_ref_",ref_aligner,"_dmpfold.png"),dpi="retina")
  
  p = ggplot(delta_df,aes(x=delta_nirmsd,y=delta_tcs,color=Family)) + geom_point(shape=1) + stat_quadrant_counts() + xlim(-35,35) + ylim(-100,100) + theme_light() + xlab(nirmsd_title) + ylab(tcs_title) + theme(axis.title = element_text(size = 9)) #+ xlim(-100,100)
  ggsave(paste0("delta_nirmsd_vs_delta_tcs_",included_struct,"_ref_",ref_aligner,"_dmpfold.png"),dpi="retina")
  p = ggplot(delta_df,aes(x=delta_nirmsd,y=delta_gdt_ts,color=Family)) + geom_point(shape=1) + stat_quadrant_counts() + xlim(-35,35) + ylim(-100,100) + theme_light() + xlab(nirmsd_title) + ylab(gdt_title) + theme(axis.title = element_text(size = 9)) #+ xlim(-100,100)
  ggsave(paste0("delta_nirmsd_vs_delta_gdt_ts_",included_struct,"_ref_",ref_aligner,"_dmpfold.png"),dpi="retina")

  p = ggplot(delta_df,aes(x=delta_gdt_ts,y=delta_tcs,color=Family)) + geom_point(shape=1) + stat_quadrant_counts() + xlim(-100,100) + ylim(-100,100) + theme_light() + xlab(gdt_title) + ylab(tcs_title) + theme(axis.title = element_text(size = 9))
  ggsave(paste0("delta_gdt_ts_vs_delta_tcs_",included_struct,"_ref_",ref_aligner,"_dmpfold.png"),dpi="retina")

  return(delta_df)

  #
  # Deltas all against all
  #

  delta_tcs_all_vs_all = intra_element_diff(merged_df$TCS)
  delta_sop_all_vs_all = intra_element_diff(merged_df$SoP)
  delta_gdt_ts_all_vs_all = intra_element_diff(merged_df$GDT_TS)
  delta_nirmsd_all_vs_all = intra_element_diff(merged_df$niRMSD)*(-1)
  delta_df_all_vs_all = data.frame(delta_tcs_all_vs_all,delta_sop_all_vs_all,delta_gdt_ts_all_vs_all,delta_nirmsd_all_vs_all)

  sop_all_vs_all_title = paste0("Δ Sum-of-Pairs scores between ", included_struct," pairs - ",ref_aligner,"_DMPFOLD vs ",ref_aligner,"_NAT")
  tcs_all_vs_all_title = paste0("Δ TCS scores between ",included_struct, " pairs - ",ref_aligner,"_DMPFOLD")
  gdt_all_vs_all_title = paste0("Δ GDT_TS scores between ",included_struct, " pairs")
  nirmsd_all_vs_all_title = paste0("Δ niRMSD scores between ", included_struct, " pairs - ",ref_aligner,"_DMPFOLD with DMPFOLD structures")

  p = ggplot(delta_df_all_vs_all,aes(x=delta_sop_all_vs_all,y=delta_tcs_all_vs_all)) + geom_point(color="darkblue",shape=1) + stat_quadrant_counts() + xlim(-100,100) + ylim(-100,100) + theme_light() + xlab(sop_all_vs_all_title) + ylab(tcs_all_vs_all_title) + theme(axis.title = element_text(size = 9))
  #ggsave(paste0("delta_sop_vs_delta_tcs_",included_struct,"_ref_",ref_aligner,"_all_vs_all_dmpfold.png"),dpi="retina")
  p = ggplot(delta_df_all_vs_all,aes(x=delta_sop_all_vs_all,y=delta_gdt_ts_all_vs_all)) + geom_point(color="darkblue",shape=1) + stat_quadrant_counts() + xlim(-100,100) + ylim(-100,100) + theme_light() + xlab(sop_all_vs_all_title) + ylab(gdt_all_vs_all_title) + theme(axis.title = element_text(size = 9))
  #ggsave(paste0("delta_sop_vs_delta_gdt_ts_",included_struct,"_ref_",ref_aligner,"_all_vs_all_dmpfold.png"),dpi="retina")

  p = ggplot(delta_df_all_vs_all,aes(x=delta_nirmsd_all_vs_all,y=delta_tcs_all_vs_all)) + geom_point(color="darkblue",shape=1) + stat_quadrant_counts() + xlim(-35,35) + ylim(-100,100) + theme_light() + xlab(nirmsd_all_vs_all_title) + ylab(tcs_all_vs_all_title) + theme(axis.title = element_text(size = 9))
  #ggsave(paste0("delta_nirmsd_vs_delta_tcs_",included_struct,"_ref_",ref_aligner,"_all_vs_all_dmpfold.png"),dpi="retina")
  p = ggplot(delta_df_all_vs_all,aes(x=delta_nirmsd_all_vs_all,y=delta_gdt_ts_all_vs_all)) + geom_point(color="darkblue",shape=1) + stat_quadrant_counts() + xlim(-35,35) + ylim(-100,100) + theme_light() + xlab(nirmsd_all_vs_all_title) + ylab(gdt_all_vs_all_title) + theme(axis.title = element_text(size = 9))
  #ggsave(paste0("delta_nirmsd_vs_delta_gdt_ts_",included_struct,"_ref_",ref_aligner,"_all_vs_all_dmpfold.png"),dpi="retina")

  p = ggplot(delta_df_all_vs_all,aes(x=delta_gdt_ts_all_vs_all,y=delta_tcs_all_vs_all)) + geom_point(color="darkblue",shape=1) + stat_quadrant_counts() + xlim(-100,100) + ylim(-100,100) + theme_light() + xlab(gdt_all_vs_all_title) + ylab(tcs_all_vs_all_title) + theme(axis.title = element_text(size = 9))
  #ggsave(paste0("delta_gdt_ts_vs_delta_tcs_",included_struct,"_ref_",ref_aligner,"_all_vs_all_dmpfold.png"),dpi="retina")
}

#
# Input niRMSD scores
#

nirmsd_ref = read.table("./selected_comparisons_nirmsd.txt",header = F,stringsAsFactors = F)
colnames(nirmsd_ref)=c("Family","Famsa","Ginsi","MSAProbs","TCoffee","PSIcoffee","3DCoffee_NAT","3DCoffee_TMalign_NAT","mTMalign_NAT","3DCoffee_DMPFOLD_REF_NAT","3DCoffee_TMalign_DMPFOLD_REF_NAT","mTMalign_DMPFOLD_REF_NAT","3DCoffee_DMPFOLD_REF_DMPFOLD","3DCoffee_TMalign_DMPFOLD_REF_DMPFOLD","mTMalign_DMPFOLD_REF_DMPFOLD","3DCoffee_AF2_REF_NAT","3DCoffee_TMalign_AF2_REF_NAT","mTMalign_AF2_REF_NAT","3DCoffee_AF2_REF_AF2","3DCoffee_TMalign_AF2_REF_AF2","mTMalign_AF2_REF_AF2") # ,"Deepblast"
nirmsd_ref = rbind(nirmsd_ref,c("Average",colMeans(nirmsd_ref[,-1])))
row.names(nirmsd_ref) = make.unique(nirmsd_ref$Family)
nirmsd_ref[colnames(nirmsd_ref)[-1]] <- sapply(nirmsd_ref[colnames(nirmsd_ref)[-1]],as.numeric)
#pheatmap(nirmsd_ref[,c(2:9)],display_numbers = T,angle_col = 45,fontsize = 10,cluster_cols = F,cluster_rows = F,number_color="black",number_format="%.3f",fontsize_number=12,filename = "nirmsd_ref_NAT.png")

nirmsd_ref_avg = read.table("./selected_comparisons_nirmsd_avg.txt",header = F,stringsAsFactors = F)
nirmsd_ref_avg = nirmsd_ref_avg[!duplicated(as.list(nirmsd_ref_avg))]
colnames(nirmsd_ref_avg)=c("Sequence","Family","Famsa","Ginsi","MSAProbs","TCoffee","PSIcoffee","3DCoffee_NAT","3DCoffee_TMalign_NAT","mTMalign_NAT","3DCoffee_DMPFOLD_REF_NAT","3DCoffee_TMalign_DMPFOLD_REF_NAT","mTMalign_DMPFOLD_REF_NAT","3DCoffee_DMPFOLD_REF_DMPFOLD","3DCoffee_TMalign_DMPFOLD_REF_DMPFOLD","mTMalign_DMPFOLD_REF_DMPFOLD","3DCoffee_AF2_REF_NAT","3DCoffee_TMalign_AF2_REF_NAT","mTMalign_AF2_REF_NAT","3DCoffee_AF2_REF_AF2","3DCoffee_TMalign_AF2_REF_AF2","mTMalign_AF2_REF_AF2") # ,"Deepblast"

nirmsd_ref_pair = read.table("./selected_comparisons_nirmsd_pair.txt",header = F,stringsAsFactors = F)
nirmsd_ref_pair = nirmsd_ref_pair[!duplicated(as.list(nirmsd_ref_pair))]
colnames(nirmsd_ref_pair)=c("Sequence_1","Sequence_2","Family","Famsa","Ginsi","MSAProbs","TCoffee","PSIcoffee","3DCoffee_NAT","3DCoffee_TMalign_NAT","mTMalign_NAT","3DCoffee_DMPFOLD_REF_NAT","3DCoffee_TMalign_DMPFOLD_REF_NAT","mTMalign_DMPFOLD_REF_NAT","3DCoffee_DMPFOLD_REF_DMPFOLD","3DCoffee_TMalign_DMPFOLD_REF_DMPFOLD","mTMalign_DMPFOLD_REF_DMPFOLD","3DCoffee_AF2_REF_NAT","3DCoffee_TMalign_AF2_REF_NAT","mTMalign_AF2_REF_NAT","3DCoffee_AF2_REF_AF2","3DCoffee_TMalign_AF2_REF_AF2","mTMalign_AF2_REF_AF2") # ,"Deepblast"

#
# Input pid scores
#

pid_ref = read.table("./selected_comparisons_pid.txt",header = F,stringsAsFactors = F)
pid_ref = pid_ref[!duplicated(as.list(pid_ref))]
colnames(pid_ref)=c("Family","Famsa","Ginsi","MSAProbs","TCoffee","PSIcoffee","3DCoffee_NAT","3DCoffee_TMalign_NAT","mTMalign_NAT","3DCoffee_DMPFOLD","3DCoffee_TMalign_DMPFOLD","mTMalign_DMPFOLD","3DCoffee_AF2","3DCoffee_TMalign_AF2","mTMalign_AF2")

pid_ref_avg = read.table("./selected_comparisons_pid.avg.txt",header = F,stringsAsFactors = F)
pid_ref_avg = pid_ref_avg[!duplicated(as.list(pid_ref_avg))]
colnames(pid_ref_avg)=c("Sequence","Famsa","Family","Ginsi","MSAProbs","TCoffee","PSIcoffee","3DCoffee_NAT","3DCoffee_TMalign_NAT","mTMalign_NAT","3DCoffee_DMPFOLD","3DCoffee_TMalign_DMPFOLD","mTMalign_DMPFOLD","3DCoffee_AF2","3DCoffee_TMalign_AF2","mTMalign_AF2")

pid_ref_pair = read.table("./selected_comparisons_pid.pair.txt",header = F,stringsAsFactors = F)
pid_ref_pair = pid_ref_pair[!duplicated(as.list(pid_ref_pair))]
colnames(pid_ref_pair)=c("Sequence_1","Sequence_2","Famsa","Family","Ginsi","MSAProbs","TCoffee","PSIcoffee","3DCoffee_NAT","3DCoffee_TMalign_NAT","mTMalign_NAT","3DCoffee_DMPFOLD","3DCoffee_TMalign_DMPFOLD","mTMalign_DMPFOLD","3DCoffee_AF2","3DCoffee_TMalign_AF2","mTMalign_AF2")



benchfam_colors = colorRampPalette(c("green","blue"))(13)
pdb_colors = colorRampPalette(c("yellow2","red"))(12)
color_num = nrow(nirmsd_ref)-1
if (color_num > 15) {
  color_codes = c(benchfam_colors,pdb_colors)
} else {
  color_codes <- colorRampPalette(brewer.pal(8, "Set1"))(color_num)
}
color_per_fam = data.frame(row.names(nirmsd_ref)[-nrow(nirmsd_ref)],color_codes)
colnames(color_per_fam) = c("Family","Color")
color_codes_for_circ = as.character(color_per_fam[order(color_per_fam$Family),]$Color)

#
# Input SoP scores
#

ref_3dcoffee_sp = read.table("./selected_comparisons_ref_3dcoffee_sp.txt",header = F,stringsAsFactors = F)
ref_3dcoffee_sp = ref_3dcoffee_sp[!duplicated(as.list(ref_3dcoffee_sp))]
colnames(ref_3dcoffee_sp)=c("Family","Famsa","Ginsi","MSAProbs","TCoffee","PSIcoffee","3DCoffee_NAT","3DCoffee_TMalign_NAT","mTMalign_NAT","3DCoffee_DMPFOLD","3DCoffee_TMalign_DMPFOLD","mTMalign_DMPFOLD","3DCoffee_AF2","3DCoffee_TMalign_AF2","mTMalign_AF2") #,"Deepblast_vs_REF"
write.table(ref_3dcoffee_sp, file="selected_comparisons_ref_3dcoffee_sp.tsv",quote=FALSE, sep="\t",row.names=FALSE)

ref_3dcoffee_avg_sp = read.table("./selected_comparisons_ref_3dcoffee_avg_sp.txt",header = F,stringsAsFactors = F)
ref_3dcoffee_avg_sp = ref_3dcoffee_avg_sp[!duplicated(as.list(ref_3dcoffee_avg_sp))]
colnames(ref_3dcoffee_avg_sp)= c("Sequence","Family","Famsa","Ginsi","MSAProbs","TCoffee","PSIcoffee","3DCoffee_NAT","3DCoffee_TMalign_NAT","mTMalign_NAT","3DCoffee_DMPFOLD","3DCoffee_TMalign_DMPFOLD","mTMalign_DMPFOLD","3DCoffee_AF2","3DCoffee_TMalign_AF2","mTMalign_AF2")
write.table(ref_3dcoffee_avg_sp, file="selected_comparisons_ref_3dcoffee_avg_sp.tsv",quote=FALSE, sep="\t",row.names=FALSE)

ref_3dcoffee_pair_sp = read.table("./selected_comparisons_ref_3dcoffee_pair_sp.txt",header = F,stringsAsFactors = F)
ref_3dcoffee_pair_sp = ref_3dcoffee_pair_sp[!duplicated(as.list(ref_3dcoffee_pair_sp))]
colnames(ref_3dcoffee_pair_sp)=c("Sequence_1","Sequence_2","Family","Famsa","Ginsi","MSAProbs","TCoffee","PSIcoffee","3DCoffee_NAT","3DCoffee_TMalign_NAT","mTMalign_NAT","3DCoffee_DMPFOLD","3DCoffee_TMalign_DMPFOLD","mTMalign_DMPFOLD","3DCoffee_AF2","3DCoffee_TMalign_AF2","mTMalign_AF2")

  #
  # SoP without loops
  #
ref_3dcoffee_sp_without_loops = read.table("./selected_comparisons_ref_3dcoffee_sp.without_loops.txt",header = F,stringsAsFactors = F)
ref_3dcoffee_sp_without_loops = ref_3dcoffee_sp_without_loops[!duplicated(as.list(ref_3dcoffee_sp_without_loops))]
colnames(ref_3dcoffee_sp_without_loops)=c("Family","Famsa","Ginsi","MSAProbs","TCoffee","PSIcoffee","3DCoffee_NAT","3DCoffee_TMalign_NAT","mTMalign_NAT","3DCoffee_DMPFOLD","3DCoffee_TMalign_DMPFOLD","mTMalign_DMPFOLD","3DCoffee_AF2","3DCoffee_TMalign_AF2","mTMalign_AF2") #,"Deepblast_vs_REF"
write.table(ref_3dcoffee_sp_without_loops, file="selected_comparisons_ref_3dcoffee_sp.without_loops.tsv",quote=FALSE, sep="\t",row.names=FALSE)

ref_3dcoffee_avg_sp_without_loops = read.table("./selected_comparisons_ref_3dcoffee_avg_sp.without_loops.txt",header = F,stringsAsFactors = F)
ref_3dcoffee_avg_sp_without_loops = ref_3dcoffee_avg_sp_without_loops[!duplicated(as.list(ref_3dcoffee_avg_sp_without_loops))]
colnames(ref_3dcoffee_avg_sp_without_loops)= c("Sequence","Family","Famsa","Ginsi","MSAProbs","TCoffee","PSIcoffee","3DCoffee_NAT","3DCoffee_TMalign_NAT","mTMalign_NAT","3DCoffee_DMPFOLD","3DCoffee_TMalign_DMPFOLD","mTMalign_DMPFOLD","3DCoffee_AF2","3DCoffee_TMalign_AF2","mTMalign_AF2")
write.table(ref_3dcoffee_avg_sp_without_loops, file="selected_comparisons_ref_3dcoffee_avg_sp_without_loops.tsv",quote=FALSE, sep="\t",row.names=FALSE)

ref_3dcoffee_pair_sp_without_loops = read.table("./selected_comparisons_ref_3dcoffee_pair_sp.without_loops.txt",header = F,stringsAsFactors = F)
ref_3dcoffee_pair_sp_without_loops = ref_3dcoffee_pair_sp_without_loops[!duplicated(as.list(ref_3dcoffee_pair_sp_without_loops))]
colnames(ref_3dcoffee_pair_sp_without_loops)=c("Sequence_1","Sequence_2","Family","Famsa","Ginsi","MSAProbs","TCoffee","PSIcoffee","3DCoffee_NAT","3DCoffee_TMalign_NAT","mTMalign_NAT","3DCoffee_DMPFOLD","3DCoffee_TMalign_DMPFOLD","mTMalign_DMPFOLD","3DCoffee_AF2","3DCoffee_TMalign_AF2","mTMalign_AF2")

  #
  # SoP same state
  #
ref_3dcoffee_sp_same_state = read.table("./selected_comparisons_ref_3dcoffee_sp.same_state.txt",header = F,stringsAsFactors = F)
ref_3dcoffee_sp_same_state = ref_3dcoffee_sp_same_state[!duplicated(as.list(ref_3dcoffee_sp_same_state))]
colnames(ref_3dcoffee_sp_same_state)=c("Family","Famsa","Ginsi","MSAProbs","TCoffee","PSIcoffee","3DCoffee_NAT","3DCoffee_TMalign_NAT","mTMalign_NAT","3DCoffee_DMPFOLD","3DCoffee_TMalign_DMPFOLD","mTMalign_DMPFOLD","3DCoffee_AF2","3DCoffee_TMalign_AF2","mTMalign_AF2") #,"Deepblast_vs_REF"
write.table(ref_3dcoffee_sp_same_state, file="selected_comparisons_ref_3dcoffee_sp.same_state.tsv",quote=FALSE, sep="\t",row.names=FALSE)

ref_3dcoffee_avg_sp_same_state = read.table("./selected_comparisons_ref_3dcoffee_avg_sp.same_state.txt",header = F,stringsAsFactors = F)
ref_3dcoffee_avg_sp_same_state = ref_3dcoffee_avg_sp_same_state[!duplicated(as.list(ref_3dcoffee_avg_sp_same_state))]
colnames(ref_3dcoffee_avg_sp_same_state)= c("Sequence","Family","Famsa","Ginsi","MSAProbs","TCoffee","PSIcoffee","3DCoffee_NAT","3DCoffee_TMalign_NAT","mTMalign_NAT","3DCoffee_DMPFOLD","3DCoffee_TMalign_DMPFOLD","mTMalign_DMPFOLD","3DCoffee_AF2","3DCoffee_TMalign_AF2","mTMalign_AF2")
write.table(ref_3dcoffee_avg_sp_same_state, file="selected_comparisons_ref_3dcoffee_avg_sp.same_state.tsv",quote=FALSE, sep="\t",row.names=FALSE)

ref_3dcoffee_pair_sp_same_state = read.table("./selected_comparisons_ref_3dcoffee_pair_sp.same_state.txt",header = F,stringsAsFactors = F)
ref_3dcoffee_pair_sp_same_state = ref_3dcoffee_pair_sp_same_state[!duplicated(as.list(ref_3dcoffee_pair_sp_same_state))]
colnames(ref_3dcoffee_pair_sp_same_state)=c("Sequence_1","Sequence_2","Family","Famsa","Ginsi","MSAProbs","TCoffee","PSIcoffee","3DCoffee_NAT","3DCoffee_TMalign_NAT","mTMalign_NAT","3DCoffee_DMPFOLD","3DCoffee_TMalign_DMPFOLD","mTMalign_DMPFOLD","3DCoffee_AF2","3DCoffee_TMalign_AF2","mTMalign_AF2")


ref_3dcoffee_pair_sp$Family = factor(ref_3dcoffee_pair_sp$Family,levels=row.names(nirmsd_ref)[-nrow(nirmsd_ref)])
o=ggplot(ref_3dcoffee_pair_sp,aes(x=`3DCoffee_DMPFOLD`,y=PSIcoffee,color=Family)) + theme_light() + theme(legend.position="top") +geom_point()+xlim(0,100)+ylim(0,100)+geom_abline(intercept =0 , slope = 1,lwd=0.2)+xlab("Sum-of-Pairs score between pairs of sequences - 3DCoffee_DMPFOLD vs 3DCoffee_NAT (%)")+ylab("Sum-of-Pairs score between pairs of sequences - PSICoffee vs 3DCoffee_NAT (%)") + scale_color_manual(values=color_codes) + ggtitle("Sum-of-Pairs score between pairs of sequences with TCS > 0")
#ggsave("ref_3dcoffee_pair_sp_score_per_family_dmpfold_vs_psicoffee.png",dpi = "retina")
o=ggplot(ref_3dcoffee_pair_sp,aes(x=`3DCoffee_AF2`,y=PSIcoffee,color=Family)) + theme_light() + theme(legend.position="top") +geom_point()+xlim(0,100)+ylim(0,100)+geom_abline(intercept =0 , slope = 1,lwd=0.2)+xlab("Sum-of-Pairs score between pairs of sequences - 3DCoffee_AF2 vs 3DCoffee_NAT (%)")+ylab("Sum-of-Pairs score between pairs of sequences - PSICoffee vs 3DCoffee_NAT (%)") + scale_color_manual(values=color_codes) + ggtitle("Sum-of-Pairs score between pairs of sequences with TCS > 0")
#ggsave("ref_3dcoffee_pair_sp_score_per_family_alphafold_vs_psicoffee.png",dpi = "retina")
o=ggplot(ref_3dcoffee_pair_sp,aes(x=`3DCoffee_DMPFOLD`,y=Ginsi,color=Family)) + theme_light() + theme(legend.position="top") +geom_point()+xlim(0,100)+ylim(0,100)+geom_abline(intercept =0 , slope = 1,lwd=0.2)+xlab("Sum-of-Pairs score between pairs of sequences - 3DCoffee_DMPFOLD vs 3DCoffee_NAT (%)")+ylab("Sum-of-Pairs score between pairs of sequences - Ginsi vs 3DCoffee_NAT (%)") + scale_color_manual(values=color_codes) + ggtitle("Sum-of-Pairs score between pairs of sequences with TCS > 0")
#ggsave("ref_3dcoffee_pair_sp_score_per_family_dmpfold_vs_ginsi.png",dpi = "retina")
o=ggplot(ref_3dcoffee_pair_sp,aes(x=`3DCoffee_AF2`,y=Ginsi,color=Family)) + theme_light() + theme(legend.position="top") +geom_point()+xlim(0,100)+ylim(0,100)+geom_abline(intercept =0 , slope = 1,lwd=0.2)+xlab("Sum-of-Pairs score between pairs of sequences - 3DCoffee_AF2 vs 3DCoffee_NAT (%)")+ylab("Sum-of-Pairs score between pairs of sequences - Ginsi vs 3DCoffee_NAT (%)") + scale_color_manual(values=color_codes) + ggtitle("Sum-of-Pairs score between pairs of sequences with TCS > 0")
#ggsave("ref_3dcoffee_pair_sp_score_per_family_alphafold_vs_ginsi.png",dpi = "retina")

ref_3dcoffee_pair_sp_same_state$Family = factor(ref_3dcoffee_pair_sp_same_state$Family,levels=row.names(nirmsd_ref)[-nrow(nirmsd_ref)])
o=ggplot(ref_3dcoffee_pair_sp_same_state,aes(x=`3DCoffee_DMPFOLD`,y=PSIcoffee,color=Family)) + theme_light() + theme(legend.position="top") +geom_point()+xlim(0,100)+ylim(0,100)+geom_abline(intercept =0 , slope = 1,lwd=0.2)+xlab("Sum-of-Pairs score between pairs of sequences - 3DCoffee_DMPFOLD vs 3DCoffee_NAT (%)")+ylab("Sum-of-Pairs score between pairs of sequences - PSICoffee vs 3DCoffee_NAT (%)") + scale_color_manual(values=color_codes) + ggtitle("Sum-of-Pairs score between pairs of sequences with TCS > 0")
ggsave("ref_3dcoffee_pair_sp_score_per_family_dmpfold_vs_psicoffee.png",dpi = "retina")
o=ggplot(ref_3dcoffee_pair_sp_same_state,aes(x=`3DCoffee_AF2`,y=PSIcoffee,color=Family)) + theme_light() + theme(legend.position="top") +geom_point()+xlim(0,100)+ylim(0,100)+geom_abline(intercept =0 , slope = 1,lwd=0.2)+xlab("Sum-of-Pairs score between pairs of sequences - 3DCoffee_AF2 vs 3DCoffee_NAT (%)")+ylab("Sum-of-Pairs score between pairs of sequences - PSICoffee vs 3DCoffee_NAT (%)") + scale_color_manual(values=color_codes) + ggtitle("Sum-of-Pairs score between pairs of sequences with TCS > 0")
ggsave("ref_3dcoffee_pair_sp_score_per_family_alphafold_vs_psicoffee.png",dpi = "retina")
o=ggplot(ref_3dcoffee_pair_sp_same_state,aes(x=`3DCoffee_DMPFOLD`,y=Ginsi,color=Family)) + theme_light() + theme(legend.position="top") +geom_point()+xlim(0,100)+ylim(0,100)+geom_abline(intercept =0 , slope = 1,lwd=0.2)+xlab("Sum-of-Pairs score between pairs of sequences - 3DCoffee_DMPFOLD vs 3DCoffee_NAT (%)")+ylab("Sum-of-Pairs score between pairs of sequences - Ginsi vs 3DCoffee_NAT (%)") + scale_color_manual(values=color_codes) + ggtitle("Sum-of-Pairs score between pairs of sequences with TCS > 0")
ggsave("ref_3dcoffee_pair_sp_score_per_family_dmpfold_vs_ginsi.png",dpi = "retina")
o=ggplot(ref_3dcoffee_pair_sp_same_state,aes(x=`3DCoffee_AF2`,y=Ginsi,color=Family)) + theme_light() + theme(legend.position="bottom") + geom_point(alpha=0.6)+xlim(0,100)+ylim(0,100)+geom_abline(intercept =0 , slope = 1,lwd=0.2)+xlab("SoP score on pairs of sequences - MSA-AF2 vs MSA-PDB (%)")+ylab("SoP score on pairs of sequences - MSA-Seq vs MSA-PDB (%)") + scale_color_manual(values=color_codes) #+ ggtitle("Sum-of-Pairs score between pairs of sequences with TCS > 0")
ggsave("ref_3dcoffee_pair_sp_score_per_family_alphafold_vs_ginsi.png",ggMarginal(o, type = "density",fill = "darkblue"),dpi = 700)

perc_of_pairs_msa_af2_superior_to_msa_seq = round(length(which(ref_3dcoffee_pair_sp_same_state$`3DCoffee_AF2` > ref_3dcoffee_pair_sp_same_state$Ginsi))/nrow(ref_3dcoffee_pair_sp_same_state),digits=2)


ref_mtmalign_sp = read.table("./selected_comparisons_ref_mtmalign_sp.txt",header = F,stringsAsFactors = F)
ref_mtmalign_sp = ref_mtmalign_sp[!duplicated(as.list(ref_mtmalign_sp))]
colnames(ref_mtmalign_sp)=c("Family","Famsa","Ginsi","MSAProbs","TCoffee","PSIcoffee","3DCoffee_NAT","3DCoffee_TMalign_NAT","mTMalign_NAT","3DCoffee_DMPFOLD","3DCoffee_TMalign_DMPFOLD","mTMalign_DMPFOLD","3DCoffee_AF2","3DCoffee_TMalign_AF2","mTMalign_AF2") #,"Deepblast_vs_REF"
write.table(ref_mtmalign_sp, file="selected_comparisons_ref_mtmalign_sp.tsv",quote=FALSE, sep="\t",row.names=FALSE)

ref_mtmalign_avg_sp = read.table("./selected_comparisons_ref_mtmalign_avg_sp.txt",header = F,stringsAsFactors = F)
ref_mtmalign_avg_sp = ref_mtmalign_avg_sp[!duplicated(as.list(ref_mtmalign_avg_sp))]
colnames(ref_mtmalign_avg_sp)= c("Sequence","Family","Famsa","Ginsi","MSAProbs","TCoffee","PSIcoffee","3DCoffee_NAT","3DCoffee_TMalign_NAT","mTMalign_NAT","3DCoffee_DMPFOLD","3DCoffee_TMalign_DMPFOLD","mTMalign_DMPFOLD","3DCoffee_AF2","3DCoffee_TMalign_AF2","mTMalign_AF2")
write.table(ref_mtmalign_avg_sp, file="selected_comparisons_ref_mtmalign_avg_sp.tsv",quote=FALSE, sep="\t",row.names=FALSE)

ref_mtmalign_pair_sp = read.table("./selected_comparisons_ref_mtmalign_pair_sp.txt",header = F,stringsAsFactors = F)
ref_mtmalign_pair_sp = ref_mtmalign_pair_sp[!duplicated(as.list(ref_mtmalign_pair_sp))]
colnames(ref_mtmalign_pair_sp)=c("Sequence_1","Sequence_2","Family","Famsa","Ginsi","MSAProbs","TCoffee","PSIcoffee","3DCoffee_NAT","3DCoffee_TMalign_NAT","mTMalign_NAT","3DCoffee_DMPFOLD","3DCoffee_TMalign_DMPFOLD","mTMalign_DMPFOLD","3DCoffee_AF2","3DCoffee_TMalign_AF2","mTMalign_AF2")

  #
  # SoP without loops
  #
ref_mtmalign_sp_without_loops = read.table("./selected_comparisons_ref_mtmalign_sp.without_loops.txt",header = F,stringsAsFactors = F)
ref_mtmalign_sp_without_loops = ref_mtmalign_sp_without_loops[!duplicated(as.list(ref_mtmalign_sp_without_loops))]
colnames(ref_mtmalign_sp_without_loops)=c("Family","Famsa","Ginsi","MSAProbs","TCoffee","PSIcoffee","3DCoffee_NAT","3DCoffee_TMalign_NAT","mTMalign_NAT","3DCoffee_DMPFOLD","3DCoffee_TMalign_DMPFOLD","mTMalign_DMPFOLD","3DCoffee_AF2","3DCoffee_TMalign_AF2","mTMalign_AF2") #,"Deepblast_vs_REF"
write.table(ref_mtmalign_sp_without_loops, file="selected_comparisons_ref_mtmalign_sp.without_loops.tsv",quote=FALSE, sep="\t",row.names=FALSE)

ref_mtmalign_avg_sp_without_loops = read.table("./selected_comparisons_ref_mtmalign_avg_sp.without_loops.txt",header = F,stringsAsFactors = F)
ref_mtmalign_avg_sp_without_loops = ref_mtmalign_avg_sp_without_loops[!duplicated(as.list(ref_mtmalign_avg_sp_without_loops))]
colnames(ref_mtmalign_avg_sp_without_loops)= c("Sequence","Family","Famsa","Ginsi","MSAProbs","TCoffee","PSIcoffee","3DCoffee_NAT","3DCoffee_TMalign_NAT","mTMalign_NAT","3DCoffee_DMPFOLD","3DCoffee_TMalign_DMPFOLD","mTMalign_DMPFOLD","3DCoffee_AF2","3DCoffee_TMalign_AF2","mTMalign_AF2")
write.table(ref_mtmalign_avg_sp_without_loops, file="selected_comparisons_ref_mtmalign_avg_sp_without_loops.tsv",quote=FALSE, sep="\t",row.names=FALSE)

ref_mtmalign_pair_sp_without_loops = read.table("./selected_comparisons_ref_mtmalign_pair_sp.without_loops.txt",header = F,stringsAsFactors = F)
ref_mtmalign_pair_sp_without_loops = ref_mtmalign_pair_sp_without_loops[!duplicated(as.list(ref_mtmalign_pair_sp_without_loops))]
colnames(ref_mtmalign_pair_sp_without_loops)=c("Sequence_1","Sequence_2","Family","Famsa","Ginsi","MSAProbs","TCoffee","PSIcoffee","3DCoffee_NAT","3DCoffee_TMalign_NAT","mTMalign_NAT","3DCoffee_DMPFOLD","3DCoffee_TMalign_DMPFOLD","mTMalign_DMPFOLD","3DCoffee_AF2","3DCoffee_TMalign_AF2","mTMalign_AF2")

  #
  # SoP same state
  #
ref_mtmalign_sp_same_state = read.table("./selected_comparisons_ref_mtmalign_sp.same_state.txt",header = F,stringsAsFactors = F)
ref_mtmalign_sp_same_state = ref_mtmalign_sp_same_state[!duplicated(as.list(ref_mtmalign_sp_same_state))]
colnames(ref_mtmalign_sp_same_state)=c("Family","Famsa","Ginsi","MSAProbs","TCoffee","PSIcoffee","3DCoffee_NAT","3DCoffee_TMalign_NAT","mTMalign_NAT","3DCoffee_DMPFOLD","3DCoffee_TMalign_DMPFOLD","mTMalign_DMPFOLD","3DCoffee_AF2","3DCoffee_TMalign_AF2","mTMalign_AF2") #,"Deepblast_vs_REF"
write.table(ref_mtmalign_sp_same_state, file="selected_comparisons_ref_mtmalign_sp.same_state.tsv",quote=FALSE, sep="\t",row.names=FALSE)

ref_mtmalign_avg_sp_same_state = read.table("./selected_comparisons_ref_mtmalign_avg_sp.same_state.txt",header = F,stringsAsFactors = F)
ref_mtmalign_avg_sp_same_state = ref_mtmalign_avg_sp_same_state[!duplicated(as.list(ref_mtmalign_avg_sp_same_state))]
colnames(ref_mtmalign_avg_sp_same_state)= c("Sequence","Family","Famsa","Ginsi","MSAProbs","TCoffee","PSIcoffee","3DCoffee_NAT","3DCoffee_TMalign_NAT","mTMalign_NAT","3DCoffee_DMPFOLD","3DCoffee_TMalign_DMPFOLD","mTMalign_DMPFOLD","3DCoffee_AF2","3DCoffee_TMalign_AF2","mTMalign_AF2")
write.table(ref_mtmalign_avg_sp_same_state, file="selected_comparisons_ref_mtmalign_avg_sp.same_state.tsv",quote=FALSE, sep="\t",row.names=FALSE)

ref_mtmalign_pair_sp_same_state = read.table("./selected_comparisons_ref_mtmalign_pair_sp.same_state.txt",header = F,stringsAsFactors = F)
ref_mtmalign_pair_sp_same_state = ref_mtmalign_pair_sp_same_state[!duplicated(as.list(ref_mtmalign_pair_sp_same_state))]
colnames(ref_mtmalign_pair_sp_same_state)=c("Sequence_1","Sequence_2","Family","Famsa","Ginsi","MSAProbs","TCoffee","PSIcoffee","3DCoffee_NAT","3DCoffee_TMalign_NAT","mTMalign_NAT","3DCoffee_DMPFOLD","3DCoffee_TMalign_DMPFOLD","mTMalign_DMPFOLD","3DCoffee_AF2","3DCoffee_TMalign_AF2","mTMalign_AF2")


ref_mtmalign_pair_sp$Family = factor(ref_mtmalign_pair_sp$Family,levels=row.names(nirmsd_ref)[-nrow(nirmsd_ref)])
o=ggplot(ref_mtmalign_pair_sp,aes(x=`mTMalign_DMPFOLD`,y=PSIcoffee,color=Family)) + theme_light() + theme(legend.position="top") +geom_point()+xlim(0,100)+ylim(0,100)+geom_abline(intercept =0 , slope = 1,lwd=0.2)+xlab("Sum-of-Pairs score between pairs of sequences - mTMalign_DMPFOLD vs mTMalign_NAT (%)")+ylab("Sum-of-Pairs score between pairs of sequences - PSICoffee vs mTMalign_NAT (%)") + scale_color_manual(values=color_codes) + ggtitle("Sum-of-Pairs score between pairs of sequences with TCS > 0")
#ggsave("ref_mtmalign_pair_sp_score_per_family_dmpfold_vs_psicoffee.png",dpi = "retina")
o=ggplot(ref_mtmalign_pair_sp,aes(x=`mTMalign_AF2`,y=PSIcoffee,color=Family)) + theme_light() + theme(legend.position="top") +geom_point()+xlim(0,100)+ylim(0,100)+geom_abline(intercept =0 , slope = 1,lwd=0.2)+xlab("Sum-of-Pairs score between pairs of sequences - mTMalign_AF2 vs mTMalign_NAT (%)")+ylab("Sum-of-Pairs score between pairs of sequences - PSICoffee vs mTMalign_NAT (%)") + scale_color_manual(values=color_codes) + ggtitle("Sum-of-Pairs score between pairs of sequences with TCS > 0")
#ggsave("ref_mtmalign_pair_sp_score_per_family_alphafold_vs_psicoffee.png",dpi = "retina")
o=ggplot(ref_mtmalign_pair_sp,aes(x=`mTMalign_DMPFOLD`,y=Ginsi,color=Family)) + theme_light() + theme(legend.position="top") +geom_point()+xlim(0,100)+ylim(0,100)+geom_abline(intercept =0 , slope = 1,lwd=0.2)+xlab("Sum-of-Pairs score between pairs of sequences - mTMalign_DMPFOLD vs mTMalign_NAT (%)")+ylab("Sum-of-Pairs score between pairs of sequences - Ginsi vs mTMalign_NAT (%)") + scale_color_manual(values=color_codes) + ggtitle("Sum-of-Pairs score between pairs of sequences with TCS > 0")
#ggsave("ref_mtmalign_pair_sp_score_per_family_dmpfold_vs_ginsi.png",dpi = "retina")
o=ggplot(ref_mtmalign_pair_sp,aes(x=`mTMalign_AF2`,y=Ginsi,color=Family)) + theme_light() + theme(legend.position="top") +geom_point()+xlim(0,100)+ylim(0,100)+geom_abline(intercept =0 , slope = 1,lwd=0.2)+xlab("Sum-of-Pairs score between pairs of sequences - mTMalign_AF2 vs mTMalign_NAT (%)")+ylab("Sum-of-Pairs score between pairs of sequences - Ginsi vs mTMalign_NAT (%)") + scale_color_manual(values=color_codes) + ggtitle("Sum-of-Pairs score between pairs of sequences with TCS > 0")
#ggsave("ref_mtmalign_pair_sp_score_per_family_alphafold_vs_ginsi.png",dpi = "retina")


#
# Input the TCS scores
#

tcs_score = read.table("./selected_comparisons_tcs.txt",header=F,stringsAsFactors = F)
colnames(tcs_score) = c("Family","Ginsi","TCoffee","PSIcoffee","3DCoffee_NAT","3DCoffee_TMalign_NAT","mTMalign_NAT","3DCoffee_DMPFOLD","3DCoffee_TMalign_DMPFOLD","mTMalign_DMPFOLD","3DCoffee_AF2","3DCoffee_TMalign_AF2","mTMalign_AF2") #"Deepblast"

tcs_score_melted = melt(tcs_score,measure.vars = 2:13)
m=ggplot(tcs_score_melted,aes(x=variable,y=value,fill=variable)) + geom_boxplot() + stat_summary(fun.y=mean, geom="point", shape=18, size=4)+ylab("TCS score")+xlab("Modes")+theme(axis.text.x=element_text(angle = 45, hjust = 1))
#ggsave("tcs_score_per_mode.png",dpi = "retina")
x=ggplot(tcs_score_melted,aes(x=Family,y=value,fill=variable))+geom_bar(position = "dodge",stat = "identity",width = 0.9,size=0.05,color="navy")+theme(axis.text.x=element_text(angle = 45, hjust = 1))+ylab("TCS score")+labs(fill='Comparisons')
#ggsave("tcs_score_per_family.png",dpi = "retina")

tcs_score_per_seq = read.table("./selected_comparisons_tcs_avg.txt",header=F,stringsAsFactors = F)
tcs_score_per_seq = tcs_score_per_seq[!duplicated(as.list(tcs_score_per_seq))]
colnames(tcs_score_per_seq) = c("Sequence","Family","Ginsi","TCoffee","PSIcoffee","3DCoffee_NAT","3DCoffee_TMalign_NAT","mTMalign_NAT","3DCoffee_DMPFOLD","3DCoffee_TMalign_DMPFOLD","mTMalign_DMPFOLD","3DCoffee_AF2","3DCoffee_TMalign_AF2","mTMalign_AF2") #,"Deepblast"
tcs_score_per_seq[,-c(1,2)]=tcs_score_per_seq[,-c(1,2)]*10


#
# PRED vs NAT structure comparisons
#

dmpfold = read.table("./dmpfold_vs_ref_pdb_comparison_selected.tsv",header = F,stringsAsFactors = F)
alphafold = read.table("./alphafold_vs_ref_pdb_comparison_selected.tsv",header = F,stringsAsFactors = F)
#dmpfold_4 = read.table("./dmpfold_4_vs_ref_pdb_comparison_selected.tsv",header = F,stringsAsFactors = F)
#dmpfold_5 = read.table("./dmpfold_5_vs_ref_pdb_comparison_selected.tsv",header = F,stringsAsFactors = F)
colnames(dmpfold)=c("Sequence","RMSD","TMscore","GDT_TS","Family")
colnames(alphafold)=c("Sequence","RMSD","TMscore","GDT_TS","Family")
#colnames(dmpfold_4)=c("Sequence","RMSD","TMscore","GDT_TS","Family")
#colnames(dmpfold_5)=c("Sequence","RMSD","TMscore","GDT_TS","Family")

complete_merged=merge(dmpfold,alphafold,by=c("Sequence"))
complete_merged = complete_merged[!duplicated(as.list(complete_merged))]
colnames(complete_merged) = c("Sequence","RMSD_DMPFOLD","TMscore_DMPFOLD","GDT_TS_DMPFOLD","Family","RMSD_AF2","TMscore_AF2","GDT_TS_AF2")
p = ggplot(complete_merged,aes(x=RMSD_AF2,y=GDT_TS_AF2,color=Family)) + geom_point(shape=1) + stat_cor(inherit.aes = FALSE,data=complete_merged,aes(x=RMSD_AF2,y=GDT_TS_AF2),method = "pearson",hjust = -1) + geom_smooth(method='lm',inherit.aes = FALSE,data=complete_merged,aes(x=RMSD_AF2,y=GDT_TS_AF2)) + theme_light()
#ggsave("RMSD_vs_GDT_TS_AF2.png",dpi="retina")
p = ggplot(complete_merged,aes(x=RMSD_DMPFOLD,y=GDT_TS_DMPFOLD,color=Family)) + geom_point(shape=1) + stat_cor(inherit.aes = FALSE,data=complete_merged,aes(x=RMSD_DMPFOLD,y=GDT_TS_DMPFOLD),method = "pearson",hjust = -1) + geom_smooth(method='lm',inherit.aes = FALSE,data=complete_merged,aes(x=RMSD_DMPFOLD,y=GDT_TS_DMPFOLD)) + theme_light()
#ggsave("RMSD_vs_GDT_TS_AF2.png",dpi="retina")


rmsd_complete=complete_merged[,c(2,6)]
colnames(rmsd_complete)=c("dmpfold","alphafold")
#png(filename = "rmsd_heatmap.png",width = 1200,height = 900,res = 200)
#superheat(as.matrix(rmsd_complete),col.dendrogram = F,row.dendrogram = F,pretty.order.rows = T,pretty.order.cols = F,heat.na.col = "white",left.label.text.size = 2,bottom.label.text.size = 3,heat.col.scheme = "viridis",legend.text.size = 7,title = "RMSD between dmpfold and experimental structures")
#dev.off()

tmscore_complete=complete_merged[,c(3,7)]
colnames(tmscore_complete)=c("dmpfold","alphafold")

tmscore_complete_with_seqs = complete_merged[,c(1,3,7,length(colnames(complete_merged)))]
colnames(tmscore_complete_with_seqs) = c("Sequence","dmpfold","alphafold","Family")
tmscore_complete_with_seqs = tmscore_complete_with_seqs[order(tmscore_complete_with_seqs$Family),]
#png(filename = "tmscore_complex_heatmap.png",width = 1200,height = 900,res = 200)
#Heatmap(tmscore_complete_with_seqs[,2:4],name="TMscore",cluster_rows=F,cluster_columns=F,show_row_names=F)+Heatmap(tmscore_complete_with_seqs$Family, name = "Family", width = unit(5, "mm"),col = colorRampPalette(brewer.pal(n=12,name="Paired"))(length(levels(as.factor(tmscore_complete_with_seqs$Family)))))
#dev.off()

gdt_ts_complete=complete_merged[,c(4,8)]
colnames(gdt_ts_complete)=c("dmpfold","alphafold")

gdt_ts_complete_with_seqs = complete_merged[,c("Sequence","GDT_TS_DMPFOLD","GDT_TS_AF2","Family")]
colnames(gdt_ts_complete_with_seqs) = c("Sequence","dmpfold","alphafold","Family")
gdt_ts_complete_with_seqs = gdt_ts_complete_with_seqs[order(gdt_ts_complete_with_seqs$Family),]
#png(filename = "gdt_ts_complex_heatmap.png",width = 1200,height = 900,res = 200)
#Heatmap(gdt_ts_complete_with_seqs[,2:4],name="GDT_TS",cluster_rows=F,cluster_columns=F,show_row_names=F) + Heatmap(gdt_ts_complete_with_seqs$Family, name = "Family", width = unit(5, "mm"),col = colorRampPalette(brewer.pal(n=12,name="Paired"))(length(levels(as.factor(gdt_ts_complete_with_seqs$Family)))))
#dev.off()

gdt_ts_complete_with_seqs_rescale = gdt_ts_complete_with_seqs
gdt_ts_complete_with_seqs_rescale$alphafold = gdt_ts_complete_with_seqs_rescale$alphafold * 100
p = ggplot(gdt_ts_complete_with_seqs_rescale, aes(x=alphafold,fill=Family)) + geom_histogram(color="black",position="stack") + theme_light() + xlab("AF2 GDT_TS")
#ggsave("gdt_ts_alphafold_histogram.png",dpi="retina")

gdt_ts_complete_with_seqs_rescale$dmpfold = gdt_ts_complete_with_seqs_rescale$dmpfold * 100
p = ggplot(gdt_ts_complete_with_seqs_rescale, aes(x=dmpfold)) + geom_histogram(color="black",fill="darkblue") + theme_light() + xlab("DMPFOLD GDT_TS")
#ggsave("gdt_ts_dmpfold_histogram.png",dpi="retina")

alphafold_plddts = read.table("./alphafold_plddts_selected.tsv",header = F, stringsAsFactors = F)
colnames(alphafold_plddts)=c("Sequence","pLDDT","Family")
agg_plddts = aggregate(x=alphafold_plddts[,-c(1,3)],by=list(alphafold_plddts$Family),FUN=mean)
colnames(agg_plddts) = c("Family","pLDDT TCS > 0")

agg_gdts = aggregate(x=gdt_ts_complete_with_seqs[,c(3)],by=list(gdt_ts_complete_with_seqs$Family),FUN=mean)
colnames(agg_gdts) = c("Family","GDT_TS TCS > 0")

for (cutoff in seq(50, 95, by=5)) {
  tmp_seqs = read.table(paste0("all_fams.list_tcs_above_",cutoff,"_ref_3dcoffee_alphafold"),header=F,stringsAsFactors=F)
  colnames(tmp_seqs) = c("Sequence")
  tmp_kept_seqs = merge(alphafold_plddts,tmp_seqs,by="Sequence")
  agg_kept_seqs = aggregate(x=tmp_kept_seqs[,-c(1,3)],by=list(tmp_kept_seqs$Family),FUN=mean)
  colnames(agg_kept_seqs) = c("Family",paste0("pLDDT TCS > ",cutoff*10))
  agg_plddts = merge(agg_plddts,agg_kept_seqs,all=T,by="Family")

  tmp_kept_seqs_gdt = merge(gdt_ts_complete_with_seqs,tmp_seqs,by="Sequence")
  agg_kept_seqs_gdt = aggregate(x=tmp_kept_seqs_gdt[,c(3)],by=list(tmp_kept_seqs_gdt$Family),FUN=mean)
  colnames(agg_kept_seqs_gdt) = c("Family",paste0("GDT_TS TCS > ",cutoff*10))
  agg_gdts = merge(agg_gdts,agg_kept_seqs_gdt,all=T,by="Family")
}
agg_plddts = agg_plddts[,-c(1)]
agg_gdts = agg_gdts[,-c(1)]*100

#
# Deltas analysis
#


#################
# Ref mtmalign
#################

#
# AF2
#

  #
  # mTMalign
  #

#ref_mtmalign_avg_sp_mTMalign_AF2 = ref_mtmalign_avg_sp[,c("Sequence","Family","mTMalign_AF2")]
ref_mtmalign_avg_sp_mTMalign_AF2 = ref_mtmalign_avg_sp_same_state[,c("Sequence","Family","mTMalign_AF2")]
colnames(ref_mtmalign_avg_sp_mTMalign_AF2)[3] = "SoP"
tcs_score_per_seq_mTMalign_AF2 = tcs_score_per_seq[,c("Sequence","Family","mTMalign_AF2")]
colnames(tcs_score_per_seq_mTMalign_AF2)[3] = "TCS"
nirmsd_ref_avg_mTMalign_AF2_ref_NAT = nirmsd_ref_avg[c("Sequence","Family","mTMalign_AF2_REF_NAT")]
#nirmsd_ref_avg_mTMalign_AF2_ref_AF2 = nirmsd_ref_avg[c("Sequence","Family","mTMalign_AF2_REF_AF2")]
colnames(nirmsd_ref_avg_mTMalign_AF2_ref_NAT)[3] = "niRMSD"
delta_analysis_AF2(tcs_score_per_seq_mTMalign_AF2,gdt_ts_complete_with_seqs[,-c(2)],alphafold_plddts,ref_mtmalign_avg_sp_mTMalign_AF2,nirmsd_ref_avg_mTMalign_AF2_ref_NAT,"mTMalign","all")

GDT_AF2_lower_quart = quantile(gdt_ts_complete_with_seqs$alphafold)[2]
#delta_analysis_AF2(tcs_score_per_seq_mTMalign_AF2,gdt_ts_complete_with_seqs[which(gdt_ts_complete_with_seqs$alphafold <= GDT_AF2_lower_quart),-c(2)],alphafold_plddts,ref_mtmalign_avg_sp_mTMalign_AF2,nirmsd_ref_avg_mTMalign_AF2_ref_NAT,"mTMalign","lower_quartile")

#ref_mtmalign_pair_sp_mTMalign_AF2 = ref_mtmalign_pair_sp[,c("Sequence_1","Sequence_2","Family","mTMalign_AF2")]
ref_mtmalign_pair_sp_mTMalign_AF2 = ref_mtmalign_pair_sp_same_state[,c("Sequence_1","Sequence_2","Family","mTMalign_AF2")]
colnames(ref_mtmalign_pair_sp_mTMalign_AF2)[4] = "SoP"
pairwise_plots(gdt_ts_complete_with_seqs[,-c(2)],ref_mtmalign_pair_sp_mTMalign_AF2,tcs_score_per_seq_mTMalign_AF2,alphafold_plddts,"mTMalign","all")

  #
  # 3DCoffee
  #

#ref_3dcoffee_avg_sp_3DCoffee_AF2 = ref_3dcoffee_avg_sp[,c("Sequence","Family","3DCoffee_AF2")]
ref_3dcoffee_avg_sp_3DCoffee_AF2 = ref_3dcoffee_avg_sp_same_state[,c("Sequence","Family","3DCoffee_AF2")]
colnames(ref_3dcoffee_avg_sp_3DCoffee_AF2)[3] = "SoP"
tcs_score_per_seq_3DCoffee_AF2 = tcs_score_per_seq[,c("Sequence","Family","3DCoffee_AF2")]
colnames(tcs_score_per_seq_3DCoffee_AF2)[3] = "TCS"
nirmsd_ref_avg_3DCoffee_AF2_ref_NAT = nirmsd_ref_avg[c("Sequence","Family","3DCoffee_AF2_REF_NAT")]
#nirmsd_ref_avg_3DCoffee_AF2_ref_AF2 = nirmsd_ref_avg[c("Sequence","Family","3DCoffee_AF2_REF_AF2")]
colnames(nirmsd_ref_avg_3DCoffee_AF2_ref_NAT)[3] = "niRMSD"
ref_3dcoffee_AF2_avg_df = delta_analysis_AF2(tcs_score_per_seq_3DCoffee_AF2,gdt_ts_complete_with_seqs[,-c(2)],alphafold_plddts,ref_3dcoffee_avg_sp_3DCoffee_AF2,nirmsd_ref_avg_3DCoffee_AF2_ref_NAT,"3DCoffee","all")

ref_3dcoffee_AF2_avg_df_sop_vs_tcs = table(ref_3dcoffee_AF2_avg_df$delta_sop >= 0, ref_3dcoffee_AF2_avg_df$delta_tcs >= 0)
ref_3dcoffee_AF2_avg_df_sop_vs_tcs = 100 * (ref_3dcoffee_AF2_avg_df_sop_vs_tcs[1] + ref_3dcoffee_AF2_avg_df_sop_vs_tcs[4]) / sum(ref_3dcoffee_AF2_avg_df_sop_vs_tcs)
ref_3dcoffee_AF2_avg_df_sop_vs_plddt = table(ref_3dcoffee_AF2_avg_df$delta_sop >= 0, ref_3dcoffee_AF2_avg_df$delta_plddt >= 0)
ref_3dcoffee_AF2_avg_df_sop_vs_plddt = 100 * (ref_3dcoffee_AF2_avg_df_sop_vs_plddt[1] + ref_3dcoffee_AF2_avg_df_sop_vs_plddt[4]) / sum(ref_3dcoffee_AF2_avg_df_sop_vs_plddt)
ref_3dcoffee_AF2_avg_df_sop_vs_gdt = table(ref_3dcoffee_AF2_avg_df$delta_sop >= 0, ref_3dcoffee_AF2_avg_df$delta_gdt_ts >= 0)
ref_3dcoffee_AF2_avg_df_sop_vs_gdt = 100 * (ref_3dcoffee_AF2_avg_df_sop_vs_gdt[1] + ref_3dcoffee_AF2_avg_df_sop_vs_gdt[4]) / sum(ref_3dcoffee_AF2_avg_df_sop_vs_gdt)
ref_3dcoffee_AF2_avg_df_delta_perc = c(ref_3dcoffee_AF2_avg_df_sop_vs_tcs,ref_3dcoffee_AF2_avg_df_sop_vs_plddt,ref_3dcoffee_AF2_avg_df_sop_vs_gdt)


ref_3dcoffee_avg_sp_ginsi = ref_3dcoffee_avg_sp_same_state[,c("Sequence","Family","Ginsi")]
colnames(ref_3dcoffee_avg_sp_ginsi)[3] = "SoP"
tcs_score_per_seq_ginsi = tcs_score_per_seq[,c("Sequence","Family","Ginsi")]
colnames(tcs_score_per_seq_ginsi)[3] = "TCS"
nirmsd_ref_avg_ginsi = nirmsd_ref_avg[c("Sequence","Family","Ginsi")]
colnames(nirmsd_ref_avg_ginsi)[3] = "niRMSD"
ref_ginsi_avg_df = delta_analysis_AF2(tcs_score_per_seq_ginsi,gdt_ts_complete_with_seqs[,-c(2)],alphafold_plddts,ref_3dcoffee_avg_sp_ginsi,nirmsd_ref_avg_ginsi,"3DCoffee","ginsi")

ref_ginsi_avg_df_sop_vs_tcs = table(ref_ginsi_avg_df$delta_sop >= 0, ref_ginsi_avg_df$delta_tcs >= 0)
ref_ginsi_avg_df_sop_vs_tcs = 100 * (ref_ginsi_avg_df_sop_vs_tcs[1] + ref_ginsi_avg_df_sop_vs_tcs[4]) / sum(ref_ginsi_avg_df_sop_vs_tcs)
ref_ginsi_avg_df_delta_perc = c(ref_ginsi_avg_df_sop_vs_tcs,NA,NA)


ref_3dcoffee_avg_sp_psicoffee = ref_3dcoffee_avg_sp_same_state[,c("Sequence","Family","PSIcoffee")]
colnames(ref_3dcoffee_avg_sp_psicoffee)[3] = "SoP"
tcs_score_per_seq_psicoffee = tcs_score_per_seq[,c("Sequence","Family","PSIcoffee")]
colnames(tcs_score_per_seq_psicoffee)[3] = "TCS"
nirmsd_ref_avg_psicoffee = nirmsd_ref_avg[c("Sequence","Family","PSIcoffee")]
colnames(nirmsd_ref_avg_psicoffee)[3] = "niRMSD"
ref_psicoffee_avg_df = delta_analysis_AF2(tcs_score_per_seq_psicoffee,gdt_ts_complete_with_seqs[,-c(2)],alphafold_plddts,ref_3dcoffee_avg_sp_psicoffee,nirmsd_ref_avg_psicoffee,"3DCoffee","psicoffee")

ref_psicoffee_avg_df_sop_vs_tcs = table(ref_psicoffee_avg_df$delta_sop >= 0, ref_psicoffee_avg_df$delta_tcs >= 0)
ref_psicoffee_avg_df_sop_vs_tcs = 100 * (ref_psicoffee_avg_df_sop_vs_tcs[1] + ref_psicoffee_avg_df_sop_vs_tcs[4]) / sum(ref_psicoffee_avg_df_sop_vs_tcs)
ref_psicoffee_avg_df_delta_perc = c(ref_psicoffee_avg_df_sop_vs_tcs,NA,NA)



ref_3dcoffee_avg_sp_3DCoffee_NAT = ref_3dcoffee_avg_sp_same_state[,c("Sequence","Family","3DCoffee_NAT")]
colnames(ref_3dcoffee_avg_sp_3DCoffee_NAT)[3] = "SoP"
tcs_score_per_seq_3DCoffee_NAT = tcs_score_per_seq[,c("Sequence","Family","3DCoffee_NAT")]
colnames(tcs_score_per_seq_3DCoffee_NAT)[3] = "TCS"
nirmsd_ref_avg_3DCoffee_NAT = nirmsd_ref_avg[c("Sequence","Family","3DCoffee_NAT")]
colnames(nirmsd_ref_avg_3DCoffee_NAT)[3] = "niRMSD"
ref_3dcoffee_NAT_avg_df = delta_analysis_AF2(tcs_score_per_seq_3DCoffee_NAT,gdt_ts_complete_with_seqs[,-c(2)],alphafold_plddts,ref_3dcoffee_avg_sp_3DCoffee_NAT,nirmsd_ref_avg_3DCoffee_NAT,"3DCoffee","3DCoffee_NAT")

ref_3dcoffee_NAT_avg_df_sop_vs_tcs = table(ref_3dcoffee_NAT_avg_df$delta_sop >= 0, ref_3dcoffee_NAT_avg_df$delta_tcs >= 0)
ref_3dcoffee_NAT_avg_df_sop_vs_tcs = 100 * (ref_3dcoffee_NAT_avg_df_sop_vs_tcs[2]) / sum(ref_3dcoffee_NAT_avg_df_sop_vs_tcs)
ref_3dcoffee_NAT_avg_df_delta_perc = c(ref_3dcoffee_NAT_avg_df_sop_vs_tcs,NA,NA)





#delta_analysis_AF2(tcs_score_per_seq_3DCoffee_AF2,gdt_ts_complete_with_seqs[which(gdt_ts_complete_with_seqs$alphafold <= GDT_AF2_lower_quart),-c(2)],alphafold_plddts,ref_3dcoffee_avg_sp_3DCoffee_AF2,nirmsd_ref_avg_3DCoffee_AF2_ref_NAT,"3DCoffee","lower_quartile")

#ref_3dcoffee_pair_sp_3dcoffee_AF2 = ref_3dcoffee_pair_sp[,c("Sequence_1","Sequence_2","Family","3DCoffee_AF2")]
ref_3dcoffee_pair_sp_3dcoffee_AF2 = ref_3dcoffee_pair_sp_same_state[,c("Sequence_1","Sequence_2","Family","3DCoffee_AF2")]
colnames(ref_3dcoffee_pair_sp_3dcoffee_AF2)[4] = "SoP"
ref_3dcoffee_AF2_pair_df = pairwise_plots(gdt_ts_complete_with_seqs[,-c(2)],ref_3dcoffee_pair_sp_3dcoffee_AF2,tcs_score_per_seq_3DCoffee_AF2,alphafold_plddts,"3DCoffee","all")

ref_3dcoffee_AF2_pair_df_gm_gdt_below_75 = ref_3dcoffee_AF2_pair_df[which(ref_3dcoffee_AF2_pair_df$GM_GDT_TS < 75),]
ref_3dcoffee_AF2_with_pairs_gm_gdt_below_75_sop_above_75 = length(which(ref_3dcoffee_AF2_pair_df_gm_gdt_below_75$SoP > 75))/nrow(ref_3dcoffee_AF2_pair_df_gm_gdt_below_75)
ref_3dcoffee_AF2_with_pairs_gm_gdt_below_75_sop_above_75_df = as.data.frame(t(as.data.frame(table(ref_3dcoffee_AF2_pair_df_gm_gdt_below_75$SoP > 75))))

ref_3dcoffee_AF2_with_pairs_gm_gdt_below_75_sop_above_75_df = ref_3dcoffee_AF2_with_pairs_gm_gdt_below_75_sop_above_75_df[2,]
row.names(ref_3dcoffee_AF2_with_pairs_gm_gdt_below_75_sop_above_75_df) = c("AF2")
colnames(ref_3dcoffee_AF2_with_pairs_gm_gdt_below_75_sop_above_75_df) = c("SoP < 75", "SoP > 75")


#
# DMPFOLD
#

  #
  # mTMalign
  #

#ref_mtmalign_avg_sp_mTMalign_DMPFOLD = ref_mtmalign_avg_sp[,c("Sequence","Family","mTMalign_DMPFOLD")]
ref_mtmalign_avg_sp_mTMalign_DMPFOLD = ref_mtmalign_avg_sp_same_state[,c("Sequence","Family","mTMalign_DMPFOLD")]
colnames(ref_mtmalign_avg_sp_mTMalign_DMPFOLD)[3] = "SoP"
tcs_score_per_seq_mTMalign_DMPFOLD = tcs_score_per_seq[,c("Sequence","Family","mTMalign_DMPFOLD")]
colnames(tcs_score_per_seq_mTMalign_DMPFOLD)[3] = "TCS"
nirmsd_ref_avg_mTMalign_DMPFOLD_ref_NAT = nirmsd_ref_avg[c("Sequence","Family","mTMalign_DMPFOLD_REF_NAT")]
#nirmsd_ref_avg_mTMalign_DMPFOLD_ref_DMPFOLD = nirmsd_ref_avg[c("Sequence","Family","mTMalign_DMPFOLD_REF_DMPFOLD")]
colnames(nirmsd_ref_avg_mTMalign_DMPFOLD_ref_NAT)[3] = "niRMSD"
#delta_analysis_DMPFOLD(tcs_score_per_seq_mTMalign_DMPFOLD,gdt_ts_complete_with_seqs[,-c(3)],ref_mtmalign_avg_sp_mTMalign_DMPFOLD,nirmsd_ref_avg_mTMalign_DMPFOLD_ref_NAT,"mTMalign","all")

GDT_DMPFOLD_lower_quart = quantile(gdt_ts_complete_with_seqs$dmpfold)[2]
#delta_analysis_DMPFOLD(tcs_score_per_seq_mTMalign_DMPFOLD,gdt_ts_complete_with_seqs[which(gdt_ts_complete_with_seqs$dmpfold <= GDT_DMPFOLD_lower_quart),-c(3)],ref_mtmalign_avg_sp_mTMalign_DMPFOLD,nirmsd_ref_avg_mTMalign_DMPFOLD_ref_NAT,"mTMalign","lower_quartile")


#ref_mtmalign_pair_sp_mTMalign_DMPFOLD = ref_mtmalign_pair_sp[,c("Sequence_1","Sequence_2","Family","mTMalign_DMPFOLD")]
ref_mtmalign_pair_sp_mTMalign_DMPFOLD = ref_mtmalign_pair_sp_same_state[,c("Sequence_1","Sequence_2","Family","mTMalign_DMPFOLD")]
colnames(ref_mtmalign_pair_sp_mTMalign_DMPFOLD)[4] = "SoP"
pairwise_plots_dmpfold(gdt_ts_complete_with_seqs[,-c(3)],ref_mtmalign_pair_sp_mTMalign_DMPFOLD,tcs_score_per_seq_mTMalign_DMPFOLD,"mTMalign","all")

  #
  # 3DCoffee
  #

#ref_3dcoffee_avg_sp_3DCoffee_DMPFOLD = ref_3dcoffee_avg_sp[,c("Sequence","Family","3DCoffee_DMPFOLD")]
ref_3dcoffee_avg_sp_3DCoffee_DMPFOLD = ref_3dcoffee_avg_sp_same_state[,c("Sequence","Family","3DCoffee_DMPFOLD")]
colnames(ref_3dcoffee_avg_sp_3DCoffee_DMPFOLD)[3] = "SoP"
tcs_score_per_seq_3DCoffee_DMPFOLD = tcs_score_per_seq[,c("Sequence","Family","3DCoffee_DMPFOLD")]
colnames(tcs_score_per_seq_3DCoffee_DMPFOLD)[3] = "TCS"
nirmsd_ref_avg_3DCoffee_DMPFOLD_ref_NAT = nirmsd_ref_avg[c("Sequence","Family","3DCoffee_DMPFOLD_REF_NAT")]
#nirmsd_ref_avg_3DCoffee_DMPFOLD_ref_DMPFOLD = nirmsd_ref_avg[c("Sequence","Family","3DCoffee_DMPFOLD_REF_DMPFOLD")]
colnames(nirmsd_ref_avg_3DCoffee_DMPFOLD_ref_NAT)[3] = "niRMSD"
ref_3dcoffee_DMPFOLD_avg_df = delta_analysis_DMPFOLD(tcs_score_per_seq_3DCoffee_DMPFOLD,gdt_ts_complete_with_seqs[,-c(3)],ref_3dcoffee_avg_sp_3DCoffee_DMPFOLD,nirmsd_ref_avg_3DCoffee_DMPFOLD_ref_NAT,"3DCoffee","all")
#delta_analysis_DMPFOLD(tcs_score_per_seq_3DCoffee_DMPFOLD,gdt_ts_complete_with_seqs[which(gdt_ts_complete_with_seqs$dmpfold <= GDT_DMPFOLD_lower_quart),-c(3)],ref_3dcoffee_avg_sp_3DCoffee_DMPFOLD,nirmsd_ref_avg_3DCoffee_DMPFOLD_ref_NAT,"3DCoffee","lower_quartile")

ref_3dcoffee_DMPFOLD_avg_df_sop_vs_tcs = table(ref_3dcoffee_DMPFOLD_avg_df$delta_sop >= 0, ref_3dcoffee_DMPFOLD_avg_df$delta_tcs >= 0)
ref_3dcoffee_DMPFOLD_avg_df_sop_vs_tcs = 100 * (ref_3dcoffee_DMPFOLD_avg_df_sop_vs_tcs[1] + ref_3dcoffee_DMPFOLD_avg_df_sop_vs_tcs[4]) / sum(ref_3dcoffee_DMPFOLD_avg_df_sop_vs_tcs)
ref_3dcoffee_DMPFOLD_avg_df_sop_vs_gdt = table(ref_3dcoffee_DMPFOLD_avg_df$delta_sop >= 0, ref_3dcoffee_DMPFOLD_avg_df$delta_gdt_ts >= 0)
ref_3dcoffee_DMPFOLD_avg_df_sop_vs_gdt = 100 * (ref_3dcoffee_DMPFOLD_avg_df_sop_vs_gdt[1] + ref_3dcoffee_DMPFOLD_avg_df_sop_vs_gdt[4]) / sum(ref_3dcoffee_DMPFOLD_avg_df_sop_vs_gdt)
ref_3dcoffee_DMPFOLD_avg_df_delta_perc = c(ref_3dcoffee_DMPFOLD_avg_df_sop_vs_tcs,NA,ref_3dcoffee_DMPFOLD_avg_df_sop_vs_gdt)


#ref_3dcoffee_pair_sp_3dcoffee_DMPFOLD = ref_3dcoffee_pair_sp[,c("Sequence_1","Sequence_2","Family","3DCoffee_DMPFOLD")]
ref_3dcoffee_pair_sp_3dcoffee_DMPFOLD = ref_3dcoffee_pair_sp_same_state[,c("Sequence_1","Sequence_2","Family","3DCoffee_DMPFOLD")]
colnames(ref_3dcoffee_pair_sp_3dcoffee_DMPFOLD)[4] = "SoP"
ref_3dcoffee_DMPFOLD_pair_df = pairwise_plots_dmpfold(gdt_ts_complete_with_seqs[,-c(3)],ref_3dcoffee_pair_sp_3dcoffee_DMPFOLD,tcs_score_per_seq_3DCoffee_DMPFOLD,"3DCoffee","all")

ref_3dcoffee_DMPFOLD_pair_df_gm_gdt_below_75 = ref_3dcoffee_DMPFOLD_pair_df[which(ref_3dcoffee_DMPFOLD_pair_df$GM_GDT_TS < 75),]
ref_3dcoffee_DMPFOLD_with_pairs_gm_gdt_below_75_sop_above_75 = length(which(ref_3dcoffee_DMPFOLD_pair_df_gm_gdt_below_75$SoP > 75))/nrow(ref_3dcoffee_DMPFOLD_pair_df_gm_gdt_below_75)
ref_3dcoffee_DMPFOLD_with_pairs_gm_gdt_below_75_sop_above_75_df = as.data.frame(t(as.data.frame(table(ref_3dcoffee_DMPFOLD_pair_df_gm_gdt_below_75$SoP > 75))))

ref_3dcoffee_DMPFOLD_with_pairs_gm_gdt_below_75_sop_above_75_df = ref_3dcoffee_DMPFOLD_with_pairs_gm_gdt_below_75_sop_above_75_df[2,]
row.names(ref_3dcoffee_DMPFOLD_with_pairs_gm_gdt_below_75_sop_above_75_df) = c("DMPFOLD")
colnames(ref_3dcoffee_DMPFOLD_with_pairs_gm_gdt_below_75_sop_above_75_df) = c("SoP < 75", "SoP > 75")


#
# CPU time vs seq length
#

seq_length = read.table("./seq_lengths",header = F,stringsAsFactors = F)
colnames(seq_length) = c("Sequence","Family","Length")
seq_realtime = read.table("./seq_realtime.txt",header = F,stringsAsFactors = F)
colnames(seq_realtime) = c("Proccess","Sequence","Family","Time")
seq_realtime_only_3 = seq_realtime[which(seq_realtime$Proccess == 'run_seq2maps' | seq_realtime$Proccess == 'run_dmpfold'),]
seq_realtime = seq_realtime[order(seq_realtime[,2]),]
seq_realtime_only_3 = seq_realtime_only_3[order(seq_realtime_only_3[,2]),]
sum_seq_realtime_only_3 = aggregate(seq_realtime_only_3$Time, by=list(Sequence = seq_realtime_only_3$Sequence),FUN=sum)
colnames(sum_seq_realtime_only_3)[2] = "Time"
sum_seq_realtime_only_3$Time = round(sum_seq_realtime_only_3$Time / (1000*60),digits=2)
length_and_realtime = merge(seq_length,sum_seq_realtime_only_3,by="Sequence")
length_and_realtime$Family = factor(length_and_realtime$Family,levels=row.names(nirmsd_ref)[-nrow(nirmsd_ref)])

p = ggscatter(length_and_realtime,x="Length",y="Time",add = "reg.line",conf.int = TRUE,color="Family",palette=color_codes, add.params = list(color = "blue",fill = "lightgray")) + stat_cor(method = "pearson") + xlab("Number of residues") + ylab("CPU time (minutes)")
ggsave("length_vs_cpu_time_dmpfold.png",dpi="retina")

seq_realtime_alphafold = read.table("./seq_realtime_alphafold.txt",header = F,stringsAsFactors = F)
colnames(seq_realtime_alphafold) = c("Proccess","Sequence","Family","Time")
seq_realtime_alphafold$Time = round(seq_realtime_alphafold$Time / (1000*60),digits=2)
length_and_realtime_alphafold = merge(seq_length,seq_realtime_alphafold,by="Sequence")
colnames(length_and_realtime_alphafold)[2] = "Family"
length_and_realtime_alphafold$Family = factor(length_and_realtime_alphafold$Family,levels=row.names(nirmsd_ref)[-nrow(nirmsd_ref)])

p = ggscatter(length_and_realtime_alphafold,x="Length",y="Time",add = "reg.line",conf.int = TRUE,color="Family",palette=color_codes, add.params = list(color = "blue",fill = "lightgray")) + stat_cor(method = "pearson") + xlab("Number of residues") + ylab("CPU time (minutes)")
ggsave("length_vs_cpu_time_alphafold.png",dpi="retina")

#
# GDT_TS vs AVG_SoP and TCS vs AVG_SoP
#

#
# Ref mtmalign
#
if (FALSE) {
  ref_mtmalign_avg_sp_and_gdt_ts_complete_with_seqs = merge(gdt_ts_complete_with_seqs,ref_mtmalign_avg_sp,by="Sequence")
colnames(ref_mtmalign_avg_sp_and_gdt_ts_complete_with_seqs)[5] = "Family"
ref_mtmalign_avg_sp_and_gdt_ts_complete_with_seqs$Family = factor(ref_mtmalign_avg_sp_and_gdt_ts_complete_with_seqs$Family,levels=row.names(nirmsd_ref)[-nrow(nirmsd_ref)])
p = ggscatter(ref_mtmalign_avg_sp_and_gdt_ts_complete_with_seqs,x="mTMalign_DMPFOLD",y="dmpfold",add = "reg.line",conf.int = TRUE,color="Family",palette=color_codes, add.params = list(color = "blue",fill = "lightgray")) + stat_cor(method = "pearson") + xlab("Sum-of-Pairs score per sequence - mTMalign_DMPFOLD vs mTMalign_NAT (%)") + ylab("GDT_TS DMPFOLD vs NAT structures") + font("xlab", size = 12) + font("ylab", size = 12)
ggsave("ref_mtmalign_avg_sp_score_per_seq_per_family_vs_gdt_score_dmpfold.png",dpi="retina")
p = ggscatter(ref_mtmalign_avg_sp_and_gdt_ts_complete_with_seqs,x="mTMalign_AF2",y="alphafold",add = "reg.line",conf.int = TRUE,color="Family",palette=color_codes, add.params = list(color = "blue",fill = "lightgray")) + stat_cor(method = "pearson") + xlab("Sum-of-Pairs score per sequence - mTMalign_AF2 vs mTMalign_NAT (%)") + ylab("GDT_TS AF2 vs NAT structures") + font("xlab", size = 12) + font("ylab", size = 12)
ggsave("ref_mtmalign_avg_sp_score_per_seq_per_family_vs_gdt_score_alphafold.png",dpi="retina")

ref_mtmalign_avg_sp_with_tcs_score_per_seq = merge(ref_mtmalign_avg_sp,tcs_score_per_seq,by="Sequence")
colnames(ref_mtmalign_avg_sp_with_tcs_score_per_seq)[2] = "Family"
ref_mtmalign_avg_sp_with_tcs_score_per_seq$mTMalign_AF2.y = ref_mtmalign_avg_sp_with_tcs_score_per_seq$mTMalign_AF2.y/10
ref_mtmalign_avg_sp_with_tcs_score_per_seq$Family = factor(ref_mtmalign_avg_sp_with_tcs_score_per_seq$Family,levels=row.names(nirmsd_ref)[-nrow(nirmsd_ref)])
p = ggscatter(ref_mtmalign_avg_sp_with_tcs_score_per_seq,x="mTMalign_DMPFOLD.x",y="mTMalign_DMPFOLD.y",add = "reg.line",conf.int = TRUE,color="Family",palette=color_codes, add.params = list(color = "blue",fill = "lightgray")) + stat_cor(method = "pearson") + xlab("Sum-of-Pairs score per sequence - mTMalign_DMPFOLD vs mTMalign_NAT (%)") + ylab("TCS score per sequence - mTMalign_DMPFOLD") + font("xlab", size = 12) + font("ylab", size = 12)
ggsave("ref_mtmalign_avg_sp_score_per_seq_per_family_vs_tcs_score_dmpfold.png",dpi="retina")
#p = ggscatter(ref_mtmalign_avg_sp_with_tcs_score_per_seq,x="mTMalign_AF2.x",y="mTMalign_AF2.y",add = "reg.line",conf.int = TRUE,color="Family",palette=color_codes, add.params = list(color = "blue",fill = "lightgray")) + stat_cor(method = "pearson") + xlab("Sum-of-Pairs score per sequence - mTMalign_AF2 vs mTMalign_NAT (%)") + ylab("TCS score per sequence - mTMalign_AF2") + font("xlab", size = 8) + font("ylab", size = 8) + xlim(c(0,100)) + ylim(c(0,100))
p = ggplot(ref_mtmalign_avg_sp_with_tcs_score_per_seq,aes(x=mTMalign_AF2.x,y=mTMalign_AF2.y,color=Family)) + geom_point() + stat_cor(inherit.aes = FALSE,data=ref_mtmalign_avg_sp_with_tcs_score_per_seq,aes(x=mTMalign_AF2.x,y=mTMalign_AF2.y),method = "pearson") + geom_abline(slope = 1, intercept = 0,lwd=0.2) + xlab("Sum-of-Pairs score per sequence - mTMalign_AF2 vs mTMalign_NAT (%)") + ylab("TCS score per sequence - mTMalign_AF2") + theme_light() + xlim(0,100) + ylim(0,100)
ggsave("ref_mtmalign_avg_sp_score_per_seq_per_family_vs_tcs_score_alphafold.png",dpi="retina")

ref_mtmalign_avg_tcs_with_gdt_ts_complete_score_per_seq = merge(tcs_score_per_seq,gdt_ts_complete_with_seqs,by="Sequence")
colnames(ref_mtmalign_avg_tcs_with_gdt_ts_complete_score_per_seq)[2] = "Family"
p = ggscatter(ref_mtmalign_avg_tcs_with_gdt_ts_complete_score_per_seq,x="mTMalign_DMPFOLD",y="dmpfold",add = "reg.line",conf.int = TRUE,color="Family",palette=color_codes, add.params = list(color = "blue",fill = "lightgray")) + stat_cor(method = "pearson") + xlab("TCS score per sequence - mTMalign_DMPFOLD") + ylab("GDT_TS DMPFOLD vs NAT structures") + font("xlab", size = 12) + font("ylab", size = 12)
ggsave("ref_mtmalign_avg_tcs_score_per_seq_per_family_vs_gdt_ts_score_dmpfold.png",dpi="retina")
p = ggscatter(ref_mtmalign_avg_tcs_with_gdt_ts_complete_score_per_seq,x="mTMalign_AF2",y="alphafold",add = "reg.line",conf.int = TRUE,color="Family",palette=color_codes, add.params = list(color = "blue",fill = "lightgray")) + stat_cor(method = "pearson") + xlab("TCS score per sequence - mTMalign_AF2") + ylab("GDT_TS AF2 vs NAT structures") + font("xlab", size = 12) + font("ylab", size = 12)
ggsave("ref_mtmalign_avg_tcs_score_per_seq_per_family_vs_gdt_ts_score_alphafold.png",dpi="retina")
}


ref_mtmalign_dmpfold_3_tcs_cutoff_list = read.table("./all_tcs_cutoffs_all_fams_selected.list_ref_mtmalign_dmpfold",header = F,stringsAsFactors = F)
colnames(ref_mtmalign_dmpfold_3_tcs_cutoff_list) = c("Family","TCS > 0","TCS > 50","TCS > 55","TCS > 60","TCS > 65","TCS > 70","TCS > 75","TCS > 80","TCS > 85","TCS > 90","TCS > 95")
weights = sapply(ref_mtmalign_dmpfold_3_tcs_cutoff_list[,-c(1)], function(x) x/sum(x))

ref_mtmalign_dmpfold_3_tcs_cutoff_sp = read.table("./all_tcs_cutoffs_all_fams_selected_ref_ginsi_psicoffee_and_mtmalign_dmpfold_vs_ref_mtmalign.sp",header = F,stringsAsFactors = F)
ref_mtmalign_dmpfold_3_tcs_cutoff_sp = ref_mtmalign_dmpfold_3_tcs_cutoff_sp[,which(sapply(ref_mtmalign_dmpfold_3_tcs_cutoff_sp,class) != "character")]
#ref_3dcoffee_dmpfold_3_tcs_cutoff_sp = ref_3dcoffee_dmpfold_3_tcs_cutoff_sp[!duplicated(as.list(ref_3dcoffee_dmpfold_3_tcs_cutoff_sp))]

col_vector = c()
for (method in c("Gins1 TCS > ", "PSICoffee TCS > ", "mTMalign_DMPFOLD TCS > ")){
  for (cutoff in seq(50, 95, by=5)) {
      col_vector = c(col_vector,paste0(method,cutoff))
  }
}
colnames(ref_mtmalign_dmpfold_3_tcs_cutoff_sp) = col_vector

ref_mtmalign_dmpfold_3_tcs_cutoff_sp$`Gins1 TCS > 0` = ref_mtmalign_sp$Ginsi
ref_mtmalign_dmpfold_3_tcs_cutoff_sp$`PSICoffee TCS > 0` = ref_mtmalign_sp$PSIcoffee
ref_mtmalign_dmpfold_3_tcs_cutoff_sp$`mTMalign_DMPFOLD TCS > 0` = ref_mtmalign_sp$`mTMalign_DMPFOLD`
ref_mtmalign_dmpfold_3_tcs_cutoff_sp = ref_mtmalign_dmpfold_3_tcs_cutoff_sp[,c(31,1:10,32,11:20,33,21:30)]

#ref_3dcoffee_dmpfold_3_tcs_cutoff_list$ref = ref_3dcoffee_dmpfold_3_tcs_cutoff_list$`TCS > 0`
#for (i in colnames(ref_3dcoffee_dmpfold_3_tcs_cutoff_list)[2:8]) {
#    ref_3dcoffee_dmpfold_3_tcs_cutoff_list[[i]] = (ref_3dcoffee_dmpfold_3_tcs_cutoff_list[[i]] / ref_3dcoffee_dmpfold_3_tcs_cutoff_list[["ref"]]) *100
#}
#ref_3dcoffee_dmpfold_3_tcs_cutoff_list = ref_3dcoffee_dmpfold_3_tcs_cutoff_list[,-c(9)]

sp_mean_seq = c()
sp_sd_seq = c()
sp_mean_psicoffee = c()
sp_sd_psicoffee = c()
sp_mean = c()
sp_sd = c()

for (mod in 1:ncol(weights)) {
  sp_mean_seq = c(sp_mean_seq,wt.mean(x=ref_mtmalign_dmpfold_3_tcs_cutoff_sp[,mod],wt=weights[,mod]))
  sp_sd_seq = c(sp_sd_seq,wt.sd(x=ref_mtmalign_dmpfold_3_tcs_cutoff_sp[,mod],wt=weights[,mod]))
  
  sp_mean_psicoffee = c(sp_mean_psicoffee,wt.mean(x=ref_mtmalign_dmpfold_3_tcs_cutoff_sp[,mod+11],wt=weights[,mod]))
  sp_sd_psicoffee = c(sp_sd_psicoffee,wt.sd(x=ref_mtmalign_dmpfold_3_tcs_cutoff_sp[,mod+11],wt=weights[,mod]))
  
  sp_mean = c(sp_mean,wt.mean(x=ref_mtmalign_dmpfold_3_tcs_cutoff_sp[,mod+22],wt=weights[,mod]))
  sp_sd = c(sp_sd,wt.sd(x=ref_mtmalign_dmpfold_3_tcs_cutoff_sp[,mod+22],wt=weights[,mod]))
}

#sp_mean_seq = t(ref_3dcoffee_dmpfold_3_tcs_cutoff_sp[1:7] %>% summarise_if(is.numeric, mean, na.rm = TRUE))[,1]
#sp_sd_seq = t(ref_3dcoffee_dmpfold_3_tcs_cutoff_sp[1:7] %>% summarise_if(is.numeric, sd, na.rm = TRUE))[,1]
#sp_mean_psicoffee = t(ref_3dcoffee_dmpfold_3_tcs_cutoff_sp[8:14] %>% summarise_if(is.numeric, mean, na.rm = TRUE))[,1]
#sp_sd_psicoffee = t(ref_3dcoffee_dmpfold_3_tcs_cutoff_sp[8:14] %>% summarise_if(is.numeric, sd, na.rm = TRUE))[,1]
#sp_mean = t(ref_3dcoffee_dmpfold_3_tcs_cutoff_sp[15:21] %>% summarise_if(is.numeric, mean, na.rm = TRUE))[,1]
#sp_sd = t(ref_3dcoffee_dmpfold_3_tcs_cutoff_sp[15:21] %>% summarise_if(is.numeric, sd, na.rm = TRUE))[,1]

perc_of_seq_mean = sapply(ref_mtmalign_dmpfold_3_tcs_cutoff_list[,-c(1)], function(x) (sum(x)/sum(ref_mtmalign_dmpfold_3_tcs_cutoff_list[,2]))*100)

ref_mtmalign_dmpfold_3_tcs_cutoff_list_and_sp_means_and_sds = data.frame(sp_mean_seq,sp_sd_seq,sp_mean_psicoffee,sp_sd_psicoffee,sp_mean, sp_sd, perc_of_seq_mean)
row.names(ref_mtmalign_dmpfold_3_tcs_cutoff_list_and_sp_means_and_sds) = colnames(ref_mtmalign_dmpfold_3_tcs_cutoff_list)[-1]
ref_mtmalign_dmpfold_3_tcs_cutoff_list_and_sp_means_and_sds$TCS_cutoff = row.names(ref_mtmalign_dmpfold_3_tcs_cutoff_list_and_sp_means_and_sds)

x = ggplot(ref_mtmalign_dmpfold_3_tcs_cutoff_list_and_sp_means_and_sds, aes(x=TCS_cutoff,group=1)) + scale_x_discrete(limits=ref_mtmalign_dmpfold_3_tcs_cutoff_list_and_sp_means_and_sds$TCS_cutoff) + geom_ribbon(aes(ymin=sp_mean_seq-sp_sd_seq, ymax=sp_mean_seq+sp_sd_seq),fill="lightgreen",alpha=0.3) + geom_ribbon(aes(ymin=sp_mean_psicoffee-sp_sd_psicoffee, ymax=sp_mean_psicoffee+sp_sd_psicoffee),fill="darkorange",alpha=0.3) + geom_ribbon(aes(ymin=sp_mean-sp_sd, ymax=sp_mean+sp_sd),fill="lightblue",alpha=0.3) + geom_line(aes(y=sp_mean_seq,col="sp_mean_seq"))+ geom_line(aes(y=sp_mean_psicoffee,col="sp_mean_psicoffee")) + geom_line(aes(y=sp_mean,col="sp_mean")) + geom_line(aes(y=perc_of_seq_mean,col="perc_of_seq_mean"),linetype="dashed")+ scale_color_manual("", values = c("sp_mean_seq" = "darkolivegreen", "sp_mean_psicoffee" = "darkorange4", "sp_mean" = "blue", "perc_of_seq_mean" = "red"),labels = c("sp_mean_seq" = "Ginsi vs mTMalign_NAT","sp_mean_psicoffee" = "PSICoffee vs mTMalign_NAT", "sp_mean" = "mTMalign_DMPFOLD vs mTMalign_NAT", "perc_of_seq_mean" = "Percent of Sequences")) + scale_y_continuous( name = "Sum-of-Pairs score (%)", sec.axis = sec_axis(~.,name="Percent of sequences included (%)"),limits=c(0,100)) + xlab("TCS cutoff") + theme_light() + theme(axis.text.x=element_text(angle = 45, hjust = 1),legend.position="top") + guides(colour = guide_legend(nrow = 2))
ggsave("ref_mtmalign_sp_score_and_perc_of_seqs_dmpfold_vs_tcs_cutoff.png",dpi="retina")

x = ggplot(ref_mtmalign_dmpfold_3_tcs_cutoff_list_and_sp_means_and_sds, aes(x=TCS_cutoff,group=1)) + 
  scale_x_discrete(limits=ref_mtmalign_dmpfold_3_tcs_cutoff_list_and_sp_means_and_sds$TCS_cutoff) + 
  geom_ribbon(aes(ymin=sp_mean_seq-sp_sd_seq, ymax=sp_mean_seq+sp_sd_seq),fill="lightgreen",alpha=0.3) +  
  geom_ribbon(aes(ymin=sp_mean-sp_sd, ymax=sp_mean+sp_sd),fill="lightblue",alpha=0.3) + 
  geom_line(aes(y=sp_mean_seq,col="sp_mean_seq"))+ 
  geom_line(aes(y=sp_mean,col="sp_mean")) +  
  geom_line(aes(y=perc_of_seq_mean,col="perc_of_seq_mean"),linetype="dashed") + 
  scale_color_manual("", values = c("sp_mean_seq" = "darkolivegreen", "sp_mean" = "blue", "perc_of_seq_mean" = "red"),labels = c("sp_mean_seq" = "SoP Ginsi vs mTMalign_NAT","sp_mean" = "SoP mTMalign_DMPFOLD vs mTMalign_NAT", "perc_of_seq_mean" = "Percent of Sequences")) + 
  scale_y_continuous( name = "Percent (%)", sec.axis = sec_axis(~.,name="Percent of sequences included (%)"),limits=c(0,100)) + xlab("TCS cutoff") + 
  theme_light() + theme(axis.text.x=element_text(angle = 45, hjust = 1),legend.position="top") + 
  guides(colour = guide_legend(nrow = 2)) + geom_vline(xintercept = "TCS > 0") + geom_vline(xintercept = "TCS > 75")
ggsave("ref_mtmalign_sp_score_and_perc_of_seqs_dmpfold_vs_tcs_cutoff_simple.png",dpi="retina")


ref_mtmalign_dmpfold_3_tcs_cutoff_list$ref = ref_mtmalign_dmpfold_3_tcs_cutoff_list[,2]
for (mode in 2:ncol(ref_mtmalign_dmpfold_3_tcs_cutoff_list)) {
  ref_mtmalign_dmpfold_3_tcs_cutoff_list[,mode] = (ref_mtmalign_dmpfold_3_tcs_cutoff_list[,mode]/ref_mtmalign_dmpfold_3_tcs_cutoff_list$ref)*100
}

ref_mtmalign_dmpfold_3_tcs_cutoff_list = ref_mtmalign_dmpfold_3_tcs_cutoff_list[,-c(13)]
ref_mtmalign_dmpfold_3_tcs_cutoff_list_melted = melt(ref_mtmalign_dmpfold_3_tcs_cutoff_list,variable.name="TCS_mode",value.name="percent_of_seq")
ref_mtmalign_dmpfold_3_tcs_cutoff_sp_melted = melt(ref_mtmalign_dmpfold_3_tcs_cutoff_sp[,-c(12:22)],variable.name="TCS_mode",value.name="SP")
ref_mtmalign_dmpfold_3_tcs_cutoff_list_and_sp = cbind(ref_mtmalign_dmpfold_3_tcs_cutoff_list_melted,ref_mtmalign_dmpfold_3_tcs_cutoff_sp_melted)
colnames(ref_mtmalign_dmpfold_3_tcs_cutoff_list_and_sp) = make.unique(colnames(ref_mtmalign_dmpfold_3_tcs_cutoff_list_and_sp))
ref_mtmalign_dmpfold_3_tcs_cutoff_list_and_sp$group = c(rep("SoP Ginsi vs mTMalign_NAT",nrow(ref_mtmalign_dmpfold_3_tcs_cutoff_list_and_sp)/2),rep("SoP mTMalign_DMPFOLD vs mTMalign_NAT",nrow(ref_mtmalign_dmpfold_3_tcs_cutoff_list_and_sp)/2))
ref_mtmalign_dmpfold_3_tcs_cutoff_list_and_sp$Family = factor(ref_mtmalign_dmpfold_3_tcs_cutoff_list_and_sp$Family,levels=row.names(nirmsd_ref)[-nrow(nirmsd_ref)])
ref_mtmalign_dmpfold_3_tcs_cutoff_list_and_sp = ref_mtmalign_dmpfold_3_tcs_cutoff_list_and_sp[order(ref_mtmalign_dmpfold_3_tcs_cutoff_list_and_sp$Family,ref_mtmalign_dmpfold_3_tcs_cutoff_list_and_sp$TCS_mode),]
circular_plot(ref_mtmalign_dmpfold_3_tcs_cutoff_list_and_sp,'SoP Ginsi vs mTMalign_NAT','SoP mTMalign_DMPFOLD vs mTMalign_NAT',"ref_mtmalign_sp_score_per_family_dmpfold_with_tcs_cutoff_circular.png",color_codes_for_circ)

for (threshold in c("75","80","85","90","95")) {
  ref_mtmalign_pair_sp_tcs_above_75 = read.table(paste0("./selected_comparisons_ref_mtmalign_pair_sp.txt_tcs_above_",threshold,"_ref_mtmalign_dmpfold"),header = F,stringsAsFactors = F)
  ref_mtmalign_pair_sp_tcs_above_75 = ref_mtmalign_pair_sp_tcs_above_75[!duplicated(as.list(ref_mtmalign_pair_sp_tcs_above_75))]
  ref_mtmalign_pair_sp_tcs_above_75 = ref_mtmalign_pair_sp_tcs_above_75[,c(1:3,5,8,9)]
  colnames(ref_mtmalign_pair_sp_tcs_above_75) = c("Sequence_1","Sequence_2","Family","Ginsi","PSIcoffee","mTMalign_DMPFOLD")
  #ref_mtmalign_pair_sp_tcs_above_75$Family = factor(ref_mtmalign_pair_sp_tcs_above_75$Family,levels=row.names(nirmsd_ref)[-nrow(nirmsd_ref)])
  ref_mtmalign_pair_sp_tcs_above_75$Family = factor(ref_mtmalign_pair_sp_tcs_above_75$Family)


  selected_fams = levels(as.factor(ref_mtmalign_pair_sp_tcs_above_75$Family))
  selected_colors = color_per_fam[color_per_fam$Family %in% selected_fams,c(2)]

  ο=ggplot(ref_mtmalign_pair_sp_tcs_above_75,aes(x=`mTMalign_DMPFOLD`,y=PSIcoffee,color=Family)) + theme_light() + theme(legend.position="top") +geom_point()+xlim(0,100)+ylim(0,100)+geom_abline(intercept =0 , slope = 1,lwd=0.2)+xlab("Sum-of-Pairs score between pairs of sequences - mTMalign_DMPFOLD vs mTMalign_NAT (%)")+ylab("Sum-of-Pairs score between pairs of sequences - PSICoffee vs mTMalign_NAT (%)") + scale_color_manual(values=as.character(selected_colors)) + ggtitle(paste0("Sum-of-Pairs score between pairs of sequences with TCS > ",threshold))
  ggsave(paste0("ref_mtmalign_pair_sp_score_per_family_dmpfold_vs_psicoffee_tcs_above_",threshold,".png"),dpi = "retina")
  ο=ggplot(ref_mtmalign_pair_sp_tcs_above_75,aes(x=`mTMalign_DMPFOLD`,y=Ginsi,color=Family)) + theme_light() + theme(legend.position="top") +geom_point()+xlim(0,100)+ylim(0,100)+geom_abline(intercept =0 , slope = 1,lwd=0.2)+xlab("Sum-of-Pairs score between pairs of sequences - mTMalign_DMPFOLD vs mTMalign_NAT (%)")+ylab("Sum-of-Pairs score between pairs of sequences - Ginsi vs mTMalign_NAT (%)") + scale_color_manual(values=as.character(selected_colors)) + ggtitle(paste0("Sum-of-Pairs score between pairs of sequences with TCS > ",threshold))
  ggsave(paste0("ref_mtmalign_pair_sp_score_per_family_dmpfold_vs_ginsi_tcs_above_",threshold,".png"),dpi = "retina")
}



ref_mtmalign_alphafold_tcs_cutoff_list = read.table("./all_tcs_cutoffs_all_fams_selected.list_ref_mtmalign_alphafold",header = F,stringsAsFactors = F)
colnames(ref_mtmalign_alphafold_tcs_cutoff_list) = c("Family","TCS > 0","TCS > 50","TCS > 55","TCS > 60","TCS > 65","TCS > 70","TCS > 75","TCS > 80","TCS > 85","TCS > 90","TCS > 95")
weights = sapply(ref_mtmalign_alphafold_tcs_cutoff_list[,-c(1)], function(x) x/sum(x))

ref_mtmalign_alphafold_tcs_cutoff_sp = read.table("./all_tcs_cutoffs_all_fams_selected_ref_ginsi_psicoffee_and_mtmalign_alphafold_vs_ref_mtmalign.sp",header = F,stringsAsFactors = F)
ref_mtmalign_alphafold_tcs_cutoff_sp = ref_mtmalign_alphafold_tcs_cutoff_sp[,which(sapply(ref_mtmalign_alphafold_tcs_cutoff_sp,class) != "character")]
#ref_mtmalign_dmpfold_3_tcs_cutoff_sp = ref_mtmalign_dmpfold_3_tcs_cutoff_sp[!duplicated(as.list(ref_mtmalign_dmpfold_3_tcs_cutoff_sp))]

col_vector = c()
for (method in c("Gins1 TCS > ", "PSICoffee TCS > ", "mTMalign_AF2 TCS > ")){
  for (cutoff in seq(50, 95, by=5)) {
      col_vector = c(col_vector,paste0(method,cutoff))
  }
}
colnames(ref_mtmalign_alphafold_tcs_cutoff_sp) = col_vector

ref_mtmalign_alphafold_tcs_cutoff_sp$`Gins1 TCS > 0` = ref_mtmalign_sp$Ginsi
ref_mtmalign_alphafold_tcs_cutoff_sp$`PSICoffee TCS > 0` = ref_mtmalign_sp$PSIcoffee
ref_mtmalign_alphafold_tcs_cutoff_sp$`mTMalign_AF2 TCS > 0` = ref_mtmalign_sp$`mTMalign_AF2`
ref_mtmalign_alphafold_tcs_cutoff_sp = ref_mtmalign_alphafold_tcs_cutoff_sp[,c(31,1:10,32,11:20,33,21:30)]

sp_mean_seq = c()
sp_sd_seq = c()
sp_mean_psicoffee = c()
sp_sd_psicoffee = c()
sp_mean = c()
sp_sd = c()
gdt_mean = c()
gdt_sd = c()
plddt_mean = c()
plddt_sd = c()

for (mod in 1:ncol(weights)) {
  sp_mean_seq = c(sp_mean_seq,wt.mean(x=ref_mtmalign_alphafold_tcs_cutoff_sp[,mod],wt=weights[,mod]))
  sp_sd_seq = c(sp_sd_seq,wt.sd(x=ref_mtmalign_alphafold_tcs_cutoff_sp[,mod],wt=weights[,mod]))
  
  sp_mean_psicoffee = c(sp_mean_psicoffee,wt.mean(x=ref_mtmalign_alphafold_tcs_cutoff_sp[,mod+11],wt=weights[,mod]))
  sp_sd_psicoffee = c(sp_sd_psicoffee,wt.sd(x=ref_mtmalign_alphafold_tcs_cutoff_sp[,mod+11],wt=weights[,mod]))
  
  sp_mean = c(sp_mean,wt.mean(x=ref_mtmalign_alphafold_tcs_cutoff_sp[,mod+22],wt=weights[,mod]))
  sp_sd = c(sp_sd,wt.sd(x=ref_mtmalign_alphafold_tcs_cutoff_sp[,mod+22],wt=weights[,mod]))

  plddt_mean = c(plddt_mean,wt.mean(x=agg_plddts[,mod],wt=weights[,mod]))
  plddt_sd = c(plddt_sd,wt.sd(x=agg_plddts[,mod],wt=weights[,mod]))

  gdt_mean = c(gdt_mean,wt.mean(x=agg_gdts[,mod],wt=weights[,mod]))
  gdt_sd = c(gdt_sd,wt.sd(x=agg_gdts[,mod],wt=weights[,mod]))
}

perc_of_seq_mean = sapply(ref_mtmalign_alphafold_tcs_cutoff_list[,-c(1)], function(x) (sum(x)/sum(ref_mtmalign_alphafold_tcs_cutoff_list[,2]))*100)

ref_mtmalign_alphafold_tcs_cutoff_list_and_sp_means_and_sds = data.frame(sp_mean_seq,sp_sd_seq,sp_mean_psicoffee,sp_sd_psicoffee,sp_mean, sp_sd,plddt_mean,plddt_sd,gdt_mean,gdt_sd,perc_of_seq_mean)
row.names(ref_mtmalign_alphafold_tcs_cutoff_list_and_sp_means_and_sds) = colnames(ref_mtmalign_alphafold_tcs_cutoff_list)[-1]
ref_mtmalign_alphafold_tcs_cutoff_list_and_sp_means_and_sds$TCS_cutoff = row.names(ref_mtmalign_alphafold_tcs_cutoff_list_and_sp_means_and_sds)

x = ggplot(ref_mtmalign_alphafold_tcs_cutoff_list_and_sp_means_and_sds, aes(x=TCS_cutoff,group=1)) + scale_x_discrete(limits=ref_mtmalign_alphafold_tcs_cutoff_list_and_sp_means_and_sds$TCS_cutoff) + geom_ribbon(aes(ymin=sp_mean_seq-sp_sd_seq, ymax=sp_mean_seq+sp_sd_seq),fill="lightgreen",alpha=0.3) + geom_ribbon(aes(ymin=sp_mean_psicoffee-sp_sd_psicoffee, ymax=sp_mean_psicoffee+sp_sd_psicoffee),fill="darkorange",alpha=0.3) + geom_ribbon(aes(ymin=sp_mean-sp_sd, ymax=sp_mean+sp_sd),fill="lightblue",alpha=0.3) + geom_line(aes(y=sp_mean_seq,col="sp_mean_seq"))+ geom_line(aes(y=sp_mean_psicoffee,col="sp_mean_psicoffee")) + geom_line(aes(y=sp_mean,col="sp_mean")) + geom_line(aes(y=plddt_mean,col="plddt_mean")) + geom_line(aes(y=gdt_mean,col="gdt_mean")) + geom_line(aes(y=perc_of_seq_mean,col="perc_of_seq_mean"),linetype="dashed")+ scale_color_manual("", values = c("sp_mean_seq" = "darkolivegreen", "sp_mean_psicoffee" = "darkorange4", "sp_mean" = "blue", "plddt_mean" = "black", "gdt_mean" = "magenta", "perc_of_seq_mean" = "red"),labels = c("sp_mean_seq" = "SoP Ginsi vs mTMalign_NAT","sp_mean_psicoffee" = "SoP PSICoffee vs mTMalign_NAT", "sp_mean" = "SoP mTMalign_AF2 vs mTMalign_NAT", "plddt_mean" = "AF2 pLDDT", "gdt_mean" = "AF2 GDT_TS", "perc_of_seq_mean" = "Percent of Sequences")) + scale_y_continuous( name = "Percent (%)", sec.axis = sec_axis(~.,name="Percent of sequences included (%)"),limits=c(0,100)) + xlab("TCS cutoff") + theme_minimal() + theme(axis.text.x=element_text(angle = 45, hjust = 1),legend.position="top") + guides(colour = guide_legend(nrow = 2)) # + geom_ribbon(aes(ymin=plddt_mean-plddt_sd, ymax=plddt_mean+plddt_sd),alpha=0.3)+ geom_ribbon(aes(ymin=gdt_mean-gdt_sd, ymax=gdt_mean+gdt_sd),alpha=0.3)
ggsave("ref_mtmalign_sp_score_and_perc_of_seqs_alphafold_vs_tcs_cutoff.png",dpi="retina")

x = ggplot(ref_mtmalign_alphafold_tcs_cutoff_list_and_sp_means_and_sds, aes(x=TCS_cutoff,group=1)) + 
  scale_x_discrete(limits=ref_mtmalign_alphafold_tcs_cutoff_list_and_sp_means_and_sds$TCS_cutoff) + 
  geom_ribbon(aes(ymin=sp_mean_seq-sp_sd_seq, ymax=sp_mean_seq+sp_sd_seq),fill="lightgreen",alpha=0.3) +  
  geom_ribbon(aes(ymin=sp_mean-sp_sd, ymax=sp_mean+sp_sd),fill="lightblue",alpha=0.3) + 
  geom_line(aes(y=sp_mean_seq,col="sp_mean_seq"))+ 
  geom_line(aes(y=sp_mean,col="sp_mean")) +  
  geom_line(aes(y=perc_of_seq_mean,col="perc_of_seq_mean"),linetype="dashed") + 
  scale_color_manual("", values = c("sp_mean_seq" = "darkolivegreen", "sp_mean" = "blue", "perc_of_seq_mean" = "red"),labels = c("sp_mean_seq" = "SoP Ginsi vs mTMalign_NAT","sp_mean" = "SoP mTMalign_AF2 vs mTMalign_NAT", "perc_of_seq_mean" = "Percent of Sequences")) + 
  scale_y_continuous( name = "Percent (%)", sec.axis = sec_axis(~.,name="Percent of sequences included (%)"),limits=c(0,100)) + xlab("TCS cutoff") + 
  theme_light() + theme(axis.text.x=element_text(angle = 45, hjust = 1),legend.position="top") + 
  guides(colour = guide_legend(nrow = 2)) + geom_vline(xintercept = "TCS > 0") + geom_vline(xintercept = "TCS > 90")
ggsave("ref_mtmalign_sp_score_and_perc_of_seqs_alphafold_vs_tcs_cutoff_simple.png",dpi="retina")


x = ggplot(ref_mtmalign_alphafold_tcs_cutoff_list_and_sp_means_and_sds, aes(x=perc_of_seq_mean,group=1)) + geom_ribbon(aes(ymin=sp_mean_seq-sp_sd_seq, ymax=sp_mean_seq+sp_sd_seq),fill="lightgreen",alpha=0.3) + geom_ribbon(aes(ymin=sp_mean_psicoffee-sp_sd_psicoffee, ymax=sp_mean_psicoffee+sp_sd_psicoffee),fill="darkorange",alpha=0.3) + geom_ribbon(aes(ymin=sp_mean-sp_sd, ymax=sp_mean+sp_sd),fill="lightblue",alpha=0.3) + geom_line(aes(y=sp_mean_seq,col="sp_mean_seq"))+ geom_line(aes(y=sp_mean_psicoffee,col="sp_mean_psicoffee")) + geom_line(aes(y=sp_mean,col="sp_mean")) + geom_line(aes(y=plddt_mean,col="plddt_mean")) + geom_line(aes(y=gdt_mean,col="gdt_mean")) + scale_color_manual("", values = c("sp_mean_seq" = "darkolivegreen", "sp_mean_psicoffee" = "darkorange4", "sp_mean" = "blue", "plddt_mean" = "black", "gdt_mean" = "magenta"),labels = c("sp_mean_seq" = "SoP Ginsi vs mTMalign_NAT","sp_mean_psicoffee" = "SoP PSICoffee vs mTMalign_NAT", "sp_mean" = "SoP mTMalign_AF2 vs mTMalign_NAT", "plddt_mean" = "AF2 pLDDT", "gdt_mean" = "AF2 GDT_TS")) + scale_x_reverse() + scale_y_continuous( name = "Percent (%)",limits=c(0,100)) + xlab("Percent of sequences included (%)") + theme_minimal() + theme(axis.text.x=element_text(angle = 45, hjust = 1),legend.position="top") + guides(colour = guide_legend(nrow = 2)) # + geom_ribbon(aes(ymin=plddt_mean-plddt_sd, ymax=plddt_mean+plddt_sd),alpha=0.3)+ geom_ribbon(aes(ymin=gdt_mean-gdt_sd, ymax=gdt_mean+gdt_sd),alpha=0.3)
ggsave("ref_mtmalign_sp_scores_gdt_plddt_vs_perc_of_seqs_alphafold.png",dpi="retina")

ref_mtmalign_alphafold_tcs_cutoff_list$ref = ref_mtmalign_alphafold_tcs_cutoff_list[,2]
for (mode in 2:ncol(ref_mtmalign_alphafold_tcs_cutoff_list)) {
  ref_mtmalign_alphafold_tcs_cutoff_list[,mode] = (ref_mtmalign_alphafold_tcs_cutoff_list[,mode]/ref_mtmalign_alphafold_tcs_cutoff_list$ref)*100
}

ref_mtmalign_alphafold_tcs_cutoff_list = ref_mtmalign_alphafold_tcs_cutoff_list[,-c(13)]
ref_mtmalign_alphafold_tcs_cutoff_list_melted = melt(ref_mtmalign_alphafold_tcs_cutoff_list,variable.name="TCS_mode",value.name="percent_of_seq")
ref_mtmalign_alphafold_tcs_cutoff_sp_melted = melt(ref_mtmalign_alphafold_tcs_cutoff_sp[,-c(12:22)],variable.name="TCS_mode",value.name="SP")
ref_mtmalign_alphafold_tcs_cutoff_list_and_sp = cbind(ref_mtmalign_alphafold_tcs_cutoff_list_melted,ref_mtmalign_alphafold_tcs_cutoff_sp_melted)
colnames(ref_mtmalign_alphafold_tcs_cutoff_list_and_sp) = make.unique(colnames(ref_mtmalign_alphafold_tcs_cutoff_list_and_sp))
ref_mtmalign_alphafold_tcs_cutoff_list_and_sp$group = c(rep("SoP Ginsi vs mTMalign_NAT",nrow(ref_mtmalign_alphafold_tcs_cutoff_list_and_sp)/2),rep("SoP mTMalign_AF2 vs mTMalign_NAT",nrow(ref_mtmalign_alphafold_tcs_cutoff_list_and_sp)/2))
ref_mtmalign_alphafold_tcs_cutoff_list_and_sp$Family = factor(ref_mtmalign_alphafold_tcs_cutoff_list_and_sp$Family,levels=row.names(nirmsd_ref)[-nrow(nirmsd_ref)])
ref_mtmalign_alphafold_tcs_cutoff_list_and_sp = ref_mtmalign_alphafold_tcs_cutoff_list_and_sp[order(ref_mtmalign_alphafold_tcs_cutoff_list_and_sp$Family,ref_mtmalign_alphafold_tcs_cutoff_list_and_sp$TCS_mode),]
circular_plot(ref_mtmalign_alphafold_tcs_cutoff_list_and_sp,'SoP Ginsi vs mTMalign_NAT','SoP mTMalign_AF2 vs mTMalign_NAT',"ref_mtmalign_sp_score_per_family_alphafold_with_tcs_cutoff_circular.png",color_codes_for_circ)


for (threshold in c("75","80","85","90","95")) {
  ref_mtmalign_pair_sp_tcs_above_75 = read.table(paste0("./selected_comparisons_ref_mtmalign_pair_sp.txt_tcs_above_",threshold,"_ref_mtmalign_alphafold"),header = F,stringsAsFactors = F)
  ref_mtmalign_pair_sp_tcs_above_75 = ref_mtmalign_pair_sp_tcs_above_75[!duplicated(as.list(ref_mtmalign_pair_sp_tcs_above_75))]
  ref_mtmalign_pair_sp_tcs_above_75 = ref_mtmalign_pair_sp_tcs_above_75[,c(1:3,5,8,12)]
  colnames(ref_mtmalign_pair_sp_tcs_above_75) = c("Sequence_1","Sequence_2","Family","Ginsi","PSIcoffee","mTMalign_AF2")
  #ref_mtmalign_pair_sp_tcs_above_75$Family = factor(ref_mtmalign_pair_sp_tcs_above_75$Family,levels=row.names(nirmsd_ref)[-nrow(nirmsd_ref)])
  ref_mtmalign_pair_sp_tcs_above_75$Family = factor(ref_mtmalign_pair_sp_tcs_above_75$Family)

  selected_fams = levels(as.factor(ref_mtmalign_pair_sp_tcs_above_75$Family))
  selected_colors = color_per_fam[color_per_fam$Family %in% selected_fams,c(2)]

  ο=ggplot(ref_mtmalign_pair_sp_tcs_above_75,aes(x=`mTMalign_AF2`,y=PSIcoffee,color=Family)) + theme_light() + theme(legend.position="top") +geom_point()+xlim(0,100)+ylim(0,100)+geom_abline(intercept =0 , slope = 1,lwd=0.2)+xlab("Sum-of-Pairs score between pairs of sequences - mTMalign_AF2 vs mTMalign_NAT (%)")+ylab("Sum-of-Pairs score between pairs of sequences - PSICoffee vs mTMalign_NAT (%)") + scale_color_manual(values=as.character(selected_colors)) + ggtitle(paste0("Sum-of-Pairs score between pairs of sequences with TCS > ",threshold))
  ggsave(paste0("ref_mtmalign_pair_sp_score_per_family_alphafold_vs_psicoffee_tcs_above_",threshold,".png"),dpi = "retina")
  ο=ggplot(ref_mtmalign_pair_sp_tcs_above_75,aes(x=`mTMalign_AF2`,y=Ginsi,color=Family)) + theme_light() + theme(legend.position="top") +geom_point()+xlim(0,100)+ylim(0,100)+geom_abline(intercept =0 , slope = 1,lwd=0.2)+xlab("Sum-of-Pairs score between pairs of sequences - mTMalign_AF2 vs mTMalign_NAT (%)")+ylab("Sum-of-Pairs score between pairs of sequences - Ginsi vs mTMalign_NAT (%)") + scale_color_manual(values=as.character(selected_colors)) + ggtitle(paste0("Sum-of-Pairs score between pairs of sequences with TCS > ",threshold))
  ggsave(paste0("ref_mtmalign_pair_sp_score_per_family_alphafold_vs_ginsi_tcs_above_",threshold,".png"),dpi = "retina")
}


#
# Ref 3dcoffee
#

ref_3dcoffee_avg_sp_and_gdt_ts_complete_with_seqs = merge(gdt_ts_complete_with_seqs,ref_3dcoffee_avg_sp,by="Sequence")
colnames(ref_3dcoffee_avg_sp_and_gdt_ts_complete_with_seqs)[5] = "Family"
ref_3dcoffee_avg_sp_and_gdt_ts_complete_with_seqs$Family = factor(ref_3dcoffee_avg_sp_and_gdt_ts_complete_with_seqs$Family,levels=row.names(nirmsd_ref)[-nrow(nirmsd_ref)])
p = ggscatter(ref_3dcoffee_avg_sp_and_gdt_ts_complete_with_seqs,x="3DCoffee_DMPFOLD",y="dmpfold",add = "reg.line",conf.int = TRUE,color="Family",palette=color_codes, add.params = list(color = "blue",fill = "lightgray")) + stat_cor(method = "pearson") + xlab("Sum-of-Pairs score per sequence - 3DCoffee_DMPFOLD vs 3DCoffee_NAT (%)") + ylab("GDT_TS DMPFOLD vs NAT structures") + font("xlab", size = 12) + font("ylab", size = 12)
#ggsave("ref_3dcoffee_avg_sp_score_per_seq_per_family_vs_gdt_score_dmpfold.png",dpi="retina")
p = ggscatter(ref_3dcoffee_avg_sp_and_gdt_ts_complete_with_seqs,x="3DCoffee_AF2",y="alphafold",add = "reg.line",conf.int = TRUE,color="Family",palette=color_codes, add.params = list(color = "blue",fill = "lightgray")) + stat_cor(method = "pearson") + xlab("Sum-of-Pairs score per sequence - 3DCoffee_AF2 vs 3DCoffee_NAT (%)") + ylab("GDT_TS AF2 vs NAT structures") + font("xlab", size = 12) + font("ylab", size = 12)
#ggsave("ref_3dcoffee_avg_sp_score_per_seq_per_family_vs_gdt_score_alphafold.png",dpi="retina")

ref_3dcoffee_avg_sp_with_tcs_score_per_seq = merge(ref_3dcoffee_avg_sp,tcs_score_per_seq,by="Sequence")
colnames(ref_3dcoffee_avg_sp_with_tcs_score_per_seq)[2] = "Family"
ref_3dcoffee_avg_sp_with_tcs_score_per_seq$Family = factor(ref_3dcoffee_avg_sp_with_tcs_score_per_seq$Family,levels=row.names(nirmsd_ref)[-nrow(nirmsd_ref)])
p = ggscatter(ref_3dcoffee_avg_sp_with_tcs_score_per_seq,x="3DCoffee_DMPFOLD.x",y="3DCoffee_DMPFOLD.y",add = "reg.line",conf.int = TRUE,color="Family",palette=color_codes, add.params = list(color = "blue",fill = "lightgray")) + stat_cor(method = "pearson") + xlab("Sum-of-Pairs score per sequence - 3DCoffee_DMPFOLD vs 3DCoffee_NAT (%)") + ylab("TCS score per sequence - 3DCoffee_DMPFOLD") + font("xlab", size = 12) + font("ylab", size = 12)
#ggsave("ref_3dcoffee_avg_sp_score_per_seq_per_family_vs_tcs_score_dmpfold.png",dpi="retina")
p = ggscatter(ref_3dcoffee_avg_sp_with_tcs_score_per_seq,x="3DCoffee_AF2.x",y="3DCoffee_AF2.y",add = "reg.line",conf.int = TRUE,color="Family",palette=color_codes, add.params = list(color = "blue",fill = "lightgray")) + stat_cor(method = "pearson") + xlab("Sum-of-Pairs score per sequence - 3DCoffee_AF2 vs 3DCoffee_NAT (%)") + ylab("TCS score per sequence - 3DCoffee_AF2") + font("xlab", size = 12) + font("ylab", size = 12)
#ggsave("ref_3dcoffee_avg_sp_score_per_seq_per_family_vs_tcs_score_alphafold.png",dpi="retina")

ref_3dcoffee_avg_tcs_with_gdt_ts_complete_score_per_seq = merge(tcs_score_per_seq,gdt_ts_complete_with_seqs,by="Sequence")
colnames(ref_3dcoffee_avg_tcs_with_gdt_ts_complete_score_per_seq)[2] = "Family"
p = ggscatter(ref_3dcoffee_avg_tcs_with_gdt_ts_complete_score_per_seq,x="3DCoffee_DMPFOLD",y="dmpfold",add = "reg.line",conf.int = TRUE,color="Family",palette=color_codes, add.params = list(color = "blue",fill = "lightgray")) + stat_cor(method = "pearson") + xlab("TCS score per sequence - 3DCoffee_DMPFOLD") + ylab("GDT_TS DMPFOLD vs NAT structures") + font("xlab", size = 12) + font("ylab", size = 12)
#ggsave("ref_3dcoffee_avg_tcs_score_per_seq_per_family_vs_gdt_ts_score_dmpfold.png",dpi="retina")
p = ggscatter(ref_3dcoffee_avg_tcs_with_gdt_ts_complete_score_per_seq,x="3DCoffee_AF2",y="alphafold",add = "reg.line",conf.int = TRUE,color="Family",palette=color_codes, add.params = list(color = "blue",fill = "lightgray")) + stat_cor(method = "pearson") + xlab("TCS score per sequence - 3DCoffee_AF2") + ylab("GDT_TS AF2 vs NAT structures") + font("xlab", size = 12) + font("ylab", size = 12)
#ggsave("ref_3dcoffee_avg_tcs_score_per_seq_per_family_vs_gdt_ts_score_alphafold.png",dpi="retina")

ref_3dcoffee_dmpfold_3_tcs_cutoff_list = read.table("./all_tcs_cutoffs_all_fams_selected.list_ref_3dcoffee_dmpfold",header = F,stringsAsFactors = F)
colnames(ref_3dcoffee_dmpfold_3_tcs_cutoff_list) = c("Family","TCS > 0","TCS > 50","TCS > 55","TCS > 60","TCS > 65","TCS > 70","TCS > 75","TCS > 80","TCS > 85","TCS > 90","TCS > 95")
weights = sapply(ref_3dcoffee_dmpfold_3_tcs_cutoff_list[,-c(1)], function(x) x/sum(x))

ref_3dcoffee_dmpfold_3_tcs_cutoff_sp = read.table("./all_tcs_cutoffs_all_fams_selected_ref_ginsi_psicoffee_and_3dcoffee_dmpfold_vs_ref_3dcoffee.sp",header = F,stringsAsFactors = F)
ref_3dcoffee_dmpfold_3_tcs_cutoff_sp = ref_3dcoffee_dmpfold_3_tcs_cutoff_sp[,which(sapply(ref_3dcoffee_dmpfold_3_tcs_cutoff_sp,class) != "character")]
#ref_3dcoffee_dmpfold_3_tcs_cutoff_sp = ref_3dcoffee_dmpfold_3_tcs_cutoff_sp[!duplicated(as.list(ref_3dcoffee_dmpfold_3_tcs_cutoff_sp))]

col_vector = c()
for (method in c("Gins1 TCS > ", "PSICoffee TCS > ", "3DCoffee_DMPFOLD TCS > ")){
  for (cutoff in seq(50, 95, by=5)) {
      col_vector = c(col_vector,paste0(method,cutoff))
  }
}
colnames(ref_3dcoffee_dmpfold_3_tcs_cutoff_sp) = col_vector

ref_3dcoffee_dmpfold_3_tcs_cutoff_sp$`Gins1 TCS > 0` = ref_3dcoffee_sp$Ginsi
ref_3dcoffee_dmpfold_3_tcs_cutoff_sp$`PSICoffee TCS > 0` = ref_3dcoffee_sp$PSIcoffee
ref_3dcoffee_dmpfold_3_tcs_cutoff_sp$`3DCoffee_DMPFOLD TCS > 0` = ref_3dcoffee_sp$`3DCoffee_DMPFOLD`
ref_3dcoffee_dmpfold_3_tcs_cutoff_sp = ref_3dcoffee_dmpfold_3_tcs_cutoff_sp[,c(31,1:10,32,11:20,33,21:30)]

#ref_3dcoffee_dmpfold_3_tcs_cutoff_list$ref = ref_3dcoffee_dmpfold_3_tcs_cutoff_list$`TCS > 0`
#for (i in colnames(ref_3dcoffee_dmpfold_3_tcs_cutoff_list)[2:8]) {
#    ref_3dcoffee_dmpfold_3_tcs_cutoff_list[[i]] = (ref_3dcoffee_dmpfold_3_tcs_cutoff_list[[i]] / ref_3dcoffee_dmpfold_3_tcs_cutoff_list[["ref"]]) *100
#}
#ref_3dcoffee_dmpfold_3_tcs_cutoff_list = ref_3dcoffee_dmpfold_3_tcs_cutoff_list[,-c(9)]

sp_mean_seq = c()
sp_sd_seq = c()
sp_mean_psicoffee = c()
sp_sd_psicoffee = c()
sp_mean = c()
sp_sd = c()

for (mod in 1:ncol(weights)) {
  sp_mean_seq = c(sp_mean_seq,wt.mean(x=ref_3dcoffee_dmpfold_3_tcs_cutoff_sp[,mod],wt=weights[,mod]))
  sp_sd_seq = c(sp_sd_seq,wt.sd(x=ref_3dcoffee_dmpfold_3_tcs_cutoff_sp[,mod],wt=weights[,mod]))
  
  sp_mean_psicoffee = c(sp_mean_psicoffee,wt.mean(x=ref_3dcoffee_dmpfold_3_tcs_cutoff_sp[,mod+11],wt=weights[,mod]))
  sp_sd_psicoffee = c(sp_sd_psicoffee,wt.sd(x=ref_3dcoffee_dmpfold_3_tcs_cutoff_sp[,mod+11],wt=weights[,mod]))
  
  sp_mean = c(sp_mean,wt.mean(x=ref_3dcoffee_dmpfold_3_tcs_cutoff_sp[,mod+22],wt=weights[,mod]))
  sp_sd = c(sp_sd,wt.sd(x=ref_3dcoffee_dmpfold_3_tcs_cutoff_sp[,mod+22],wt=weights[,mod]))
}

#sp_mean_seq = t(ref_3dcoffee_dmpfold_3_tcs_cutoff_sp[1:7] %>% summarise_if(is.numeric, mean, na.rm = TRUE))[,1]
#sp_sd_seq = t(ref_3dcoffee_dmpfold_3_tcs_cutoff_sp[1:7] %>% summarise_if(is.numeric, sd, na.rm = TRUE))[,1]
#sp_mean_psicoffee = t(ref_3dcoffee_dmpfold_3_tcs_cutoff_sp[8:14] %>% summarise_if(is.numeric, mean, na.rm = TRUE))[,1]
#sp_sd_psicoffee = t(ref_3dcoffee_dmpfold_3_tcs_cutoff_sp[8:14] %>% summarise_if(is.numeric, sd, na.rm = TRUE))[,1]
#sp_mean = t(ref_3dcoffee_dmpfold_3_tcs_cutoff_sp[15:21] %>% summarise_if(is.numeric, mean, na.rm = TRUE))[,1]
#sp_sd = t(ref_3dcoffee_dmpfold_3_tcs_cutoff_sp[15:21] %>% summarise_if(is.numeric, sd, na.rm = TRUE))[,1]

perc_of_seq_mean = sapply(ref_3dcoffee_dmpfold_3_tcs_cutoff_list[,-c(1)], function(x) (sum(x)/sum(ref_3dcoffee_dmpfold_3_tcs_cutoff_list[,2]))*100)

ref_3dcoffee_dmpfold_3_tcs_cutoff_list_and_sp_means_and_sds = data.frame(sp_mean_seq,sp_sd_seq,sp_mean_psicoffee,sp_sd_psicoffee,sp_mean, sp_sd, perc_of_seq_mean)
row.names(ref_3dcoffee_dmpfold_3_tcs_cutoff_list_and_sp_means_and_sds) = colnames(ref_3dcoffee_dmpfold_3_tcs_cutoff_list)[-1]
ref_3dcoffee_dmpfold_3_tcs_cutoff_list_and_sp_means_and_sds$TCS_cutoff = row.names(ref_3dcoffee_dmpfold_3_tcs_cutoff_list_and_sp_means_and_sds)

x = ggplot(ref_3dcoffee_dmpfold_3_tcs_cutoff_list_and_sp_means_and_sds, aes(x=TCS_cutoff,group=1)) + scale_x_discrete(limits=ref_3dcoffee_dmpfold_3_tcs_cutoff_list_and_sp_means_and_sds$TCS_cutoff) + geom_ribbon(aes(ymin=sp_mean_seq-sp_sd_seq, ymax=sp_mean_seq+sp_sd_seq),fill="lightgreen",alpha=0.3) + geom_ribbon(aes(ymin=sp_mean_psicoffee-sp_sd_psicoffee, ymax=sp_mean_psicoffee+sp_sd_psicoffee),fill="darkorange",alpha=0.3) + geom_ribbon(aes(ymin=sp_mean-sp_sd, ymax=sp_mean+sp_sd),fill="lightblue",alpha=0.3) + geom_line(aes(y=sp_mean_seq,col="sp_mean_seq"))+ geom_line(aes(y=sp_mean_psicoffee,col="sp_mean_psicoffee")) + geom_line(aes(y=sp_mean,col="sp_mean")) + geom_line(aes(y=perc_of_seq_mean,col="perc_of_seq_mean"),linetype="dashed")+ scale_color_manual("", values = c("sp_mean_seq" = "darkolivegreen", "sp_mean_psicoffee" = "darkorange4", "sp_mean" = "blue", "perc_of_seq_mean" = "red"),labels = c("sp_mean_seq" = "Ginsi vs 3DCoffee_NAT","sp_mean_psicoffee" = "PSICoffee vs 3DCoffee_NAT", "sp_mean" = "3DCoffee_DMPFOLD vs 3DCoffee_NAT", "perc_of_seq_mean" = "Percent of Sequences")) + scale_y_continuous( name = "Sum-of-Pairs score (%)", sec.axis = sec_axis(~.,name="Percent of sequences included (%)"),limits=c(0,100)) + xlab("TCS cutoff") + theme_minimal() + theme(axis.text.x=element_text(angle = 45, hjust = 1),legend.position="top") + guides(colour = guide_legend(nrow = 2))
ggsave("ref_3dcoffee_sp_score_and_perc_of_seqs_dmpfold_vs_tcs_cutoff.png",dpi="retina")


x = ggplot(ref_3dcoffee_dmpfold_3_tcs_cutoff_list_and_sp_means_and_sds, aes(x=TCS_cutoff,group=1)) + 
  scale_x_discrete(limits=ref_3dcoffee_dmpfold_3_tcs_cutoff_list_and_sp_means_and_sds$TCS_cutoff) + 
  geom_ribbon(aes(ymin=sp_mean_seq-sp_sd_seq, ymax=sp_mean_seq+sp_sd_seq),fill="lightgreen",alpha=0.3) +  
  geom_ribbon(aes(ymin=sp_mean-sp_sd, ymax=sp_mean+sp_sd),fill="lightblue",alpha=0.3) + 
  geom_line(aes(y=sp_mean_seq,col="sp_mean_seq"))+ 
  geom_line(aes(y=sp_mean,col="sp_mean")) +  
  geom_line(aes(y=perc_of_seq_mean,col="perc_of_seq_mean"),linetype="dashed") + 
  scale_color_manual("", values = c("sp_mean_seq" = "darkolivegreen", "sp_mean" = "blue", "perc_of_seq_mean" = "red"),labels = c("sp_mean_seq" = "SoP Ginsi vs 3DCoffee_NAT","sp_mean" = "SoP 3DCoffee_DMPFOLD vs 3DCoffee_NAT", "perc_of_seq_mean" = "Percent of Sequences")) + 
  scale_y_continuous( name = "Percent (%)", sec.axis = sec_axis(~.,name="Percent of sequences included (%)"),limits=c(0,100)) + xlab("TCS cutoff") + 
  theme_light() + theme(axis.text.x=element_text(angle = 45, hjust = 1),legend.position="top") + 
  guides(colour = guide_legend(nrow = 2)) + geom_vline(xintercept = "TCS > 0") + geom_vline(xintercept = "TCS > 75")
ggsave("ref_3dcoffee_sp_score_and_perc_of_seqs_dmpfold_vs_tcs_cutoff_simple.png",dpi="retina")



ref_3dcoffee_dmpfold_3_tcs_cutoff_list$ref = ref_3dcoffee_dmpfold_3_tcs_cutoff_list[,2]
for (mode in 2:ncol(ref_3dcoffee_dmpfold_3_tcs_cutoff_list)) {
  ref_3dcoffee_dmpfold_3_tcs_cutoff_list[,mode] = (ref_3dcoffee_dmpfold_3_tcs_cutoff_list[,mode]/ref_3dcoffee_dmpfold_3_tcs_cutoff_list$ref)*100
}

ref_3dcoffee_dmpfold_3_tcs_cutoff_list = ref_3dcoffee_dmpfold_3_tcs_cutoff_list[,-c(13)]
ref_3dcoffee_dmpfold_3_tcs_cutoff_list_melted = melt(ref_3dcoffee_dmpfold_3_tcs_cutoff_list,variable.name="TCS_mode",value.name="percent_of_seq")
ref_3dcoffee_dmpfold_3_tcs_cutoff_sp_melted = melt(ref_3dcoffee_dmpfold_3_tcs_cutoff_sp[,-c(12:22)],variable.name="TCS_mode",value.name="SP",id.vars=NULL)
ref_3dcoffee_dmpfold_3_tcs_cutoff_list_and_sp = cbind(ref_3dcoffee_dmpfold_3_tcs_cutoff_list_melted,ref_3dcoffee_dmpfold_3_tcs_cutoff_sp_melted)
colnames(ref_3dcoffee_dmpfold_3_tcs_cutoff_list_and_sp) = make.unique(colnames(ref_3dcoffee_dmpfold_3_tcs_cutoff_list_and_sp))
ref_3dcoffee_dmpfold_3_tcs_cutoff_list_and_sp$group = c(rep("SoP Ginsi vs 3DCoffee_NAT",nrow(ref_3dcoffee_dmpfold_3_tcs_cutoff_list_and_sp)/2),rep("SoP 3DCoffee_DMPFOLD vs 3DCoffee_NAT",nrow(ref_3dcoffee_dmpfold_3_tcs_cutoff_list_and_sp)/2))
ref_3dcoffee_dmpfold_3_tcs_cutoff_list_and_sp$Family = factor(ref_3dcoffee_dmpfold_3_tcs_cutoff_list_and_sp$Family,levels=row.names(nirmsd_ref)[-nrow(nirmsd_ref)])
ref_3dcoffee_dmpfold_3_tcs_cutoff_list_and_sp = ref_3dcoffee_dmpfold_3_tcs_cutoff_list_and_sp[order(ref_3dcoffee_dmpfold_3_tcs_cutoff_list_and_sp$Family,ref_3dcoffee_dmpfold_3_tcs_cutoff_list_and_sp$TCS_mode),]
circular_plot(ref_3dcoffee_dmpfold_3_tcs_cutoff_list_and_sp,'SoP Ginsi vs 3DCoffee_NAT','SoP 3DCoffee_DMPFOLD vs 3DCoffee_NAT',"ref_3dcoffee_sp_score_per_family_dmpfold_with_tcs_cutoff_circular.png",color_codes_for_circ)

for (threshold in c("75","80","85","90","95")) {
  ref_3dcoffee_pair_sp_tcs_above_75 = read.table(paste0("./selected_comparisons_ref_3dcoffee_pair_sp.txt_tcs_above_",threshold,"_ref_3dcoffee_dmpfold"),header = F,stringsAsFactors = F)
  ref_3dcoffee_pair_sp_tcs_above_75 = ref_3dcoffee_pair_sp_tcs_above_75[!duplicated(as.list(ref_3dcoffee_pair_sp_tcs_above_75))]
  ref_3dcoffee_pair_sp_tcs_above_75 = ref_3dcoffee_pair_sp_tcs_above_75[,c(1:3,5,8,9)]
  colnames(ref_3dcoffee_pair_sp_tcs_above_75) = c("Sequence_1","Sequence_2","Family","Ginsi","PSIcoffee","3DCoffee_DMPFOLD")
  #ref_3dcoffee_pair_sp_tcs_above_75$Family = factor(ref_3dcoffee_pair_sp_tcs_above_75$Family,levels=row.names(nirmsd_ref)[-nrow(nirmsd_ref)])
  ref_3dcoffee_pair_sp_tcs_above_75$Family = factor(ref_3dcoffee_pair_sp_tcs_above_75$Family)


  selected_fams = levels(as.factor(ref_3dcoffee_pair_sp_tcs_above_75$Family))
  selected_colors = color_per_fam[color_per_fam$Family %in% selected_fams,c(2)]

  ο=ggplot(ref_3dcoffee_pair_sp_tcs_above_75,aes(x=`3DCoffee_DMPFOLD`,y=PSIcoffee,color=Family)) + theme_light() + theme(legend.position="top") +geom_point()+xlim(0,100)+ylim(0,100)+geom_abline(intercept =0 , slope = 1,lwd=0.2)+xlab("Sum-of-Pairs score between pairs of sequences - 3DCoffee_DMPFOLD vs 3DCoffee_NAT (%)")+ylab("Sum-of-Pairs score between pairs of sequences - PSICoffee vs 3DCoffee_NAT (%)") + scale_color_manual(values=as.character(selected_colors)) + ggtitle(paste0("Sum-of-Pairs score between pairs of sequences with TCS > ",threshold))
  ggsave(paste0("ref_3dcoffee_pair_sp_score_per_family_dmpfold_vs_psicoffee_tcs_above_",threshold,".png"),dpi = "retina")
  ο=ggplot(ref_3dcoffee_pair_sp_tcs_above_75,aes(x=`3DCoffee_DMPFOLD`,y=Ginsi,color=Family)) + theme_light() + theme(legend.position="top") +geom_point()+xlim(0,100)+ylim(0,100)+geom_abline(intercept =0 , slope = 1,lwd=0.2)+xlab("Sum-of-Pairs score between pairs of sequences - 3DCoffee_DMPFOLD vs 3DCoffee_NAT (%)")+ylab("Sum-of-Pairs score between pairs of sequences - Ginsi vs 3DCoffee_NAT (%)") + scale_color_manual(values=as.character(selected_colors)) + ggtitle(paste0("Sum-of-Pairs score between pairs of sequences with TCS > ",threshold))
  ggsave(paste0("ref_3dcoffee_pair_sp_score_per_family_dmpfold_vs_ginsi_tcs_above_",threshold,".png"),dpi = "retina")
}



ref_3dcoffee_alphafold_tcs_cutoff_list = read.table("./all_tcs_cutoffs_all_fams_selected.list_ref_3dcoffee_alphafold",header = F,stringsAsFactors = F)
colnames(ref_3dcoffee_alphafold_tcs_cutoff_list) = c("Family","TCS > 0","TCS > 50","TCS > 55","TCS > 60","TCS > 65","TCS > 70","TCS > 75","TCS > 80","TCS > 85","TCS > 90","TCS > 95")
weights = sapply(ref_3dcoffee_alphafold_tcs_cutoff_list[,-c(1)], function(x) x/sum(x))

ref_3dcoffee_alphafold_tcs_cutoff_sp = read.table("./all_tcs_cutoffs_all_fams_selected_ref_ginsi_psicoffee_and_3dcoffee_alphafold_vs_ref_3dcoffee.sp",header = F,stringsAsFactors = F)
ref_3dcoffee_alphafold_tcs_cutoff_sp = ref_3dcoffee_alphafold_tcs_cutoff_sp[,which(sapply(ref_3dcoffee_alphafold_tcs_cutoff_sp,class) != "character")]
#ref_3dcoffee_dmpfold_3_tcs_cutoff_sp = ref_3dcoffee_dmpfold_3_tcs_cutoff_sp[!duplicated(as.list(ref_3dcoffee_dmpfold_3_tcs_cutoff_sp))]

col_vector = c()
for (method in c("Gins1 TCS > ", "PSICoffee TCS > ", "3DCoffee_AF2 TCS > ")){
  for (cutoff in seq(50, 95, by=5)) {
      col_vector = c(col_vector,paste0(method,cutoff))
  }
}
colnames(ref_3dcoffee_alphafold_tcs_cutoff_sp) = col_vector

ref_3dcoffee_alphafold_tcs_cutoff_sp$`Gins1 TCS > 0` = ref_3dcoffee_sp$Ginsi
ref_3dcoffee_alphafold_tcs_cutoff_sp$`PSICoffee TCS > 0` = ref_3dcoffee_sp$PSIcoffee
ref_3dcoffee_alphafold_tcs_cutoff_sp$`3DCoffee_AF2 TCS > 0` = ref_3dcoffee_sp$`3DCoffee_AF2`
ref_3dcoffee_alphafold_tcs_cutoff_sp = ref_3dcoffee_alphafold_tcs_cutoff_sp[,c(31,1:10,32,11:20,33,21:30)]

#ref_3dcoffee_dmpfold_3_tcs_cutoff_list$ref = ref_3dcoffee_dmpfold_3_tcs_cutoff_list$`TCS > 0`
#for (i in colnames(ref_3dcoffee_dmpfold_3_tcs_cutoff_list)[2:8]) {
#    ref_3dcoffee_dmpfold_3_tcs_cutoff_list[[i]] = (ref_3dcoffee_dmpfold_3_tcs_cutoff_list[[i]] / ref_3dcoffee_dmpfold_3_tcs_cutoff_list[["ref"]]) *100
#}
#ref_3dcoffee_dmpfold_3_tcs_cutoff_list = ref_3dcoffee_dmpfold_3_tcs_cutoff_list[,-c(9)]

sp_mean_seq = c()
sp_sd_seq = c()
sp_mean_psicoffee = c()
sp_sd_psicoffee = c()
sp_mean = c()
sp_sd = c()
gdt_mean = c()
gdt_sd = c()
plddt_mean = c()
plddt_sd = c()

for (mod in 1:ncol(weights)) {
  sp_mean_seq = c(sp_mean_seq,wt.mean(x=ref_3dcoffee_alphafold_tcs_cutoff_sp[,mod],wt=weights[,mod]))
  sp_sd_seq = c(sp_sd_seq,wt.sd(x=ref_3dcoffee_alphafold_tcs_cutoff_sp[,mod],wt=weights[,mod]))
  
  sp_mean_psicoffee = c(sp_mean_psicoffee,wt.mean(x=ref_3dcoffee_alphafold_tcs_cutoff_sp[,mod+11],wt=weights[,mod]))
  sp_sd_psicoffee = c(sp_sd_psicoffee,wt.sd(x=ref_3dcoffee_alphafold_tcs_cutoff_sp[,mod+11],wt=weights[,mod]))
  
  sp_mean = c(sp_mean,wt.mean(x=ref_3dcoffee_alphafold_tcs_cutoff_sp[,mod+22],wt=weights[,mod]))
  sp_sd = c(sp_sd,wt.sd(x=ref_3dcoffee_alphafold_tcs_cutoff_sp[,mod+22],wt=weights[,mod]))

  plddt_mean = c(plddt_mean,wt.mean(x=agg_plddts[,mod],wt=weights[,mod]))
  plddt_sd = c(plddt_sd,wt.sd(x=agg_plddts[,mod],wt=weights[,mod]))

  gdt_mean = c(gdt_mean,wt.mean(x=agg_gdts[,mod],wt=weights[,mod]))
  gdt_sd = c(gdt_sd,wt.sd(x=agg_gdts[,mod],wt=weights[,mod]))
}

perc_of_seq_mean = sapply(ref_3dcoffee_alphafold_tcs_cutoff_list[,-c(1)], function(x) (sum(x)/sum(ref_3dcoffee_alphafold_tcs_cutoff_list[,2]))*100)

ref_3dcoffee_alphafold_tcs_cutoff_list_and_sp_means_and_sds = data.frame(sp_mean_seq,sp_sd_seq,sp_mean_psicoffee,sp_sd_psicoffee,sp_mean, sp_sd,plddt_mean,plddt_sd,gdt_mean,gdt_sd,perc_of_seq_mean)
row.names(ref_3dcoffee_alphafold_tcs_cutoff_list_and_sp_means_and_sds) = colnames(ref_3dcoffee_alphafold_tcs_cutoff_list)[-1]
ref_3dcoffee_alphafold_tcs_cutoff_list_and_sp_means_and_sds$TCS_cutoff = row.names(ref_3dcoffee_alphafold_tcs_cutoff_list_and_sp_means_and_sds)

x = ggplot(ref_3dcoffee_alphafold_tcs_cutoff_list_and_sp_means_and_sds, aes(x=TCS_cutoff,group=1)) + scale_x_discrete(limits=ref_3dcoffee_alphafold_tcs_cutoff_list_and_sp_means_and_sds$TCS_cutoff) + geom_ribbon(aes(ymin=sp_mean_seq-sp_sd_seq, ymax=sp_mean_seq+sp_sd_seq),fill="lightgreen",alpha=0.3) + geom_ribbon(aes(ymin=sp_mean_psicoffee-sp_sd_psicoffee, ymax=sp_mean_psicoffee+sp_sd_psicoffee),fill="darkorange",alpha=0.3) + geom_ribbon(aes(ymin=sp_mean-sp_sd, ymax=sp_mean+sp_sd),fill="lightblue",alpha=0.3) + geom_line(aes(y=sp_mean_seq,col="sp_mean_seq"))+ geom_line(aes(y=sp_mean_psicoffee,col="sp_mean_psicoffee")) + geom_line(aes(y=sp_mean,col="sp_mean")) + geom_line(aes(y=plddt_mean,col="plddt_mean")) + geom_line(aes(y=gdt_mean,col="gdt_mean")) + geom_line(aes(y=perc_of_seq_mean,col="perc_of_seq_mean"),linetype="dashed")+ scale_color_manual("", values = c("sp_mean_seq" = "darkolivegreen", "sp_mean_psicoffee" = "darkorange4", "sp_mean" = "blue", "plddt_mean" = "black", "gdt_mean" = "magenta", "perc_of_seq_mean" = "red"),labels = c("sp_mean_seq" = "SoP Ginsi vs 3DCoffee_NAT","sp_mean_psicoffee" = "SoP PSICoffee vs 3DCoffee_NAT", "sp_mean" = "SoP 3DCoffee_AF2 vs 3DCoffee_NAT", "plddt_mean" = "AF2 pLDDT", "gdt_mean" = "AF2 GDT_TS", "perc_of_seq_mean" = "Percent of Sequences")) + scale_y_continuous( name = "Percent (%)", sec.axis = sec_axis(~.,name="Percent of sequences included (%)"),limits=c(0,100)) + xlab("TCS cutoff") + theme_minimal() + theme(axis.text.x=element_text(angle = 45, hjust = 1),legend.position="top") + guides(colour = guide_legend(nrow = 2)) # + geom_ribbon(aes(ymin=plddt_mean-plddt_sd, ymax=plddt_mean+plddt_sd),alpha=0.3)+ geom_ribbon(aes(ymin=gdt_mean-gdt_sd, ymax=gdt_mean+gdt_sd),alpha=0.3)
ggsave("ref_3dcoffee_sp_score_and_perc_of_seqs_alphafold_vs_tcs_cutoff.png",dpi="retina")


x = ggplot(ref_3dcoffee_alphafold_tcs_cutoff_list_and_sp_means_and_sds, aes(x=TCS_cutoff,group=1)) + 
  scale_x_discrete(limits=ref_3dcoffee_alphafold_tcs_cutoff_list_and_sp_means_and_sds$TCS_cutoff) + 
  geom_ribbon(aes(ymin=sp_mean_seq-sp_sd_seq, ymax=sp_mean_seq+sp_sd_seq),fill="lightgreen",alpha=0.3) +  
  geom_ribbon(aes(ymin=sp_mean-sp_sd, ymax=sp_mean+sp_sd),fill="lightblue",alpha=0.3) + 
  geom_line(aes(y=sp_mean_seq,col="sp_mean_seq"))+ 
  geom_line(aes(y=sp_mean,col="sp_mean")) +  
  geom_line(aes(y=perc_of_seq_mean,col="perc_of_seq_mean"),linetype="dashed") + 
  scale_color_manual("", values = c("sp_mean_seq" = "darkolivegreen", "sp_mean" = "blue", "perc_of_seq_mean" = "red"),labels = c("sp_mean_seq" = "SoP Ginsi vs 3DCoffee_NAT","sp_mean" = "SoP 3DCoffee_AF2 vs 3DCoffee_NAT", "perc_of_seq_mean" = "Percent of Sequences")) + 
  scale_y_continuous( name = "Percent (%)", sec.axis = sec_axis(~.,name="Percent of sequences included (%)"),limits=c(0,100)) + xlab("TCS cutoff") + 
  theme_light() + theme(axis.text.x=element_text(angle = 45, hjust = 1),legend.position="top") + 
  guides(colour = guide_legend(nrow = 2)) + geom_vline(xintercept = "TCS > 0") + geom_vline(xintercept = "TCS > 75")
ggsave("ref_3dcoffee_sp_score_and_perc_of_seqs_alphafold_vs_tcs_cutoff_simple.png",dpi="retina")



x = ggplot(ref_3dcoffee_alphafold_tcs_cutoff_list_and_sp_means_and_sds, aes(x=perc_of_seq_mean,group=1)) + geom_ribbon(aes(ymin=sp_mean_seq-sp_sd_seq, ymax=sp_mean_seq+sp_sd_seq),fill="lightgreen",alpha=0.3) + geom_ribbon(aes(ymin=sp_mean_psicoffee-sp_sd_psicoffee, ymax=sp_mean_psicoffee+sp_sd_psicoffee),fill="darkorange",alpha=0.3) + geom_ribbon(aes(ymin=sp_mean-sp_sd, ymax=sp_mean+sp_sd),fill="lightblue",alpha=0.3) + geom_line(aes(y=sp_mean_seq,col="sp_mean_seq"))+ geom_line(aes(y=sp_mean_psicoffee,col="sp_mean_psicoffee")) + geom_line(aes(y=sp_mean,col="sp_mean")) + geom_line(aes(y=plddt_mean,col="plddt_mean")) + geom_line(aes(y=gdt_mean,col="gdt_mean")) + scale_color_manual("", values = c("sp_mean_seq" = "darkolivegreen", "sp_mean_psicoffee" = "darkorange4", "sp_mean" = "blue", "plddt_mean" = "black", "gdt_mean" = "magenta"),labels = c("sp_mean_seq" = "SoP Ginsi vs 3DCoffee_NAT","sp_mean_psicoffee" = "SoP PSICoffee vs 3DCoffee_NAT", "sp_mean" = "SoP 3DCoffee_AF2 vs 3DCoffee_NAT", "plddt_mean" = "AF2 pLDDT", "gdt_mean" = "AF2 GDT_TS")) + scale_x_reverse() + scale_y_continuous( name = "Percent (%)",limits=c(0,100)) + xlab("Percent of sequences included (%)") + theme_minimal() + theme(axis.text.x=element_text(angle = 45, hjust = 1),legend.position="top") + guides(colour = guide_legend(nrow = 2)) # + geom_ribbon(aes(ymin=plddt_mean-plddt_sd, ymax=plddt_mean+plddt_sd),alpha=0.3)+ geom_ribbon(aes(ymin=gdt_mean-gdt_sd, ymax=gdt_mean+gdt_sd),alpha=0.3)
#ggsave("ref_3dcoffee_sp_scores_gdt_plddt_vs_perc_of_seqs_alphafold.png",dpi="retina")

ref_3dcoffee_alphafold_tcs_cutoff_list$ref = ref_3dcoffee_alphafold_tcs_cutoff_list[,2]
for (mode in 2:ncol(ref_3dcoffee_alphafold_tcs_cutoff_list)) {
  ref_3dcoffee_alphafold_tcs_cutoff_list[,mode] = (ref_3dcoffee_alphafold_tcs_cutoff_list[,mode]/ref_3dcoffee_alphafold_tcs_cutoff_list$ref)*100
}

ref_3dcoffee_alphafold_tcs_cutoff_list = ref_3dcoffee_alphafold_tcs_cutoff_list[,-c(13)]
ref_3dcoffee_alphafold_tcs_cutoff_list_melted = melt(ref_3dcoffee_alphafold_tcs_cutoff_list,variable.name="TCS_mode",value.name="percent_of_seq")
ref_3dcoffee_alphafold_tcs_cutoff_sp_melted = melt(ref_3dcoffee_alphafold_tcs_cutoff_sp[,-c(12:22)],variable.name="TCS_mode",value.name="SP")
ref_3dcoffee_alphafold_tcs_cutoff_list_and_sp = cbind(ref_3dcoffee_alphafold_tcs_cutoff_list_melted,ref_3dcoffee_alphafold_tcs_cutoff_sp_melted)
colnames(ref_3dcoffee_alphafold_tcs_cutoff_list_and_sp) = make.unique(colnames(ref_3dcoffee_alphafold_tcs_cutoff_list_and_sp))
ref_3dcoffee_alphafold_tcs_cutoff_list_and_sp$group = c(rep("SoP Ginsi vs 3DCoffee_NAT",nrow(ref_3dcoffee_alphafold_tcs_cutoff_list_and_sp)/2),rep("SoP 3DCoffee_AF2 vs 3DCoffee_NAT",nrow(ref_3dcoffee_alphafold_tcs_cutoff_list_and_sp)/2))
ref_3dcoffee_alphafold_tcs_cutoff_list_and_sp$Family = factor(ref_3dcoffee_alphafold_tcs_cutoff_list_and_sp$Family,levels=row.names(nirmsd_ref)[-nrow(nirmsd_ref)])
ref_3dcoffee_alphafold_tcs_cutoff_list_and_sp = ref_3dcoffee_alphafold_tcs_cutoff_list_and_sp[order(ref_3dcoffee_alphafold_tcs_cutoff_list_and_sp$Family,ref_3dcoffee_alphafold_tcs_cutoff_list_and_sp$TCS_mode),]
circular_plot(ref_3dcoffee_alphafold_tcs_cutoff_list_and_sp,'SoP Ginsi vs 3DCoffee_NAT','SoP 3DCoffee_AF2 vs 3DCoffee_NAT',"ref_3dcoffee_sp_score_per_family_alphafold_with_tcs_cutoff_circular.png",color_codes_for_circ)


for (threshold in c("75","80","85","90","95")) {
  ref_3dcoffee_pair_sp_tcs_above_75 = read.table(paste0("./selected_comparisons_ref_3dcoffee_pair_sp.txt_tcs_above_",threshold,"_ref_3dcoffee_alphafold"),header = F,stringsAsFactors = F)
  ref_3dcoffee_pair_sp_tcs_above_75 = ref_3dcoffee_pair_sp_tcs_above_75[!duplicated(as.list(ref_3dcoffee_pair_sp_tcs_above_75))]
  ref_3dcoffee_pair_sp_tcs_above_75 = ref_3dcoffee_pair_sp_tcs_above_75[,c(1:3,5,8,12)]
  colnames(ref_3dcoffee_pair_sp_tcs_above_75) = c("Sequence_1","Sequence_2","Family","Ginsi","PSIcoffee","3DCoffee_AF2")
  #ref_3dcoffee_pair_sp_tcs_above_75$Family = factor(ref_3dcoffee_pair_sp_tcs_above_75$Family,levels=row.names(nirmsd_ref)[-nrow(nirmsd_ref)])
  ref_3dcoffee_pair_sp_tcs_above_75$Family = factor(ref_3dcoffee_pair_sp_tcs_above_75$Family)

  selected_fams = levels(as.factor(ref_3dcoffee_pair_sp_tcs_above_75$Family))
  selected_colors = color_per_fam[color_per_fam$Family %in% selected_fams,c(2)]

  ο=ggplot(ref_3dcoffee_pair_sp_tcs_above_75,aes(x=`3DCoffee_AF2`,y=PSIcoffee,color=Family)) + theme_light() + theme(legend.position="top") +geom_point()+xlim(0,100)+ylim(0,100)+geom_abline(intercept =0 , slope = 1,lwd=0.2)+xlab("Sum-of-Pairs score between pairs of sequences - 3DCoffee_AF2 vs 3DCoffee_NAT (%)")+ylab("Sum-of-Pairs score between pairs of sequences - PSICoffee vs 3DCoffee_NAT (%)") + scale_color_manual(values=as.character(selected_colors)) + ggtitle(paste0("Sum-of-Pairs score between pairs of sequences with TCS > ",threshold))
  ggsave(paste0("ref_3dcoffee_pair_sp_score_per_family_alphafold_vs_psicoffee_tcs_above_",threshold,".png"),dpi = "retina")
  ο=ggplot(ref_3dcoffee_pair_sp_tcs_above_75,aes(x=`3DCoffee_AF2`,y=Ginsi,color=Family)) + theme_light() + theme(legend.position="top") +geom_point()+xlim(0,100)+ylim(0,100)+geom_abline(intercept =0 , slope = 1,lwd=0.2)+xlab("Sum-of-Pairs score between pairs of sequences - 3DCoffee_AF2 vs 3DCoffee_NAT (%)")+ylab("Sum-of-Pairs score between pairs of sequences - Ginsi vs 3DCoffee_NAT (%)") + scale_color_manual(values=as.character(selected_colors)) + ggtitle(paste0("Sum-of-Pairs score between pairs of sequences with TCS > ",threshold))
  ggsave(paste0("ref_3dcoffee_pair_sp_score_per_family_alphafold_vs_ginsi_tcs_above_",threshold,".png"),dpi = "retina")
}


#
# Summary table
#

#
# Ref mtmalign
#
if (FALSE) {
sop_sum = as.data.frame(colMeans(ref_mtmalign_sp[,c("Ginsi","TCoffee","PSIcoffee","mTMalign_AF2","mTMalign_DMPFOLD")]))
colnames(sop_sum)[1] = "SoP"
nirmsd_sum = as.data.frame(t(nirmsd_ref[13,c("Ginsi","TCoffee","PSICoffee","mTMalign_AF2_REF_NAT","mTMalign_DMPFOLD_REF_NAT")]))
colnames(nirmsd_sum)[1] = "NiRMSD"
tcs_sum = as.data.frame(colMeans(tcs_score[,c("Ginsi","TCoffee","PSIcoffee","mTMalign_AF2","mTMalign_DMPFOLD")]))
colnames(tcs_sum)[1] = "TCS"
sum_table = merge(sop_sum,tcs_sum,by=0,all=T)
row.names(sum_table) = sum_table$Row.names
row.names(nirmsd_sum)[3] = "PSIcoffee"
row.names(nirmsd_sum)[4] = "mTMalign_AF2"
row.names(nirmsd_sum)[5] = "mTMalign_DMPFOLD"
sum_table = merge(sum_table,nirmsd_sum,by=0,all=T)
row.names(sum_table) = sum_table$Row.names
sum_table = sum_table[,-c(1,2)]
sum_table$TCS = sum_table$TCS/10
sum_table = round(sum_table,digits=2)
row.names(sum_table)[1] = "Gins1"
sum_table = sum_table[c(5,1,4,2,3),]
write.table(sum_table,"Table1_summary.tsv",quote=F,row.names=T,col.names=T,sep="\t")
}

#
# Ref 3dcoffee
#

quadrants_df = data.frame(t(data.frame(ref_ginsi_avg_df_delta_perc,ref_psicoffee_avg_df_delta_perc,ref_3dcoffee_AF2_avg_df_delta_perc,ref_3dcoffee_DMPFOLD_avg_df_delta_perc,ref_3dcoffee_NAT_avg_df_delta_perc)))
row.names(quadrants_df) = c("MSA-Seq","MSA-PSI","MSA-AF2","MSA-DMP","MSA-PDB")
colnames(quadrants_df) = c("SoP vs TCS","SoP vs pLDDT","SoP vs GDT-TS")

sop_sum = as.data.frame(colMeans(ref_3dcoffee_sp_same_state[,c("Ginsi","PSIcoffee","3DCoffee_NAT","3DCoffee_AF2","3DCoffee_DMPFOLD")]))
colnames(sop_sum)[1] = "SoP"
nirmsd_sum = as.data.frame(t(nirmsd_ref[13,c("Ginsi","PSIcoffee","3DCoffee_NAT","3DCoffee_AF2_REF_NAT","3DCoffee_DMPFOLD_REF_NAT")]))
colnames(nirmsd_sum)[1] = "NiRMSD"
tcs_sum = as.data.frame(colMeans(tcs_score[,c("Ginsi","PSIcoffee","3DCoffee_NAT","3DCoffee_AF2","3DCoffee_DMPFOLD")]))
colnames(tcs_sum)[1] = "TCS"
sum_table = merge(sop_sum,tcs_sum,by=0,all=T)
row.names(sum_table) = sum_table$Row.names
row.names(nirmsd_sum)[2] = "PSIcoffee"
row.names(nirmsd_sum)[4] = "3DCoffee_AF2"
row.names(nirmsd_sum)[5] = "3DCoffee_DMPFOLD"
sum_table = merge(sum_table,nirmsd_sum,by=0,all=T)
row.names(sum_table) = sum_table$Row.names
sum_table$TCS = sum_table$TCS/10
sum_table = sum_table[,-c(1,2)]
row.names(sum_table) = c("MSA-AF2","MSA-DMP","MSA-PDB","MSA-Seq","MSA-PSI")
sum_table = merge(sum_table,quadrants_df,by=0)
row.names(sum_table) = sum_table$Row.names
sum_table = sum_table[,-c(1)]
sum_table = round(sum_table,digits=2)
sum_table = sum_table[c(5,4,1,2,3),]
write.table(sum_table,"Table1_summary.tsv",quote=F,row.names=T,col.names=T,sep="\t")



#
# Big table
#

  #
  # Ref 3dcoffee
  #
tmp_nirmsd_ref_pair = nirmsd_ref_pair
colnames(tmp_nirmsd_ref_pair)[1:2] = c("Sequence_2","Sequence_1")

msa_seq = Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "Family", all.x = TRUE),
list(ref_3dcoffee_sp_same_state[,c("Family","Ginsi")],tcs_score[,c("Family","Ginsi")],nirmsd_ref[-13,c("Family","Ginsi")],pid_ref[,c("Family","Ginsi")]))
colnames(msa_seq) = c("Family","SoP","TCS","NiRMSD","PID")
write.table(msa_seq,"msa_seq.tsv",quote=FALSE,row.names=FALSE,col.names=T,sep="\t")

msa_seq_avg = Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = c("Family","Sequence"), all.x = TRUE),
list(ref_3dcoffee_avg_sp_same_state[,c("Family","Sequence","Ginsi")],tcs_score_per_seq[,c("Family","Sequence","Ginsi")],nirmsd_ref_avg[,c("Family","Sequence","Ginsi")],pid_ref_avg[,c("Family","Sequence","Ginsi")],gdt_ts_complete_with_seqs[,c("Family","Sequence","alphafold")],alphafold_plddts))
colnames(msa_seq_avg) = c("Family","Sequence","SoP","TCS","NiRMSD","PID","GDT_TS","pLDDT")
msa_seq_avg$GDT_TS = round(msa_seq_avg$GDT_TS*100,digits=2)
msa_seq_avg$pLDDT = round(msa_seq_avg$pLDDT,digits=2)
write.table(msa_seq_avg,"msa_seq_avg.tsv",quote=FALSE,row.names=FALSE,col.names=T,sep="\t")


msa_seq_pair = Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = c("Family","Sequence_1","Sequence_2"), all.x = TRUE),
list(ref_3dcoffee_pair_sp_same_state[,c("Family","Sequence_1","Sequence_2","Ginsi")],tmp_nirmsd_ref_pair[,c("Family","Sequence_1","Sequence_2","Ginsi")],pid_ref_pair[,c("Family","Sequence_1","Sequence_2","Ginsi")]))
msa_seq_pair = msa_seq_pair %>% 
  left_join(tcs_score_per_seq[,c("Family","Sequence","Ginsi")],by = c("Sequence_1" = "Sequence", "Family" = "Family")) %>%
  left_join(tcs_score_per_seq[,c("Family","Sequence","Ginsi")],by = c("Sequence_2" = "Sequence", "Family" = "Family")) %>%
  left_join(gdt_ts_complete_with_seqs[,c("Family","Sequence","alphafold")],by = c("Sequence_1" = "Sequence", "Family" = "Family")) %>%
  left_join(gdt_ts_complete_with_seqs[,c("Family","Sequence","alphafold")],by = c("Sequence_2" = "Sequence", "Family" = "Family")) %>% 
  left_join(alphafold_plddts,by = c("Sequence_1" = "Sequence", "Family" = "Family")) %>%
  left_join(alphafold_plddts,by = c("Sequence_2" = "Sequence", "Family" = "Family"))
colnames(msa_seq_pair) = c("Family","Sequence_1","Sequence_2","SoP","NiRMSD","PID","TCS_1","TCS_2","GDT_TS_1","GDT_TS_2","pLDDT_1","pLDDT_2")
msa_seq_pair$GDT_TS_1 = round(msa_seq_pair$GDT_TS_1*100,digits=2)
msa_seq_pair$GDT_TS_2 = round(msa_seq_pair$GDT_TS_2*100,digits=2)
msa_seq_pair$pLDDT_1 = round(msa_seq_pair$pLDDT_1,digits=2)
msa_seq_pair$pLDDT_2 = round(msa_seq_pair$pLDDT_2,digits=2)  
write.table(msa_seq_pair,"msa_seq_pair.tsv",quote=FALSE,row.names=FALSE,col.names=T,sep="\t")



msa_psi = Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "Family", all.x = TRUE),
list(ref_3dcoffee_sp_same_state[,c("Family","PSIcoffee")],tcs_score[,c("Family","PSIcoffee")],nirmsd_ref[-13,c("Family","PSIcoffee")],pid_ref[,c("Family","PSIcoffee")]))
colnames(msa_psi) = c("Family","SoP","TCS","NiRMSD","PID")
write.table(msa_psi,"msa_psi.tsv",quote=FALSE,row.names=FALSE,col.names=T,sep="\t")

msa_psi_avg = Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = c("Family","Sequence"), all.x = TRUE),
list(ref_3dcoffee_avg_sp_same_state[,c("Family","Sequence","PSIcoffee")],tcs_score_per_seq[,c("Family","Sequence","PSIcoffee")],nirmsd_ref_avg[,c("Family","Sequence","PSIcoffee")],pid_ref_avg[,c("Family","Sequence","PSIcoffee")],gdt_ts_complete_with_seqs[,c("Family","Sequence","alphafold")],alphafold_plddts))
colnames(msa_psi_avg) = c("Family","Sequence","SoP","TCS","NiRMSD","PID","GDT_TS","pLDDT")
msa_psi_avg$GDT_TS = round(msa_psi_avg$GDT_TS*100,digits=2)
msa_psi_avg$pLDDT = round(msa_psi_avg$pLDDT,digits=2)
write.table(msa_psi_avg,"msa_psi_avg.tsv",quote=FALSE,row.names=FALSE,col.names=T,sep="\t")


msa_psi_pair = Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = c("Family","Sequence_1","Sequence_2"), all.x = TRUE),
list(ref_3dcoffee_pair_sp_same_state[,c("Family","Sequence_1","Sequence_2","PSIcoffee")],tmp_nirmsd_ref_pair[,c("Family","Sequence_1","Sequence_2","PSIcoffee")],pid_ref_pair[,c("Family","Sequence_1","Sequence_2","PSIcoffee")]))
msa_psi_pair = msa_psi_pair %>% 
  left_join(tcs_score_per_seq[,c("Family","Sequence","PSIcoffee")],by = c("Sequence_1" = "Sequence", "Family" = "Family")) %>%
  left_join(tcs_score_per_seq[,c("Family","Sequence","PSIcoffee")],by = c("Sequence_2" = "Sequence", "Family" = "Family")) %>%
  left_join(gdt_ts_complete_with_seqs[,c("Family","Sequence","alphafold")],by = c("Sequence_1" = "Sequence", "Family" = "Family")) %>%
  left_join(gdt_ts_complete_with_seqs[,c("Family","Sequence","alphafold")],by = c("Sequence_2" = "Sequence", "Family" = "Family")) %>% 
  left_join(alphafold_plddts,by = c("Sequence_1" = "Sequence", "Family" = "Family")) %>%
  left_join(alphafold_plddts,by = c("Sequence_2" = "Sequence", "Family" = "Family"))
colnames(msa_psi_pair) = c("Family","Sequence_1","Sequence_2","SoP","NiRMSD","PID","TCS_1","TCS_2","GDT_TS_1","GDT_TS_2","pLDDT_1","pLDDT_2")
msa_psi_pair$GDT_TS_1 = round(msa_psi_pair$GDT_TS_1*100,digits=2)
msa_psi_pair$GDT_TS_2 = round(msa_psi_pair$GDT_TS_2*100,digits=2)
msa_psi_pair$pLDDT_1 = round(msa_psi_pair$pLDDT_1,digits=2)
msa_psi_pair$pLDDT_2 = round(msa_psi_pair$pLDDT_2,digits=2)  
write.table(msa_psi_pair,"msa_psi_pair.tsv",quote=FALSE,row.names=FALSE,col.names=T,sep="\t")




msa_af2 = Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "Family", all.x = TRUE),
list(ref_3dcoffee_sp_same_state[,c("Family","3DCoffee_AF2")],tcs_score[,c("Family","3DCoffee_AF2")],nirmsd_ref[-13,c("Family","3DCoffee_AF2_REF_NAT")],pid_ref[,c("Family","3DCoffee_AF2")]))
colnames(msa_af2) = c("Family","SoP","TCS","NiRMSD","PID")
write.table(msa_af2,"msa_af2.tsv",quote=FALSE,row.names=FALSE,col.names=T,sep="\t")


msa_af2_avg = Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = c("Family","Sequence"), all.x = TRUE),
list(ref_3dcoffee_avg_sp_same_state[,c("Family","Sequence","3DCoffee_AF2")],tcs_score_per_seq[,c("Family","Sequence","3DCoffee_AF2")],nirmsd_ref_avg[,c("Family","Sequence","3DCoffee_AF2_REF_NAT")],pid_ref_avg[,c("Family","Sequence","3DCoffee_AF2")],gdt_ts_complete_with_seqs[,c("Family","Sequence","alphafold")],alphafold_plddts))
colnames(msa_af2_avg) = c("Family","Sequence","SoP","TCS","NiRMSD","PID","GDT_TS","pLDDT")
msa_af2_avg$GDT_TS = round(msa_af2_avg$GDT_TS*100,digits=2)
msa_af2_avg$pLDDT = round(msa_af2_avg$pLDDT,digits=2)
write.table(msa_af2_avg,"msa_af2_avg.tsv",quote=FALSE,row.names=FALSE,col.names=T,sep="\t")


msa_af2_pair = Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = c("Family","Sequence_1","Sequence_2"), all.x = TRUE),
list(ref_3dcoffee_pair_sp_same_state[,c("Family","Sequence_1","Sequence_2","3DCoffee_AF2")],tmp_nirmsd_ref_pair[,c("Family","Sequence_1","Sequence_2","3DCoffee_AF2_REF_NAT")],pid_ref_pair[,c("Family","Sequence_1","Sequence_2","3DCoffee_AF2")]))
msa_af2_pair = msa_af2_pair %>% 
  left_join(tcs_score_per_seq[,c("Family","Sequence","3DCoffee_AF2")],by = c("Sequence_1" = "Sequence", "Family" = "Family")) %>%
  left_join(tcs_score_per_seq[,c("Family","Sequence","3DCoffee_AF2")],by = c("Sequence_2" = "Sequence", "Family" = "Family")) %>%
  left_join(gdt_ts_complete_with_seqs[,c("Family","Sequence","alphafold")],by = c("Sequence_1" = "Sequence", "Family" = "Family")) %>%
  left_join(gdt_ts_complete_with_seqs[,c("Family","Sequence","alphafold")],by = c("Sequence_2" = "Sequence", "Family" = "Family")) %>% 
  left_join(alphafold_plddts,by = c("Sequence_1" = "Sequence", "Family" = "Family")) %>%
  left_join(alphafold_plddts,by = c("Sequence_2" = "Sequence", "Family" = "Family"))
colnames(msa_af2_pair) = c("Family","Sequence_1","Sequence_2","SoP","NiRMSD","PID","TCS_1","TCS_2","GDT_TS_1","GDT_TS_2","pLDDT_1","pLDDT_2")
msa_af2_pair$GDT_TS_1 = round(msa_af2_pair$GDT_TS_1*100,digits=2)
msa_af2_pair$GDT_TS_2 = round(msa_af2_pair$GDT_TS_2*100,digits=2)
msa_af2_pair$pLDDT_1 = round(msa_af2_pair$pLDDT_1,digits=2)
msa_af2_pair$pLDDT_2 = round(msa_af2_pair$pLDDT_2,digits=2)  
write.table(msa_af2_pair,"msa_af2_pair.tsv",quote=FALSE,row.names=FALSE,col.names=T,sep="\t")


msa_pdb = Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "Family", all.x = TRUE),
list(ref_3dcoffee_sp_same_state[,c("Family","3DCoffee_NAT")],tcs_score[,c("Family","3DCoffee_NAT")],nirmsd_ref[-13,c("Family","3DCoffee_NAT")],pid_ref[,c("Family","3DCoffee_NAT")]))
colnames(msa_pdb) = c("Family","SoP","TCS","NiRMSD","PID")
write.table(msa_pdb,"msa_pdb.tsv",quote=FALSE,row.names=FALSE,col.names=T,sep="\t")


msa_pdb_avg = Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = c("Family","Sequence"), all.x = TRUE),
list(ref_3dcoffee_avg_sp_same_state[,c("Family","Sequence","3DCoffee_NAT")],tcs_score_per_seq[,c("Family","Sequence","3DCoffee_NAT")],nirmsd_ref_avg[,c("Family","Sequence","3DCoffee_NAT")],pid_ref_avg[,c("Family","Sequence","3DCoffee_NAT")],gdt_ts_complete_with_seqs[,c("Family","Sequence","alphafold")],alphafold_plddts))
colnames(msa_pdb_avg) = c("Family","Sequence","SoP","TCS","NiRMSD","PID","GDT_TS","pLDDT")
msa_pdb_avg$GDT_TS = round(msa_pdb_avg$GDT_TS*100,digits=2)
msa_pdb_avg$pLDDT = round(msa_pdb_avg$pLDDT,digits=2)
write.table(msa_pdb_avg,"msa_pdb_avg.tsv",quote=FALSE,row.names=FALSE,col.names=T,sep="\t")


msa_pdb_pair = Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = c("Family","Sequence_1","Sequence_2"), all.x = TRUE),
list(ref_3dcoffee_pair_sp_same_state[,c("Family","Sequence_1","Sequence_2","3DCoffee_NAT")],tmp_nirmsd_ref_pair[,c("Family","Sequence_1","Sequence_2","3DCoffee_NAT")],pid_ref_pair[,c("Family","Sequence_1","Sequence_2","3DCoffee_NAT")]))
msa_pdb_pair = msa_pdb_pair %>% 
  left_join(tcs_score_per_seq[,c("Family","Sequence","3DCoffee_NAT")],by = c("Sequence_1" = "Sequence", "Family" = "Family")) %>%
  left_join(tcs_score_per_seq[,c("Family","Sequence","3DCoffee_NAT")],by = c("Sequence_2" = "Sequence", "Family" = "Family")) %>%
  left_join(gdt_ts_complete_with_seqs[,c("Family","Sequence","alphafold")],by = c("Sequence_1" = "Sequence", "Family" = "Family")) %>%
  left_join(gdt_ts_complete_with_seqs[,c("Family","Sequence","alphafold")],by = c("Sequence_2" = "Sequence", "Family" = "Family")) %>% 
  left_join(alphafold_plddts,by = c("Sequence_1" = "Sequence", "Family" = "Family")) %>%
  left_join(alphafold_plddts,by = c("Sequence_2" = "Sequence", "Family" = "Family"))
colnames(msa_pdb_pair) = c("Family","Sequence_1","Sequence_2","SoP","NiRMSD","PID","TCS_1","TCS_2","GDT_TS_1","GDT_TS_2","pLDDT_1","pLDDT_2")
msa_pdb_pair$GDT_TS_1 = round(msa_pdb_pair$GDT_TS_1*100,digits=2)
msa_pdb_pair$GDT_TS_2 = round(msa_pdb_pair$GDT_TS_2*100,digits=2)
msa_pdb_pair$pLDDT_1 = round(msa_pdb_pair$pLDDT_1,digits=2)
msa_pdb_pair$pLDDT_2 = round(msa_pdb_pair$pLDDT_2,digits=2)
write.table(msa_pdb_pair,"msa_pdb_pair.tsv",quote=FALSE,row.names=FALSE,col.names=T,sep="\t")





  #
  # Ref mtmalign
  #

if (FALSE) {
tmp_nirmsd_ref_pair = nirmsd_ref_pair
colnames(tmp_nirmsd_ref_pair)[1:2] = c("Sequence_2","Sequence_1")

msa_seq = Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "Family", all.x = TRUE),
        list(ref_mtmalign_sp[,c("Family","Ginsi")],tcs_score[,c("Family","Ginsi")],nirmsd_ref[-13,c("Family","Ginsi")],pid_ref[,c("Family","Ginsi")]))
colnames(msa_seq) = c("Family","SoP","TCS","NiRMSD","PID")
write.table(msa_seq,"msa_seq.tsv",quote=FALSE,row.names=FALSE,col.names=T,sep="\t")

msa_seq_avg = Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = c("Family","Sequence"), all.x = TRUE),
        list(ref_mtmalign_avg_sp[,c("Family","Sequence","Ginsi")],tcs_score_per_seq[,c("Family","Sequence","Ginsi")],nirmsd_ref_avg[,c("Family","Sequence","Ginsi")],pid_ref_avg[,c("Family","Sequence","Ginsi")],gdt_ts_complete_with_seqs[,c("Family","Sequence","alphafold")],alphafold_plddts))
colnames(msa_seq_avg) = c("Family","Sequence","SoP","TCS","NiRMSD","PID","GDT_TS","pLDDT")
msa_seq_avg$GDT_TS = round(msa_seq_avg$GDT_TS*100,digits=2)
msa_seq_avg$pLDDT = round(msa_seq_avg$pLDDT,digits=2)
write.table(msa_seq_avg,"msa_seq_avg.tsv",quote=FALSE,row.names=FALSE,col.names=T,sep="\t")


msa_seq_pair = Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = c("Family","Sequence_1","Sequence_2"), all.x = TRUE),
        list(ref_mtmalign_pair_sp[,c("Family","Sequence_1","Sequence_2","Ginsi")],tmp_nirmsd_ref_pair[,c("Family","Sequence_1","Sequence_2","Ginsi")],pid_ref_pair[,c("Family","Sequence_1","Sequence_2","Ginsi")]))
msa_seq_pair = msa_seq_pair %>% 
  left_join(tcs_score_per_seq[,c("Family","Sequence","Ginsi")],by = c("Sequence_1" = "Sequence", "Family" = "Family")) %>%
  left_join(tcs_score_per_seq[,c("Family","Sequence","Ginsi")],by = c("Sequence_2" = "Sequence", "Family" = "Family")) %>%
  left_join(gdt_ts_complete_with_seqs[,c("Family","Sequence","alphafold")],by = c("Sequence_1" = "Sequence", "Family" = "Family")) %>%
  left_join(gdt_ts_complete_with_seqs[,c("Family","Sequence","alphafold")],by = c("Sequence_2" = "Sequence", "Family" = "Family")) %>% 
  left_join(alphafold_plddts,by = c("Sequence_1" = "Sequence", "Family" = "Family")) %>%
  left_join(alphafold_plddts,by = c("Sequence_2" = "Sequence", "Family" = "Family"))
colnames(msa_seq_pair) = c("Family","Sequence_1","Sequence_2","SoP","NiRMSD","PID","TCS_1","TCS_2","GDT_TS_1","GDT_TS_2","pLDDT_1","pLDDT_2")
msa_seq_pair$GDT_TS_1 = round(msa_seq_pair$GDT_TS_1*100,digits=2)
msa_seq_pair$GDT_TS_2 = round(msa_seq_pair$GDT_TS_2*100,digits=2)
msa_seq_pair$pLDDT_1 = round(msa_seq_pair$pLDDT_1,digits=2)
msa_seq_pair$pLDDT_2 = round(msa_seq_pair$pLDDT_2,digits=2)  
write.table(msa_seq_pair,"msa_seq_pair.tsv",quote=FALSE,row.names=FALSE,col.names=T,sep="\t")



msa_af2 = Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "Family", all.x = TRUE),
        list(ref_mtmalign_sp[,c("Family","mTMalign_AF2")],tcs_score[,c("Family","mTMalign_AF2")],nirmsd_ref[-13,c("Family","mTMalign_AF2_REF_NAT")],pid_ref[,c("Family","mTMalign_AF2")]))
colnames(msa_af2) = c("Family","SoP","TCS","NiRMSD","PID")
write.table(msa_af2,"msa_af2.tsv",quote=FALSE,row.names=FALSE,col.names=T,sep="\t")


msa_af2_avg = Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = c("Family","Sequence"), all.x = TRUE),
        list(ref_mtmalign_avg_sp[,c("Family","Sequence","mTMalign_AF2")],tcs_score_per_seq[,c("Family","Sequence","mTMalign_AF2")],nirmsd_ref_avg[,c("Family","Sequence","mTMalign_AF2_REF_NAT")],pid_ref_avg[,c("Family","Sequence","mTMalign_AF2")],gdt_ts_complete_with_seqs[,c("Family","Sequence","alphafold")],alphafold_plddts))
colnames(msa_af2_avg) = c("Family","Sequence","SoP","TCS","NiRMSD","PID","GDT_TS","pLDDT")
msa_af2_avg$GDT_TS = round(msa_af2_avg$GDT_TS*100,digits=2)
msa_af2_avg$pLDDT = round(msa_af2_avg$pLDDT,digits=2)
write.table(msa_af2_avg,"msa_af2_avg.tsv",quote=FALSE,row.names=FALSE,col.names=T,sep="\t")


msa_af2_pair = Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = c("Family","Sequence_1","Sequence_2"), all.x = TRUE),
        list(ref_mtmalign_pair_sp[,c("Family","Sequence_1","Sequence_2","mTMalign_AF2")],tmp_nirmsd_ref_pair[,c("Family","Sequence_1","Sequence_2","mTMalign_AF2_REF_NAT")],pid_ref_pair[,c("Family","Sequence_1","Sequence_2","mTMalign_AF2")]))
msa_af2_pair = msa_af2_pair %>% 
  left_join(tcs_score_per_seq[,c("Family","Sequence","mTMalign_AF2")],by = c("Sequence_1" = "Sequence", "Family" = "Family")) %>%
  left_join(tcs_score_per_seq[,c("Family","Sequence","mTMalign_AF2")],by = c("Sequence_2" = "Sequence", "Family" = "Family")) %>%
  left_join(gdt_ts_complete_with_seqs[,c("Family","Sequence","alphafold")],by = c("Sequence_1" = "Sequence", "Family" = "Family")) %>%
  left_join(gdt_ts_complete_with_seqs[,c("Family","Sequence","alphafold")],by = c("Sequence_2" = "Sequence", "Family" = "Family")) %>% 
  left_join(alphafold_plddts,by = c("Sequence_1" = "Sequence", "Family" = "Family")) %>%
  left_join(alphafold_plddts,by = c("Sequence_2" = "Sequence", "Family" = "Family"))
colnames(msa_af2_pair) = c("Family","Sequence_1","Sequence_2","SoP","NiRMSD","PID","TCS_1","TCS_2","GDT_TS_1","GDT_TS_2","pLDDT_1","pLDDT_2")
msa_af2_pair$GDT_TS_1 = round(msa_af2_pair$GDT_TS_1*100,digits=2)
msa_af2_pair$GDT_TS_2 = round(msa_af2_pair$GDT_TS_2*100,digits=2)
msa_af2_pair$pLDDT_1 = round(msa_af2_pair$pLDDT_1,digits=2)
msa_af2_pair$pLDDT_2 = round(msa_af2_pair$pLDDT_2,digits=2)  
write.table(msa_af2_pair,"msa_af2_pair.tsv",quote=FALSE,row.names=FALSE,col.names=T,sep="\t")


msa_pdb = Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "Family", all.x = TRUE),
        list(ref_mtmalign_sp[,c("Family","mTMalign_NAT")],tcs_score[,c("Family","mTMalign_NAT")],nirmsd_ref[-13,c("Family","mTMalign_NAT")],pid_ref[,c("Family","mTMalign_NAT")]))
colnames(msa_pdb) = c("Family","SoP","TCS","NiRMSD","PID")
write.table(msa_pdb,"msa_pdb.tsv",quote=FALSE,row.names=FALSE,col.names=T,sep="\t")


msa_pdb_avg = Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = c("Family","Sequence"), all.x = TRUE),
        list(ref_mtmalign_avg_sp[,c("Family","Sequence","mTMalign_NAT")],tcs_score_per_seq[,c("Family","Sequence","mTMalign_NAT")],nirmsd_ref_avg[,c("Family","Sequence","mTMalign_NAT")],pid_ref_avg[,c("Family","Sequence","mTMalign_NAT")],gdt_ts_complete_with_seqs[,c("Family","Sequence","alphafold")],alphafold_plddts))
colnames(msa_pdb_avg) = c("Family","Sequence","SoP","TCS","NiRMSD","PID","GDT_TS","pLDDT")
msa_pdb_avg$GDT_TS = round(msa_pdb_avg$GDT_TS*100,digits=2)
msa_pdb_avg$pLDDT = round(msa_pdb_avg$pLDDT,digits=2)
write.table(msa_pdb_avg,"msa_pdb_avg.tsv",quote=FALSE,row.names=FALSE,col.names=T,sep="\t")


msa_pdb_pair = Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = c("Family","Sequence_1","Sequence_2"), all.x = TRUE),
        list(ref_mtmalign_pair_sp[,c("Family","Sequence_1","Sequence_2","mTMalign_NAT")],tmp_nirmsd_ref_pair[,c("Family","Sequence_1","Sequence_2","mTMalign_NAT")],pid_ref_pair[,c("Family","Sequence_1","Sequence_2","mTMalign_NAT")]))
msa_pdb_pair = msa_pdb_pair %>% 
  left_join(tcs_score_per_seq[,c("Family","Sequence","mTMalign_NAT")],by = c("Sequence_1" = "Sequence", "Family" = "Family")) %>%
  left_join(tcs_score_per_seq[,c("Family","Sequence","mTMalign_NAT")],by = c("Sequence_2" = "Sequence", "Family" = "Family")) %>%
  left_join(gdt_ts_complete_with_seqs[,c("Family","Sequence","alphafold")],by = c("Sequence_1" = "Sequence", "Family" = "Family")) %>%
  left_join(gdt_ts_complete_with_seqs[,c("Family","Sequence","alphafold")],by = c("Sequence_2" = "Sequence", "Family" = "Family")) %>% 
  left_join(alphafold_plddts,by = c("Sequence_1" = "Sequence", "Family" = "Family")) %>%
  left_join(alphafold_plddts,by = c("Sequence_2" = "Sequence", "Family" = "Family"))
colnames(msa_pdb_pair) = c("Family","Sequence_1","Sequence_2","SoP","NiRMSD","PID","TCS_1","TCS_2","GDT_TS_1","GDT_TS_2","pLDDT_1","pLDDT_2")
msa_pdb_pair$GDT_TS_1 = round(msa_pdb_pair$GDT_TS_1*100,digits=2)
msa_pdb_pair$GDT_TS_2 = round(msa_pdb_pair$GDT_TS_2*100,digits=2)
msa_pdb_pair$pLDDT_1 = round(msa_pdb_pair$pLDDT_1,digits=2)
msa_pdb_pair$pLDDT_2 = round(msa_pdb_pair$pLDDT_2,digits=2)
write.table(msa_pdb_pair,"msa_pdb_pair.tsv",quote=FALSE,row.names=FALSE,col.names=T,sep="\t")

}

#
# Supp Table 2
#

supp_table_2 = rbind(ref_3dcoffee_AF2_with_pairs_gm_gdt_below_75_sop_above_75_df,ref_3dcoffee_DMPFOLD_with_pairs_gm_gdt_below_75_sop_above_75_df)
write.table(supp_table_2,"Supp_table_2.tsv",quote=FALSE,row.names=TRUE,col.names=TRUE,sep="\t")



#
# Select the best aligners based on the lowest niRMSD scores 
#


if (FALSE) {
ref_3dcoffee_tmalign_sp = read.table("./selected_comparisons_ref_3dcoffee_TMalign_sp.txt",header = F,stringsAsFactors = F)
ref_3dcoffee_tmalign_sp = ref_3dcoffee_tmalign_sp[!duplicated(as.list(ref_3dcoffee_tmalign_sp))]
colnames(ref_3dcoffee_tmalign_sp)=c("Family","Famsa","Ginsi","MSAProbs","TCoffee","PSIcoffee","3DCoffee_NAT","3DCoffee_TMalign_NAT","mTMalign_NAT","3DCoffee_PRED_3_ITER","3DCoffee_TMalign_PRED_3_ITER","mTMalign_PRED_3_ITER","3DCoffee_PRED_4_ITER","3DCoffee_TMalign_PRED_4_ITER","mTMalign_PRED_4_ITER","3DCoffee_PRED_5_ITER","3DCoffee_TMalign_PRED_5_ITER","mTMalign_PRED_5_ITER") #,"Deepblast_vs_REF"
write.table(ref_3dcoffee_tmalign_sp, file="selected_comparisons_ref_3dcoffee_TMalign_sp.tsv",quote=FALSE, sep="\t",row.names=FALSE)

ref_3dcoffee_tmalign_avg_sp = read.table("./selected_comparisons_ref_3dcoffee_TMalign_avg_sp.txt",header = F,stringsAsFactors = F)
ref_3dcoffee_tmalign_avg_sp = ref_3dcoffee_tmalign_avg_sp[!duplicated(as.list(ref_3dcoffee_tmalign_avg_sp))]
colnames(ref_3dcoffee_tmalign_avg_sp) = c("Sequence","Family","Famsa","Ginsi","MSAProbs","TCoffee","PSIcoffee","3DCoffee_NAT","3DCoffee_TMalign_NAT","mTMalign_NAT","3DCoffee_PRED_3_ITER","3DCoffee_TMalign_PRED_3_ITER","mTMalign_PRED_3_ITER","3DCoffee_PRED_4_ITER","3DCoffee_TMalign_PRED_4_ITER","mTMalign_PRED_4_ITER","3DCoffee_PRED_5_ITER","3DCoffee_TMalign_PRED_5_ITER","mTMalign_PRED_5_ITER")
write.table(ref_3dcoffee_tmalign_avg_sp, file="selected_comparisons_ref_3dcoffee_TMalign_avg_sp.tsv",quote=FALSE, sep="\t",row.names=FALSE)

ref_3dcoffee_tmalign_pair_sp = read.table("./selected_comparisons_ref_3dcoffee_TMalign_pair_sp.txt",header = F,stringsAsFactors = F)
ref_3dcoffee_tmalign_pair_sp = ref_3dcoffee_tmalign_pair_sp[!duplicated(as.list(ref_3dcoffee_tmalign_pair_sp))]
colnames(ref_3dcoffee_tmalign_pair_sp)=c("Sequence_1","Sequence_2","Family","Famsa","Ginsi","MSAProbs","TCoffee","PSIcoffee","3DCoffee_NAT","3DCoffee_TMalign_NAT","mTMalign_NAT","3DCoffee_PRED_3_ITER","3DCoffee_TMalign_PRED_3_ITER","mTMalign_PRED_3_ITER","3DCoffee_PRED_4_ITER","3DCoffee_TMalign_PRED_4_ITER","mTMalign_PRED_4_ITER","3DCoffee_PRED_5_ITER","3DCoffee_TMalign_PRED_5_ITER","mTMalign_PRED_5_ITER")


for(fam in best_aln$Family){
  my_vec = which(best_aln$Family == fam)
  my_best_nat = as.character(best_aln[["best_nat_aln"]][my_vec])
  my_best_pred = as.character(best_aln[["best_pred_aln"]][my_vec])
  my_best_seq = as.character(best_aln[["best_seq_aln"]][my_vec])

  if (my_best_nat == "3DCoffee_NAT") {
    tmp_df = ref_3dcoffee_sp[my_vec,c(1,which(colnames(ref_3dcoffee_sp) == my_best_seq),6,which(colnames(ref_3dcoffee_sp) == my_best_pred))]
    colnames(tmp_df)[4] = "Best_PRED_3_ITER"
    colnames(tmp_df)[2] = "Best_SEQ"

    my_subset = which(ref_3dcoffee_avg_sp$Family == fam)
    tmp_avg_df = ref_3dcoffee_avg_sp[my_subset,c(1,2,which(colnames(ref_3dcoffee_avg_sp) == my_best_seq),7,which(colnames(ref_3dcoffee_avg_sp) == my_best_pred))]
    colnames(tmp_avg_df)[5] = "Best_PRED_3_ITER"
    colnames(tmp_avg_df)[3] = "Best_SEQ"

    my_subset_pair = which(ref_3dcoffee_pair_sp$Family == fam)
    tmp_pair_df = ref_3dcoffee_pair_sp[my_subset_pair,c(1:3,which(colnames(ref_3dcoffee_pair_sp) == my_best_seq),8,which(colnames(ref_3dcoffee_pair_sp) == my_best_pred))]
    colnames(tmp_pair_df)[6] = "Best_PRED_3_ITER"
    colnames(tmp_pair_df)[4] = "Best_SEQ"

    if (my_vec == 1){
      new_df = tmp_df
      new_avg_df = tmp_avg_df
      new_pair_df = tmp_pair_df
    } else {
      new_df = rbind(new_df,tmp_df)
      new_avg_df = rbind(new_avg_df,tmp_avg_df)
      new_pair_df = rbind(new_pair_df,tmp_pair_df)
    }
  } else if (my_best_nat == "3DCoffee_TMalign_NAT") {
    tmp_df = ref_3dcoffee_tmalign_sp[my_vec,c(1,which(colnames(ref_3dcoffee_tmalign_sp) == my_best_seq),6,which(colnames(ref_3dcoffee_tmalign_sp) == my_best_pred))]
    colnames(tmp_df)[4] = "Best_PRED_3_ITER"
    colnames(tmp_df)[2] = "Best_SEQ"

    my_subset = which(ref_3dcoffee_tmalign_avg_sp$Family == fam)
    tmp_avg_df = ref_3dcoffee_tmalign_avg_sp[my_subset,c(1,2,which(colnames(ref_3dcoffee_tmalign_avg_sp) == my_best_seq),7,which(colnames(ref_3dcoffee_tmalign_avg_sp) == my_best_pred))]
    colnames(tmp_avg_df)[5] = "Best_PRED_3_ITER"
    colnames(tmp_avg_df)[3] = "Best_SEQ"

    my_subset_pair = which(ref_3dcoffee_tmalign_pair_sp$Family == fam)
    tmp_pair_df = ref_3dcoffee_tmalign_pair_sp[my_subset_pair,c(1:3,which(colnames(ref_3dcoffee_tmalign_pair_sp) == my_best_seq),8,which(colnames(ref_3dcoffee_tmalign_pair_sp) == my_best_pred))]
    colnames(tmp_pair_df)[6] = "Best_PRED_3_ITER"
    colnames(tmp_pair_df)[4] = "Best_SEQ"

    if (my_vec == 1){
      new_df = tmp_df
      new_avg_df = tmp_avg_df
      new_pair_df = tmp_pair_df
    } else {
      new_df = rbind(new_df,tmp_df)
      new_avg_df = rbind(new_avg_df,tmp_avg_df)
      new_pair_df = rbind(new_pair_df,tmp_pair_df)      
    }
  } else {
    tmp_df = ref_mtmalign_sp[my_vec,c(1,which(colnames(ref_mtmalign_sp) == my_best_seq),6,which(colnames(ref_mtmalign_sp) == my_best_pred))]
    colnames(tmp_df)[4] = "Best_PRED_3_ITER"
    colnames(tmp_df)[2] = "Best_SEQ"

    my_subset = which(ref_mtmalign_avg_sp$Family == fam)
    tmp_avg_df = ref_mtmalign_avg_sp[my_subset,c(1,2,which(colnames(ref_mtmalign_avg_sp) == my_best_seq),7,which(colnames(ref_mtmalign_avg_sp) == my_best_pred))]
    colnames(tmp_avg_df)[5] = "Best_PRED_3_ITER"
    colnames(tmp_avg_df)[3] = "Best_SEQ"

    my_subset_pair = which(ref_mtmalign_pair_sp$Family == fam)
    tmp_pair_df = ref_mtmalign_pair_sp[my_subset_pair,c(1:3,which(colnames(ref_mtmalign_pair_sp) == my_best_seq),8,which(colnames(ref_mtmalign_pair_sp) == my_best_pred))]
    colnames(tmp_pair_df)[6] = "Best_PRED_3_ITER"
    colnames(tmp_pair_df)[4] = "Best_SEQ"

    if (my_vec == 1){
      new_df = tmp_df
      new_avg_df = tmp_avg_df
      new_pair_df = tmp_pair_df
    } else {
      new_df = rbind(new_df,tmp_df)
      new_avg_df = rbind(new_avg_df,tmp_avg_df)
      new_pair_df = rbind(new_pair_df,tmp_pair_df)     
    }
  }

}

write.table(new_df, file="selected_comparisons_best_aln_sp.tsv",quote=FALSE, sep="\t",row.names=FALSE)
write.table(new_df, file="selected_comparisons_best_aln_avg_sp.tsv",quote=FALSE, sep="\t",row.names=FALSE)
ο=ggplot(new_pair_df,aes(x=Best_PRED_3_ITER,y=PSIcoffee,color=Family))+geom_point()+xlim(0,100)+ylim(0,100)+geom_abline(intercept =0 , slope = 1,lwd=0.2)+xlab("Sum-of-Pairs score between pairs of sequences - Best_PRED vs Best_NAT (%)")+ylab("Sum-of-Pairs score between pairs of sequences - PSICoffee vs Best_NAT (%)")+ theme_minimal() + scale_color_manual(values=color_codes) #+ggtitle("Pairwise SP scores")
ggsave("best_pred_pair_sp_score_per_family_3_iter_vs_psicoffee.png",dpi = "retina")
ο=ggplot(new_pair_df,aes(x=Best_PRED_3_ITER,y=Best_SEQ,color=Family))+geom_point()+xlim(0,100)+ylim(0,100)+geom_abline(intercept =0 , slope = 1,lwd=0.2)+xlab("Sum-of-Pairs score between pairs of sequences - Best_PRED vs Best_NAT (%)")+ylab("Sum-of-Pairs score between pairs of sequences - Best_SEQ vs Best_NAT (%)")+ theme_minimal() + scale_color_manual(values=color_codes) #+ggtitle("Pairwise SP scores")
ggsave("best_pred_pair_sp_score_per_family_3_iter_vs_best_seq.png",dpi = "retina")
}





if (FALSE) {
  nirmsd_ref = read.table("./selected_comparisons_nirmsd.txt",header = F,stringsAsFactors = F)
colnames(nirmsd_ref)=c("Family","Famsa","Ginsi","MSAProbs","TCoffee","PSICoffee","3DCoffee_NAT","3DCoffee_TMalign_NAT","mTMalign_NAT","3DCoffee_PRED_3_ITER","3DCoffee_TMalign_PRED_3_ITER","mTMalign_PRED_3_ITER","3DCoffee_PRED_4_ITER","3DCoffee_TMalign_PRED_4_ITER","mTMalign_PRED_4_ITER","3DCoffee_PRED_5_ITER","3DCoffee_TMalign_PRED_5_ITER","mTMalign_PRED_5_ITER") # ,"Deepblast"
nirmsd_ref = rbind(nirmsd_ref,c("Average",colMeans(nirmsd_ref[,-1])))
row.names(nirmsd_ref) = nirmsd_ref$Family
nirmsd_ref[colnames(nirmsd_ref)[-1]] <- sapply(nirmsd_ref[colnames(nirmsd_ref)[-1]],as.numeric)
pheatmap(nirmsd_ref[,-c(1,13:length(nirmsd_ref))],display_numbers = T,angle_col = 45,fontsize = 10,cluster_cols = F,cluster_rows = F,number_color="black",fontsize_number=12,filename = "nirmsd_ref_NAT.png")

  test=apply(nirmsd_ref[-nrow(nirmsd_ref),7:9], 1, function(x) which(x == min(x))[1])
best_nat_aln = colnames(nirmsd_ref[,7:9])[test]

test=apply(nirmsd_ref[-nrow(nirmsd_ref),10:12], 1, function(x) which(x == min(x))[1])
best_pred_aln = colnames(nirmsd_ref[,10:12])[test]

test=apply(nirmsd_ref[-nrow(nirmsd_ref),2:5], 1, function(x) which(x == min(x))[1])
best_seq_aln = colnames(nirmsd_ref[,2:5])[test]

best_aln = data.frame(nirmsd_ref[-nrow(nirmsd_ref),c(1,6)],best_seq_aln,best_nat_aln,best_pred_aln)
colnames(best_aln)[1] = "Family"
}

#
# Analysis based on TCS-cutoff with best aligners
#

if (FALSE){
  for(fam in best_aln$Family){
  my_vec = which(best_aln$Family == fam)
  my_best_nat = as.character(best_aln[["best_nat_aln"]][my_vec])
  my_best_pred = as.character(best_aln[["best_pred_aln"]][my_vec])
  my_best_seq = as.character(best_aln[["best_seq_aln"]][my_vec])

  test_aln = tolower(my_best_seq)
  tcs_ref = gsub("_tmalign","_TMalign",tolower(gsub("_PRED_3_ITER","_dmpfold",my_best_pred)))
  ref_aln = gsub("_tmalign","_TMalign",tolower(gsub("_NAT","",my_best_nat)))

  myfilename = paste0("all_tcs_cutoffs_all_fams_selected_ref_",test_aln,"_vs_ref_",ref_aln,".sp_ref_",tcs_ref)
  all_tcs_cutoffs_all_fams_best_seq = read.table(myfilename,header = F,stringsAsFactors = F)
  all_tcs_cutoffs_all_fams_best_seq = all_tcs_cutoffs_all_fams_best_seq[my_vec,]
  myfilename = paste0("all_tcs_cutoffs_all_fams_selected_ref_psicoffee_vs_ref_",ref_aln,".sp_ref_",tcs_ref)
  all_tcs_cutoffs_all_fams_psicoffee = read.table(myfilename,header = F,stringsAsFactors = F)
  all_tcs_cutoffs_all_fams_psicoffee = all_tcs_cutoffs_all_fams_psicoffee[my_vec,]
  myfilename = paste0("all_tcs_cutoffs_all_fams_selected_ref_",tcs_ref,"_vs_ref_",ref_aln,".sp_ref_",tcs_ref)
  all_tcs_cutoffs_all_fams_best_pred = read.table(myfilename,header = F,stringsAsFactors = F)
  all_tcs_cutoffs_all_fams_best_pred = all_tcs_cutoffs_all_fams_best_pred[my_vec,]
  all_tcs_cutoffs_all_fams_collect = cbind(all_tcs_cutoffs_all_fams_best_seq,all_tcs_cutoffs_all_fams_psicoffee,all_tcs_cutoffs_all_fams_best_pred)
  colnames(all_tcs_cutoffs_all_fams_collect) = c("Best_SEQ TCS > 500","Best_SEQ TCS > 550","Best_SEQ TCS > 600","Best_SEQ TCS > 650","Best_SEQ TCS > 700","Best_SEQ TCS > 750","PSICoffee TCS > 500","PSICoffee TCS > 550","PSICoffee TCS > 600","PSICoffee TCS > 650","PSICoffee TCS > 700","PSICoffee TCS > 750","Best_PRED TCS > 500","Best_PRED TCS > 550","Best_PRED TCS > 600","Best_PRED TCS > 650","Best_PRED TCS > 700","Best_PRED TCS > 750")

  myfilename = paste0("all_tcs_cutoffs_all_fams_selected.list_ref_",tcs_ref)
  tcs_cutoff_list = read.table(myfilename,header = F,stringsAsFactors = F)
  tcs_cutoff_list = tcs_cutoff_list[my_vec,]
  colnames(tcs_cutoff_list) = c("Family","TCS > 0","TCS > 500","TCS > 550","TCS > 600","TCS > 650","TCS > 700","TCS > 750")

  myfilename = paste0(fam,"/",fam,"_selected_ref_",test_aln,"_vs_ref_",ref_aln,".sp.pair_tcs_above_75_ref_",tcs_ref)
  if (file.info(myfilename)$size != 1) {
    best_seq_pair_tcs_above_75 = read.table(myfilename,header = F,stringsAsFactors = F)
    myfilename = paste0(fam,"/",fam,"_selected_ref_psicoffee_vs_ref_",ref_aln,".sp.pair_tcs_above_75_ref_",tcs_ref)
    psicoffee_pair_tcs_above_75 = read.table(myfilename,header = F,stringsAsFactors = F)
    myfilename = paste0(fam,"/",fam,"_selected_ref_",tcs_ref,"_vs_ref_",ref_aln,".sp.pair_tcs_above_75_ref_",tcs_ref)
    best_pred_pair_tcs_above_75 = read.table(myfilename,header = F,stringsAsFactors = F)
    best_aln_tcs_above_75 = cbind(best_seq_pair_tcs_above_75,psicoffee_pair_tcs_above_75,best_pred_pair_tcs_above_75)
    best_aln_tcs_above_75 = best_aln_tcs_above_75[!duplicated(as.list(best_aln_tcs_above_75))]
    colnames(best_aln_tcs_above_75) = c("Sequence_1","Sequence_2","Family","Best_SEQ","PSIcoffee","Best_PRED")
  
    if (my_vec == 1){
    best_aln_tcs_above_75_df = best_aln_tcs_above_75
    } else {
    best_aln_tcs_above_75_df = rbind(best_aln_tcs_above_75_df,best_aln_tcs_above_75)
    }
  }
  if (my_vec == 1){
    tcs_df = all_tcs_cutoffs_all_fams_collect
    tcs_list = tcs_cutoff_list
  } else {
    tcs_df = rbind(tcs_df,all_tcs_cutoffs_all_fams_collect)
    tcs_list = rbind(tcs_list,tcs_cutoff_list)
  }
}

weights = sapply(tcs_list[,-c(1)], function(x) x/sum(x))
tcs_df$`Best_SEQ TCS > 0` = new_df$Best_SEQ
tcs_df$`PSIcoffee TCS > 0` = new_df$PSIcoffee
tcs_df$`Best_PRED TCS > 0` = new_df$`Best_PRED_3_ITER`
tcs_df = tcs_df[,c(19,1:6,20,7:12,21,13:18)]

sp_mean_seq = c()
sp_sd_seq = c()
sp_mean_psicoffee = c()
sp_sd_psicoffee = c()
sp_mean = c()
sp_sd = c()

for (mod in 1:ncol(weights)) {
  sp_mean_seq = c(sp_mean_seq,wt.mean(x=tcs_df[,mod],wt=weights[,mod]))
  sp_sd_seq = c(sp_sd_seq,wt.sd(x=tcs_df[,mod],wt=weights[,mod]))
  
  sp_mean_psicoffee = c(sp_mean_psicoffee,wt.mean(x=tcs_df[,mod+7],wt=weights[,mod]))
  sp_sd_psicoffee = c(sp_sd_psicoffee,wt.sd(x=tcs_df[,mod+7],wt=weights[,mod]))
  
  sp_mean = c(sp_mean,wt.mean(x=tcs_df[,mod+14],wt=weights[,mod]))
  sp_sd = c(sp_sd,wt.sd(x=tcs_df[,mod+14],wt=weights[,mod]))
}

perc_of_seq_mean = sapply(tcs_list[,-c(1)], function(x) (sum(x)/sum(tcs_list[,2]))*100)

tcs_list_and_sp_means_and_sds = data.frame(sp_mean_seq,sp_sd_seq,sp_mean_psicoffee,sp_sd_psicoffee,sp_mean, sp_sd, perc_of_seq_mean)
row.names(tcs_list_and_sp_means_and_sds) = colnames(tcs_list)[-1]
tcs_list_and_sp_means_and_sds$TCS_cutoff = row.names(tcs_list_and_sp_means_and_sds)

x = ggplot(tcs_list_and_sp_means_and_sds, aes(x=TCS_cutoff,group=1)) + scale_x_discrete(limits=tcs_list_and_sp_means_and_sds$TCS_cutoff) + geom_ribbon(aes(ymin=sp_mean_seq-sp_sd_seq, ymax=sp_mean_seq+sp_sd_seq),fill="lightgreen",alpha=0.3) + geom_ribbon(aes(ymin=sp_mean_psicoffee-sp_sd_psicoffee, ymax=sp_mean_psicoffee+sp_sd_psicoffee),fill="darkorange",alpha=0.3) + geom_ribbon(aes(ymin=sp_mean-sp_sd, ymax=sp_mean+sp_sd),fill="lightblue",alpha=0.3) + geom_line(aes(y=sp_mean_seq,col="sp_mean_seq"))+ geom_line(aes(y=sp_mean_psicoffee,col="sp_mean_psicoffee")) + geom_line(aes(y=sp_mean,col="sp_mean")) + geom_line(aes(y=perc_of_seq_mean,col="perc_of_seq_mean"),linetype="dashed")+ scale_color_manual("", values = c("sp_mean_seq" = "darkolivegreen", "sp_mean_psicoffee" = "darkorange4", "sp_mean" = "blue", "perc_of_seq_mean" = "red"),labels = c("sp_mean_seq" = "Best_SEQ vs Best_NAT","sp_mean_psicoffee" = "PSICoffee vs Best_NAT", "sp_mean" = "Best_PRED vs Best_NAT", "perc_of_seq_mean" = "Percent of Sequences")) + scale_y_continuous( name = "Sum-of-Pairs score (%)", sec.axis = sec_axis(~.,name="Percent of sequences included (%)"),limits=c(0,100)) + xlab("TCS cutoff") + theme_minimal() + theme(axis.text.x=element_text(angle = 45, hjust = 1),legend.position="top") + guides(colour = guide_legend(nrow = 2))
ggsave("best_aln_sp_score_and_perc_of_seqs_3_iter_vs_tcs_cutoff.png",dpi="retina")

tcs_list$ref = tcs_list[,2]
for (mode in 2:ncol(tcs_list)) {
  tcs_list[,mode] = (tcs_list[,mode]/tcs_list$ref)*100
}

tcs_list = tcs_list[,-c(9)]
tcs_list_melted = melt(tcs_list,variable.name="TCS_mode",value.name="percent_of_seq")
tcs_df_melted = melt(tcs_df[,-c(8:14)],variable.name="TCS_mode",value.name="SP")
tcs_list_and_df = cbind(tcs_list_melted,tcs_df_melted)
colnames(tcs_list_and_df) = make.unique(colnames(tcs_list_and_df))
tcs_list_and_df$group = c(rep("SoP Best_SEQ vs Best_NAT",nrow(tcs_list_and_df)/2),rep("SoP Best_PRED vs Best_NAT",nrow(tcs_list_and_df)/2))
tcs_list_and_df = tcs_list_and_df[order(tcs_list_and_df$Family,tcs_list_and_df$TCS_mode),]
circular_plot(tcs_list_and_df,"SoP Best_SEQ vs Best_NAT","SoP Best_PRED vs Best_NAT","best_aln_sp_score_per_family_3_iter_with_tcs_cutoff_circular.png")

selected_fams = levels(as.factor(best_aln_tcs_above_75_df$Family))
selected_colors = color_per_fam[color_per_fam$Family %in% selected_fams,c(2)]

ο=ggplot(best_aln_tcs_above_75_df,aes(x=`Best_PRED`,y=PSIcoffee,color=Family))+geom_point()+xlim(0,100)+ylim(0,100)+geom_abline(intercept =0 , slope = 1,lwd=0.2)+xlab("Sum-of-Pairs score between pairs of sequences - Best_PRED vs Best_NAT (%)")+ylab("Sum-of-Pairs score between pairs of sequences - PSICoffee vs Best_NAT (%)") + theme_minimal() + scale_color_manual(values=as.character(selected_colors)) #+ggtitle("Sum-of-Pairs score between pairs of sequences with TCS > 750")
ggsave("best_pred_pair_sp_score_per_family_3_iter_vs_psicoffee_tcs_above_75.png",dpi = "retina")
ο=ggplot(best_aln_tcs_above_75_df,aes(x=`Best_PRED`,y=Best_SEQ,color=Family))+geom_point()+xlim(0,100)+ylim(0,100)+geom_abline(intercept =0 , slope = 1,lwd=0.2)+xlab("Sum-of-Pairs score between pairs of sequences - Best_PRED vs Best_NAT (%)")+ylab("Sum-of-Pairs score between pairs of sequences - Best_SEQ vs Best_NAT (%)") + theme_minimal() + scale_color_manual(values=as.character(selected_colors)) #+ggtitle("Sum-of-Pairs score between pairs of sequences with TCS > 750")
ggsave("best_pred_pair_sp_score_per_family_3_iter_vs_best_seq_tcs_above_75.png",dpi = "retina")
}


#mtmalign_perc_of_seq = c()
#for (fam in levels(as.factor(tcs_score_per_seq$Family))) {
#    subset = which(tcs_score_per_seq$Family == fam)
#    tot_num = length(subset)
#    num = length(which(tcs_score_per_seq[subset,]$mTMalign_3_ITER > 700))
#    mtmalign_perc_of_seq = c(mtmalign_perc_of_seq,round( (num/tot_num*100), digits=2))
#}

#p_3dcoffee_perc_of_seq = c()
#for (fam in levels(as.factor(tcs_score_per_seq$Family))) {
#    subset = which(tcs_score_per_seq$Family == fam)
#    tot_num = length(subset)
#    num = length(which(tcs_score_per_seq[subset,]$'3dcoffee_3_ITER' > 700))
#    p_3dcoffee_perc_of_seq = c(p_3dcoffee_perc_of_seq,round( (num/tot_num*100), digits=2))
#}

#mydf_ref_3dcoffee_with_tcs_score = merge(mydf_ref_3dcoffee,tcs_score,by="Family")
#mydf_ref_3dcoffee_with_tcs_score$perc_of_seq = p_3dcoffee_perc_of_seq
#o=ggplot(mydf_ref_3dcoffee_with_tcs_score,aes(x=`3DCoffee_3_ITER_vs_REF`,y=Psicoffee_vs_REF,color=Family))+geom_point()+xlim(0,100)+ylim(0,100)+geom_abline(intercept =0 , slope = 1, lwd=0.2)+xlab("3DCoffee_3_ITER_vs_REF")+ylab("Psicoffee_vs_REF")+ggtitle("SP score per family")+labs(size = "TCS 3DCoffee_3_ITER")
#ο=ggplot(mydf_ref_3dcoffee,aes(x=`3DCoffee_3_ITER_vs_REF`,y=Psicoffee_vs_REF,color=Family))+geom_point()+xlim(0,100)+ylim(0,100)+geom_abline(intercept =0 , slope = 1,lwd=0.2)+xlab("3DCoffee_3_ITER_vs_REF")+ylab("Psicoffee_vs_REF")+ggtitle("SP score per family")
#ggsave("ref_3dcoffee_sp_score_per_family_3_iter_vs_psicoffee.png",dpi = "retina")
#ο=ggplot(mydf_ref_3dcoffee_with_tcs_score,aes(x=`3DCoffee_3_ITER_vs_REF`,y=MSAprobs_vs_REF,color=Family))+geom_point()+xlim(0,100)+ylim(0,100)+geom_abline(intercept =0 , slope = 1,lwd=0.2)+xlab("3DCoffee_3_ITER_vs_REF")+ylab("MSAprobs_vs_REF")+ggtitle("SP score per family")+labs(size = "TCS 3DCoffee_3_ITER")
#ggsave("ref_3dcoffee_sp_score_per_family_3_iter_vs_msaprobs.png",dpi = "retina")

#k = ggplot(mydf_ref_3dcoffee_with_tcs_score,aes(x=`3DCoffee_3_ITER_vs_REF`,y=`3DCoffee_3_ITER`,color=Family,size=perc_of_seq))+geom_point()+xlim(0,100)+geom_abline(intercept =0 , slope = 1,lwd=0.2)+xlab("3DCoffee_3_ITER_vs_REF")+ylab("TCS 3DCoffee_3_ITER")+ggtitle("SP score vs TCS per family - TCS seq cutoff 700")+labs(size = "Percent of sequences included")
#ggsave("ref_3dcoffee_sp_score_per_family_3_iter_vs_tcs_score_3_iter.png",dpi = "retina")



#mydf_ref_mtmalign_with_tcs_score = merge(mydf_ref_mtmalign,tcs_score,by="Family")

#mydf_ref_mtmalign_with_tcs_score$perc_of_seq = mtmalign_perc_of_seq

#ο=ggplot(mydf_ref_mtmalign_with_tcs_score,aes(x=`mTMalign_3_ITER_vs_REF`,y=Psicoffee_vs_REF,color=Family))+geom_point()+xlim(0,100)+ylim(0,100)+geom_abline(intercept =0 , slope = 1,lwd=0.2)+xlab("mTMalign_3_ITER_vs_REF")+ylab("Psicoffee_vs_REF")+ggtitle("SP score per family")+labs(size = "TCS mTMalign_3_ITER")
#ggsave("ref_mtmalign_sp_score_per_family_3_iter_vs_psicoffee.png",dpi = "retina")
#ο=ggplot(mydf_ref_mtmalign_with_tcs_score,aes(x=`mTMalign_3_ITER_vs_REF`,y=MSAprobs_vs_REF,color=Family))+geom_point()+xlim(0,100)+ylim(0,100)+geom_abline(intercept =0 , slope = 1,lwd=0.2)+xlab("mTMalign_3_ITER_vs_REF")+ylab("MSAprobs_vs_REF")+ggtitle("SP score per family")+labs(size = "TCS mTMalign_3_ITER")
#ggsave("ref_mtmalign_sp_score_per_family_3_iter_vs_msaprobs.png",dpi = "retina")

#k=ggplot(mydf_ref_mtmalign_with_tcs_score,aes(x=`mTMalign_3_ITER_vs_REF`,y=`mTMalign_3_ITER`,color=Family,size=perc_of_seq))+geom_point()+xlim(0,100)+geom_abline(intercept =0 , slope = 1,lwd=0.2)+xlab("mTMalign_3_ITER_vs_REF")+ylab("TCS mTMalign_3_ITER")+ggtitle("SP score vs TCS per family - TCS seq cutoff 700")+labs(size = "Percent of sequences included")
#ggsave("ref_mtmalign_sp_score_per_family_3_iter_vs_tcs_score_3_iter.png",dpi = "retina")




#mydf_ref_3dcoffee_melted = melt(mydf_ref_3dcoffee,measure.vars = 2:length(colnames(mydf_ref_3dcoffee)))
#m=ggplot(mydf_ref_3dcoffee_melted, aes(x=variable,y=value,fill=variable)) + geom_boxplot() + stat_summary(fun.y=mean, geom="point", shape=18, size=4)+ylab("SP score")+xlab("Modes")+theme(axis.text.x=element_text(angle = 45, hjust = 1))
#ggsave("ref_3dcoffee_sp_score_per_mode.png")
#x=ggplot(mydf_ref_3dcoffee_melted,aes(x=Family,y=value,fill=variable))+geom_bar(position = "dodge",stat = "identity",width = 0.87)+theme(axis.text.x=element_text(angle = 45, hjust = 1),legend.text = element_text(size = 6))+ylab("SP score")+labs(fill='Comparisons')
#ggsave("ref_3dcoffee_sp_score_per_family.png",dpi = "retina")

#mydf_ref_mtmalign_melted = melt(mydf_ref_mtmalign,measure.vars = 2:length(colnames(mydf_ref_mtmalign)))
#m=ggplot(mydf_ref_mtmalign_melted, aes(x=variable,y=value,fill=variable)) + geom_boxplot() + stat_summary(fun.y=mean, geom="point", shape=18, size=4)+ylab("SP score")+xlab("Modes")+theme(axis.text.x=element_text(angle = 45, hjust = 1))
#ggsave("ref_mtmalign_sp_score_per_mode.png")
#x=ggplot(mydf_ref_mtmalign_melted,aes(x=Family,y=value,fill=variable))+geom_bar(position = "dodge",stat = "identity",width = 0.87)+theme(axis.text.x=element_text(angle = 45, hjust = 1),legend.text = element_text(size = 6))+ylab("SP score")+labs(fill='Comparisons')
#ggsave("ref_mtmalign_sp_score_per_family.png",dpi = "retina")

#mydf_ref_3dcoffee_TMalign_melted = melt(mydf_ref_3dcoffee_TMalign,measure.vars = 2:length(colnames(mydf_ref_3dcoffee_TMalign)))
#m=ggplot(mydf_ref_3dcoffee_TMalign_melted, aes(x=variable,y=value,fill=variable)) + geom_boxplot() + stat_summary(fun.y=mean, geom="point", shape=18, size=4)+ylab("SP score")+xlab("Modes")+theme(axis.text.x=element_text(angle = 45, hjust = 1))
#ggsave("ref_3dcoffee_TMalign_sp_score_per_mode.png")
#x=ggplot(mydf_ref_3dcoffee_TMalign_melted,aes(x=Family,y=value,fill=variable))+geom_bar(position = "dodge",stat = "identity",width = 0.87)+theme(axis.text.x=element_text(angle = 45, hjust = 1),legend.text = element_text(size = 6))+ylab("SP score")+labs(fill='Comparisons')
#ggsave("ref_3dcoffee_TMalign_sp_score_per_family.png",dpi = "retina")

#ο=ggplot(ref_3dcoffee_avg_sp_with_tcs_score_per_seq,aes(x=`3DCoffee_3_ITER_vs_REF`,y=Psicoffee_vs_REF,color=Family))+geom_point()+xlim(0,100)+ylim(0,100)+geom_abline(intercept =0 , slope = 1,lwd=0.2)+xlab("3DCoffee_3_ITER_vs_REF")+ylab("Psicoffee_vs_REF")+ggtitle("Average SP scores per sequence")+labs(size = "TCS 3DCoffee_3_ITER")
#ggsave("ref_3dcoffee_avg_sp_score_per_seq_per_family_3_iter_vs_psicoffee.png",dpi = "retina")
#ο=ggplot(ref_3dcoffee_avg_sp_with_tcs_score_per_seq,aes(x=`3DCoffee_3_ITER_vs_REF`,y=MSAprobs_vs_REF,color=Family))+geom_point()+xlim(0,100)+ylim(0,100)+geom_abline(intercept =0 , slope = 1,lwd=0.2)+xlab("3DCoffee_3_ITER_vs_REF")+ylab("MSAprobs_vs_REF")+ggtitle("Average SP scores per sequence")+labs(size = "TCS 3DCoffee_3_ITER")
#ggsave("ref_3dcoffee_avg_sp_score_per_seq_per_family_3_iter_vs_msaprobs.png",dpi = "retina")

#per_included = round(length(ref_3dcoffee_avg_sp_with_tcs_score_per_seq[which(ref_3dcoffee_avg_sp_with_tcs_score_per_seq$`3DCoffee_3_ITER` > 700),]$`3DCoffee_3_ITER`)/length(ref_3dcoffee_avg_sp_with_tcs_score_per_seq$`3DCoffee_3_ITER`)*100,digits = 2)
#ο=ggplot(ref_3dcoffee_avg_sp_with_tcs_score_per_seq[which(ref_3dcoffee_avg_sp_with_tcs_score_per_seq$`3DCoffee_3_ITER` > 700),],aes(x=`3DCoffee_3_ITER_vs_REF`,y=`3DCoffee_3_ITER`,color=Family))+geom_point()+xlim(0,100)+xlab("SP 3DCoffee_3_ITER_vs_REF")+ylab("TCS 3DCoffee_3_ITER")+ggtitle(paste0("Average SP vs TCS scores per sequence - TCS cutoff 700 - ",per_included))
#ggsave("ref_3dcoffee_avg_sp_vs_tcs_score_per_seq_per_family_3_iter.png",dpi = "retina")

#ref_mtmalign_avg_sp_with_tcs_score_per_seq = merge(ref_mtmalign_avg_sp,tcs_score_per_seq,by="Sequence")
#colnames(ref_mtmalign_avg_sp_with_tcs_score_per_seq)[2] = "Family"

#ο=ggplot(ref_mtmalign_avg_sp_with_tcs_score_per_seq,aes(x=`mTMalign_3_ITER_vs_REF`,y=Psicoffee_vs_REF,color=Family))+geom_point()+xlim(0,100)+ylim(0,100)+geom_abline(intercept =0 , slope = 1,lwd=0.2)+xlab("mTMalign_3_ITER_vs_REF")+ylab("Psicoffee_vs_REF")+ggtitle("Average SP scores per sequence")+labs(size = "TCS mTMalign_3_ITER")
#ggsave("ref_mtmalign_avg_sp_score_per_seq_per_family_3_iter_vs_psicoffee.png",dpi = "retina")
#ο=ggplot(ref_mtmalign_avg_sp_with_tcs_score_per_seq,aes(x=`mTMalign_3_ITER_vs_REF`,y=MSAprobs_vs_REF,color=Family))+geom_point()+xlim(0,100)+ylim(0,100)+geom_abline(intercept =0 , slope = 1,lwd=0.2)+xlab("mTMalign_3_ITER_vs_REF")+ylab("MSAprobs_vs_REF")+ggtitle("Average SP scores per sequence")+labs(size = "TCS mTMalign_3_ITER")
#ggsave("ref_mtmalign_avg_sp_score_per_seq_per_family_3_iter_vs_msaprobs.png",dpi = "retina")

#per_included = round(length(ref_mtmalign_avg_sp_with_tcs_score_per_seq[which(ref_mtmalign_avg_sp_with_tcs_score_per_seq$`mTMalign_3_ITER` > 700),]$`mTMalign_3_ITER`)/length(ref_mtmalign_avg_sp_with_tcs_score_per_seq$`mTMalign_3_ITER`)*100,digits = 2)
#ο=ggplot(ref_mtmalign_avg_sp_with_tcs_score_per_seq[which(ref_mtmalign_avg_sp_with_tcs_score_per_seq$`mTMalign_3_ITER` > 700),],aes(x=`mTMalign_3_ITER_vs_REF`,y=`mTMalign_3_ITER`,color=Family))+geom_point()+xlim(0,100)+xlab("SP mTMalign_3_ITER_vs_REF")+ylab("TCS mTMalign_3_ITER")+ggtitle(paste0("Average SP vs TCS scores per sequence - TCS cutoff 700 - ",per_included))
#ggsave("ref_mtmalign_avg_sp_vs_tcs_score_per_seq_per_family_3_iter.png",dpi = "retina")



#nirmsd_ref_dmpfold = read.table("./selected_comparisons_nirmsd_ref_dmpfold.txt",header = F,stringsAsFactors = F)
#nirmsd_ref_dmpfold = cbind(nirmsd_ref$Family[1:(length(nirmsd_ref$Family)-1)],nirmsd_ref_dmpfold)
#colnames(nirmsd_ref_dmpfold)=c("Family","Famsa","Ginsi","MSAprobs","Tcoffee","Psicoffee","3DCoffee_NAT","3DCoffee_TMalign_NAT","mTMalign_NAT","3DCoffee_3_ITER","3DCoffee_TMalign_3_ITER","mTMalign_3_ITER","3DCoffee_4_ITER","3DCoffee_TMalign_4_ITER","mTMalign_4_ITER","3DCoffee_5_ITER","3DCoffee_TMalign_5_ITER","mTMalign_5_ITER")
#nirmsd_ref_dmpfold$Family = as.character(nirmsd_ref_dmpfold$Family)
#nirmsd_ref_dmpfold = rbind(nirmsd_ref_dmpfold,c("Average",colMeans(nirmsd_ref_dmpfold[,-1])))
#row.names(nirmsd_ref_dmpfold) = nirmsd_ref_dmpfold$Family
#nirmsd_ref_dmpfold[colnames(nirmsd_ref_dmpfold)[-1]] <- sapply(nirmsd_ref_dmpfold[colnames(nirmsd_ref_dmpfold)[-1]],as.numeric)
#pheatmap(nirmsd_ref_dmpfold[,-1],display_numbers = T,angle_col = 45,fontsize = 6,cluster_cols = F,filename = "nirmsd_ref_dmpfold.png")



#mydf_ref_3dcoffee_for_tcs_melted=melt(mydf_ref_3dcoffee[,c(1,5:length(colnames(mydf_ref_3dcoffee)))])
#tcs_score_for_sp_melted=tcs_score_melted[which(tcs_score_melted$variable!='3DCoffee_NAT'),]

#ref_3dcoffee_tcs_vs_sp = cbind(tcs_score_for_sp_melted,mydf_ref_3dcoffee_for_tcs_melted)
#colnames(ref_3dcoffee_tcs_vs_sp) = c("Family","Mode","TCS","Family.1","Mode.1","SP")
#p = ggscatter(ref_3dcoffee_tcs_vs_sp, x = "SP", y = "TCS",color = "Mode", size = 1, conf.int = FALSE) + stat_cor()
#ggpar(p, legend = "right") %>%
#ggexport(filename="ref_3dcoffee_sp_vs_tcs.png")

#boot_struct = read.table("./output.txt",header = F,stringsAsFactors = F)
#colnames(boot_struct) = c("families","avg_SP","rep","structures")
#collect_avg_sp = c()
#for (struct in levels(as.factor(boot_struct$structures))) {
#    collect_avg_sp = c(collect_avg_sp,mean(boot_struct[which(boot_struct$structures == struct),]$avg_SP))
#}
#boot_struct_avg = data.frame(levels(as.factor(boot_struct$structures)),collect_avg_sp)
#colnames(boot_struct_avg) = c("structures","avg_SP")
#b = ggplot(boot_struct_avg, aes(x=as.integer(structures),y=avg_SP)) + geom_point() + geom_line(color="navy") + geom_hline(yintercept=mean(mydf_ref_mtmalign$mTMalign_3_ITER_vs_REF)) + ylab("Average SP for kept families") + xlab("Number of structures per family included")
#ggsave("bootstrap_struct_avg_barplot_ref_mtmalign.png", dpi = "retina")
#b = ggplot(boot_struct_avg, aes(x=structures,y=avg_SP)) + geom_point() + geom_hline(yintercept=mean(mydf_ref_3dcoffee$`3DCoffee_3_ITER_vs_REF`)) + ylab("Average SP for kept families") + xlab("Number of structures per family included")
#ggsave("bootstrap_struct_avg_barplot_ref_3dcoffee.png", dpi = "retina")
#o = ggplot(boot_struct, aes(x=avg_SP,y=structures,color=families)) + geom_point() + scale_color_gradient(low="blue", high="red") + xlab("Average SP for kept families") + ylab("Number of structures per family included") + labs(color = "Number of families kept")
#ggsave("bootstrap_struct.png", dpi = "retina")
