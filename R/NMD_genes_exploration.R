NMD_all <- reactome_gmt %>%
  filter(grepl("EJC",ont)) %>%
  pull(gene) %>% unique()

RPLs <- c(unique(NMD_all[grepl("RPL",NMD_all)]),"UBA52")
RPSs <- c(unique(NMD_all[grepl("RPS",NMD_all)]),"FAU")

NMD_gene_list <- list(mRNA_EJC_complex=c("ETF1","NCBP2","NCBP1","PABPC1",
                                         "EIF4G1","UPF1","UPF2","RNPS1","RBM8A",
                                         "MAGOH","MAGOHB","UPF3A","UPF3B",
                                         "EIF4A3","CASC3","GSPT1","GSPT2"),
                  SMG_complex=c("SMG8","SMG9","SMG1"),
                  RPLs=RPLs,
                  RPSs=RPSs,
                  NMD_factors=c("SMG6","SMG7","SMG5","PNRC2","DCP1A"),
                  PP2A_complex=c("PPP2R1A","PPP2R2A","PPP2CA")
                  )

for(i in 1:length(NMD_gene_list)){
  vec <- NMD_gene_list[[i]]
  df <- data.frame(NMD_step=names(NMD_gene_list)[i],
                   gene=vec)
  if(i==1){
    NMD_gene_df <- df
  }else{
    NMD_gene_df <- rbind(NMD_gene_df,df)
  }
}


NMD_DETs <- DETs_all_pc %>%
  filter(series=="ACE2-MOI2.0",
         external_gene_name %in% NMD_gene_df$gene) %>%
  left_join(NMD_gene_df,by=c("external_gene_name"="gene"))

NMD_DETs_prep <- NMD_DETs %>%
  dplyr::select(external_gene_name,expr,log2FoldChangeShrink,NMD_step) %>%
  dplyr::rename(gene=1,log2FC=3)

DTPs_all_prep <- DTPs_all %>%
  dplyr::filter(time %in% c("10h","24h"),
                gene %in% NMD_all) %>%
  dplyr::mutate(expr=paste0("translatome_",time)) %>%
  dplyr::select(gene,expr,log2FC) %>%
  left_join(NMD_gene_df,by="gene")

#Remove RPLs and RPSs
plot_df <- rbind(NMD_DETs_prep,DTPs_all_prep) %>%
  filter(!NMD_step %in% c("RPLs","RPSs"))

order <- plot_df %>%
  dplyr::filter(expr=="unproductive") %>%
  dplyr::group_by(NMD_step) %>%
  dplyr::arrange(desc(log2FC),.by_group = TRUE) %>%
  dplyr::ungroup()

plot_df_2 <- plot_df %>%
  filter(gene %in% order$gene) %>%
  mutate(gene=factor(x = gene,levels = order$gene))

plot_mx <- reshape2::dcast(plot_df_2,formula = gene~expr,value.var = "log2FC") %>%
  dplyr::select(1,2,5,3,4) %>%
  column_to_rownames("gene")

plot_mx[is.na(plot_mx)] <- 0

annot <- NMD_gene_df %>%
  filter(!NMD_step %in% c("RPSs","RPLs")) %>%
  column_to_rownames("gene")

paletteLength <- 1000
myColor <- colorRampPalette(c("blue", "white", "red3"))(paletteLength)
myBreaks <- c(seq(min(plot_mx), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(plot_mx)/paletteLength, max(plot_mx), length.out=floor(paletteLength/2)))

pheatmap::pheatmap(plot_mx,cluster_rows = F,cluster_cols = F,
                   annotation_row = annot,
                   gaps_col = 2,gaps_row = c(16,20,23),
                   breaks = myBreaks,color = myColor,
                   show_rownames = T)

#Only RPLs and RPSs
plot_df <- rbind(NMD_DETs_prep,DTPs_all_prep) %>%
  filter(NMD_step %in% c("RPLs","RPSs"))

order <- plot_df %>%
  dplyr::filter(expr=="unproductive") %>%
  dplyr::group_by(NMD_step) %>%
  dplyr::arrange(desc(log2FC),.by_group = TRUE) %>%
  dplyr::ungroup()

plot_df_2 <- plot_df %>%
  filter(gene %in% order$gene) %>%
  mutate(gene=factor(x = gene,levels = order$gene))

plot_mx <- reshape2::dcast(plot_df_2,formula = gene~expr,value.var = "log2FC") %>%
  dplyr::select(1,2,5,3,4) %>%
  column_to_rownames("gene")

plot_mx[is.na(plot_mx)] <- 0

annot <- NMD_gene_df %>%
  filter(NMD_step %in% c("RPSs","RPLs")) %>%
  column_to_rownames("gene")

paletteLength <- 1000
myColor <- colorRampPalette(c("blue", "white", "red3"))(paletteLength)
myBreaks <- c(seq(min(plot_mx), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(plot_mx)/paletteLength, max(plot_mx), length.out=floor(paletteLength/2)))

pheatmap::pheatmap(plot_mx,cluster_rows = F,cluster_cols = F,
                   annotation_row = annot,
                   gaps_col = 2,gaps_row = 49,
                   breaks = myBreaks,color = myColor,
                   show_rownames = T)
