library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(igraph)
library(fgsea)
library(clusterProfiler)
library(IsoformSwitchAnalyzeR)
library(ggrastr)
library(ggpubr)
library(ggrepel)
library(openxlsx)
options(stringsAsFactors = F)


#id dictionary
load("data/id_dictionary.RData")
load("data/gene_type_dictionary.RData")

#load custom gmt object
reactome_gmt <- readRDS("data/all_level_reactome.RDS")

reactome_gmt_tx <- reactome_gmt %>%
  left_join(id_dictionary[,c(6,5)],by=c("gene"="external_gene_name")) %>%
  filter(!is.na(external_transcript_name))

reactome_gmt_fgsea <- split(reactome_gmt$gene,reactome_gmt$ont)
reactome_gmt_tx_fgsea <- split(reactome_gmt_tx$external_transcript_name,
                               reactome_gmt_tx$ont)

#Create separate terms for each step of the NDM reactome pathway
NMD_all <- reactome_gmt %>%
  filter(grepl("EJC",ont)) %>%
  pull(gene) %>% unique()

RPLs <- c(unique(NMD_all[grepl("RPL",NMD_all)]),"UBA52")
RPSs <- c(unique(NMD_all[grepl("RPS",NMD_all)]),"FAU")

NMD_gene_list <- list(`mRNA-EJC complex`=c("ETF1","NCBP2","NCBP1","PABPC1",
                                         "EIF4G1","UPF1","UPF2","RNPS1","RBM8A",
                                         "MAGOH","MAGOHB","UPF3A","UPF3B",
                                         "EIF4A3","CASC3","GSPT1","GSPT2"),
                      `SMG complex`=c("SMG8","SMG9","SMG1"),
                      `NMD factors`=c("SMG6","SMG7","SMG5","PNRC2","DCP1A"),
                      `PP2A complex`=c("PPP2R1A","PPP2R2A","PPP2CA"),
                      RPLs=RPLs,
                      RPSs=RPSs)

for(i in 1:length(NMD_gene_list)){
  vec <- NMD_gene_list[[i]]
  df <- data.frame(term=names(NMD_gene_list)[i],
                   genes=vec)
  if(i==1){
    NMD_gene_df <- df
  }else{
    NMD_gene_df <- rbind(NMD_gene_df,df)
  }
}

#Include NMD steps in custom gmt
custom_gmt_final <- readRDS("data/custom_gmt_final.RDS")
custom_gmt_final <- rbind(custom_gmt_final,NMD_gene_df)
custom_gmt_final_fgsea <- split(custom_gmt_final$genes,
                                custom_gmt_final$term)

term_order <- unique(custom_gmt_final$term)
term_class <- data.frame(term=term_order,
                         term_class=c(rep("Innate immune response",6),
                                      rep("Class I MHC",3),
                                      rep("Splicing",9),
                                      rep("mRNA nuclear export",1),
                                      rep("NMD",7)))

term_class_order <- unique(term_class$term_class)

#series to label
series <- c("WT-MOI0.2","WT-MOI2.0","ACE2-MOI0.2","ACE2-MOI2.0")
non_coding_types <- c("retained_intron","nonsense_mediated_decay",
                      "processed_transcript")

all_types <- c("protein_coding",non_coding_types)

series_to_label <- data.frame(series=paste0("Series",c(2,5,6,16)),
                              expr=series)

series_expr <- c(paste0(series,c("_productive")),
                 paste0(series,c("_unproductive")))

series_type <- apply(expand.grid(series, all_types), 1, paste0, collapse="_")

series_all <- c(series_expr,series_type)

####Run DESeq2 for each series summing counts of productive and unproductive transcripts####
# expression_files <- list.files("data",pattern = "all_expr",full.names = T)[c(3,4,5,2)]
# for(i in 1:length(expression_files)){
#   load(expression_files[i])
# 
#   #load count tables for genes and transcripts
#   expr_tx <- as.data.frame(all_expr_3DRNAseq$trans_count) %>%
#     rownames_to_column("ensembl_transcript_id_version") %>%
#     dplyr::left_join(id_dictionary,by = "ensembl_transcript_id_version") %>%
#     dplyr::select(-c(1,8,9,10,13)) %>%
#     dplyr::select(c(7,8,10,9),everything())
# 
#   expr_gene <- as.data.frame(all_expr_3DRNAseq$genes_counts) %>%
#     rownames_to_column("ensembl_gene_id_version") %>%
#     left_join(gene_type_dictionary,by="ensembl_gene_id_version") %>%
#     dplyr::select(9,2:7) %>%
#     dplyr::filter(!duplicated(external_gene_name)) %>%
#     column_to_rownames("external_gene_name") %>%
#     round() %>%
#     as.matrix()
# 
#   #run DESeq2 for all genes
#   samples <- colnames(expr_gene)
#   metadata <- data.frame(sample_name=samples,
#                          class=ifelse(grepl("Mock",samples),
#                                       "Mock","SARS_CoV_2"))
# 
#   dds_gene <- DESeqDataSetFromMatrix(countData = expr_gene,
#                                      colData = metadata,
#                                      design = ~class,tidy = F)
# 
#   dds_gene <- DESeq(dds_gene)
#   res_gene <- results(dds_gene,
#                       contrast = c("class","SARS_CoV_2","Mock"))
# 
#   res_gene_out <- as.data.frame(res_gene) %>%
#     rownames_to_column("external_gene_name") %>%
#     left_join(gene_type_dictionary[,-1],by = "external_gene_name")
# 
#   res_gene_shrink <- lfcShrink(dds = dds_gene,
#                                contrast = c("class","SARS_CoV_2","Mock"),
#                                res = res_gene,
#                                type = "normal")
# 
#   res_gene_shrink_out <- as.data.frame(res_gene_shrink) %>%
#     rownames_to_column("external_gene_name") %>%
#     left_join(gene_type_dictionary[,-1],by = "external_gene_name") %>%
#     dplyr::select(external_gene_name,log2FoldChange) %>%
#     dplyr::rename(log2FoldChangeShrink=log2FoldChange)
# 
#   res_gene_out_def <- res_gene_out %>%
#     left_join(res_gene_shrink_out,by="external_gene_name") %>%
#     dplyr::select(1:3,log2FoldChangeShrink,everything()) %>%
#     dplyr::mutate(expr="DEGs",series=series[i])
# 
#   #run DESeq2 for all transcripts
#   expr_tx_all <- expr_tx %>%
#     dplyr::select(-c(2,3,4)) %>%
#     column_to_rownames("external_transcript_name") %>%
#     as.matrix() %>%
#     round()
# 
#   dds_tx_all <- DESeqDataSetFromMatrix(countData = expr_tx_all,
#                                 colData = metadata,
#                                 design = ~class,tidy = F)
#   dds_tx_all <- DESeq(dds_tx_all)
#   res_tx_all <- results(dds_tx_all,contrast = c("class","SARS_CoV_2","Mock"))
# 
#   res_tx_all_out <- as.data.frame(res_tx_all) %>%
#     rownames_to_column("external_transcript_name") %>%
#     left_join(id_dictionary[,c(5,6,9,8)],by="external_transcript_name")
# 
#   res_tx_all_shrink <- lfcShrink(dds = dds_tx_all,
#                                  res = res_tx_all,
#                                  contrast = c("class","SARS_CoV_2","Mock"),
#                                  type = "normal")
# 
#   res_tx_all_shrink_out <- as.data.frame(res_tx_all_shrink) %>%
#     rownames_to_column("external_transcript_name") %>%
#     left_join(id_dictionary[,c(5,6,9,8)],by="external_transcript_name") %>%
#     dplyr::select(external_transcript_name,log2FoldChange) %>%
#     dplyr::rename(log2FoldChangeShrink=log2FoldChange)
# 
#   res_tx_all_out_def <- res_tx_all_out %>%
#     left_join(res_tx_all_shrink_out,by="external_transcript_name") %>%
#     dplyr::select(1:3,log2FoldChangeShrink,everything()) %>%
#     dplyr::mutate(expr="DETs_all",series=series[i])
# 
#   #run DESeq2 with the count sum of transcripts of the same category as new
#   #transcripts (protein coding isoforms are summed to make the "productive component"
#   #of the gene. Non-protein coding isoforms are summed to make the "unproductive
#   #component" of the gene). This is done only for protein coding genes, and the
#   #isoforms used to account for the gene expression are only those annotated as
#   #protein_coding, retained_intron, nonsense_mediated_decay, and processed_transcript
# 
#   expr_tx_pc <- expr_tx %>%
#     dplyr::filter(gene_type=="protein_coding") %>%
#     dplyr::mutate(transcript_class=ifelse(transcript_type=="protein_coding",
#                                    "productive","unproductive")) %>%
#     dplyr::group_by(external_gene_name,transcript_class) %>%
#     dplyr::summarise_if(.predicate = is.numeric,.funs = sum) %>%
#     dplyr::ungroup() %>%
#     dplyr::mutate(feature_id=paste0(external_gene_name,"_",transcript_class)) %>%
#     dplyr::select(-c(external_gene_name,transcript_class)) %>%
#     dplyr::select(feature_id,everything()) %>%
#     column_to_rownames("feature_id") %>%
#     as.matrix() %>%
#     round()
# 
#   dds_tx_pc <- DESeqDataSetFromMatrix(countData = expr_tx_pc,
#                                 colData = metadata,
#                                 design = ~class,tidy = F)
#   dds_tx_pc <- DESeq(dds_tx_pc)
#   res_tx_pc <- results(dds_tx_pc,contrast = c("class","SARS_CoV_2","Mock"))
# 
#   res_tx_pc_out <- as.data.frame(res_tx_pc) %>%
#     rownames_to_column("external_gene_name")
# 
#   res_tx_pc_shrink <- lfcShrink(dds = dds_tx_pc,
#                                 res = res_tx_pc,
#                                 contrast = c("class","SARS_CoV_2","Mock"),
#                                 type = "normal")
# 
#   res_tx_pc_shrink_out <- as.data.frame(res_tx_pc_shrink) %>%
#     rownames_to_column("external_gene_name") %>%
#     dplyr::select(external_gene_name,log2FoldChange) %>%
#     dplyr::rename(log2FoldChangeShrink=log2FoldChange)
# 
#   res_tx_pc_out_def <- res_tx_pc_out %>%
#     left_join(res_tx_pc_shrink_out,by="external_gene_name") %>%
#     separate(external_gene_name,sep = "_",into = c("external_gene_name","expr")) %>%
#     dplyr::select(1:4,log2FoldChangeShrink,everything()) %>%
#     dplyr::mutate(series=series[i])
# 
#   if(i==1){
#     DEGs_all_new <- res_gene_out_def
#     DETs_all_new <- res_tx_all_out_def
#     DETs_pc_new <- res_tx_pc_out_def
#   }else{
#     DEGs_all_new <- rbind(DEGs_all_new,res_gene_out_def)
#     DETs_all_new <- rbind(DETs_all_new,res_tx_all_out_def)
#     DETs_pc_new <- rbind(DETs_pc_new,res_tx_pc_out_def)
#   }
# }
#  
# new_DE_objects <- list(DEGs_all_new,DETs_all_new,DETs_pc_new)

# saveRDS(new_DE_objects,file = "series_separated/data/new_DE_objects.RDS")

####Process DE analysis results####
new_DE_objects <- readRDS("data/new_DE_objects.RDS")

DEGs_all_new <- new_DE_objects[[1]]
DETs_all_new <- new_DE_objects[[2]]
DETs_pc_new <- new_DE_objects[[3]]

#create an equivalent DETs file for the "productive" and "unproductive" components
#of protein coding genes considering the component to be the mean log2FC of all
#protein coding isoforms for the productive component and the mean log2FC of all
#non-protein coding isoforms (retained_intron, NMD, etc) for the unproductive 
#component.
DETs_all_collapsed <- DETs_all_new %>%
  filter(gene_type=="protein_coding") %>%
  mutate(expr=ifelse(transcript_type=="protein_coding",
                                 "productive","unproductive")) %>%
  group_by(external_gene_name,expr,series) %>%
  summarise(log2FoldChange_mean=mean(log2FoldChange),
            log2FoldChangeShrink_mean=mean(log2FoldChangeShrink)) %>%
  ungroup() %>%
  left_join(DETs_pc_new,by=c("external_gene_name","expr","series"))

saveRDS(DETs_all_collapsed,file = "data/DETs_all_collapsed.RDS")

#plot number of DEGs and DETs
DEGs_plot <- DEGs_all_new %>%
  dplyr::mutate(external_transcript_name=external_gene_name,
                transcript_type="gene") %>%
  dplyr::select(external_gene_name,external_transcript_name,
                gene_type,transcript_type,
                log2FoldChangeShrink,padj,
                series,expr) %>%
  dplyr::rename(log2FC=5,FDR=6) %>%
  dplyr::mutate(expr="DEGs") %>%
  dplyr::filter(!is.na(FDR)) %>%
  dplyr::mutate(is.sig=ifelse(abs(log2FC) > 2 & FDR < 0.05,
                              yes=1,no=0),
                dir=ifelse(log2FC>0,"up","down"))

DETs_plot <- DETs_all_new %>%
  dplyr::select(external_gene_name,external_transcript_name,
                gene_type,transcript_type,
                log2FoldChangeShrink,padj,
                series,expr) %>%
  dplyr::rename(log2FC=5,FDR=6) %>%
  dplyr::mutate(expr="DETs") %>%
  dplyr::filter(!is.na(FDR)) %>%
  dplyr::mutate(is.sig=ifelse(abs(log2FC) > 2 & FDR < 0.05,
                              yes=1,no=0),
                dir=ifelse(log2FC>0,"up","down"))

DEGs_DETs_all <- rbind(DEGs_plot,DETs_plot)

p1 <- DEGs_DETs_all %>%
  dplyr::filter(expr=="DEGs") %>%
  dplyr::filter(!is.na(transcript_type)) %>%
  dplyr::mutate(is.sig=ifelse(dir=="up",is.sig,-is.sig)) %>%
  dplyr::mutate(transcript_type=ifelse(gene_type %in% c("protein_coding","lncRNA"),
                                       gene_type,"other\nnon-coding genes")) %>%
  dplyr::group_by(expr,series,transcript_type,dir) %>%
  dplyr::summarise(sigs=sum(is.sig)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(series=factor(series,levels = series_order)) %>%
  dplyr::mutate(dir=factor(dir,levels = c("up","down"))) %>%
  dplyr::mutate(transcript_type=factor(transcript_type,
                                       levels = rev(c("protein_coding","lncRNA",
                                                      "other\nnon-coding genes")))) %>%
  ggplot(aes(x=series,y=sigs,fill=dir))+
  geom_col(alpha=0.8,color="black",size=0.2)+
  geom_text(aes(label=abs(sigs),y=(abs(sigs)+150)*sign(sigs)),size=2)+
  scale_y_continuous(limits = c(-2000,2000),
                     labels = c(2000,1000,0,1000,2000),
                     name = "")+
  scale_x_discrete(name="")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+
  facet_grid(facets = expr~transcript_type)

p2 <- DEGs_DETs_all %>%
  dplyr::filter(expr != "DEGs") %>%
  dplyr::filter(!is.na(transcript_type)) %>%
  dplyr::mutate(is.sig=ifelse(dir=="up",is.sig,-is.sig)) %>%
  dplyr::group_by(expr,series,transcript_type,dir) %>%
  dplyr::summarise(sigs=sum(is.sig)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(series=factor(series,levels = series_order)) %>%
  dplyr::filter(transcript_type %in% c("protein_coding",
                                       "retained_intron",
                                       "nonsense_mediated_decay",
                                       "processed_transcript")) %>%
  dplyr::mutate(transcript_type=factor(transcript_type,
                                       levels = rev(c("protein_coding",
                                                      "retained_intron",
                                                      "nonsense_mediated_decay",
                                                      "processed_transcript")))) %>%
  dplyr::mutate(dir=factor(dir,levels = c("down","up"))) %>%
  ggplot(aes(x=series,y=sigs,fill=dir))+
  geom_col(alpha=0.8,color="black",size=0.2)+
  geom_text(aes(label=abs(sigs),y=(abs(sigs)+300)*sign(sigs)),size=2)+
  theme_bw()+
  scale_fill_manual(values = c("steelblue","tomato3"))+
  scale_y_continuous(limits = c(-2500,2500),breaks = c(-2000,-1000,0,1000,2000),
                     labels = c(2000,1000,0,1000,2000),
                     name = "")+
  scale_x_discrete(name="")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+
  facet_grid(facets = expr~transcript_type)

p <- ggarrange(plotlist = list(p1,p2),
               ncol = 2,nrow = 1,widths = c(2,3),
               common.legend = F,legend = "none")
pdf("figures/figure_01/DEGs_DETs_per_series.pdf",
    width = 12.53,height = 3.15)
print(p)
dev.off()

tableS1_list <- list(DEGs=DEGs_all_new,DETs=DETs_all_new)
write.xlsx(x = tableS1_list,file = "tables/table_S1.xlsx")

#Plot ratio between number of productive over unproductive isoforms in each series
prod_unpro_ratio <- DETs_plot %>%
  dplyr::filter(gene_type=="protein_coding",
                is.sig==1) %>%
  dplyr::filter(transcript_type %in% c("protein_coding",
                                       "retained_intron",
                                       "nonsense_mediated_decay",
                                       "processed_transcript")) %>%
  dplyr::mutate(number=1,
                tx_class=ifelse(transcript_type=="protein_coding",
                                "productive","unproductive")) %>%
  dplyr::group_by(tx_class,series,dir) %>%
  dplyr::summarise(tx_number=sum(number)) %>%
  dplyr::ungroup() %>%
  reshape2::dcast(formula = series~tx_class+dir) %>%
  mutate(up_ratio=productive_up/(unproductive_up+productive_up),
         down_ratio=productive_down/(unproductive_down+productive_down))

prod_unpro_ratio[is.na(prod_unpro_ratio)] <- 1

prod_unpro_ratio <- reshape2::melt(prod_unpro_ratio[,c(1,6,7)])

p <- prod_unpro_ratio %>%
  mutate(series=factor(series,levels = series_order)) %>%
  ggplot(aes(x=series,y=value,group=variable,color=variable))+
  geom_point()+
  geom_line()+
  scale_y_continuous(limits = c(0,1),name = "% productive DETs")+
  scale_color_manual(values = c("tomato3","steelblue"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
        legend.position = "none")+
  facet_wrap(facets = ~variable,nrow = 2)
pdf(file = "figures/figure_01/percent_productive_DETs.pdf",
    width = 1.75,height = 5)
print(p)
dev.off()

#Run functional enrichment analysis using all genesets available:
#all_terms_pathways and custom_gmt (20 pathways custom selected)
#Will use log2FC shrink values and DETs_pc_all instead of DETs_all_collapsed
custom_gmt_final_tx <- custom_gmt_final %>%
  left_join(id_dictionary,by=c("genes"="external_gene_name"))

custom_gmt_final_tx_fgsea <- split(custom_gmt_final_tx$external_transcript_name,
                                   custom_gmt_final_tx$term)

###Run DEGs and DETs (DETs all and DETs_pc) enrichment against custom gmt#####
# for(i in 1:length(series)){
#   srs <- series[i]
# 
#   #DEGs
#   DEGs_df <- DEGs_all_new %>%
#     filter(series==srs,
#            !is.na(log2FoldChangeShrink))
# 
#   ranks <- DEGs_df$log2FoldChangeShrink
#   names(ranks) <- DEGs_df$external_gene_name
# 
#   fgsea_DEGs_custom <- fgseaSimple(pathways = custom_gmt_final_fgsea,
#                                    stats = ranks,nperm = 10000)
# 
#   print(paste0("fgsea custom DEGs ",srs," done"))
# 
#   fgsea_DEGs <- fgsea_DEGs_custom
#   fgsea_DEGs$expr <- "DEGs"
#   fgsea_DEGs$series <- srs
# 
#   #DETs_all
#   DETs_df <- DETs_all_new %>%
#     filter(series==srs,
#            !is.na(log2FoldChangeShrink))
# 
#   prod_df <- DETs_df %>%
#     filter(gene_type=="protein_coding",
#            transcript_type=="protein_coding")
# 
#   ranks <- prod_df$log2FoldChangeShrink
#   names(ranks) <- prod_df$external_transcript_name
# 
#   fgsea_DETs_prod_custom <- fgseaSimple(pathways = custom_gmt_final_tx_fgsea,
#                                         stats = ranks,nperm = 10000)
# 
#   print(paste0("fgsea custom DETs prod ",srs," done"))
# 
#   fgsea_DETs_prod <- fgsea_DETs_prod_custom
# 
#   fgsea_DETs_prod$expr <- "DETs_prod"
#   fgsea_DETs_prod$series <- srs
# 
#   unprod_df <- DETs_df %>%
#     filter(gene_type=="protein_coding",
#            transcript_type!="protein_coding")
# 
#   ranks <- unprod_df$log2FoldChangeShrink
#   names(ranks) <- unprod_df$external_transcript_name
# 
#   fgsea_DETs_unprod_custom <- fgseaSimple(pathways = custom_gmt_final_tx_fgsea,
#                                         stats = ranks,nperm = 10000)
# 
#   print(paste0("fgsea custom DETs unprod ",srs," done"))
# 
#   fgsea_DETs_unprod <- fgsea_DETs_unprod_custom
# 
#   fgsea_DETs_unprod$expr <- "DETs_unprod"
#   fgsea_DETs_unprod$series <- srs
# 
#   fgsea_DETs <- rbind(fgsea_DETs_prod,fgsea_DETs_unprod)
# 
#   #DETs_pc
#   DETs_pc_df <- DETs_pc_new %>%
#     filter(!is.na(log2FoldChangeShrink),
#            series==srs)
# 
#   prod_df <- DETs_pc_df %>%
#     filter(expr=="productive")
# 
#   ranks <- prod_df$log2FoldChangeShrink
#   names(ranks) <- prod_df$external_gene_name
# 
#   fgsea_DETs_pc_prod_custom <- fgseaSimple(pathways = custom_gmt_final_fgsea,
#                                         stats = ranks,nperm = 10000)
# 
#   print(paste0("fgsea custom DETs_pc prod ",srs," done"))
# 
#   fgsea_DETs_pc_prod <- fgsea_DETs_pc_prod_custom
# 
#   fgsea_DETs_pc_prod$expr <- "DETs_pc_prod"
#   fgsea_DETs_pc_prod$series <- srs
# 
#   unprod_df <- DETs_pc_df %>%
#     filter(expr=="unproductive")
# 
#   ranks <- unprod_df$log2FoldChangeShrink
#   names(ranks) <- unprod_df$external_gene_name
# 
#   fgsea_DETs_pc_unprod_custom <- fgseaSimple(pathways = custom_gmt_final_fgsea,
#                                           stats = ranks,nperm = 10000)
# 
#   print(paste0("fgsea custom DETs_pc unprod ",srs," done"))
# 
#   fgsea_DETs_pc_unprod <- fgsea_DETs_pc_unprod_custom
# 
#   fgsea_DETs_pc_unprod$expr <- "DETs_pc_unprod"
#   fgsea_DETs_pc_unprod$series <- srs
# 
#   fgsea_DETs_pc <- rbind(fgsea_DETs_pc_prod,fgsea_DETs_pc_unprod)
# 
#   if(i==1){
#     fgsea_DEGs_all <- fgsea_DEGs
#     fgsea_DETs_all <- fgsea_DETs
#     fgsea_DETs_pc_all <- fgsea_DETs_pc
#   }else{
#     fgsea_DEGs_all <- rbind(fgsea_DEGs_all,fgsea_DEGs)
#     fgsea_DETs_all <- rbind(fgsea_DETs_all,fgsea_DETs)
#     fgsea_DETs_pc_all <- rbind(fgsea_DETs_pc_all,fgsea_DETs_pc)
#   }
#   print(paste0(srs," done"))
# }
# 
# fgsea_res_list <- list(fgsea_DEGs_all,fgsea_DETs_all,fgsea_DETs_pc_all)
# saveRDS(fgsea_res_list,file = "data/fgsea_res_DESeq2_all.RDS")


####Visualize and process fgsea results####
fgsea_res_list <- readRDS("series_separated/data/fgsea_res_DESeq2_all.RDS")

fgsea_DEGs_all <- fgsea_res_list[[1]]
fgsea_DETs_all <- fgsea_res_list[[2]]
fgsea_DETs_pc_all <- fgsea_res_list[[3]]

series_order <- series

p <- fgsea_DETs_pc_all %>%
  filter(pathway %in% term_order) %>%
  mutate(is.sig=ifelse(pval<0.05,"yes","no"),
         class=ifelse(expr=="DETs_pc_prod",
                      "productive","unproductive")) %>%
  left_join(term_class,by = c("pathway"="term")) %>%
  mutate(series=factor(series,levels = series_order),
         term_class=factor(term_class,levels = term_class_order),
         pathway=factor(pathway,levels = rev(term_order))) %>%
  ggplot(aes(x=series,y=pathway,size=-log(pval),color=NES,alpha=is.sig))+
  geom_point()+
  scale_color_gradient2(low = "blue",mid = "white",high = "red3")+
  scale_alpha_manual(values = c(0.1,1))+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  facet_grid(rows = vars(term_class),
             cols = vars(class),
             scales = "free", space = "free_y")
p

p <- fgsea_DETs_all %>%
  filter(pathway %in% term_order) %>%
  mutate(is.sig=ifelse(pval<0.05,"yes","no"),
         class=ifelse(expr=="DETs_prod",
                      "productive","unproductive")) %>%
  left_join(term_class,by = c("pathway"="term")) %>%
  mutate(series=factor(series,levels = series_order),
         term_class=factor(term_class,levels = term_class_order),
         pathway=factor(pathway,levels = rev(term_order))) %>%
  ggplot(aes(x=series,y=pathway,size=-log(pval),color=NES,alpha=is.sig))+
  geom_point()+
  scale_color_gradient2(low = "blue",mid = "white",high = "red3")+
  scale_alpha_manual(values = c(0.1,1))+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  facet_grid(rows = vars(term_class),
             cols = vars(class),
             scales = "free", space = "free_y")
p

fgsea_DEGs_DETs <- rbind(fgsea_DEGs_all,fgsea_DETs_all)

p <- fgsea_DEGs_DETs %>%
  filter(pathway %in% term_order[-c(8,10,15,16,18,20,22:24)]) %>%
  mutate(series=factor(series,levels = series_order)) %>%
  mutate(is.sig=ifelse(pval<0.05,"yes","no"),
         experiment=ifelse(expr=="DEGs","DEGs","DETs"),
         class=ifelse(expr=="DETs_prod",
                      "Productive\nisoforms",
                      ifelse(expr=="DETs_unprod",
                             "Unproductive\nisoforms","Gene\nlevel"))) %>%
  left_join(term_class,by = c("pathway"="term")) %>%
  mutate(term_class=factor(term_class,levels = term_class_order),
         pathway=factor(pathway,levels = rev(term_order))) %>%
  ggplot(aes(x=series,y=pathway,size=-log(pval),color=NES,alpha=is.sig))+
  geom_point()+
  scale_color_gradient2(low = "blue",mid = "white",high = "red3")+
  scale_alpha_manual(values = c(0.1,1))+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  facet_grid(rows = vars(term_class),
             cols = vars(class),
             scales = "free", space = "free_y")
pdf(file = "figures/figure_02/enrichment_DEGs_DETs_new.pdf",
    width = 7.3,height = 6.77)
print(p)
dev.off()

#fgsea enrichment only for ACE2 cells (first) and WT cells (next)
p <- fgsea_DEGs_DETs %>%
  filter(pathway %in% term_order[-c(8,10,15,16,18,20,22:24)]) %>%
  filter(series %in% c("ACE2-MOI0.2","ACE2-MOI2.0")) %>%
  mutate(series=factor(series,levels = series_order)) %>%
  mutate(is.sig=ifelse(pval<0.05,"yes","no"),
         experiment=ifelse(expr=="DEGs","DEGs","DETs"),
         class=ifelse(expr=="DETs_prod",
                      "Productive\nisoforms",
                      ifelse(expr=="DETs_unprod",
                             "Unproductive\nisoforms","Gene\nlevel"))) %>%
  left_join(term_class,by = c("pathway"="term")) %>%
  mutate(term_class=factor(term_class,levels = term_class_order),
         pathway=factor(pathway,levels = rev(term_order))) %>%
  ggplot(aes(x=series,y=pathway,size=-log(pval),color=NES,alpha=is.sig))+
  geom_point()+
  scale_size_continuous(range = c(2,12))+
  scale_color_gradient2(low = "blue",mid = "white",high = "red3")+
  scale_alpha_manual(values = c(0.1,1))+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,vjust = 0.5,angle = 90,size = 10),
        axis.text.y = element_text(size = 10),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10))+
  facet_grid(rows = vars(term_class),
             cols = vars(class),
             scales = "free", space = "free_y")
pdf(file = "figures/figure_02/enrichment_DEGs_DETs_new_ACE2.pdf",
    width = 9.02,height = 9.64)
print(p)
dev.off()

p <- fgsea_DEGs_DETs %>%
  filter(pathway %in% term_order[-c(8,10,15,16,18,20,22:24)]) %>%
  filter(series %in% c("WT-MOI0.2","WT-MOI2.0")) %>%
  mutate(series=factor(series,levels = series_order)) %>%
  mutate(is.sig=ifelse(pval<0.05,"yes","no"),
         experiment=ifelse(expr=="DEGs","DEGs","DETs"),
         class=ifelse(expr=="DETs_prod",
                      "Productive\nisoforms",
                      ifelse(expr=="DETs_unprod",
                             "Unproductive\nisoforms","Gene\nlevel"))) %>%
  left_join(term_class,by = c("pathway"="term")) %>%
  mutate(term_class=factor(term_class,levels = term_class_order),
         pathway=factor(pathway,levels = rev(term_order))) %>%
  ggplot(aes(x=series,y=pathway,size=-log(pval),color=NES,alpha=is.sig))+
  geom_point()+
  scale_color_gradient2(low = "blue",mid = "white",high = "red3")+
  scale_alpha_manual(values = c(0.1,1))+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,vjust = 0.5,angle = 90,size = 12),
        axis.text.y = element_text(size = 12))+
  facet_grid(rows = vars(term_class),
             cols = vars(class),
             scales = "free", space = "free_y")
pdf(file = "figures/figure_02/enrichment_DEGs_DETs_new_WT.pdf",
    width = 7,height = 8)
print(p)
dev.off()

p <- fgsea_DEGs_DETs %>%
  filter(pathway %in% term_order[-c(8,10,15,16,18,20,22:24)]) %>%
  filter(series == "ACE2-MOI2.0") %>%
  filter(expr!="DEGs") %>%
  mutate(series=factor(series,levels = series_order)) %>%
  mutate(is.sig=ifelse(pval<0.05,"yes","no"),
         experiment=ifelse(expr=="DEGs","DEGs","DETs"),
         class=ifelse(expr=="DETs_prod",
                      "Productive\nisoforms",
                      ifelse(expr=="DETs_unprod",
                             "Unproductive\nisoforms","Gene\nlevel"))) %>%
  left_join(term_class,by = c("pathway"="term")) %>%
  mutate(term_class=factor(term_class,levels = term_class_order),
         pathway=factor(pathway,levels = rev(term_order))) %>%
  ggplot(aes(x=series,y=pathway,size=-log(pval),color=NES,alpha=is.sig))+
  geom_point()+
  scale_size_continuous(range = c(2,12))+
  scale_color_gradient2(low = "blue",mid = "white",high = "red3")+
  scale_alpha_manual(values = c(0.1,1))+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,vjust = 0.5,angle = 90,size = 10),
        axis.text.y = element_text(size = 10),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10))+
  facet_grid(rows = vars(term_class),
             cols = vars(class),
             scales = "free", space = "free_y")
pdf(file = "figures/figure_02/enrichment_DEGs_DETs_new_ACE2-MOI2.0.pdf",
    width = 7.5,height = 9.64)
print(p)
dev.off()

#find genes not DE with DE isoforms in each term
#find genes DE with only unproductive DE isoforms in each term
#plot log2FC plots for these cases

#plot log2FC unproductive vs log2FC productive
genes_to_label <- readRDS("data/DETs_DTPs_common_genes.RDS")

innate_immune_genes_label <- c("OAS2","CXCL2","CXCL8","MX1","IFIH1",
                               "CCL2","RELB","TANK","NFKB1","IL6","REL",
                               "RSAD2","OAS1","IRF7","TRIM22","OAS3","MYD88",
                               "RELA","CHUK","STAT2","NFKB2")

class_I_MHC_genes_label <- c("HLA-B","HLA-A","B2M","UBA7","PSMB5","PSMB3",
                             "PSMB6","UBA1","TAP2","NEDD4L","RBBP6","HERC5",
                             "HERC6","PSMB9","UBE2V2","RNF213")

high_unprod_genes <- c("RPS27L","U2AF1","EIF4E","POM121C")

NMD_genes <- c("UPF1","UPF2","ETF1","EIF4G1")

genes_to_label <- unique(c(genes_to_label,innate_immune_genes_label,
                           class_I_MHC_genes_label,high_unprod_genes,NMD_genes))

saveRDS(genes_to_label,file = "series_separated/data/DETs_DTPs_common_genes_2.RDS")

prod <- DETs_all_collapsed %>%
  dplyr::filter(expr=="productive") %>%
  dplyr::select(external_gene_name,log2FoldChangeShrink,series) %>%
  dplyr::rename(log2FCproductive=2)

prod_unprod_plot <- DETs_all_collapsed %>%
  dplyr::filter(expr=="unproductive") %>%
  dplyr::select(external_gene_name,log2FoldChangeShrink,series) %>%
  dplyr::rename(log2FCunproductive=2) %>%
  dplyr::left_join(prod,by=c("external_gene_name","series"))

df_plot <- prod_unprod_plot %>%
  mutate(unprod_to_prod=abs(log2FCunproductive - log2FCproductive),
         is.higher=ifelse(log2FCunproductive > log2FCproductive,
                          yes = "unprod",no="prod")) %>%
  left_join(custom_gmt_final,by=c("external_gene_name"="genes")) %>%
  left_join(DEGs_all_new,by=c("external_gene_name","series")) %>%
  #filter(abs(log2FCunprod)>1.5 | abs(log2FCprod) > 1.5) %>%
  mutate(series=factor(series,levels = series_order),
         term=factor(term,levels = term_order),
         is.term=ifelse(is.na(term),"no","yes"),
         is.sig.gene=ifelse(padj<0.05,"yes","no"),
         prod_to_gene=log2FCproductive/log2FoldChangeShrink,
         unprod_to_gene=log2FCunproductive/log2FoldChangeShrink) %>%
  filter(series=="ACE2-MOI2.0") %>%
  filter(external_gene_name %in% custom_gmt_final$genes) %>%
  filter(!is.na(is.higher)) %>%
  mutate(labs=ifelse(external_gene_name %in% genes_to_label,
                     external_gene_name,"")) %>%
  # mutate(labs=ifelse((abs(log2FCunproductive)>2 | abs(log2FCproductive)>2) & is.sig.gene=="yes",
  #                    yes = external_gene_name,"")) %>%
  filter(!is.na(is.sig.gene)) %>%
  left_join(term_class,by="term") %>%
  mutate(term_class=factor(term_class,levels = term_class_order)) %>%
  filter(!term %in% c("Nonsense-mediated decay","RPLs","RPSs"))

p <- df_plot %>%
  ggplot(aes(x=log2FCunproductive,y = log2FCproductive,
             color=is.higher))+
  geom_point(aes(alpha=is.sig.gene),size=2.5)+
  geom_abline(linetype="dashed",size=.2)+
  geom_text_repel(data = subset(df_plot,labs!=""),
                  aes(label=labs),color="black",max.overlaps = 1000,size=4,
                  force = T,force_pull = -0.01,
                  segment.linetype=2,segment.size=.2,
                  point.padding = unit(1, "lines"))+
  scale_alpha_manual(values = c(0.1,1))+
  scale_y_continuous(limits = c(-6,6))+
  scale_x_continuous(limits = c(-6,6))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        strip.text.x = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  facet_wrap(facets = ~term_class,scales = "fixed")
pdf(file = "series_separated/figures/revised_figures/log2FCunprod_vs_log2FCprod.pdf",
    width = 11,height = 7.45)
print(p)
dev.off()


####IZA RODAR A PARTIR DAQUI!!!####
#plot log2FC plots for selected genes
DETs_DEGs_plot <- rbind(DEGs_plot,DETs_plot) #DEGs_plot was created in line 527, DETS_plot in line 541

tx_types <- as.data.frame(table(DETs_plot$transcript_type)) %>%
  arrange(desc(Freq)) %>%
  pull(Var1) %>%
  as.character()

tx_types_order <- c("gene",tx_types)

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

cols = c(ggplotColours(n = 5),rep("grey",length(tx_types_order)-5))
cols <- cols[c(1,4,3,2,5,6:length(cols))]
names(cols) <- tx_types_order

#genes2plot <- genes_to_label
genes2plot <- DETs_all_new %>%
  filter(series=="ACE2-MOI2.0") %>%
  filter(gene_type=="protein_coding") %>%
  filter(abs(log2FoldChangeShrink) > 2,padj < 0.05) %>%
  pull(external_gene_name) %>%
  unique()

for(i in 1:length(genes2plot)){
  gn <- genes2plot[i]
  log2FC_plot_df <- DETs_DEGs_plot %>%
    filter(external_gene_name==gn) %>%
    #filter(series=="ACE2-MOI2.0") %>%
    mutate(series=factor(series,levels = series_order),
           transcript_type=factor(transcript_type,levels = tx_types_order),
           is.sig=ifelse(FDR<0.05,"yes","no")) %>%
    dplyr::group_by(transcript_type) %>%
    dplyr::arrange(desc(log2FC),.by_group = T) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(external_transcript_name=factor(external_transcript_name,levels = unique(external_transcript_name)))
  if(length(log2FC_plot_df$external_gene_name)==0){
    next()
  }
  cols_plot <- cols[which(names(cols) %in% as.character(unique(log2FC_plot_df$transcript_type)))]
  ntxs <- length(unique(log2FC_plot_df$external_transcript_name))
  if(ntxs < 5){
    ntxs <- 10
  }
  if(ntxs > 15){
    ntxs <- 15
  }
  fil <- paste0("series_separated/figures/revised_figures/log2FC_plots/",
                gn,"_lfc_plot_ACE2_MOI2.0.pdf")
  p <- log2FC_plot_df %>%
    ggplot(aes(x=external_transcript_name,y=log2FC,
               fill=transcript_type,alpha=is.sig))+
    geom_col()+
    geom_hline(yintercept = 0)+
    scale_alpha_manual(values = c(0.2,1))+
    scale_fill_manual(values = cols_plot)+
    theme_bw()+
    theme(axis.text.x = element_text(vjust = .5,hjust = 1,angle = 90))+
    facet_wrap(facets = ~series)
  pdf(file = fil,width = ntxs,height = 4.7)
  print(p)
  dev.off()

  print(i/length(genes2plot)*100)
}

#Plot pannel with selected genes
selected_genes <- c("IRF7","IL6","OAS3","HLA-B","RBM5","RPL29")

log2FC_plot_df <- DETs_DEGs_plot %>%
  filter(external_gene_name %in% selected_genes) %>%
  filter(series=="ACE2-MOI2.0") %>%
  mutate(series=factor(series,levels = series_order),
         transcript_type=factor(transcript_type,levels = tx_types_order),
         is.sig=ifelse(FDR<0.05,"yes","no")) %>%
  dplyr::group_by(transcript_type) %>%
  dplyr::arrange(desc(log2FC),.by_group = T) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(external_transcript_name=factor(external_transcript_name,levels = unique(external_transcript_name))) %>%
  dplyr::mutate(short_name=gsub(external_transcript_name,
                         pattern = paste0(paste0(selected_genes,"-"),collapse = "|"),
                         replacement = ""))
cols_plot <- cols[which(names(cols) %in% as.character(unique(log2FC_plot_df$transcript_type)))]

plot_list <- list()
for(i in 1:length(selected_genes)){
  p <- log2FC_plot_df %>%
    filter(external_gene_name==selected_genes[i]) %>%
    ggplot(aes(x=external_transcript_name,y=log2FC,
               fill=transcript_type,alpha=is.sig))+
    geom_col()+
    geom_hline(yintercept = 0)+
    scale_alpha_manual(values = c(0.2,1))+
    scale_fill_manual(values = cols_plot)+
    scale_y_continuous(limits = c(-8,8))+
    scale_x_discrete(labels = log2FC_plot_df %>% filter(external_gene_name==selected_genes[i]) %>%
                         pull(short_name))+
    theme_bw()+
    theme(axis.text.x = element_text(vjust = .5,hjust = 1,angle = 90,size = 8),
          axis.title.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    facet_wrap(facets = ~external_gene_name,scales = "free_x",nrow = 1)
  plot_list[[i]] <- p
}

p <- ggarrange(plotlist = plot_list,
               ncol = 6,common.legend = T,
               widths = c(1,0.9,.7,1,1.3,0.9),
               legend = "bottom",align = "hv")
pdf(file = "series_separated/figures/revised_figures/log2FC_plots_selected_genes.pdf",
    width = 10,height = 2.2)
print(p)
dev.off()


p <- fgsea_DETs_all %>%
  mutate(is.sig=ifelse(padj<0.05,"yes","no")) %>%
  ggplot(aes(x=NES,y=-log(padj),color=is.sig))+
  geom_point()
p  
