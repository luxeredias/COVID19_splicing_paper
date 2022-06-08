library(dplyr)
library(ggplot2)
library(tidyverse)
library(data.table)
library(clusterProfiler)
options(stringsAsFactors = F)

#set working directory
wd <- "/home/thomaz/projects/covid19/covid19_expression/"
setwd(wd)

#id dictionary
load("data/id_dictionary.RData")

#load custom gmt
custom_gmt_final <- readRDS("series_separated/data/custom_gmt_final.RDS")
custom_gmt_final_fgsea <- split(custom_gmt_final$genes,
                                custom_gmt_final$term)
reactome_gmt <- readRDS("series_separated/data/all_level_reactome.RDS")
reactome_gmt_fgsea <- split(reactome_gmt$gene,
                            reactome_gmt$ont)
term_order=unique(custom_gmt_final$term)
term_class <- data.frame(term=term_order,
                         term_class=c(rep("Innate immune response",6),
                                      rep("Class I MHC",3),
                                      rep("Splicing",9),
                                      rep("mRNA export\nand NMD",2)))
term_class_order <- unique(term_class$term_class)

#load DE files (DETs, DTPs)
DE_files_txtome <- readRDS("series_separated/data/new_DE_objects.RDS")
DETs_all_pc <- DE_files_txtome[[3]]
DETs_all_new <- DE_files_txtome[[2]]
DETs_all_collapsed <- readRDS("series_separated/data/DETs_all_collapsed.RDS")

translatome <- fread("data/differentially_translated_proteins_edited.tsv",
                     dec = ",",check.names = T)

DTPs <- translatome %>%
  #filter(P.value.24h < 0.05,abs(Ratio.24h) > 0.5) %>%
  filter(Species.Names01=="Homo sapiens OX=9606") %>%
  filter(Gene.Symbol01 != "") %>%
  separate_rows(Gene.Symbol01,sep = ";") %>%
  mutate(Gene.Symbol01=gsub(" ","",Gene.Symbol01))

col2time <- data.frame(colA=paste0("Ratio.",c("2h","6h","10h","24h")),
                       colB=paste0("P.value.",c("2h","6h","10h","24h")),
                       time=c("2h","6h","10h","24h"))

DTPs_melt_log2FC <- reshape2::melt(DTPs[,c(3,28:31)]) %>%
  dplyr::rename(gene=1,col=2,log2FC=3) %>%
  dplyr::left_join(col2time,by=c("col"="colA")) %>%
  dplyr::select(gene,time,log2FC)

DTPs_melt_pval <- reshape2::melt(DTPs[,c(3,32:35)]) %>%
  dplyr::rename(gene=1,col=2,pval=3) %>%
  dplyr::left_join(col2time,by=c("col"="colB")) %>%
  dplyr::select(gene,time,pval)

DTPs_all <- DTPs_melt_log2FC %>%
  left_join(DTPs_melt_pval,by=c("gene","time"))

DTPs_all_sig <- DTPs_all %>%
  filter(pval < 0.1)

#plot DTPs
p <- DTPs_all_sig %>%
  mutate(number=1,
         dir=ifelse(log2FC<0,"down","up")) %>%
  group_by(time,dir) %>%
  summarise(DTPs=sum(number)) %>%
  ungroup() %>%
  mutate(DTPs=ifelse(dir=="up",DTPs,DTPs*-1),
         time=factor(time,levels = c("2h","6h","10h","24h")),
         dir=factor(dir,levels = c("up","down"))) %>%
  ggplot(aes(x=time,y=DTPs,fill=dir))+
  geom_col()+
  geom_hline(yintercept = 0)+
  geom_text(aes(label=abs(DTPs),y=(abs(DTPs)+50)*sign(DTPs)),size=4)+
  scale_y_continuous(limits = c(-500,500),
                     labels = c(500,250,0,250,500),
                     name = "")+
  theme_bw()
pdf(file = "series_separated/figures/revised_figures/DTPs_numbers.pdf",
    width = 3.62,height = 3)
print(p)
dev.off()

#Venn diagram of intersection of genes with upregulated unprodcutive isoforms (in ACE2-MOI2.0 cell)
#and downregulated peptides (at any time)
DETs_ACE2MOI2_sig <- DETs_all_new %>%
  filter(series=="ACE2-MOI2.0",
         padj < 0.05,abs(log2FoldChangeShrink) > 2)

genes_with_sig_DETs <- unique(DETs_ACE2MOI2_sig$external_gene_name) #3989 genes

genes_with_sig_DTPs <- unique(DTPs_all_sig$gene) #606 genes

genes_with_sig_DETs_and_DTPs <- intersect(genes_with_sig_DETs,
                                          genes_with_sig_DTPs) #237 genes

translatome_down <- DTPs_all_sig %>%
  filter(log2FC < 0) %>%
  pull(gene) %>% unique()

translatome_up <- DTPs_all_sig %>%
  filter(log2FC > 0) %>%
  pull(gene) %>% unique()

venn_df <- DETs_ACE2MOI2_sig %>%
  filter(external_gene_name %in% genes_with_sig_DETs)
  mutate(dir=ifelse(log2FoldChangeShrink<0,"down","up"),
         group_dir=paste0(expr,"_",dir)) %>%
  filter(!is.na(dir))

DETs_list <- split(venn_df$external_gene_name,venn_df$group_dir)

#Plot heatmap with selected genes with high unprod and low protein (selected terms)
DETs_all_plot <- DETs_all_pc %>%
  dplyr::filter(grepl("ACE2-MOI2.0",series)) %>%
  reshape2::dcast(formula = external_gene_name~series+expr,value.var = "log2FoldChange")

colnames(DETs_all_plot) <- gsub(x = colnames(DETs_all_plot),pattern = "-","_")

genes_with_sig_DETs_df <- DETs_all_plot %>%
  mutate(ACE2_MOI2.0_productive=ifelse(is.na(ACE2_MOI2.0_productive),
                                       yes = 0,no = ACE2_MOI2.0_productive),
         ACE2_MOI2.0_unproductive=ifelse(is.na(ACE2_MOI2.0_unproductive),
                                         yes = 0,no = ACE2_MOI2.0_unproductive),
         dir=ifelse(ACE2_MOI2.0_unproductive>ACE2_MOI2.0_productive,
                    yes = "unproductive",
                    no = "productive")) %>%
  filter(external_gene_name %in% genes_with_sig_DETs) %>%
  mutate(dir=ifelse(is.na(dir),yes = "no_difference",no = dir))

p <- genes_with_sig_DETs_df %>%
  filter(external_gene_name %in% genes_with_sig_DETs_and_DTPs) %>%
  mutate(dir=ifelse(ACE2_MOI2.0_unproductive>ACE2_MOI2.0_productive,
                    yes = "unproductive",
                    no = "productive"),
         title="Genes with DETs and DTPs") %>%
  ggplot(aes(x = ACE2_MOI2.0_unproductive,y = ACE2_MOI2.0_productive,
             color=dir))+
  geom_point()+
  geom_abline(intercept = 0,
              slope = 1)+
  scale_x_continuous(limits = c(-6,6))+
  scale_y_continuous(limits = c(-6,6))+
  facet_wrap(facets = ~title)+
  theme_bw()
pdf(file = "series_separated/figures/revised_figures/genes_DETs_DTPs_prod_unprod.pdf",
    width = 4.7,height = 3.6)
print(p)
dev.off()

p <- genes_with_sig_DETs_df %>%
  filter(external_gene_name %in% genes_with_sig_DETs_and_DTPs) %>%
  left_join(DTPs_all_sig %>% filter(time=="24h"),by=c("external_gene_name"="gene")) %>%
  mutate(dir=ifelse(ACE2_MOI2.0_unproductive>ACE2_MOI2.0_productive,
                    yes = "unproductive",
                    no = "productive"),
         title="Genes with DETs and DTPs",
         dir_translatome=ifelse(log2FC>0,"up","down"),
         dir_translatome=ifelse(is.na(dir_translatome),"not significant",
                                dir_translatome),
         number=1) %>%
  group_by(dir,dir_translatome) %>%
  summarise(value=sum(number)) %>%
  ungroup() %>%
  mutate(title="expression in translatome 24h") %>%
  ggplot(aes(x = dir,y = value,fill = dir_translatome))+
  geom_col()+
  scale_fill_manual(values = c("blue","grey90","red3"))+
  theme_bw()+
  facet_wrap(facets = ~title)
pdf(file = "series_separated/figures/revised_figures/genes_DETs_DTPs_bar_plot.pdf",
    width = 3.44,height = 3.63)
print(p)
dev.off()

DETs_list <- split(genes_with_sig_DETs_df$external_gene_name,
                   genes_with_sig_DETs_df$dir)

translatome_list <- list(translatome_up=translatome_up,
                         translatome_down=translatome_down)

enr <- lapply(DETs_list,FUN = fgsea::fora,pathways=translatome_list,
              universe=unique(c(genes_with_sig_DETs,genes_with_sig_DTPs)))


for(i in 1:length(enr)){
  df <- enr[[i]]
  df$DET_part <- names(enr)[i]
  if(i==1){
    enr_all <- df
  }else{
    enr_all <- rbind(enr_all,df)
  }
}

p <- enr_all %>%
  ggplot(aes(x = pathway,y=DET_part,size=overlap))+
  geom_label(aes(label=overlap,fill=overlap))+
  scale_fill_gradient(low = "pink",high = "red3")+
  theme_bw()
p

# Plot heatmap
DTPs_all_plot <- DTPs_all_sig %>%
  #filter(time=="24h") %>%
  reshape2::dcast(formula = gene~time,value.var = "log2FC",fun.aggregate = sum) %>%
  dplyr::select(1,4,5,2,3) %>%
  dplyr::rename(external_gene_name=1,
                translatome_2h=2,
                translatome_6h=3,
                translatome_10h=4,
                translatome_24h=5)

df_plot <- DTPs_all_plot %>%
  dplyr::left_join(DETs_all_plot,by="external_gene_name") %>%
  dplyr::left_join(custom_gmt_final,by=c("external_gene_name"="genes"))

df_plot_selected <- df_plot %>%
  dplyr::filter(!is.na(term)) %>%
  dplyr::select(1,6,7,2:5,8) %>%
  dplyr::left_join(term_class,by="term") %>%
  dplyr::group_by(term_class) %>%
  dplyr::arrange(desc(ACE2_MOI2.0_unproductive),desc(translatome_24h)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(factor(term_class, levels = term_class_order))

plot_matrix <- df_plot_selected %>%
  dplyr::select(1:3,6,7) %>%
  dplyr::filter(!duplicated(external_gene_name)) %>%
  column_to_rownames("external_gene_name") %>%
  as.matrix()

plot_matrix[is.na(plot_matrix)] <- 0

annot <- df_plot_selected %>%
  filter(!duplicated(external_gene_name)) %>%
  dplyr::select(external_gene_name,term_class) %>%
  mutate(term_class=ifelse(grepl("NMD",term_class),"mRNA export and NMD",term_class)) %>%
  column_to_rownames("external_gene_name")

# DETs_DTPs_common_genes <- unique(rownames(annot))
# saveRDS(DETs_DTPs_common_genes,file = "series_separated/data/DETs_DTPs_common_genes.RDS")

paletteLength <- 1000
myColor <- colorRampPalette(c("blue", "white", "red3"))(paletteLength)
myBreaks <- c(seq(min(plot_matrix), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(plot_matrix)/paletteLength, max(plot_matrix), length.out=floor(paletteLength/2)))

pdf(file = "series_separated/figures/revised_figures/DETs_ACE2_MOI2.0_DTPs_24h_heatmap.pdf",
    width = 4.2,height = 11)
pheatmap::pheatmap(plot_matrix,gaps_col = c(2),
                   annotation_row = annot,
                   gaps_row = c(1,7,36),
                   cluster_rows = F,cluster_cols = F,
                   breaks = myBreaks,color = myColor,annotation_legend = T)
dev.off()

#Run DTPs fgsea and plot
times <- unique(DTPs_all$time)
for(i in 1:length(times)){
  df <- DTPs_all %>%
    filter(time==times[i])
  
  ranks <- df$log2FC
  names(ranks) <- df$gene
  
  enr <- fgseaSimple(pathways = custom_gmt_final_fgsea,stats = ranks,nperm = 10000)
  enr$expr <- "translatome"
  enr$series <- times[i]
  
  enr_reactome <- fgseaSimple(pathways = reactome_gmt_fgsea,stats = ranks,nperm = 10000)
  enr_reactome$expr <- "translatome"
  enr_reactome$series <- times[i]
  
  enr <- rbind(enr,enr_reactome)
  
  up <- df %>%
    filter(pval < 0.1,
           log2FC>0) %>%
    pull(gene) %>%
    unique()
  
  down <- df %>%
    filter(pval < 0.1,
           log2FC<0) %>%
    pull(gene) %>%
    unique()
  
  ora_up <- fora(pathways = custom_gmt_final_fgsea,
                 genes = up,universe = unique(custom_gmt_final$genes))
  ora_up$series <- times[i]
  ora_up$dir <- "up"
  ora_down <- fora(pathways = custom_gmt_final_fgsea,
                   genes = down,universe = unique(custom_gmt_final$genes))
  ora_down$series <- times[i]
  ora_down$dir <- "down"
  
  ora <- rbind(ora_up,ora_down)
  
  if(i==1){
    fgsea_DTPs <- enr
    ora_DTPs <- ora
  }else{
    fgsea_DTPs <- rbind(fgsea_DTPs,enr)
    ora_DTPs <- rbind(ora_DTPs,ora)
  }
}

p <- fgsea_DTPs %>%
  filter(pathway %in% term_order) %>%
  left_join(term_class,by = c("pathway"="term")) %>%
  mutate(is.sig=ifelse(pval < 0.2,"yes","no"),
         term_class=factor(term_class,levels = term_class_order),
         pathway=factor(pathway,levels = rev(term_order)),
         series=factor(series,levels = times)) %>%
  ggplot(aes(x=series,y=pathway,size=-log(pval),color=NES,alpha=is.sig))+
  geom_point()+
  scale_color_gradient2(low = "blue",mid = "white",high = "red3")+
  #scale_alpha_continuous(range = c(0.5,1))+
  scale_alpha_manual(values = c(0.1,1))+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  facet_grid(rows = vars(term_class),
             cols = vars(expr),
             scales = "free", space = "free_y")
p

p <- ora_DTPs %>%
  left_join(term_class,by = c("pathway"="term")) %>%
  mutate(is.sig=ifelse(pval < 0.05,"yes","no"),
         term_class=factor(term_class,levels = term_class_order),
         pathway=factor(pathway,levels = rev(term_order)),
         series=factor(series,levels = times)) %>%
  filter(overlap > 1) %>%
  ggplot(aes(x=series,y=pathway,size=overlap,color=dir,alpha=is.sig))+
  geom_point()+
  scale_color_manual(values = c("blue","red3"))+
  scale_size_continuous(range = c(1.5,9))+
  scale_alpha_manual(values = c(0.5,1))+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  facet_grid(rows = vars(term_class),
             cols = vars(dir),
             scales = "free_y", space = "free_y")
pdf(file = "series_separated/figures/revised_figures/DTPs_ora_up_down.pdf",
    width = 6.12,height = 5.01)
print(p)
dev.off()

# #correlation between DETs prod or unprod vs DTPs (log2FC)
# DETs_prod_unpro <- DETs_all_pc %>%
#   filter(series=="ACE2-MOI2.0") %>%
#   dcast(formula = external_gene_name~expr,value.var = "log2FoldChangeShrink") %>%
#   left_join(DTPs_all_plot,by="external_gene_name") %>%
#   filter(complete.cases(.)) %>%
#   left_join(custom_gmt_final,by=c("external_gene_name"="genes")) %>%
#   left_join(term_class,by="term") %>%
#   mutate(term=ifelse(is.na(term),"Remaining genes",term),
#          term_class=ifelse(is.na(term_class),"Remaining genes",term_class)) %>%
#   mutate(unproductive_to_productive_distance=abs(unproductive-productive))
# 
# p <- DETs_prod_unpro %>%
#   # filter(unproductive > 0,
#   #        translatome_24h < 0) %>%
#   mutate(term_class=factor(term_class,levels = c(term_class_order,"Remaining genes"))) %>%
#   #filter(!is.na(term)) %>%
#   ggplot(aes(x=unproductive,
#              y=translatome_24h,
#              color=term_class))+
#   geom_point()+
#   geom_smooth()+
#   facet_grid(rows = vars(term_class))+
#   #geom_text_repel(aes(label=external_gene_name))+
#   theme_bw()
# p

#plot venn diagram of genes/proteins in common between DETs and DTPs
DETs_all_new <- new_DE_objects[[2]]

DETs_sig <- DETs_all_new %>%
  filter(padj < 0.05,abs(log2FoldChangeShrink) > 2)

genes_with_sig_tx <- DETs_sig %>%
  filter(series=="ACE2-MOI2.0") %>%
  pull(external_gene_name) %>%
  unique()

proteins_sig_any <- DTPs_all %>%
  filter(pval < 0.1) %>%
  pull(gene) %>%
  unique()
 
gene_prot_venn_list <- list(RNAseq=genes_with_sig_tx,
                            Translatome=proteins_sig_any)

intersect(genes_with_sig_tx,proteins_sig_any)

VennDiagram::venn.diagram(x = gene_prot_venn_list,units = "px",
                          filename = "series_separated/figures/revised_figures/RNAseq_translatome_venn.svg",
                          imagetype = "svg",resolution = 300,height = 10,width = 10,hyper.test = T)

DTPs_all_sig_terms <- DTPs_all_sig %>%
  left_join(custom_gmt_final,by=c("gene"="genes")) %>%
  left_join(reactome_gmt,by="gene")
