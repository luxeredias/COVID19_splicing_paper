library(dplyr)
library(stringr)
library(data.table)
library(reshape2)
library(tidyr)
library(clusterProfiler)
library(pheatmap)
library(tibble)
library(ggplot2)
library(fgsea)
options(stringsAsFactors = F)

#set workin directory
wd <- "~/projects/covid19/covid19_expression/"
setwd(wd)

#custom gmt
custom_gmt_final <- readRDS(file = "series_separated/data/custom_gmt_final.RDS")
terms_order <- unique(custom_gmt_final$term)

term_class <- data.frame(term=terms_order,
                         term_class=c(rep("Innate immune response",6),
                                      rep("Class I MHC",3),
                                      rep("Splicing",9),
                                      rep("mRNA export\nand NMD",2)))
term_class_order <- unique(term_class$term_class)

#run enrichment for Stukalov PPI
#load Kogan gmt
kogan_gmt <- read.gmt("data/krogan_ppi_genesets.gmt")
kogan_gmt$term <- str_remove(kogan_gmt$term,pattern = ".*COVID19-")
kogan_gmt$term <- str_remove(kogan_gmt$term,pattern = " protein host PPI from Krogan")

for(i in 1:length(unique(kogan_gmt$term))){
  pt <- unique(kogan_gmt$term)[i]
  genes <- kogan_gmt %>%
    dplyr::filter(term==pt) %>%
    dplyr::pull(gene) %>% unique()
  enr_kogan <- enricher(gene = genes,pAdjustMethod = "BH",TERM2GENE = custom_gmt_final)
  if(class(enr_kogan)=="NULL"){
    next()
  }else{
    enr_kogan <- enr_kogan@result
    enr_kogan$protein <- pt
    if(i==1){
      enr_kogan_all <- enr_kogan
    }else{
      enr_kogan_all <- rbind(enr_kogan_all,enr_kogan)
    }
  }
}

save(enr_kogan_all,file = "series_separated/data/enrichment_results_Krogan_PPI.RData")

p <- enr_kogan_all %>%
  mutate(ID=factor(ID,levels = rev(terms_order))) %>%
  ggplot(aes(x=protein,y=ID,color=-log(p.adjust),size=-log(p.adjust)))+
  geom_point()+
  scale_color_gradient(low = "white",high = "red3")+
  theme_bw()+
  theme(legend.position="top")
pdf(file = "series_separated/figures/revised_figures/Kogan_enrichment_new.pdf",
    width = 8.2,height = 2.6)
print(p)
dev.off()

#run enrichment for Stukalov PPI
stukalov_ppi <- fread("series_separated/data/stukalov_ppi_v2.0.csv")

ppi_stukalov <- stukalov_ppi %>%
  mutate(host_protein=gsub(host_protein,pattern = "\\...",replacement = ""),
         viral_protein=gsub(viral_protein,pattern = "_macroD",replacement = "")) %>%
  filter(!grepl("CoV",host_protein))

viral_proteins <- unique(ppi_stukalov$viral_protein)

for(i in 1:length(viral_proteins)){
  vp <- viral_proteins[i]
  df <- stukalov_ppi %>%
    filter(viral_protein==vp)
  
  genes <- gsub(unique(df$host_protein),pattern = "\\...",replacement = "")
  
  fisher_res <- enricher(gene = genes,pAdjustMethod = "BH",TERM2GENE = custom_gmt_final)
  if(class(fisher_res)=="NULL"){
    next()
  }else{
    fisher_res <- fisher_res@result
    fisher_res$series <- vp
    fisher_res$expr <- "stukalov_ppi"
    if(i==1){
      fisher_res_all_stukalov <- fisher_res
    }else{
      fisher_res_all_stukalov <- rbind(fisher_res_all_stukalov,fisher_res)
    }
  }
}

p <- fisher_res_all_stukalov %>%
  mutate(ID=factor(ID,levels = rev(terms_order))) %>%
  mutate(is.sig=ifelse(pvalue<0.05,"yes","no")) %>%
  ggplot(aes(x=series,y=ID,size=-log(pvalue),color=-log(pvalue),alpha=is.sig))+
  geom_point()+
  scale_color_gradient2(low = "blue",mid = "white",high = "red3")+
  scale_alpha_manual(values = c(0.1,1))+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  facet_wrap(facets = ~expr)
p

#Banerjee et all protein-RNA interaction network
banerjee_prnai <- fread("series_separated/data/banerjee_prnai.csv")

banerjee_prnai_melt <- melt(banerjee_prnai,id.vars = "Gene Name") %>%
  dplyr::rename(viral_protein=1,gene=2,score=3) %>%
  dplyr::group_by(viral_protein,gene) %>%
  dplyr::summarise(score=sum(score)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(score > 0)

viral_proteins <- unique(banerjee_prnai_melt$viral_protein)
for(i in 1:length(viral_proteins)){
  vp <- viral_proteins[i]
  df <- banerjee_prnai_melt %>%
    filter(viral_protein==vp)
  
  genes <- as.character(df$gene)
  
  fisher_res <- enricher(gene = genes,pAdjustMethod = "BH",TERM2GENE = custom_gmt_final)
  if(class(fisher_res)=="NULL"){
    if(i==1){
      fisher_res_all_banerjee <- as.data.frame(matrix(nrow = 0,ncol = 11))
      colnames(fisher_res_all_banerjee) <- colnames(fisher_res_all_stukalov)
    }
    next()
  }else{
    fisher_res <- fisher_res@result
    fisher_res$series <- vp
    fisher_res$expr <- "banerjee_prnai"
    if(i==1){
      fisher_res_all_banerjee <- fisher_res
    }else{
      fisher_res_all_banerjee <- rbind(fisher_res_all_banerjee,fisher_res)
    }
  }
}

p <- fisher_res_all_banerjee %>%
  mutate(ID=factor(ID,levels = rev(terms_order))) %>%
  mutate(is.sig=ifelse(pvalue<0.05,"yes","no")) %>%
  ggplot(aes(x=series,y=ID,size=-log(pvalue),color=-log(pvalue),alpha=is.sig))+
  geom_point()+
  scale_color_gradient2(low = "blue",mid = "white",high = "red3")+
  scale_alpha_manual(values = c(0.1,1))+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  facet_wrap(facets = ~expr)
p    

#plot results for both PPIs together
colnames(enr_kogan_all)[10] <- "series"
enr_kogan_all$expr <- "krogan_ppi"

fisher_res_all <- do.call(rbind,list(enr_kogan_all,fisher_res_all_stukalov,fisher_res_all_banerjee))

p <- fisher_res_all %>%
  left_join(term_class,by=c("ID"="term")) %>%
  mutate(ID=factor(ID,levels = rev(terms_order))) %>%
  mutate(term_class=factor(term_class,levels = term_class_order)) %>%
  mutate(is.sig=ifelse(pvalue<0.05,"yes","no")) %>%
  mutate(series_expr=paste0(series,"_",expr)) %>%
  mutate(series=gsub(series,pattern = " ",replacement = "")) %>%
  mutate(series=gsub(series,pattern = "Orf",replacement = "ORF")) %>%
  mutate(series=gsub(series,pattern = "Nsp",replacement = "NSP")) %>%
  mutate(expr=ifelse(grepl("banerjee",expr),"Banerjee",
                     ifelse(grepl("krogan",expr),"Gordon",
                            "Stukalov"))) %>%
  ggplot(aes(x=series,y=ID,size=-log(pvalue),color=expr))+
  geom_point()+
  #scale_color_gradient(low = "pink",high = "red3")+
  #scale_alpha_manual(values = c(0.5,1))+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90,vjust = 0.5))+
  facet_grid(rows = vars(term_class),
             cols = vars(expr),scales = "free",space = "free")
pdf(file = "series_separated/figures/revised_figures/viral_host_enrichment.pdf",
    width = 8.5,height = 7.25)
print(p)
dev.off()

#plot networks
stukalov_nw <- stukalov_ppi %>%
  dplyr::select(viral_protein,host_protein) %>%
  dplyr::mutate(host_protein=gsub(host_protein,pattern = "\\...",replacement = "")) %>%
  dplyr::filter(host_protein %in% unique(custom_gmt_final$genes)) %>%
  dplyr::filter(!duplicated(.)) %>%
  dplyr::rename(from=1,to=2) %>%
  dplyr::mutate(origin="Stukalov")

krogan_nw <- kogan_gmt %>%
  dplyr::filter(gene %in% unique(custom_gmt_final$genes)) %>%
  dplyr::mutate(term=gsub(term,pattern = " ",replacement = ""),
                term=gsub(term,pattern = "Orf",replacement = "ORF"),
                term=gsub(term,pattern = "Nsp",replacement = "NSP")) %>%
  dplyr::rename(from=1,to=2) %>%
  dplyr::mutate(origin="Krogan")

banerjee_nw <- banerjee_prnai_melt %>%
  dplyr::filter(gene %in% unique(custom_gmt_final$genes)) %>%
  dplyr::select(1,2) %>%
  dplyr::rename(from=1,to=2) %>%
  dplyr::mutate(origin="Banerjee")

genes_in_ppi <- unique(c(stukalov_nw$to,krogan_nw$to,as.character(banerjee_nw$to)))

gmt_nw <- custom_gmt_final %>%
  filter(genes %in% genes_in_ppi) %>%
  mutate(origin="annotation") %>%
  dplyr::rename(from=1,to=2)

nw <- do.call(rbind,list(stukalov_nw,krogan_nw,banerjee_nw))

library(igraph)
library(ggraph)

g <- graph_from_data_frame(nw,directed = F)
E(g)$origin <- nw$origin
E(g)$type <- ifelse(E(g)$origin=="Banerjee",yes = "RNA",no = "protein")
V(g)$class <- c(rep("Virus",18),rep("Host",83))
V(g)$degree <- igraph::degree(g)

p <- ggraph(graph = g,layout = "kk")+
  geom_edge_link(aes(color=origin),width=1)+
  geom_node_point(aes(size=degree,color=class))+
  scale_color_manual(values = c("blue","red3"))+
  scale_size_continuous(range = c(3,12))+
  geom_node_text(aes(label=V(g)$name),repel = F,max.overlaps=100)+
  theme_void()
pdf(file = "series_separated/figures/revised_figures/viral_host_network.pdf",
    width = 11.9,height = 10)
print(p)
dev.off()

#RNA only network
g2 <- delete.edges(g, which(E(g)$type == "protein"))

nodes_elim_RNA <- as.numeric(which(degree(g2) > 0))

g2 <- subgraph(graph = g2,vids = nodes_elim_RNA)

p <- ggraph(graph = g2,layout = "kk")+
  geom_edge_link(aes(color=origin),width=1)+
  geom_node_point(aes(size=degree,color=class))+
  scale_color_manual(values = c("blue","red3"))+
  scale_size_continuous(range = c(3,12))+
  geom_node_text(aes(label=V(g2)$name),repel = F,max.overlaps=100)+
  theme_void()
pdf(file = "series_separated/figures/revised_figures/viral_host_network_RNA.pdf",
    width = 11.9,height = 10)
print(p)
dev.off()

#Protein only network
g3 <- delete.edges(g, which(E(g)$type == "RNA"))

nodes_elim_protein <- as.numeric(which(degree(g3) > 0))

g3 <- subgraph(graph = g3,vids = nodes_elim_protein)

p <- ggraph(graph = g3,layout = "kk")+
  geom_edge_link(aes(color=origin),width=1)+
  geom_node_point(aes(size=degree,color=class))+
  scale_color_manual(values = c("blue","red3"))+
  scale_size_continuous(range = c(3,12))+
  geom_node_text(aes(label=V(g3)$name),repel = F,max.overlaps=100)+
  theme_void()
pdf(file = "series_separated/figures/revised_figures/viral_host_network_protein.pdf",
    width = 11.9,height = 10)
print(p)
dev.off()

#Plot network informations
krogan_info <- data.frame(viral_proteins=length(unique(kogan_gmt$term)),
                          host_proteins=length(unique(kogan_gmt$gene)),
                          origin="Gordon")

stukalov_info <- data.frame(viral_proteins=length(unique(ppi_stukalov$viral_protein)),
                            host_proteins=length(unique(ppi_stukalov$host_protein)),
                            origin="Stukalov")

banerjee_info <- data.frame(viral_proteins=length(unique(banerjee_prnai_melt$viral_protein)),
                            host_proteins=length(unique(banerjee_prnai_melt$gene)),
                            origin="Banerjee")

all_info <- do.call("rbind",list(krogan_info,stukalov_info,banerjee_info))

#Number of target proteins by function
nw_terms <- nw %>%
  left_join(custom_gmt_final,by=c("to"="genes")) %>%
  dplyr::mutate(number=1) 

p <- nw_terms %>%
  dplyr::group_by(from,origin,term) %>%
  dplyr::summarise(number=sum(number)) %>%
  dplyr::ungroup() %>%
  left_join(term_class,by="term") %>%
  mutate(term=factor(term,levels = rev(terms_order)),
         term_class=factor(term_class,levels = term_class_order)) %>%
  ggplot(aes(x=from,y=term,color=origin,size=number))+
  geom_point()+
  scale_size_continuous(range = c(2.5,8))+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90,vjust = 0.5))+
  facet_grid(rows = vars(term_class),
             cols = vars(origin),scales = "free",space = "free")
pdf(file = "series_separated/figures/revised_figures/viral_host_overlap.pdf",
    width = 9.42,height = 7.48)
print(p)
dev.off()

p <- nw_terms %>%
  group_by(from,origin) %>%
  summarise(number=sum(number)) %>%
  ungroup() %>%
  ggplot(aes(x=from,y=number,fill=origin))+
  geom_col()+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,vjust = 0.5,angle = 90))
pdf(file = "series_separated/figures/revised_figures/viral_host_overlap_numbers.pdf",
    width = 4.25,height = 2.86)
print(p)
dev.off()

p <- nw_terms %>%
  group_by(term,origin) %>%
  summarise(number=sum(number)) %>%
  ungroup() %>%
  left_join(term_class,by="term") %>%
  mutate(term=factor(term,levels = terms_order),
         term_class=factor(term_class,levels = term_class_order)) %>%
  ggplot(aes(x=term,y=number,fill=origin))+
  geom_col()+
  theme_bw()+
  coord_flip()+
  #theme(axis.text.x = element_text(hjust = 1,vjust = 0.5,angle = 90))+
  facet_grid(rows = vars(term_class),
             scales = "free_y",space = "free")
pdf(file = "series_separated/figures/revised_figures/viral_host_overlap_term_class.pdf",
    width = 6.17,height = 7)
print(p)
dev.off()

####Run stand alone enrichment of Banerjee PRNAI against Reactome
banerjee_prnai_list <- split(as.character(banerjee_prnai_melt$gene),
                             banerjee_prnai_melt$viral_protein)

reactome_df <- readRDS("series_separated/data/all_level_reactome.RDS")
reactome_fgsea <- split(reactome_df$gene,
                        reactome_df$ont)

univ <- unique(c(reactome_df$gene,as.character(banerjee_prnai_melt$gene)))
banerjee_reactome_enrichment <- lapply(banerjee_prnai_list,FUN = fgsea::fora,
                                       pathways=reactome_fgsea,
                                       universe=univ)
