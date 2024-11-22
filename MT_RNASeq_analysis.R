#--------------------------------------------------------------------------------------------------------------------#

#This .R code file has been written to execute all the analyses related to the RNA-Seq experiment described in 
#the master's thesis, and to generate the appropriate visualizations presented in the report. 
#Author: Daniel Lapuente Hernández

#1.------------------------------------------------------------------------------------------------------------------#
### Libraries required ### 

library(biomaRt)
library(openxlsx)
library(edgeR)
library(tidyverse)
library(factoextra)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library("variancePartition")
library(relaimpo)
library(VennDiagram)
library(ggtern)

#2.------------------------------------------------------------------------------------------------------------------#
### Object preparation for RNA-Seq analysis ### 

first_time=F
if(first_time == T){
  
  #Counts table
  counts=read.xlsx("raw_counts.xlsx", sheet = 1, colNames=T)
  rownames(counts) = counts$Gene.ID
  counts = counts[,2:ncol(counts)]
  
  #Feature data table
  load(file="human_mart.Rdata")
  
  feature_data = getBM(attributes = c("ensembl_gene_id","external_gene_name","description","start_position","end_position","chromosome_name","gene_biotype"), 
                       filters = "ensembl_gene_id", values = rownames(counts) , mart = human_mart)
  
  rownames(feature_data) = feature_data$ensembl_gene_id
  feature_data = feature_data[2:ncol(feature_data)] 
  colnames(feature_data)[1] = "gene_name"
  
  counts = counts[which(rownames(counts) %in% rownames(feature_data)),]
  
  counts=counts[order(rownames(counts)),]
  feature_data=feature_data[order(rownames(feature_data)),]
  
  #Metadata table
  metadata = as.data.frame(matrix(data=NA,nrow=ncol(counts),ncol=5))
  colnames(metadata) = c("cultured_with", "harvest_time","donor","control_group","rnaseq_batch")
  rownames(metadata) = colnames(counts)
  metadata$rnaseq_batch = rep("flow_cell_1", nrow(metadata))
  for (i in 1:nrow(metadata)){
    
    aux = rownames(metadata)[i]
    
    if (grepl("CTL", aux)){
      metadata[i,"cultured_with"] = "None"
      metadata[i,"harvest_time"] = "SN3"
      metadata[i,"donor"] = strsplit(aux, "_")[[1]][2]
      metadata[i,"control_group"] = "None"
    }else{
      metadata[i,"harvest_time"] = substr(aux, start = 1, stop = 3)
      metadata[i,"donor"] = strsplit(aux, "_")[[1]][2]
      
      if (substr(aux, start = 4, stop = 4) == "a"){
        metadata[i,"cultured_with"] = "ASC"
        metadata[i,"control_group"] = "Control1"
      }
      if (substr(aux, start = 4, stop = 4) == "b"){
        metadata[i,"cultured_with"] = "Myocyte"
        metadata[i,"control_group"] = "Control2"
      }
      if (substr(aux, start = 4, stop = 4) == "c"){
        metadata[i,"cultured_with"] = "Tenocyte"
        metadata[i,"control_group"] = "Control1"
      }
      if (substr(aux, start = 4, stop = 4) == "d"){
        metadata[i,"cultured_with"] = "Schwann"
        metadata[i,"control_group"] = "Control2"
      }
    }
  }
  
  metadata$control_group[1:5] = "Control1"
  metadata$control_group[6:10] = "Control2"
  
  metadata$donor = factor(metadata$donor, levels=c("D1","D2","D3","D4","D5",
                                                   "D6","D7","D8","D9","D10"))
  metadata$harvest_time = factor(metadata$harvest_time, levels=c("SN3","SN1","SN5"))
  
  #Check if dimensions are correct
  dim(counts)
  dim(feature_data)
  dim(metadata)
  length(which(rownames(counts)==rownames(feature_data)))
  length(which(colnames(counts)==rownames(metadata)))
  
  saveRDS(counts, file="./Inputs/counts.rds")
  saveRDS(metadata, file="./Inputs/metadata.rds")
  saveRDS(feature_data, file="./Inputs/feature_data.rds")
} else {
  counts = readRDS("./Inputs/counts.rds")
  metadata = readRDS("./Inputs/metadata.rds")
  feature_data = readRDS("./Inputs/feature_data.rds")
}

#3.------------------------------------------------------------------------------------------------------------------#
### Lowly expressed genes filtering and outliers removal ### 

#Gene filtering
filter=function(df_counts = counts, median_TPM, n_groups, df_metadata = metadata){
  
  metadata = df_metadata
  counts = df_counts
  
  number_of_groups=length(unique(metadata[,"cultured_with"]))
  medians=counts[,1:number_of_groups]
  colnames(medians)=unique(metadata[,"cultured_with"])
  for(i in 1:number_of_groups)
    medians[,i]=0
  for(i in 1:number_of_groups)
  {
    set=which(metadata[,"cultured_with"]==colnames(medians)[i])
    chunk=counts[,set]
    medians[,i]=apply(chunk,1,median)
  }
  levels_above_threshold <- function(row) {
    return(summary(factor(row > median_TPM, levels=c(TRUE, FALSE)))[1])
  }
  medians$levels_above_threshold=apply(medians,1,levels_above_threshold)
  genes_to_keep=rownames(medians)[which(medians$levels_above_threshold>=n_groups)]
  
  return(genes_to_keep)
}

genes_to_keep = filter(counts, median_TPM=1, n_groups=1, metadata)

counts = counts[genes_to_keep,]
feature_data = feature_data[genes_to_keep, ]

length(which(rownames(counts)==rownames(feature_data)))
length(which(colnames(counts)==rownames(metadata)))

#Sample outlier identification and removal
alt_counts = counts
alt_metadata = metadata

dge <- DGEList(counts = alt_counts)
dge <- calcNormFactors(dge)
v=voom(dge)
alt_counts=v$E

pca <-prcomp(t(alt_counts))
sum_pca = data.frame(summary(pca)$importance[,c(1:6)])
pca_data <- data.frame(pca$x)
pca_data=cbind(alt_metadata[,c("cultured_with","harvest_time","donor")],pca_data[,1:6])

pca_data$cultured_with=factor(pca_data$cultured_with,levels=unique(pca_data$cultured_with))
pca_data$harvest_time=factor(pca_data$harvest_time,levels=c("SN1","SN3","SN5"))
pca_data$donor=factor(pca_data$donor,levels=c("D1","D2","D3","D4","D5",
                                              "D6","D7","D8","D9","D10"))

colors_cultured_with = c("grey","#6CA8FF", "#F76D7C", "#F7DA63", "#5CD276")
breaks_cultured_with = c("None","ASC", "Myocyte", "Tenocyte", "Schwann")

ggplot(pca_data)+geom_point(aes(x=PC1,y=PC2,color=cultured_with),size=3)+
  xlab(paste0("PC1: ",round(100*sum_pca$PC1[2],digits=2),"% variance"))+
  ylab(paste0("PC2: ",round(100*sum_pca$PC2[2],digits=2),"% variance"))+
  scale_color_manual(breaks=breaks_cultured_with, values=colors_cultured_with)

#4 outliers identified (CTL_D10, SN3a5_D5, "SN3d4_D9", "SN5d4_D9")
samples_to_remove = c("CTL_D10", "SN3a5_D5", "SN3d4_D9", "SN5d4_D9")
counts = counts[,which(!(colnames(counts) %in% samples_to_remove))]
metadata = metadata[which(!(rownames(metadata) %in% samples_to_remove)),]
length(which(colnames(counts)==rownames(metadata)))

#4.------------------------------------------------------------------------------------------------------------------#
### Exploratory analysis (PCA) ### 

alt_counts = counts
alt_metadata = metadata

dge <- DGEList(counts = alt_counts)
dge <- calcNormFactors(dge)
v=voom(dge)
alt_counts=v$E

#PCA 
pca <-prcomp(t(alt_counts))
sum_pca = data.frame(summary(pca)$importance[,c(1:6)])
pca_data <- data.frame(pca$x)
pca_data=cbind(alt_metadata[,c("cultured_with","harvest_time","donor")],pca_data[,1:6])

pca_data$cultured_with=factor(pca_data$cultured_with,levels=unique(pca_data$cultured_with))
pca_data$harvest_time=factor(pca_data$harvest_time,levels=c("SN1","SN3","SN5"))
pca_data$donor=factor(pca_data$donor,levels=c("D1","D2","D3","D4","D5",
                                              "D6","D7","D8","D9","D10"))

colors_cultured_with = c("grey","#6CA8FF", "#F76D7C", "#F7DA63", "#5CD276")
breaks_cultured_with = c("None","ASC", "Myocyte", "Tenocyte", "Schwann")


#Figure 4A: "PC1 vs. PC2 plot highlighting co-culture variable" 
ggplot(pca_data)+geom_point(aes(x=PC1,y=PC2,color=cultured_with),size=3)+
  xlab(paste0("PC1: ",round(100*sum_pca$PC1[2],digits=2),"% variance"))+
  ylab(paste0("PC2: ",round(100*sum_pca$PC2[2],digits=2),"% variance"))+
  scale_color_manual(breaks=breaks_cultured_with, values=colors_cultured_with)


#Figure 4B: "PC1 vs. PC2 plot highlighting harvest time variable"
ggplot(pca_data)+geom_point(aes(x=PC1,y=PC2,color=harvest_time),size=3)+
  xlab(paste0("PC1: ",round(100*sum_pca$PC1[2],digits=2),"% variance"))+
  ylab(paste0("PC2: ",round(100*sum_pca$PC2[2],digits=2),"% variance"))


#Mann-Whitney U tests data presented in Table 1: "Mann-Whitney U tests results on PC2 harvest time clusters"
pairwise.wilcox.test(pca_data$PC2, pca_data$harvest_time, p.adjust.method = "fdr")

#Regress out donor and time effects from expression data
design <- model.matrix(~donor+cultured_with+harvest_time:cultured_with,data=alt_metadata)
fit <- lmFit(alt_counts, design)
donor_effect <- fit$coefficients[, grepl("donor",colnames(fit$coefficients))]
time_effect <- fit$coefficients[, grepl("harvest_time",colnames(fit$coefficients))]
time_effect = time_effect[,!(grepl("cultured_withNone",colnames(time_effect)))]

for (donor in colnames(donor_effect)){
  col_idx = rownames(design)[which(t(design[, donor])==1)]
  alt_counts[,col_idx] = alt_counts[,col_idx] - donor_effect[,donor]
}

for (time in colnames(time_effect)){
  col_idx = rownames(design)[which(t(design[, time])==1)]
  alt_counts[,col_idx] = alt_counts[,col_idx] - time_effect[,time]
}

pca <-prcomp(t(alt_counts))
sum_pca = data.frame(summary(pca)$importance[,c(1:6)])
pca_data <- data.frame(pca$x)
pca_data=cbind(alt_metadata[,c("cultured_with","harvest_time","donor")],pca_data[,1:6])


#Figure 5: "PC1 vs. PC2 plot highlighting co-culture variable, after regressing out donor and harvest time effects from expression data"
ggplot(pca_data)+geom_point(aes(x=PC1,y=PC2,color=cultured_with),size=3)+
  xlab(paste0("PC1: ",round(100*sum_pca$PC1[2],digits=2),"% variance"))+
  ylab(paste0("PC2: ",round(100*sum_pca$PC2[2],digits=2),"% variance"))+
  scale_color_manual(breaks=breaks_cultured_with, values=colors_cultured_with)

#5.------------------------------------------------------------------------------------------------------------------#
### Differential expression analysis (DESeq2) on harvest time effects ### 
 
count_DEGs=function(table,FC_name="logFC",name="adj.P.Val",threshold=0.05){
  down=-length(which(table[,name]<threshold & table[,FC_name]<=-0.2))
  up=length(which(table[,name]<threshold & table[,FC_name]>=0.2))
  return(c(down,up))
}

id_DEGs=function(table,FC_name="logFC",name="adj.P.Val",threshold=0.05){
  down=rownames(table[which(table[,name]<threshold & table[,FC_name]<=-0.2),])
  up=rownames(table[which(table[,name]<threshold & table[,FC_name]>=0.2),])
  return(list(down,up))
}

alt_counts = counts
alt_metadata = metadata

alt_metadata$cultured_with = relevel(factor(alt_metadata$cultured_with),ref="None")
alt_metadata$harvest_time = relevel(factor(alt_metadata$harvest_time),ref="SN1")
alt_metadata$donor = relevel(factor(alt_metadata$donor),ref="D1")

#DESEQ2 (fixed effects)
design=model.matrix(~cultured_with+harvest_time+donor,data=alt_metadata)
dds <- DESeqDataSetFromMatrix(countData = round(alt_counts),colData = alt_metadata,design=design)
dds <- DESeq(dds)
dds <- dds[which(mcols(dds)$betaConv),]

SN3_vs_SN1_preliminary_DESEQ2 <- results(dds, name="harvest_timeSN3")
SN3_vs_SN1 <- lfcShrink(dds, coef="harvest_timeSN3", type="ashr")
SN3_vs_SN1$stat = SN3_vs_SN1_preliminary_DESEQ2$stat

alt_metadata$harvest_time = relevel(factor(alt_metadata$harvest_time),ref="SN3")

design=model.matrix(~cultured_with+harvest_time+donor,data=alt_metadata)
dds <- DESeqDataSetFromMatrix(countData = round(alt_counts),colData = alt_metadata,design=design)
dds <- DESeq(dds)
dds <- dds[which(mcols(dds)$betaConv),]

SN5_vs_SN3_preliminary_DESEQ2 <- results(dds, name="harvest_timeSN5")
SN5_vs_SN3 <- lfcShrink(dds, coef="harvest_timeSN5", type="ashr")
SN5_vs_SN3$stat = SN5_vs_SN3_preliminary_DESEQ2$stat

DEGs_HP_vs_phase=data.frame(
  SN3_vs_SN1=count_DEGs(table=SN3_vs_SN1,FC_name="log2FoldChange",name="padj",threshold=0.05),
  SN5_vs_SN3=count_DEGs(table=SN5_vs_SN3,FC_name="log2FoldChange",name="padj",threshold=0.05))
rownames(DEGs_HP_vs_phase)=c("Down-regulated vs prev. HP", "Up-regulated vs prev. HP")

data_plot_figureDEGs <- data.frame(
  Condition = rep(c("SN3 vs SN1", "SN5 vs SN3"), 2),
  Direction = c("Downregulated", "Downregulated", "Upregulated", "Upregulated"),
  Count = c(DEGs_HP_vs_phase[1,1], DEGs_HP_vs_phase[1,2],
            DEGs_HP_vs_phase[2,1], DEGs_HP_vs_phase[2,2])  
)

data_plot_figureDEGs$Direction = factor(data_plot_figureDEGs$Direction, levels=c("Upregulated","Downregulated"))


#Figure 6: "DEGs between sequential co-culture phases (FDR < 5% & |logFC| > 0.2)" 
ggplot(data_plot_figureDEGs, aes(x = Condition, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "identity", width = 0.6) +
  scale_y_continuous(
    name = "# DEGs",
    limits = c(-1550, 1550),
    breaks = seq(-1550, 1550, by = 250),  
    labels = abs(seq(-1500, 1500, by = 250))  
  )+
  geom_hline(yintercept = 0, color = "black", size = 0.5)+
  scale_fill_manual(
    values = c("Upregulated" = "skyblue", "Downregulated" = "salmon"),
    labels = c("Upregulated" = "More highly\nexpressed later", 
               "Downregulated" = "More highly\nexpressed earlier")  
  )+ xlab(NULL)+
  theme(legend.title = element_blank(),
        axis.text.x = element_text(size = 11, color = "black"),
        legend.key.height = unit(3, "lines"),
        legend.text = element_text(size = 10)) +
  geom_text(aes(label = abs(Count), vjust = ifelse(Count > 0, -0.5, 1.5)), size =3.5)


SN3_vs_SN1 = as.data.frame(SN3_vs_SN1)
SN5_vs_SN3 = as.data.frame(SN5_vs_SN3)

SN3_vs_SN1$DE = "NO"
SN3_vs_SN1$DE[SN3_vs_SN1$log2FoldChange >= 0.2 & SN3_vs_SN1$padj < 0.05] = "UP"
SN3_vs_SN1$DE[SN3_vs_SN1$log2FoldChange <= -0.2 & SN3_vs_SN1$padj < 0.05] = "DOWN"

SN5_vs_SN3$DE = "NO"
SN5_vs_SN3$DE[SN5_vs_SN3$log2FoldChange >= 0.2 & SN5_vs_SN3$padj < 0.05] = "UP"
SN5_vs_SN3$DE[SN5_vs_SN3$log2FoldChange <= -0.2 & SN5_vs_SN3$padj < 0.05] = "DOWN"

#Set functional enrichment for Up and Down regulated genes in each case
GO_UP_SN3_vs_SN1 <- enrichGO(
  gene          = rownames(SN3_vs_SN1)[which(SN3_vs_SN1$DE=="UP")],
  universe      = rownames(SN3_vs_SN1),
  keyType       = 'ENSEMBL',
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05)@result

GO_DOWN_SN3_vs_SN1 <- enrichGO(
  gene          = rownames(SN3_vs_SN1)[which(SN3_vs_SN1$DE=="DOWN")],
  universe      = rownames(SN3_vs_SN1),
  keyType       = 'ENSEMBL',
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05)@result

GO_UP_SN5_vs_SN3 <- enrichGO(
  gene          = rownames(SN5_vs_SN3)[which(SN5_vs_SN3$DE=="UP")],
  universe      = rownames(SN5_vs_SN3),
  keyType       = 'ENSEMBL',
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05)@result

GO_DOWN_SN5_vs_SN3 <- enrichGO(
  gene          = rownames(SN5_vs_SN3)[which(SN5_vs_SN3$DE=="DOWN")],
  universe      = rownames(SN5_vs_SN3),
  keyType       = 'ENSEMBL',
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05)@result

GO_UP_SN3_vs_SN1 = GO_UP_SN3_vs_SN1[which(GO_UP_SN3_vs_SN1$p.adjust <= 0.05),]
GO_DOWN_SN3_vs_SN1 = GO_DOWN_SN3_vs_SN1[which(GO_DOWN_SN3_vs_SN1$p.adjust <= 0.05),]
GO_UP_SN5_vs_SN3 = GO_UP_SN5_vs_SN3[which(GO_UP_SN5_vs_SN3$p.adjust <= 0.05),]
GO_DOWN_SN5_vs_SN3 = GO_DOWN_SN5_vs_SN3[which(GO_DOWN_SN5_vs_SN3$p.adjust <= 0.05),]

GO_UP_SN3_vs_SN1$logOR <- with(GO_UP_SN3_vs_SN1, {
  a <- Count  
  N <- as.numeric(sapply(strsplit(GeneRatio, "/"), `[`, 2))  
  b <- as.numeric(sapply(strsplit(BgRatio, "/"), `[`, 1))  
  M <- as.numeric(sapply(strsplit(BgRatio, "/"), `[`, 2))  
  
  log2((a / (N - a)) / (b / (M - b)))  
})

GO_DOWN_SN3_vs_SN1$logOR <- with(GO_DOWN_SN3_vs_SN1, {
  a <- Count  
  N <- as.numeric(sapply(strsplit(GeneRatio, "/"), `[`, 2)) 
  b <- as.numeric(sapply(strsplit(BgRatio, "/"), `[`, 1)) 
  M <- as.numeric(sapply(strsplit(BgRatio, "/"), `[`, 2)) 
  
  log2((a / (N - a)) / (b / (M - b)))  
})

GO_UP_SN5_vs_SN3$logOR <- with(GO_UP_SN5_vs_SN3, {
  a <- Count 
  N <- as.numeric(sapply(strsplit(GeneRatio, "/"), `[`, 2)) 
  b <- as.numeric(sapply(strsplit(BgRatio, "/"), `[`, 1)) 
  M <- as.numeric(sapply(strsplit(BgRatio, "/"), `[`, 2)) 
  
  log2((a / (N - a)) / (b / (M - b)))  
})

GO_DOWN_SN5_vs_SN3$logOR <- with(GO_DOWN_SN5_vs_SN3, {
  a <- Count 
  N <- as.numeric(sapply(strsplit(GeneRatio, "/"), `[`, 2)) 
  b <- as.numeric(sapply(strsplit(BgRatio, "/"), `[`, 1))  
  M <- as.numeric(sapply(strsplit(BgRatio, "/"), `[`, 2))  
  
  log2((a / (N - a)) / (b / (M - b)))  
})

GO_UP_SN3_vs_SN1 = GO_UP_SN3_vs_SN1[1:ifelse(nrow(GO_UP_SN3_vs_SN1)<10, nrow(GO_UP_SN3_vs_SN1), 10), ]
GO_DOWN_SN3_vs_SN1 = GO_DOWN_SN3_vs_SN1[1:ifelse(nrow(GO_DOWN_SN3_vs_SN1)<10, nrow(GO_DOWN_SN3_vs_SN1), 10), ]
GO_UP_SN5_vs_SN3 = GO_UP_SN5_vs_SN3[1:ifelse(nrow(GO_UP_SN5_vs_SN3)<10, nrow(GO_UP_SN5_vs_SN3), 10), ]
GO_DOWN_SN5_vs_SN3 = GO_DOWN_SN5_vs_SN3[1:ifelse(nrow(GO_DOWN_SN5_vs_SN3)<10, nrow(GO_DOWN_SN5_vs_SN3), 10), ]

GO_UP_SN3_vs_SN1$p.adjust = -log10(GO_UP_SN3_vs_SN1$p.adjust)
GO_DOWN_SN3_vs_SN1$p.adjust = -log10(GO_DOWN_SN3_vs_SN1$p.adjust)
GO_UP_SN5_vs_SN3$p.adjust = -log10(GO_UP_SN5_vs_SN3$p.adjust)
GO_DOWN_SN5_vs_SN3$p.adjust = -log10(GO_DOWN_SN5_vs_SN3$p.adjust)


#Figure 7A: "Upregulated biological pathways in SN3 vs.SN1 harvest times" 
ggplot(GO_UP_SN3_vs_SN1, aes(x = reorder(Description, logOR), y = logOR)) +
  geom_bar(stat = "identity", color = "black", width = 0.05) +  
  geom_point(aes(color = p.adjust, size = Count), shape = 16, stroke = 0.5) +  
  scale_color_gradientn(
    colors = c("darkblue", "green", "yellow"),
    name = expression(-log[10](FDR))
  ) +
  scale_size_continuous(range = c(3, 8), name = "# Genes") + 
  labs(
    x = NULL,  
    y = expression(log[2](OR)),
  ) +
  coord_flip() +
  theme_minimal() +
  guides(
    size = guide_legend(order = 1),   
    color = guide_colorbar(order = 2) 
  )+
  theme(
    axis.text.y = element_text(size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "right")


#Figure 7B: "Downregulated biological pathways in SN3 vs.SN1 harvest times" 
ggplot(GO_DOWN_SN3_vs_SN1, aes(x = reorder(Description, logOR), y = logOR)) +
  geom_bar(stat = "identity", color = "black", width = 0.05) +  
  geom_point(aes(color = p.adjust, size = Count), shape = 16, stroke = 0.5) +  
  scale_color_gradientn(
    colors = c("darkblue", "green", "yellow"),
    name = expression(-log[10](FDR))
  ) +
  scale_size_continuous(range = c(3, 8), name = "# Genes") + 
  labs(
    x = NULL,  
    y = expression(log[2](OR)),
  ) +
  coord_flip() +
  theme_minimal() +
  guides(
    size = guide_legend(order = 1),   
    color = guide_colorbar(order = 2) 
  )+
  theme(
    axis.text.y = element_text(size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "right")


#Figure 8A: "Upregulated biological pathways in SN5 vs.SN3 harvest times" 
ggplot(GO_UP_SN5_vs_SN3, aes(x = reorder(Description, logOR), y = logOR)) +
  geom_bar(stat = "identity", color = "black", width = 0.05) +  
  geom_point(aes(color = p.adjust, size = Count), shape = 16, stroke = 0.5) +  
  scale_color_gradientn(
    colors = c("darkblue", "green", "yellow"),
    name = expression(-log[10](FDR))
  ) +
  scale_size_continuous(range = c(3, 8), name = "# Genes") + 
  labs(
    x = NULL,  
    y = expression(log[2](OR)),
  ) +
  coord_flip() +
  theme_minimal() +
  guides(
    size = guide_legend(order = 1),   
    color = guide_colorbar(order = 2) 
  )+
  theme(
    axis.text.y = element_text(size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "right")


#Figure 8B: "Downregulated biological pathways in SN5 vs.SN3 harvest times" 
ggplot(GO_DOWN_SN5_vs_SN3, aes(x = reorder(Description, logOR), y = logOR)) +
  geom_bar(stat = "identity", color = "black", width = 0.05) +  
  geom_point(aes(color = p.adjust, size = Count), shape = 16, stroke = 0.5) +  
  scale_color_gradientn(
    colors = c("darkblue", "green", "yellow"),
    name = expression(-log[10](FDR))
  ) +
  scale_size_continuous(range = c(3, 8), name = "# Genes") + 
  labs(
    x = NULL,  
    y = expression(log[2](OR)),
  ) +
  coord_flip() +
  theme_minimal() +
  guides(
    size = guide_legend(order = 1),   
    color = guide_colorbar(order = 2) 
  )+
  theme(
    axis.text.y = element_text(size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "right")

#6.------------------------------------------------------------------------------------------------------------------#
### Differential expression analysis (DESeq2) on co-culture signaling effects ### 

alt_counts = counts
alt_metadata = metadata

ASC_vs_None_preliminary_DESEQ2 <- results(dds, name="cultured_withASC")
ASC_vs_None <- lfcShrink(dds, coef="cultured_withASC", type="ashr")
ASC_vs_None$stat = ASC_vs_None_preliminary_DESEQ2$stat

Tenocyte_vs_None_preliminary_DESEQ2 <- results(dds, name="cultured_withTenocyte")
Tenocyte_vs_None <- lfcShrink(dds, coef="cultured_withTenocyte", type="ashr")
Tenocyte_vs_None$stat = Tenocyte_vs_None_preliminary_DESEQ2$stat

Myocyte_vs_None_preliminary_DESEQ2 <- results(dds, name="cultured_withMyocyte")
Myocyte_vs_None <- lfcShrink(dds, coef="cultured_withMyocyte", type="ashr")
Myocyte_vs_None$stat = Myocyte_vs_None_preliminary_DESEQ2$stat

Schwann_vs_None_preliminary_DESEQ2 <- results(dds, name="cultured_withSchwann")
Schwann_vs_None <- lfcShrink(dds, coef="cultured_withSchwann", type="ashr")
Schwann_vs_None$stat = Schwann_vs_None_preliminary_DESEQ2$stat

ASC_vs_None = as.data.frame(ASC_vs_None)
Tenocyte_vs_None = as.data.frame(Tenocyte_vs_None)
Myocyte_vs_None = as.data.frame(Myocyte_vs_None)
Schwann_vs_None = as.data.frame(Schwann_vs_None)

DEGs_CW = data.frame(
  ASC=count_DEGs(table=ASC_vs_None,FC_name="log2FoldChange",name="padj",threshold=0.05),
  Tenocyte=count_DEGs(table=Tenocyte_vs_None,FC_name="log2FoldChange",name="padj",threshold=0.05),
  Myocyte=count_DEGs(table=Myocyte_vs_None,FC_name="log2FoldChange",name="padj",threshold=0.05),
  Schwann=count_DEGs(table=Schwann_vs_None,FC_name="log2FoldChange",name="padj",threshold=0.05))

data_plot_figureDEGs <- data.frame(
  Condition = rep(c("ASC vs CTRL", "Myocyte vs CTRL", "Tenocyte vs CTRL", "Schwann vs CTRL"), 2),
  Direction = c(rep("Downregulated",4), rep("Upregulated",4)),
  Count = c(DEGs_CW[1,1], DEGs_CW[1,3],DEGs_CW[1,2], DEGs_CW[1,4],
            DEGs_CW[2,1], DEGs_CW[2,3],DEGs_CW[2,2], DEGs_CW[2,4])  
)

data_plot_figureDEGs$Direction = factor(data_plot_figureDEGs$Direction, levels=c("Upregulated","Downregulated"))
data_plot_figureDEGs$Condition = factor(data_plot_figureDEGs$Condition, levels=c("ASC vs CTRL","Myocyte vs CTRL",
                                                                           "Tenocyte vs CTRL","Schwann vs CTRL"))


#Figure 10: "DEGs between co-cultured vs. monocultured M2 macrophages (FDR < 5% & |logFC| > 0.2)" 
ggplot(data_plot_figureDEGs, aes(x = Condition, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "identity", width = 0.6) +
  scale_y_continuous(
    name = "# DEGs",
    limits = c(-60, 180),
    breaks = seq(-60, 180, by = 25),  
    labels = abs(seq(-60, 180, by = 25))  
  )+
  geom_hline(yintercept = 0, color = "black", size = 0.5)+
  scale_fill_manual(
    values = c("Upregulated" = "skyblue", "Downregulated" = "salmon"),
    labels = c("Upregulated" = "More highly expressed in\n co-cultured M2 macs.", 
               "Downregulated" = "More highly expressed in\n monocultured M2 macs.")  
  )+ xlab(NULL)+
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle=45, size = 8, color = "black",vjust=0.6),
        legend.key.height = unit(3, "lines"),
        legend.text = element_text(size = 10)) +
  geom_text(aes(label = abs(Count), vjust = ifelse(Count > 0, -0.5, 1.5)), size = 3.5)


ASC_vs_None$DE = "NO"
ASC_vs_None$DE[ASC_vs_None$log2FoldChange >= 0.2 & ASC_vs_None$padj < 0.05] = "UP"
ASC_vs_None$DE[ASC_vs_None$log2FoldChange <= -0.2 & ASC_vs_None$padj < 0.05] = "DOWN"

Myocyte_vs_None$DE = "NO"
Myocyte_vs_None$DE[Myocyte_vs_None$log2FoldChange >= 0.2 & Myocyte_vs_None$padj < 0.05] = "UP"
Myocyte_vs_None$DE[Myocyte_vs_None$log2FoldChange <= -0.2 & Myocyte_vs_None$padj < 0.05] = "DOWN"

Tenocyte_vs_None$DE = "NO"
Tenocyte_vs_None$DE[Tenocyte_vs_None$log2FoldChange >= 0.2 & Tenocyte_vs_None$padj < 0.05] = "UP"
Tenocyte_vs_None$DE[Tenocyte_vs_None$log2FoldChange <= -0.2 & Tenocyte_vs_None$padj < 0.05] = "DOWN"

Schwann_vs_None$DE = "NO"
Schwann_vs_None$DE[Schwann_vs_None$log2FoldChange >= 0.2 & Schwann_vs_None$padj < 0.05] = "UP"
Schwann_vs_None$DE[Schwann_vs_None$log2FoldChange <= -0.2 & Schwann_vs_None$padj < 0.05] = "DOWN"


#Figure 11A: "Venn diagram of upregulated DEGs in each co-culture vs. monocultured M2 macs" 
venn.diagram(
  x = list(id_DEGs(table=ASC_vs_None,FC_name="log2FoldChange",name="padj",threshold=0.05)[[2]],
           id_DEGs(table=Myocyte_vs_None,FC_name="log2FoldChange",name="padj",threshold=0.05)[[2]],
           id_DEGs(table=Tenocyte_vs_None,FC_name="log2FoldChange",name="padj",threshold=0.05)[[2]],
           id_DEGs(table=Schwann_vs_None,FC_name="log2FoldChange",name="padj",threshold=0.05)[[2]]),
  category.names = c("ASC", "Myocyte","Tenocyte","Schwann"), filename="figure11A.jpeg",
  # Circles
  lwd = 1,lty = 'blank',fill=c("#6CA8FF", "#F76D7C", "#F7DA63", "#5CD276"),
  # Numbers
  cex = 1, fontface = "bold", fontfamily = "mono",
  cat.default.pos = "outer")


#Figure 11B: "Venn diagram of downregulated DEGs in each co-culture vs. monocultured M2 macs" 
venn.diagram(
  x = list(id_DEGs(table=ASC_vs_None,FC_name="log2FoldChange",name="padj",threshold=0.05)[[1]],
           id_DEGs(table=Myocyte_vs_None,FC_name="log2FoldChange",name="padj",threshold=0.05)[[1]],
           id_DEGs(table=Tenocyte_vs_None,FC_name="log2FoldChange",name="padj",threshold=0.05)[[1]],
           id_DEGs(table=Schwann_vs_None,FC_name="log2FoldChange",name="padj",threshold=0.05)[[1]]),
  category.names = c("ASC", "Myocyte","Tenocyte","Schwann"), filename="figure11B.jpeg",
  # Circles
  lwd = 1,lty = 'blank',fill=c("#6CA8FF", "#F76D7C", "#F7DA63", "#5CD276"),
  # Numbers
  cex = 1, fontface = "bold", fontfamily = "mono",
  cat.default.pos = "outer")
 

alt_metadata$cultured_with = relevel(factor(alt_metadata$cultured_with),ref="ASC")
alt_metadata$harvest_time = relevel(factor(alt_metadata$harvest_time),ref="SN3")
alt_metadata$donor = relevel(factor(alt_metadata$donor),ref="D1")

design=model.matrix(~cultured_with+harvest_time+donor,data=alt_metadata)
dds <- DESeqDataSetFromMatrix(countData = round(alt_counts),colData = alt_metadata,design=design)
dds <- DESeq(dds)
dds <- dds[which(mcols(dds)$betaConv),]

Tenocyte_vs_None_preliminary_DESEQ2 <- results(dds, name="cultured_withTenocyte")
Tenocyte_vs_None <- lfcShrink(dds, coef="cultured_withTenocyte", type="ashr")
Tenocyte_vs_None$stat = Tenocyte_vs_None_preliminary_DESEQ2$stat

Myocyte_vs_None_preliminary_DESEQ2 <- results(dds, name="cultured_withMyocyte")
Myocyte_vs_None <- lfcShrink(dds, coef="cultured_withMyocyte", type="ashr")
Myocyte_vs_None$stat = Myocyte_vs_None_preliminary_DESEQ2$stat

Schwann_vs_None_preliminary_DESEQ2 <- results(dds, name="cultured_withSchwann")
Schwann_vs_None <- lfcShrink(dds, coef="cultured_withSchwann", type="ashr")
Schwann_vs_None$stat = Schwann_vs_None_preliminary_DESEQ2$stat

Tenocyte_vs_None = as.data.frame(Tenocyte_vs_None)
Myocyte_vs_None = as.data.frame(Myocyte_vs_None)
Schwann_vs_None = as.data.frame(Schwann_vs_None)

DEGs_ASC = data.frame(
  Tenocyte=count_DEGs(table=Tenocyte_vs_None,FC_name="log2FoldChange",name="padj",threshold=0.05),
  Myocyte=count_DEGs(table=Myocyte_vs_None,FC_name="log2FoldChange",name="padj",threshold=0.05),
  Schwann=count_DEGs(table=Schwann_vs_None,FC_name="log2FoldChange",name="padj",threshold=0.05))

data_plot_figureDEGs<- data.frame(
  Condition = rep(c("Myocyte vs. ASC", "Tenocyte vs. ASC", "Schwann vs. ASC"), 2),
  Direction = c(rep("Downregulated",3), rep("Upregulated",3)),
  Count = c(DEGs_ASC[1,2], DEGs_ASC[1,1],DEGs_ASC[1,3],
            DEGs_ASC[2,2], DEGs_ASC[2,1],DEGs_ASC[2,3]))

data_plot_figureDEGs$Direction = factor(data_plot_figureDEGs$Direction, levels=c("Upregulated","Downregulated"))
data_plot_figureDEGs$Condition = factor(data_plot_figureDEGs$Condition, levels=c("Myocyte vs. ASC",
                                                                           "Tenocyte vs. ASC", "Schwann vs. ASC"))

#Figure 12A: "DEGs between PRS lines vs. PRS CK STORM (FDR < 5% & |logFC| > 0.2)" 
ggplot(data_plot_figureDEGs, aes(x = Condition, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "identity", width = 0.6) +
  scale_y_continuous(
    name = "# DEGs",
    limits = c(-100, 100),
    breaks = seq(-100, 100, by = 25),  
    labels = abs(seq(-100, 100, by = 25))  
  )+
  geom_hline(yintercept = 0, color = "black", size = 0.5)+
  scale_fill_manual(
    values = c("Upregulated" = "skyblue", "Downregulated" = "salmon"),
    labels = c("Upregulated" = "Upregulated \n vs. ASC co-cultures", 
               "Downregulated" = "Downregulated\n vs. ASC co-cultures")  
  )+ xlab(NULL)+
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle=45, size = 8, color = "black",vjust=0.6),
        legend.key.height = unit(3, "lines"),
        legend.text = element_text(size = 10)) +
  geom_text(aes(label = ifelse(abs(Count) > 0, abs(Count),""), vjust = ifelse(Count > 0, -0.5, 1)), size = 3.5)


alt_metadata$cultured_with = relevel(factor(alt_metadata$cultured_with),ref="Myocyte")
alt_metadata$harvest_time = relevel(factor(alt_metadata$harvest_time),ref="SN3")
alt_metadata$donor = relevel(factor(alt_metadata$donor),ref="D1")

design=model.matrix(~cultured_with+harvest_time+donor,data=alt_metadata)
dds <- DESeqDataSetFromMatrix(countData = round(alt_counts),colData = alt_metadata,design=design)
dds <- DESeq(dds)
dds <- dds[which(mcols(dds)$betaConv),]

Tenocyte_vs_None_preliminary_DESEQ2 <- results(dds, name="cultured_withTenocyte")
Tenocyte_vs_None <- lfcShrink(dds, coef="cultured_withTenocyte", type="ashr")
Tenocyte_vs_None$stat = Tenocyte_vs_None_preliminary_DESEQ2$stat

ASC_vs_None_preliminary_DESEQ2 <- results(dds, name="cultured_withASC")
ASC_vs_None <- lfcShrink(dds, coef="cultured_withASC", type="ashr")
ASC_vs_None$stat = ASC_vs_None_preliminary_DESEQ2$stat

Schwann_vs_None_preliminary_DESEQ2 <- results(dds, name="cultured_withSchwann")
Schwann_vs_None <- lfcShrink(dds, coef="cultured_withSchwann", type="ashr")
Schwann_vs_None$stat = Schwann_vs_None_preliminary_DESEQ2$stat

ASC_vs_None = as.data.frame(ASC_vs_None)
Tenocyte_vs_None = as.data.frame(Tenocyte_vs_None)
Schwann_vs_None = as.data.frame(Schwann_vs_None)

DEGs_MYO = data.frame(
  ASC = count_DEGs(table=ASC_vs_None,FC_name="log2FoldChange",name="padj",threshold=0.05),
  Tenocyte=count_DEGs(table=Tenocyte_vs_None,FC_name="log2FoldChange",name="padj",threshold=0.05),
  Schwann=count_DEGs(table=Schwann_vs_None,FC_name="log2FoldChange",name="padj",threshold=0.05))

data_plot_figureDEGs<- data.frame(
  Condition = rep(c("ASC vs. Myocyte", "Tenocyte vs. Myocyte", "Schwann vs. Myocyte"), 2),
  Direction = c(rep("Downregulated",3), rep("Upregulated",3)),
  Count = c(DEGs_MYO[1,1], DEGs_MYO[1,2],DEGs_MYO[1,3],
            DEGs_MYO[2,1], DEGs_MYO[2,2],DEGs_MYO[2,3]))

data_plot_figureDEGs$Direction = factor(data_plot_figureDEGs$Direction, levels=c("Upregulated","Downregulated"))
data_plot_figureDEGs$Condition = factor(data_plot_figureDEGs$Condition, levels=c("ASC vs. Myocyte", "Tenocyte vs. Myocyte", 
                                                                                 "Schwann vs. Myocyte"))


#Figure 12B: "DEGs between PRS lines vs. PRS MYO (FDR < 5% & |logFC| > 0.2)" 
ggplot(data_plot_figureDEGs, aes(x = Condition, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "identity", width = 0.6) +
  scale_y_continuous(
    name = "# DEGs",
    limits = c(-101, 100),
    breaks = seq(-100, 100, by = 25),  
    labels = abs(seq(-100, 100, by = 25))  
  )+
  geom_hline(yintercept = 0, color = "black", size = 0.5)+
  scale_fill_manual(
    values = c("Upregulated" = "skyblue", "Downregulated" = "salmon"),
    labels = c("Upregulated" = "Upregulated \n vs. myocyte co-cultures", 
               "Downregulated" = "Downregulated\n vs. myocyte co-cultures")  
  )+ xlab(NULL)+
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle=45, size = 8, color = "black",vjust=0.6),
        legend.key.height = unit(3, "lines"),
        legend.text = element_text(size = 10)) +
  geom_text(aes(label = ifelse(abs(Count) > 0, abs(Count),""), vjust = ifelse(Count > 0, -0.5, 1)), size = 3.5)

#7.------------------------------------------------------------------------------------------------------------------#
### Relative importance study  ### 

alt_counts = counts
alt_metadata = metadata

dge <- DGEList(counts = alt_counts)
dge <- calcNormFactors(dge)
v=voom(dge)
alt_counts=v$E

plot_df = as.data.frame(matrix(data=NA, nrow=nrow(alt_counts), ncol=3))
colnames(plot_df) = c("R2_cultured_with", "R2_harvest_time", "R2_donor")

for (i in 1:nrow(alt_counts)) {

  df <- cbind(Expression=alt_counts[rownames(alt_counts)[i], ], alt_metadata[,c("cultured_with",
                                                                      "harvest_time",
                                                                      "donor")])
  model <- lm(Expression ~ ., data = df) 
  rel_importance <- calc.relimp(model, type = "lmg", rela = TRUE)
  plot_df[i,] = rel_importance$lmg  
}

rownames(plot_df) = rownames(alt_counts)
colnames(plot_df) = c("Co-culture", "Harvest Time", "Donor")


#Figure 13: " Ternary plot of relative importance (𝑅2) for donor, harvest time and co-culture variables" 
ggtern(plot_df, aes(x = `Co-culture`, y = `Harvest Time`, z = Donor)) +
  geom_point(alpha=0.1, size = 1)+
  stat_density_tern(aes(fill = ..level..), geom = "polygon", alpha = 0.8,
                    h = c(1.5, 1.5, 1.5)) + 
  scale_fill_viridis_c(option = "D") +  
  geom_point(alpha=0.001)+
  theme_rgbw() +
  labs(
    x = "Co-culture",
    y = "Harvest Time",
    z = "Donor",
    fill = "Gene density"
  )