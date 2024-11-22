#--------------------------------------------------------------------------------------------------------------------#

#This .R code file has been written to execute all the analyses related to the multiplex experiment described in 
#the master's thesis, and to generate the appropriate visualizations presented in the report. 
#Author: Daniel Lapuente Hernández

#1.------------------------------------------------------------------------------------------------------------------#
### Libraries required ### 

library(openxlsx)
library(biomaRt)
library(limma)
library(edgeR)
library(ggplot2)

#2.------------------------------------------------------------------------------------------------------------------#
### Object preparation for multiplex data ### 

first_time = F
if (first_time == T){
  
  #General Data
  data1=read.xlsx("processed_files/Plex1.xlsx", sheet = 1, colNames=T)
  data2=read.xlsx("processed_files/Plex2.xlsx", sheet = 1, colNames=T)
  data3=read.xlsx("processed_files/Plex3.xlsx", sheet = 1, colNames=T)
  data4=read.xlsx("processed_files/Plex4.xlsx", sheet = 1, colNames=T)
  data5=read.xlsx("processed_files/Plex5.xlsx", sheet = 1, colNames=T)
  data6=read.xlsx("processed_files/Plex6.xlsx", sheet = 1, colNames=T)
  data7=read.xlsx("processed_files/Plex7.xlsx", sheet = 1, colNames=T)
  data8=read.xlsx("processed_files/Plex8.xlsx", sheet = 1, colNames=T)
  
  rownames(data1) = data1$Sample.ID
  rownames(data2) = data2$Sample.ID
  rownames(data3) = data3$Sample.ID
  rownames(data4) = data4$Sample.ID
  rownames(data5) = data5$Sample.ID
  rownames(data6) = data6$Sample.ID
  rownames(data7) = data7$Sample.ID
  rownames(data8) = data8$Sample.ID
  
  data = rbind(data1, data2)
  data = rbind(data, data3)
  data = rbind(data, data4)
  data = rbind(data, data5)
  data = rbind(data, data6)
  data = rbind(data, data7)
  data = rbind(data, data8)
  
  data = data[,2:ncol(data)]
  
  rm(data1, data2, data3, data4, data5, data6, data7, data8)
  
  #General Metadata
  metadata1=read.xlsx("processed_files/Plex1.xlsx", sheet = 2, colNames=T)
  metadata2=read.xlsx("processed_files/Plex2.xlsx", sheet = 2, colNames=T)
  metadata3=read.xlsx("processed_files/Plex3.xlsx", sheet = 2, colNames=T)
  metadata4=read.xlsx("processed_files/Plex4.xlsx", sheet = 2, colNames=T)
  metadata5=read.xlsx("processed_files/Plex5.xlsx", sheet = 2, colNames=T)
  metadata6=read.xlsx("processed_files/Plex6.xlsx", sheet = 2, colNames=T)
  metadata7=read.xlsx("processed_files/Plex7.xlsx", sheet = 2, colNames=T)
  metadata8=read.xlsx("processed_files/Plex8.xlsx", sheet = 2, colNames=T)
  
  rownames(metadata1) = metadata1$Sample.ID
  rownames(metadata2) = metadata2$Sample.ID
  rownames(metadata3) = metadata3$Sample.ID
  rownames(metadata4) = metadata4$Sample.ID
  rownames(metadata5) = metadata5$Sample.ID
  rownames(metadata6) = metadata6$Sample.ID
  rownames(metadata7) = metadata7$Sample.ID
  rownames(metadata8) = metadata8$Sample.ID
  
  metadata = rbind(metadata1, metadata2)
  metadata = rbind(metadata, metadata3)
  metadata = rbind(metadata, metadata4)
  metadata = rbind(metadata, metadata5)
  metadata = rbind(metadata, metadata6)
  metadata = rbind(metadata, metadata7)
  metadata = rbind(metadata, metadata8)
  
  metadata = metadata[,c("Sample.Name","Well","DilutionFactor","type","harvest_time","donor", "batch")]
  
  rm(metadata1, metadata2, metadata3, metadata4, metadata5, metadata6, metadata7, metadata8)
  
  #General Feature Data
  load(file="human_mart.Rdata")
  
  plex_prot_list = c(
    "ENSG00000181092", "ENSG00000172156", "ENSG00000138685", "ENSG00000164400", 
    "ENSG00000163739", "ENSG00000019991", "ENSG00000197919", "ENSG00000111537", 
    "ENSG00000115008", "ENSG00000125538", "ENSG00000136689", "ENSG00000113520", 
    "ENSG00000113525", "ENSG00000136244", "ENSG00000169429", "ENSG00000145839", 
    "ENSG00000136634", "ENSG00000168811","ENSG00000150782", "ENSG00000169245", 
    "ENSG00000128342", "ENSG00000108691", "ENSG00000277632", "ENSG00000275302", 
    "ENSG00000196611", "ENSG00000149968", "ENSG00000134259", "ENSG00000119630", 
    "ENSG00000271503", "ENSG00000049130", "ENSG00000107562", "ENSG00000102265", 
    "ENSG00000232810", "ENSG00000226979","ENSG00000112715", "ENSG00000165197"
  )
  
  feature_data = getBM(attributes = c("ensembl_gene_id","external_gene_name","description","start_position","end_position","chromosome_name","gene_biotype"), 
                       filters = "ensembl_gene_id", values = plex_prot_list , mart = human_mart)
  
  rownames(feature_data) = c("HGF", "SCF", "TIMP-1", "SDF-1.alpha", "MCP-1.(CCL2)", "IFN.gamma",
                             "VEGF-A", "IL-4", "IL-5", "IL-1.alpha", "PlGF-1", "IL-1.beta", "LIF",
                             "NGF.beta", "IL-6", "IL-10", "IL-1RA", "FGF-2" , "IL-9", "MMP-3", "IL-18",
                             "GRO.alpha.(CXCL1)", "GM-CSF", "VEGF-D", "IL-12p70", "IP-10.(CXCL10)",
                             "IL-8.(CXCL8)", "Eotaxin.(CCL11)", "Adiponectin", "MMP-1", "IFN.alpha",
                             "TNF.beta", "TNF.alpha", "RANTES.(CCL5)", "MIP-1.beta.(CCL4)", "MIP-1.alpha.(CCL3)")
  
  feature_data = feature_data[colnames(data),]
  
  feature_data$upper_OOR_subs = 3*c(825852.48,692.5,10808.62,47926.25,3695.16,23014.52,2513.35,27530.31,1432.83,7214.63,132376.42,45918.14,22998.33,44288.34,6894.22,41639.46,1916.41,27357.27,29959.91,2091.42,9824.96,11597.63,208.05,54469.06,19642.76,68205.83,20408.99,4715.53,138.4,3878.88,44650.13,142608.20,18052.12,25396.68,17258.52,7207.67)
  feature_data$lower_OOR_subs = (1/3)*c(147.35,0.69,3.33,12.38,3.02,6.28,0.70,9.75,0.51,1.79,88.22,12.65,7.71,10.60,2.18,9.77,0.50,8.06,11.88,1.56,4.27,4.28,2.25,7.25,6.95,16.10,6.46,1.51,0.80,1.01,12.97,29.21,5.75,8.00,4.92,1.87)
  
  #Check order is correct
  length(rownames(data) == rownames(metadata))
  length(colnames(data) == rownames(feature_data))
  
  #Save objects
  saveRDS(data, file="./Inputs/data.rds")
  saveRDS(feature_data, file="./Inputs/feature_data.rds")
  saveRDS(metadata, file="./Inputs/metadata.rds")
}

#3.------------------------------------------------------------------------------------------------------------------#
### Calculate mean values of duplicate wells ### 

first_time = F
if (first_time == T){
  
  data = readRDS("./Inputs/data.rds")
  feature_data = readRDS("./Inputs/feature_data.rds")
  metadata = readRDS("./Inputs/metadata.rds")
  
  data1_20  = data[which((metadata$DilutionFactor == "d1/20")),]
  metadata1_20  = metadata[which((metadata$DilutionFactor == "d1/20")),]
  
  data1_20_sample_list = c()
  for (rowname in rownames(data1_20)){
    sample = strsplit(rowname, "_1/20")[[1]][1]
    if (!(sample %in% data1_20_sample_list)){
      data1_20_sample_list = c(data1_20_sample_list, sample)
    }
  }
  
  data1_20_deduplicated = as.data.frame(matrix(data=NA, nrow=length(data1_20_sample_list), ncol=36))
  rownames(data1_20_deduplicated) = data1_20_sample_list
  colnames(data1_20_deduplicated) = colnames(data1_20)
  
  metadata1_20_deduplicated = as.data.frame(matrix(data=NA, nrow=length(data1_20_sample_list), ncol=6))
  rownames(metadata1_20_deduplicated) = data1_20_sample_list
  colnames(metadata1_20_deduplicated) = c("Sample.Name","DilutionFactor","type","harvest_time","donor","batch")
  
  #Custom function to calculate mean considering OOR< and OOR>
  calculate_mean <- function(column) {
    numeric_values <- suppressWarnings(as.numeric(column))
    if (all(column == "OOR<", na.rm = TRUE)) {
      return("OOR<")
    } else if (all(column == "OOR>", na.rm = TRUE)) {
      return("OOR>")
    } else {
      mean_value <- as.numeric(mean(numeric_values, na.rm = TRUE))
      return(mean_value)
    }
  }
  
  for (sample in data1_20_sample_list){
    
    aux = data1_20[which(grepl(sample, rownames(data1_20))),]
    aux2 = metadata1_20[which(grepl(sample, rownames(data1_20))),][,c("Sample.Name","DilutionFactor","type","harvest_time","donor","batch")]
    
    data1_20_deduplicated[sample, ] = sapply(aux, calculate_mean)
    
    metadata1_20_deduplicated[sample, c("Sample.Name","DilutionFactor","type","harvest_time","donor","batch")] = aux2[1,]
  }
  
  metadata1_20 = metadata1_20_deduplicated
  metadata1_20$type = factor(metadata1_20$type)
  metadata1_20$harvest_time = factor(metadata1_20$harvest_time)
  metadata1_20$donor = factor(metadata1_20$donor)
  metadata1_20$batch = factor(metadata1_20$batch)
  
  data1_20 = data1_20_deduplicated
  
  for (i in 1:nrow(data1_20)){
    for (j in 1:ncol(data1_20)){
      if (data1_20[i,j] == "OOR<"){
        data1_20[i,j] = feature_data$lower_OOR_subs[j]
      }else if(data1_20[i,j] == "OOR>"){
        data1_20[i,j] = feature_data$upper_OOR_subs[j]
      }
    }
  }
  
  rm(data1_20_deduplicated, metadata1_20_deduplicated, aux, aux2)
  
  saveRDS(data1_20, file="./Inputs/data1_20.rds")
  saveRDS(metadata1_20, file="./Inputs/metadata1_20.rds") 
}

#4.------------------------------------------------------------------------------------------------------------------#
### Exploratory analysis (PCA) ### 

data = readRDS("./Inputs/data1_20.rds")
data = as.data.frame(lapply(data, function(x) as.numeric(x)))
data = t(data)
feature_data = readRDS("./Inputs/feature_data.rds")
metadata = readRDS("./Inputs/metadata1_20.rds")
metadata$type = as.character(metadata$type)

metadata$type[which(metadata$type == "Mac M2")] = "M2 Mac"
metadata$type[which(metadata$type == "PRS CK Storm")] = "PRS CK STORM"
metadata$type[which(metadata$type == "PRS MIO")] = "PRS MYO"
metadata$type = factor(metadata$type)

data = log2(data + 1)

dge <- DGEList(counts = data)
dge <- calcNormFactors(dge)
v=voom(dge)
data=v$E
colnames(data) = rownames(metadata)
rownames(data) = rownames(feature_data)

alt_data = data
alt_metadata = metadata

#a. PCA
pca <-prcomp(t(alt_data))
sum_pca = data.frame(summary(pca)$importance[,c(1:6)])
pca_data <- data.frame(pca$x)
pca_data=cbind(alt_metadata[,c("type","harvest_time","donor")],pca_data[,1:6])

pca_data$type=factor(pca_data$type,levels=c("M2 Mac", "ASC", "Tenocytes", "Myocytes", "Schwann",
                                            "PRS CK STORM", "PRS TENO", "PRS MYO", "PRS NEURO"))
pca_data$donor=factor(pca_data$donor,levels=c("1","2","3","4","5","6","7","8","9","10","None"))


#Figure 14A: "PC1 vs PC2 plots highlighting type of culture variable, subsetted per PRS line" 
ggplot(pca_data[which((grepl("PRSCKStorm", rownames(pca_data)) == T)|
                        (grepl("MacM2", rownames(pca_data)) == T)|
                        (grepl("ASC", rownames(pca_data)) == T)),], aes(x=PC1,y=PC2,color=type))+geom_point(size=3)+
  xlab(paste0("PC1: ",round(100*sum_pca$PC1[2],digits=2),"% variance"))+
  ylab(paste0("PC2: ",round(100*sum_pca$PC2[2],digits=2),"% variance"))+
  scale_color_manual(breaks = c("M2 Mac", "ASC", "Tenocytes", "Myocytes", "Schwann",
                                "PRS CK STORM", "PRS TENO", "PRS MYO", "PRS NEURO"),
                     values=c("#FC97FC", "#59F9F2","yellow","#8B0000","#00FF66", 
                              "#6CA8FF", "#F7DA63", "#F76D7C", "#5CD276"))


#Figure 14B: : "PC1 vs PC2 plots highlighting type of culture variable, subsetted per PRS line" 
ggplot(pca_data[which((grepl("PRSMIO", rownames(pca_data)) == T)|
                        (grepl("MacM2", rownames(pca_data)) == T)|
                        (grepl("Myocyte", rownames(pca_data)) == T)),], aes(x=PC1,y=PC2,color=type))+geom_point(size=3)+
  xlab(paste0("PC1: ",round(100*sum_pca$PC1[2],digits=2),"% variance"))+
  ylab(paste0("PC2: ",round(100*sum_pca$PC2[2],digits=2),"% variance"))+
  scale_color_manual(breaks = c("M2 Mac", "ASC", "Tenocytes", "Myocytes", "Schwann",
                                "PRS CK STORM", "PRS TENO", "PRS MYO", "PRS NEURO"),
                     values=c("#FC97FC", "#59F9F2","yellow","#8B0000","#00FF66", 
                              "#6CA8FF", "#F7DA63", "#F76D7C", "#5CD276"))


#Figure 14C: : "PC1 vs PC2 plots highlighting type of culture variable, subsetted per PRS line" 
ggplot(pca_data[which((grepl("PRSTENO", rownames(pca_data)) == T)|
                        (grepl("MacM2", rownames(pca_data)) == T)|
                        (grepl("Tenocyte", rownames(pca_data)) == T)),], aes(x=PC1,y=PC2,color=type))+geom_point(size=3)+
  xlab(paste0("PC1: ",round(100*sum_pca$PC1[2],digits=2),"% variance"))+
  ylab(paste0("PC2: ",round(100*sum_pca$PC2[2],digits=2),"% variance"))+
  scale_color_manual(breaks = c("M2 Mac", "ASC", "Tenocytes", "Myocytes", "Schwann",
                                "PRS CK STORM", "PRS TENO", "PRS MYO", "PRS NEURO"),
                     values=c("#FC97FC", "#59F9F2","yellow","#8B0000","#00FF66", 
                              "#6CA8FF", "#F7DA63", "#F76D7C", "#5CD276"))


#Figure 14D: : "PC1 vs PC2 plots highlighting type of culture variable, subsetted per PRS line" 
ggplot(pca_data[which((grepl("PRSNEURO", rownames(pca_data)) == T)|
                        (grepl("MacM2", rownames(pca_data)) == T)|
                        (grepl("Schwann", rownames(pca_data)) == T)),], aes(x=PC1,y=PC2,color=type))+geom_point(size=3)+
  xlab(paste0("PC1: ",round(100*sum_pca$PC1[2],digits=2),"% variance"))+
  ylab(paste0("PC2: ",round(100*sum_pca$PC2[2],digits=2),"% variance"))+
  scale_color_manual(breaks = c("M2 Mac", "ASC", "Tenocytes", "Myocytes", "Schwann",
                                "PRS CK STORM", "PRS TENO", "PRS MYO", "PRS NEURO"),
                     values=c("#FC97FC", "#59F9F2","yellow","#8B0000","#00FF66", 
                              "#6CA8FF", "#F7DA63", "#F76D7C", "#5CD276"))


#Mann-Whitney U tests data presented in Table 2: "Mann-Whitney U tests results on PC2 culture clusters"
pairwise.wilcox.test(pca_data$PC2, pca_data$type, p.adjust.method = "fdr")


#Figure 15A: : "PC1 vs. PC2 plot highlighting harvest time variable, subsetted per conservation method" 
ggplot(pca_data[which((grepl("FLiq", rownames(pca_data)) == T)|
                        (grepl("FLyo", rownames(pca_data)) == T)),], aes(x=PC1,y=PC2,color=harvest_time))+geom_point(size=3)+
  xlab(paste0("PC1: ",round(100*sum_pca$PC1[2],digits=2),"% variance"))+
  ylab(paste0("PC2: ",round(100*sum_pca$PC2[2],digits=2),"% variance"))+
  scale_color_manual(breaks = c("F. Liquid", "F. Lyophilized", "SN1", "SN3", "SN5"),
                     values=c("#E76BF3", "#A3A500","#F8766D","#00BF7D","#00B0F6"))


#Figure 15B: : "PC1 vs. PC2 plot highlighting harvest time variable, subsetted per intermediate supernatant collection points" 
ggplot(pca_data[which((grepl("SN1", rownames(pca_data)) == T)|
                        (grepl("SN3", rownames(pca_data)) == T)|
                        (grepl("SN5", rownames(pca_data)) == T)),], aes(x=PC1,y=PC2,color=harvest_time))+geom_point(size=3)+
  xlab(paste0("PC1: ",round(100*sum_pca$PC1[2],digits=2),"% variance"))+
  ylab(paste0("PC2: ",round(100*sum_pca$PC2[2],digits=2),"% variance"))+
  scale_color_manual(breaks = c("F. Liquid", "F. Lyophilized", "SN1", "SN3", "SN5"),
                     values=c("#E76BF3", "#A3A500","#F8766D","#00BF7D","#00B0F6"))

#5.------------------------------------------------------------------------------------------------------------------#
### Differential expression analysis (limma) on culture type effects ### 

count_DEGs=function(table,FC_name="logFC",name="adj.P.Val",threshold=0.05){
  down=-length(which(table[,name]<threshold & table[,FC_name]<=-0.5))
  up=length(which(table[,name]<threshold & table[,FC_name]>=0.5))
  return(c(down,up))
}

id_DEGs=function(table,FC_name="logFC",name="adj.P.Val",threshold=0.05){
  down=rownames(table[which(table[,name]<threshold & table[,FC_name]<=-0.5),])
  up=rownames(table[which(table[,name]<threshold & table[,FC_name]>=0.5),])
  return(list(down,up))
}

alt_metadata = metadata
alt_data = data

alt_metadata$type = relevel(factor(alt_metadata$type),ref="M2 Mac")
alt_metadata$donor = relevel(factor(alt_metadata$donor),ref="1")
alt_metadata$harvest_time = relevel(factor(alt_metadata$harvest_time),ref="SN1")

design = model.matrix(~type+harvest_time+donor+batch, data=alt_metadata)
fit=lmFit(alt_data, design=design)
fit=eBayes(fit, robust = TRUE)

#DEGs in type variable
ASC_vs_MacM2 = topTable(fit, "typeASC", number=nrow(alt_data))
Tenocytes_vs_MacM2 = topTable(fit, "typeTenocytes", number=nrow(alt_data))
Myocytes_vs_MacM2 = topTable(fit, "typeMyocytes", number=nrow(alt_data))
Schwann_vs_MacM2 = topTable(fit, "typeSchwann", number=nrow(alt_data))
PRSCKStorm_vs_MacM2 = topTable(fit, "typePRS CK STORM", number=nrow(alt_data))
PRSTENO_vs_MacM2 = topTable(fit, "typePRS TENO", number=nrow(alt_data))
PRSMYO_vs_MacM2 = topTable(fit, "typePRS MYO", number=nrow(alt_data))
PRSNEURO_vs_MacM2 = topTable(fit, "typePRS NEURO", number=nrow(alt_data))

DEGs_type_vs_MacM2=data.frame(
  ASC=count_DEGs(table=ASC_vs_MacM2,FC_name="logFC",name="adj.P.Val",threshold=0.05),
  Tenocytes=count_DEGs(table=Tenocytes_vs_MacM2,FC_name="logFC",name="adj.P.Val",threshold=0.05),
  Myocytes=count_DEGs(table=Myocytes_vs_MacM2,FC_name="logFC",name="adj.P.Val",threshold=0.05),
  Schwann=count_DEGs(table=Schwann_vs_MacM2,FC_name="logFC",name="adj.P.Val",threshold=0.05),
  PRSCKSTORM=count_DEGs(table=PRSCKStorm_vs_MacM2,FC_name="logFC",name="adj.P.Val",threshold=0.05),
  PRSTENO=count_DEGs(table=PRSTENO_vs_MacM2,FC_name="logFC",name="adj.P.Val",threshold=0.05),
  PRSMYO=count_DEGs(table=PRSMYO_vs_MacM2,FC_name="logFC",name="adj.P.Val",threshold=0.05),
  PRSNEURO=count_DEGs(table=PRSNEURO_vs_MacM2,FC_name="logFC",name="adj.P.Val",threshold=0.05))
rownames(DEGs_type_vs_MacM2)=c("Down-regulated in type vs. Mac M2", "Up-regulated in type vs. Mac M2")

#Data presented in Table 3: "DE proteins across culture secretomes (FDR < 5% & |logFC| > 0.5)"
id_DEGs(table=ASC_vs_MacM2,FC_name="logFC",name="adj.P.Val",threshold=0.05)
id_DEGs(table=Tenocytes_vs_MacM2,FC_name="logFC",name="adj.P.Val",threshold=0.05)
id_DEGs(table=Myocytes_vs_MacM2,FC_name="logFC",name="adj.P.Val",threshold=0.05)
id_DEGs(table=Schwann_vs_MacM2,FC_name="logFC",name="adj.P.Val",threshold=0.05)
id_DEGs(table=PRSCKStorm_vs_MacM2,FC_name="logFC",name="adj.P.Val",threshold=0.05)
id_DEGs(table=PRSTENO_vs_MacM2,FC_name="logFC",name="adj.P.Val",threshold=0.05)
id_DEGs(table=PRSMYO_vs_MacM2,FC_name="logFC",name="adj.P.Val",threshold=0.05)
id_DEGs(table=PRSNEURO_vs_MacM2,FC_name="logFC",name="adj.P.Val",threshold=0.05)

#Sinergy study 
levels(alt_metadata$type) = make.names(levels(alt_metadata$type))
levels(alt_metadata$harvest_time) = make.names(levels(alt_metadata$harvest_time))

design = model.matrix(~0+type+harvest_time+donor+batch, data=alt_metadata)
fit=lmFit(alt_data, design=design)
fit=eBayes(fit, robust = TRUE)

contrast_matrix <- makeContrasts(
  PRSCKStorm_vs_Sum = `typePRS.CK.STORM` - (typeASC*1/6 + `typeM2.Mac`*5/6),
    PRSMYO_vs_Sum = `typePRS.MYO` - (typeMyocytes*1/6 + `typeM2.Mac`*5/6),
    PRSTENO_vs_Sum = `typePRS.TENO` - (typeTenocytes*1/6 + `typeM2.Mac`*5/6),
    PRSNEURO_vs_Sum = `typePRS.NEURO` - (typeSchwann*1/6 + `typeM2.Mac`*5/6),
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
PRSCKStorm_sinergy <- topTable(fit2, coef="PRSCKStorm_vs_Sum",number=nrow(alt_data))
PRSMYO_sinergy <- topTable(fit2, coef="PRSMYO_vs_Sum",number=nrow(alt_data))
PRSTENO_sinergy <- topTable(fit2, coef="PRSTENO_vs_Sum",number=nrow(alt_data))
PRSNEURO_sinergy <- topTable(fit2, coef="PRSNEURO_vs_Sum",number=nrow(alt_data))

#Data presented in Table 4: "DE proteins in cocultures vs. weighted sum of constituent monocultures (FDR < 5% & |logFC| > 0.5)"
id_DEGs(table=PRSCKStorm_sinergy,FC_name="logFC",name="adj.P.Val",threshold=0.05)
id_DEGs(table=PRSMYO_sinergy,FC_name="logFC",name="adj.P.Val",threshold=0.05)
id_DEGs(table=PRSTENO_sinergy,FC_name="logFC",name="adj.P.Val",threshold=0.05)
id_DEGs(table=PRSNEURO_sinergy,FC_name="logFC",name="adj.P.Val",threshold=0.05)

#Specifity
alt_metadata$type = relevel(factor(alt_metadata$type),ref="PRS.CK.STORM")
alt_metadata$donor = relevel(factor(alt_metadata$donor),ref="1")
alt_metadata$harvest_time = relevel(factor(alt_metadata$harvest_time),ref="SN1")

design = model.matrix(~type+harvest_time+donor+batch, data=alt_metadata)
fit=lmFit(alt_data, design=design)
fit=eBayes(fit, robust = TRUE)

PRSTENO_vs_PRSCKStorm = topTable(fit, "typePRS.TENO", number=nrow(alt_data))
PRSMIO_vs_PRSCKStorm = topTable(fit, "typePRS.MYO", number=nrow(alt_data))
PRSNEURO_vs_PRSCKStorm = topTable(fit, "typePRS.NEURO", number=nrow(alt_data))

DEGs_PRSs_vs_PRSCKStorm=data.frame(
  PRSTENO=count_DEGs(table=PRSTENO_vs_PRSCKStorm,FC_name="logFC",name="adj.P.Val",threshold=0.05),
  PRSMYO=count_DEGs(table=PRSMIO_vs_PRSCKStorm,FC_name="logFC",name="adj.P.Val",threshold=0.05),
  PRSNEURO=count_DEGs(table=PRSNEURO_vs_PRSCKStorm,FC_name="logFC",name="adj.P.Val",threshold=0.05))
rownames(DEGs_PRSs_vs_PRSCKStorm)=c("Down-regulated in PRS vs. PRS CK Storm", "Up-regulated in PRS vs. PRS CK Storm")

alt_metadata$type = relevel(factor(alt_metadata$type),ref="PRS.MYO")

design = model.matrix(~type+harvest_time+donor+batch, data=alt_metadata)
fit=lmFit(alt_data, design=design)
fit=eBayes(fit, robust = TRUE)

PRSTENO_vs_PRSMYO = topTable(fit, "typePRS.TENO", number=nrow(alt_data))
PRSNEURO_vs_PRSMYO = topTable(fit, "typePRS.NEURO", number=nrow(alt_data))

alt_metadata$type = relevel(factor(alt_metadata$type),ref="PRS.TENO")

design = model.matrix(~type+harvest_time+donor+batch, data=alt_metadata)
fit=lmFit(alt_data, design=design)
fit=eBayes(fit, robust = TRUE)

PRSNEURO_vs_PRSTENO = topTable(fit, "typePRS.NEURO", number=nrow(alt_data))

#Data presented in Table 5: "DE proteins across PRS lines (FDR < 5% & |logFC| > 0.5)"
id_DEGs(table=PRSTENO_vs_PRSCKStorm,FC_name="logFC",name="adj.P.Val",threshold=0.05)
id_DEGs(table=PRSMIO_vs_PRSCKStorm,FC_name="logFC",name="adj.P.Val",threshold=0.05)
id_DEGs(table=PRSNEURO_vs_PRSCKStorm,FC_name="logFC",name="adj.P.Val",threshold=0.05)
id_DEGs(table=PRSTENO_vs_PRSMYO,FC_name="logFC",name="adj.P.Val",threshold=0.05)
id_DEGs(table=PRSNEURO_vs_PRSMYO,FC_name="logFC",name="adj.P.Val",threshold=0.05)
id_DEGs(table=PRSNEURO_vs_PRSTENO,FC_name="logFC",name="adj.P.Val",threshold=0.05)

#6.------------------------------------------------------------------------------------------------------------------#
### Differential expression analysis (limma) on harvest time and conservation methods effects ### 

alt_metadata = metadata
alt_data = data

alt_metadata$type = relevel(factor(alt_metadata$type),ref="M2 Mac")
alt_metadata$donor = relevel(factor(alt_metadata$donor),ref="1")
alt_metadata$harvest_time = relevel(factor(alt_metadata$harvest_time),ref="SN1")

design = model.matrix(~type+harvest_time+donor+batch, data=alt_metadata)
fit=lmFit(alt_data, design=design)
fit=eBayes(fit, robust = TRUE)

SN3_vs_SN1 = topTable(fit, "harvest_timeSN3", number=nrow(alt_data))

alt_metadata$harvest_time = relevel(factor(alt_metadata$harvest_time),ref="SN3")

design = model.matrix(~type+harvest_time+donor+batch, data=alt_metadata)
fit=lmFit(alt_data, design=design)
fit=eBayes(fit, robust = TRUE)

SN5_vs_SN3 = topTable(fit, "harvest_timeSN5", number=nrow(alt_data))

alt_metadata$harvest_time = relevel(factor(alt_metadata$harvest_time),ref="F. Lyophilized")

design = model.matrix(~type+harvest_time+donor+batch, data=alt_metadata)
fit=lmFit(alt_data, design=design)
fit=eBayes(fit, robust = TRUE)

FLiq_vs_FLyo = topTable(fit, "harvest_timeF. Liquid", number=nrow(alt_data))

#Data presented in Table 6: "DE proteins attending to comparisons regarding the harvest time variable (FDR < 5% & |logFC| > 0.5)"
id_DEGs(table=SN3_vs_SN1,FC_name="logFC",name="adj.P.Val",threshold=0.05)
id_DEGs(table=SN5_vs_SN3,FC_name="logFC",name="adj.P.Val",threshold=0.05)
id_DEGs(table=FLiq_vs_FLyo,FC_name="logFC",name="adj.P.Val",threshold=0.05)
