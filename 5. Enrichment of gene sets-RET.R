##Enrichment of gene sets
library(clusterProfiler)
library(ggplot2)

#Load reference datasets
GO_file = "MSigDB_BREAST.gmt"
data=gmtPathways(GO_file)
data=data.frame(data)

#Select significant up/down genes
upregulated_genes=rownames(subset(DESeq_output2, DESeq_output2$log2FoldChange>1))
downregulated_genes=rownames(subset(DESeq_output2, DESeq_output2$log2FoldChange < -1))

#Enrichment analysis
genes <- upregulated_genes
egmt <- enricher(genes, pvalueCutoff=0.05, qvalueCutoff=0.05, pAdjustMethod = "BH", TERM2GENE=data)
file=data.frame(egmt)
rownames(file)=NULL
file2 <- file[-c(2)]
file2$GR <- sub("\\/.*","\\/", file2$BgRatio, perl=TRUE)
file2$GR = stringr::str_replace(file2$GR, "/" , "")
file2$GR = as.numeric(file2$Count)/as.numeric(file2$GR)
up = file2

genes <- downregulated_genes
egmt <- enricher(genes, pvalueCutoff=0.05, qvalueCutoff=0.05, pAdjustMethod = "BH", TERM2GENE=data)
file=data.frame(egmt)
rownames(file)=NULL
file2 <- file[-c(2)]
file2$GR <- sub("\\/.*","\\/", file2$BgRatio, perl=TRUE)
file2$GR = stringr::str_replace(file2$GR, "/" , "")
file2$GR = as.numeric(file2$Count)/as.numeric(file2$GR)
down = file2

#SUBSET BASED ON GENE INTERESTS
setz=up #run the below analysis for "down" then replace with "up" and re-run

setz=subset(setz,setz$p.adjust<0.05)
setz1=rbind(setz[grepl("E2|ESTROGEN|ESTRADIOL|ESR1|ESR2|_ESR_|ESRRA", setz$ID),]) 
setz1=setz1[order(setz1$GR),]
setz1=unique(setz1)
setz2=rbind(setz[grepl("TAMOXIFEN|RESISTANCE|SERM|FULVESTRANT|ENDOCRINE", setz$ID),])
setz2=setz2[order(setz2$GR),]
setz2=unique(setz2)
setz3=rbind(setz[grepl("MAPK|ERK1|ERK2|MAP2K|MEK", setz$ID),])#MAPK pathways
setz3=setz3[order(setz3$GR),]
setz3=unique(setz3)
setz4=rbind(setz[grepl("PI3K|AKT|MTOR", setz$ID),])#PI3K AKT pathways
setz4=setz4[order(setz4$GR),]
setz4=unique(setz4)
setz5=rbind(setz[grepl("KRAS", setz$ID),])#KRAS pathways
setz5=setz5[order(setz5$GR),]
setz5=unique(setz5)
setz6=rbind(setz[grepl("EGFR", setz$ID),])#EGFR pathways
setz6=setz6[order(setz6$GR),]
setz6=unique(setz6)
setz7=rbind(setz[grepl("VEGF", setz$ID),])#VEGF pathways
setz7=setz7[order(setz7$GR),]
setz7=unique(setz7)
setz8=rbind(setz[grepl("RET", setz$ID),])#RET pathways
setz8=setz8[order(setz8$GR),]
setz8=unique(setz8)
setz9=rbind(setz[grepl("ERBB", setz$ID),]) #ERRB2 pathways
setz9=setz9[order(setz9$GR),]
setz9=unique(setz9)
setz10=rbind(setz[grepl("MYC", setz$ID),])#MYC pathways
setz10=setz10[order(setz10$GR),]
setz10=unique(setz10)
setz11=rbind(setz[grepl("P53|APOPTOSIS|CYCLE|G2M|PHASE|DEATH|CYCLIN|CCND", setz$ID),]) # p53 AND CYCLE PATHWAYS
setz11=setz11[order(setz11$GR),]
setz11=unique(setz11)
setz12=setz[grepl("FGF", setz$ID),]
setz12=unique(setz12)
setz13=setz[grepl("LUMINAL|BASAL|NORMAL|SUBTYPE|INVASIVE", setz$ID),] ##phenotype
setz13=setz13[order(setz13$GR),]
setz13=unique(setz13)
setz14=setz[grepl("JAK|STAT|JNK|JUN", setz$ID),]
setz14=unique(setz14)
setz15=setz[grepl("TNF|NFKB|IKBKB", setz$ID),]
setz15=unique(setz15)
setz16=setz[grepl("WNT|CATENIN", setz$ID),]
setz16=unique(setz16)
setz1$category="Estrogen & ESR1 sets"
setz2$category="AE"
setz3$category="MAPK sets"
setz4$category="AKT"
setz5$category="KRAS sets"
setz6$category="EGFR"
setz7$category="VEGF sets"
setz8$category="RET sets"
setz9$category="ERBB2"
setz10$category="MYC"
setz11$category="Cell cycle sets"
setz12$category="FGFR"
setz13$category="Cell phenotype sets"
setz14$category="JAK/STAT sets"
setz15$category="NFKB"
setz16$category="WNT sets"

dffz_down=rbind(setz1,setz2,setz3,setz4,setz5,setz6,setz7,setz9,setz10,setz11,setz12,setz13,setz14,setz15,setz16) #run after "down" as setz
dffz_down$group="down"

dffz_up=rbind(setz1,setz2,setz3,setz4,setz5,setz6,setz7,setz9,setz10,setz11,setz12,setz13,setz14,setz15,setz16) #run after "up" as setz
dffz_up$group="up"

#modify directionality of plotting by allocating negative GR for downregulated sets and positive GR for upregulated sets
newdown=dffz_down
newdown$GRmod=-dffz_down$GR
newup=dffz_up
newup$GRmod=dffz_up$GR
dffz=rbind(newdown,newup)

#assign factors for grouped plotting
dffz$Enr_f = factor(dffz$category, levels=c('Estrogen & ESR1 sets','AE','MAPK sets','AKT','EGFR','VEGF sets','ERBB2','KRAS sets','MYC','Cell cycle sets','Cell phenotype sets','JNK sets', 'NFKB','WNT sets', 'FGFR'))
dffz$group=ifelse(dffz$group=="up","UPREGULATED","DOWNREGULATED")
dffz$Enr_group=factor(dffz$group, levels=c('DOWNREGULATED','UPREGULATED'))

#plot_everything
ggplot(data = dffz, aes(x = GRmod, y = ID,fill=p.adjust)) + geom_bar(stat="identity")+facet_grid(Enr_f~., scales = "free", space = "free", switch="y") + theme(axis.title.y=element_blank(),legend.position = "bottom",text=element_text(size=13))+xlab("Gene Ratio by Set")

#filter out any sets that pop-up unrelated to breast cancer
dffz=dffz[!grepl("SYNAPTIC|ARSENIC|GOBP_ENDOCRINE|GOBP_REGULATION_OF_ENDOCRINE|NEUROENDOCRINE|NFE2L2|PROSTATE|LUNG|KIDNEY|E2F3|MEISSNER|HEPATOBLASTOMA|ESOPHAGEAL|NEUROTRANSMITTER|INFANT|GLIOBLASTOMA|PLATEAU|BRUINS|GTPASE|DOXORUBICIN|MARTINEZ|INTESTINE|LTE2|PANCREAS|FETAL|DEPRIVATION|AGING|KRAS.DF|DOCETACEL|E2F1|HP_|ANATOMICAL|NCX|METHYLATED|YOSHIMURA|MTOR_UP.|AKT_UP|VEGF|P53_DN|HAMAI|KRAS.BREAST|PEREZ|JNK_DN|JAK2_DN|UP.N4.|YANG|VANTVEER",dffz$ID),]

#plot specific groups/categories by specifying keyword through grepl
ggplot(data = dffz[grepl("Estrogen|AE|MAPK|AKT",dffz$Enr_f),], aes(x = GRmod, y = ID,fill=p.adjust)) + geom_bar(stat="identity")+facet_grid(Enr_f~Enr_group, scales = "free", space = "free", switch="y") + theme(axis.title.y=element_blank(),legend.position = "right",text=element_text(size=13))+xlab("Gene Ratio by Set")+ scale_fill_continuous(name = "FDR")


