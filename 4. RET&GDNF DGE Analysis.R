#RET DGE Analysis
library(fabricatr)
library(DESeq2)
library(dplyr)
library(ComplexHeatmap)
library(EnhancedVolcano)

#RET&GDNF co-expression group assignment
RET=subset(ER, rownames(ER) == "RET")
RET=as.data.frame(t(RET))
RET=RET[rownames(ERclinical),,drop=FALSE]
RET=cbind(RET,ERclinical)
GDNF=subset(ER, rownames(ER) == "ATF2")
GDNF=as.data.frame(t(GDNF))
names(GDNF)="GDNF"
GDNF=GDNF[rownames(ERclinical),,drop=FALSE]
df=cbind(RET,GDNF)
splitRET <- split_quantile(x = df$RET, type = 4)
df$groupRET <- splitRET
splitGDNF <- split_quantile(x = df$GDNF, type = 4)
df$groupGDNF <- splitGDNF
ord=df[order(df$RET),]
groups <- ord
GRP <- ifelse(groups$groupRET == "4" & groups$groupGDNF == "4", "Above", "Below")
groups$categ <- GRP
groups=groups[c(7)]
Lower75 <- subset(groups, groups$categ=="Below")
Upper25 <- subset(groups, groups$categ=="Above")
patients <- rbind(Lower75, Upper25)

#DESeq2
genes <- subset(ER, rownames(ER)!="RET")
genes <- subset(genes, rownames(genes)!="ATF2")
completecounts <- na.omit(genes)
cc=completecounts %>% select(rownames(patients))
roundedcounts=round(cc)
patients$categ=factor(patients$categ)
dds <- DESeqDataSetFromMatrix(countData = roundedcounts,
                              colData = patients,
                              design = ~ categ)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$categ <- relevel(dds$categ, ref = "Below")
dds <- DESeq(dds)
res <- results(dds)
DESeq_output <- data.frame(res)
DESeq_output2 <- na.omit(DESeq_output)


#Selecting significantly regulated genes
df <- DESeq_output2
df.top <- df[ (df$baseMean > 50) & (df$padj < 0.05) & (abs(df$log2FoldChange) > 1),]
df.top <- df.top[order(df.top$log2FoldChange, decreasing = TRUE),]

#Normalizing matrix for heatmap visualization
rlog_out <- vst(dds, blind=FALSE) #get normalized count data from dds object
coldata <- groups
mat<-assay(rlog_out)[rownames(df.top), rownames(coldata)] #sig genes x samples
colnames(mat) <- rownames(coldata)
base_mean <- rowMeans(mat)
mat.scaled <- t(apply(mat, 1, scale)) #center and scale each column (Z-score) then transpose
colnames(mat.scaled)<-colnames(mat)
num_keep <- 40
Genez <- rbind(tail(df.top[order(df.top$log2FoldChange),], num_keep), (head(df.top[order(df.top$log2FoldChange),], num_keep)))
Genez <- Genez[order(Genez$log2FoldChange, decreasing = TRUE),]
rows_keep <- rownames(Genez)
matrixx <- mat.scaled[rows_keep,]

#Save gene signature 
num_keep <- 40
up <- rownames(tail(df.top[order(df.top$log2FoldChange),], num_keep))
down <- rownames(head(df.top[order(df.top$log2FoldChange),], num_keep))
library(fgsea)
new = list("UPREG_fullR4G4" = up,
           "DOWNREG_fullR4G4" = down)
writeGmtPathways(new, "HIGHRETGDNF_top40.gmt")


#Group heatmap color annotations
#Patient groups
RETgrps <- patients
RETgrps
RETgrps<- ifelse(RETgrps == "Above" , "High RET/GDNF", "Non-high RET/GDNF") 
RETgrps
colnames(RETgrps)=NULL
rownames(RETgrps)=NULL
RETgrps
RETT1 <- as.matrix(RETgrps)
ann <- data.frame(RETT1)
colnames(ann) <- c('RET & GDNF Co-Expression Group')
colours <- list('RET & GDNF Co-Expression Group' = c("High RET/GDNF" = "darkolivegreen4", "Non-high RET/GDNF" =  "darkolivegreen2"))
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_width = unit(c(1, 2), 'cm'),
                            gap = unit(0.5, 'mm'),
                            show_annotation_name = FALSE)

#Gene groups
Genez$regulation <- ifelse(Genez$log2FoldChange <0, "Downregulated", "Upregulated")
ReguLevelz <- subset(Genez, select = c(regulation))
colnames(ReguLevelz)=NULL
rownames(ReguLevelz)=NULL
ReguLevelz
REG <- as.matrix(ReguLevelz)
REG=data.frame(REG)
colnames(REG) <- c('Top DEG')
colours2 <- list('Top DEG' = c("Downregulated" = "cornflowerblue", "Upregulated" = "indianred"))
genecol <- HeatmapAnnotation(df = REG,
                             which = 'row',
                             col = colours2,
                             annotation_width = unit(c(1, 2), 'cm'),
                             gap = unit(0.5, 'mm'),
                             show_annotation_name = FALSE)
Heatmap(matrixx, 
        cluster_rows = FALSE,
        name="Z-score",
        cluster_columns = FALSE,
        row_names_side = "right",
        row_names_gp = gpar(fontsize = 10),
        column_title = "(n=795)",
        column_title_gp = gpar(fontsize=10),
        column_title_side = "bottom",
        show_column_names = FALSE,
        show_row_names = TRUE,
        row_names_max_width = unit(10,"in"),
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        column_order = order(as.numeric(gsub("column", "", rownames(ord)))),
        row_order = order(as.numeric(gsub("row", "", rownames(Genez)))),
        top_annotation = colAnn,
        right_annotation = genecol,
        height = unit(300, "mm"),
        width = unit(150, "mm"))


#Volcano plot
sc <- DESeq_output2
keyvals <- ifelse(sc$log2FoldChange < -1 & sc$padj < 0.05, 'blue', ifelse(sc$log2FoldChange > 1 & sc$padj < 0.05, 'red','grey'))
nrow(sc[sc$log2FoldChange > 1 & sc$padj < 0.05, ])
nrow(sc[sc$log2FoldChange < -1 & sc$padj < 0.05, ])

keyvals[is.na(keyvals)] <- 'black'
  names(keyvals)[keyvals == 'red'] <- 'Up-regulated: 195'
  names(keyvals)[keyvals == 'grey'] <- 'Intermediate'
  names(keyvals)[keyvals == 'blue'] <- 'Down-regulated: 945'

  
genez=c("FGF10,UGT2B4,UGT2B11,TAOK1,MMP9,MUC16,FOXC1,HSD17B1,KRT13,LEP,TRIM29,KLK11,SFRP1,NPY5R,SPHK1,SLPI,AREG,NTRK2,VTN")
genez=unlist(strsplit(genez,","))
  
genelist=subset(sc, rownames(sc) %in% genez)
   
EnhancedVolcano(sc,
                  lab = rownames(sc),
                  x = 'log2FoldChange',
                  y = 'padj',
                  title = NULL,
                  subtitle="",
                  subtitleLabSize = 12,
                  axisLabSize = 14,
                  pCutoff = 0.05,
                  FCcutoff = 1,
                  pointSize = 3,
                  labSize = 4,
                  colCustom=keyvals,
                  colAlpha = 0.4,
                  boxedLabels=TRUE,
                  selectLab = rownames(genelist),
                  ylab = bquote(~-Log[10] ~ FDR),
                  drawConnectors = TRUE,
                  gridlines.major = FALSE,
                  gridlines.minor = FALSE,
                  legendPosition = 'right',
                  captionLabSize = 7,
                  legendLabSize = 15,
                  legendIconSize = 3.0,
                  arrowheads = FALSE,
                  max.overlaps = 40,
                  widthConnectors = 1,
                  cutoffLineType = 'dashed',
                  cutoffLineCol = 'black',
                  cutoffLineWidth = 0.4,
                  colConnectors = 'black',
                  borderWidth = 0.7,
                  borderColour = 'black')
  
  