#RET DGE Analysis
library(fabricatr)
library(DESeq2)
library(dplyr)
library(ComplexHeatmap)
library(EnhancedVolcano)

#RET group assignment
RET=subset(ER, rownames(ER)=="RET")
RET2=data.frame(t(RET))
quantile(RET2$RET)
splitRET <- split_quantile(x = RET2$RET, type = 4)
RET2$quartile <- splitRET
ord=RET2
ord=ord[order(ord$RET),]
groups <- ord[c(2)]
GRP <- ifelse(groups$quartile == "1", "Lower75", ifelse(groups$quartile == "2", "Lower75", ifelse(groups$quartile == "3", "Lower75", "Upper25")))
groups$categ <- GRP
groups <- groups[c(2)]
Lower75 <- subset(groups, groups$categ=="Lower75")
Upper25 <- subset(groups, groups$categ=="Upper25")
patients <- rbind(Lower75, Upper25)

#DESeq2
genes <- subset(ER, rownames(ER)!="RET")
completecounts <- na.omit(genes)
cc=completecounts %>% select(rownames(patients))
roundedcounts=round(cc)
groups$categ=factor(groups$categ)
dds <- DESeqDataSetFromMatrix(countData = roundedcounts,
                              colData = groups,
                              design = ~ categ)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$categ <- relevel(dds$categ, ref = "Lower75")
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


#Group heatmap color annotations
#Patient groups
RETgrps <- patients
RETgrps
RETgrps<- ifelse(RETgrps == "Upper25" , "R4", "Non-R4") 
RETgrps
colnames(RETgrps)=NULL
rownames(RETgrps)=NULL
RETgrps
RETT1 <- as.matrix(RETgrps)
ann <- data.frame(RETT1)
colnames(ann) <- c('RET Expression Group')
colours <- list('RET Expression Group' = c("R4" = "darkolivegreen4", "Non-R4" =  "darkolivegreen2"))
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
  names(keyvals)[keyvals == 'red'] <- 'Up-regulated: 191'
  names(keyvals)[keyvals == 'grey'] <- 'Intermediate'
  names(keyvals)[keyvals == 'blue'] <- 'Down-regulated: 445'
                   
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
                boxedLabels=FALSE,
                selectLab = rownames(sc)[which(names(keyvals) %in% c('Up-regulated: 191', 'Down-regulated: 445'))],
                ylab = bquote(~-Log[10] ~ FDR),
                drawConnectors = FALSE,
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