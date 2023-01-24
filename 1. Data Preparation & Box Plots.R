### Data Preparation & Box Plots
library(dplyr)
library(tidyr)

#Input files
allgenecounts=read.csv("data_mrna_seq_v2_rsem.csv", header=TRUE)
clinical=read.csv("PAM50_ER.csv",header=TRUE,row.names=1) 
#Remove rows with empty gene symbols
allgenecounts=subset(allgenecounts, allgenecounts$Hugo_Symbol!="")
#Average duplicate genes
n_occur=data.frame(table(allgenecounts$Hugo_Symbol))
duplicates=n_occur[n_occur$Freq > 1,]
dup_genes=duplicates$Var1
duplicates=subset(allgenecounts, allgenecounts$Hugo_Symbol %in% dup_genes)
duplicates=arrange(duplicates, duplicates$Hugo_Symbol)
dup_genes=duplicates$Hugo_Symbol
duplicates=as.data.frame(t(duplicates))
duplicates=duplicates[-c(1,2),]
dup_genes
ncol(duplicates)
duplicates=duplicates%>%mutate_all(as.numeric)
duplicates$med1=colMeans(duplicates[c(1,2)],na.rm=TRUE)
duplicates$med2=colMeans(duplicates[c(3,4)],na.rm=TRUE)
duplicates$med3=colMeans(duplicates[c(5,6)],na.rm=TRUE)
duplicates$med4=colMeans(duplicates[c(7,8)],na.rm=TRUE)
duplicates$med5=colMeans(duplicates[c(9,10)],na.rm=TRUE)
duplicates$med6=colMeans(duplicates[c(11,12)],na.rm=TRUE)
duplicates$med7=colMeans(duplicates[c(13,14)],na.rm=TRUE)
duplicates=duplicates[-c(1:14)]
names(duplicates)= unique(dup_genes)
duplicates=as.data.frame(t(duplicates))
uniques=subset(allgenecounts, !(allgenecounts$Hugo_Symbol %in% dup_genes))
rownames(uniques)=uniques$Hugo_Symbol
uniques=uniques[-c(1,2)]
allgenecounts=rbind(duplicates,uniques)
names(allgenecounts) = gsub("[.]", "-", names(allgenecounts))
names(allgenecounts) = gsub("-01", "", names(allgenecounts))


#Subset RET values for all patients 
RET=subset(allgenecounts, rownames(allgenecounts) == "RET")
RET=as.data.frame(t(RET))
RET=RET[rownames(clinical),,drop=FALSE]
RET=cbind(RET,clinical)

#PAM50 plot
library(ggplot2)
library(ggpubr)
plot=ggplot(RET, aes(x = PAM50, y = RET)) + 
  geom_boxplot(fill="grey", outlier.shape = NA) + theme(legend.position = "none") + 
  xlab("PAM50 Subtype") + 
  ylab("RET")+ 
  theme_classic()+ geom_signif(comparisons = list(c("Basal", "LumA"), c("Basal", "LumB"), c("HER2E", "LumA"), c("HER2E", "LumB"), c("LumA", "LumB"), c("LumA", "normal-like"), c("LumB", "normal-like")), map_signif_level=TRUE,step_increase = 0.03, tip_length = 0.01, test="t.test")+ theme_classic()+ theme(axis.text.x=element_text(size=12))
plot
plot+ylim(c(0,4000))
bact <- compare_means(RET~PAM50, data = RET,method = 't.test')
bact

#LFC Changes for PAM50
library(limma)
data=RET[c("RET")]
data$RET=log2(data$RET)
matx=data.frame(t(apply(data, 1, as.numeric)))
names(matx)=names(data)
matx=as.matrix(matx)
RET2=RET
RET2$PAM50=replace(RET2$PAM50, RET2$PAM50=="normal-like", "normal")
design=model.matrix(~0+PAM50, RET2)
contrasts = makeContrasts(PAM50LumB-PAM50LumA, levels = design) #modify groups according to contrast interests, eg. ERLumA-ERHER2E or ERLumB-ERnormal-like, etc.
fit = lmFit(matx, design)
fit2 = contrasts.fit(fit, contrasts)
fit2 = eBayes(fit2)
top=topTable(fit2, number=Inf)

#ER plot
df=RET
RET_ER=subset(df, df$ER!="[Not Evaluated]")
RET_ER=subset(RET_ER, RET_ER$ER!="Indeterminate")
RET_ER=subset(RET_ER, RET_ER$ER!="NA")
ggplot(RG_all_ER, aes(x = ER, y = RET)) + 
  geom_boxplot(fill="grey") + theme(legend.position = "none") +
  xlab("ER Status") + 
  ylab("RET")+ 
  theme_classic() +  geom_signif(comparisons = list(c("Positive", "Negative")), map_signif_level=TRUE, step_increase = 0.05, tip_length = 0.01, test="t.test") 
bact <- compare_means(RET~ER, data = RG_all,method = 't.test')
bact

#LFC Changes for ER
data=RET_ER[c("RET")]
data$RET=log2(data$RET)
matx=data.frame(t(apply(data, 1, as.numeric)))
names(matx)=names(data)
matx=as.matrix(matx)
design=model.matrix(~0+ER, RET_ER)
contrasts = makeContrasts(ERPositive-ERNegative, levels = design)
fit = lmFit(matx, design)
fit2 = contrasts.fit(fit, contrasts)
fit2 = eBayes(fit2)
top=topTable(fit2, number=Inf)


#-----------R&G


#Subset RET and GDNF values for all patients and prepare data
RET = subset(allgenecounts, rownames(allgenecounts) == "RET")
GDNF = subset(allgenecounts, rownames(allgenecounts) == "ATF2")
rownames(GDNF)="GDNF"
RG_all=data.frame(t(rbind(RET,GDNF)))
RG_all=RG_all[rownames(clinical),,drop=FALSE]
RG_all=cbind(RG_all,clinical)

##PAM50 plot
library(ggplot2)
library(ggpubr)
plot=ggplot(RG_all, aes(x = PAM50, y = RET)) + 
  geom_boxplot(fill="grey", outlier.shape = NA) + theme(legend.position = "none") + 
  xlab("PAM50 Subtype") + 
  ylab("RET")+ 
  theme_classic()+ geom_signif(comparisons = list(c("Basal", "LumA"), c("Basal", "LumB"), c("HER2E", "LumA"), c("HER2E", "LumB"), c("LumA", "LumB"), c("LumA", "normal-like"), c("LumB", "normal-like")), map_signif_level=TRUE,step_increase = 0.03, tip_length = 0.01, test="t.test")+ theme_classic()+ theme(axis.text.x=element_text(size=12))
plot
plot+ylim(c(0,4000))
bact <- compare_means(RET~PAM50, data = RG_all,method = 't.test')
bact

#LFC Changes for PAM50
library(limma)
data=RG_all[c("RET")]
data$RET=log2(data$RET)
matx=data.frame(t(apply(data, 1, as.numeric)))
names(matx)=names(data)
matx=as.matrix(matx)
RG_all_2=RG_all
RG_all_2$PAM50=replace(RG_all$PAM50, RG_all$PAM50=="normal-like", "normal")
design=model.matrix(~0+PAM50, RG_all_2)
contrasts = makeContrasts(PAM50LumB-PAM50LumA, levels = design) #modify groups according to contrast interests, eg. ERLumA-ERHER2E or ERLumB-ERnormal-like, etc.
fit = lmFit(matx, design)
fit2 = contrasts.fit(fit, contrasts)
fit2 = eBayes(fit2)
top=topTable(fit2, number=Inf)

##ER plot
RG_all_ER=subset(RG_all, RG_all$ER!="[Not Evaluated]")
RG_all_ER=subset(RG_all_ER, RG_all_ER$ER!="Indeterminate")
RG_all_ER=subset(RG_all_ER, RG_all_ER$ER!="NA")
ggplot(RG_all_ER, aes(x = ER, y = RET)) + 
  geom_boxplot(fill="grey") + theme(legend.position = "none") +
  xlab("ER Status") + 
  ylab("RET")+ 
  theme_classic() +  geom_signif(comparisons = list(c("Positive", "Negative")), map_signif_level=TRUE, step_increase = 0.05, tip_length = 0.01, test="t.test") 
bact <- compare_means(RET~ER, data = RG_all,method = 't.test')
bact

#LFC Changes for ER
data=RG_all_ER[c("RET")]
data$RET=log2(data$RET)
matx=data.frame(t(apply(data, 1, as.numeric)))
names(matx)=names(data)
matx=as.matrix(matx)
design=model.matrix(~0+ER, RG_all_ER)
contrasts = makeContrasts(ERPositive-ERNegative, levels = design)
fit = lmFit(matx, design)
fit2 = contrasts.fit(fit, contrasts)
fit2 = eBayes(fit2)
top=topTable(fit2, number=Inf)


