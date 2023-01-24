##Endocrine therapy treated tumors database analysis
library(NOISeq)
library(ggplot2)
library(ggpubr)
library(singscore)
library(clusterProfiler)

##Data preparation
groups=read.csv("scotland_groupdata.csv",header=TRUE,row.names=1)#Samples list prepared on excel from sheet provided in supplementary data by xia et al, with duplicates averaged and annotated by _m 
groups$double=ifelse(groups$subtype=="iR","Intrinsic",ifelse(groups$subtype=="aR", "Adaptive", ifelse(groups$subtype=="Sensitive Pre-adaptive","PMP",ifelse(groups$subtype=="Sensitive","Sensitive", ""))))
groups$subtype=ifelse(groups$subtype=="iR","Intrinsic",ifelse(groups$subtype=="aR", "Adaptive", ifelse(groups$subtype=="Sensitive Pre-adaptive","Sensitive",ifelse(groups$subtype=="Sensitive","Sensitive", ""))))

matrix=read.table("meanbio_scotland.txt", header=T, sep="\t")#Gene matrix was prepared according to duplicate samples by calculating means of duplicate samples based on matrix from xia et al
names(matrix) <- gsub("[.]", "-", names(matrix))
names(matrix) <- gsub("X", "", names(matrix))

#Data filtering, normalization and log2 transformation
mat=matrix[rowSums(matrix) > 0,]
data=uqua(mat, long = 1000, lc = 0, k = 0)
my_data=data.frame(data)
dataz=log2(my_data)
names(dataz) <- gsub("[.]", "-", names(dataz))
names(dataz) <- gsub("X", "", names(dataz))


##Analysis
#Select RET and GDNF log counts
selectedRET=dataz[c("RET"),]
selectedRET=data.frame(t(selectedRET))
selectedRET$response=groups$type
selectedRET$responsetype=groups$subtype
selectedRET$double=groups$double

selectedGDNF=dataz[c("ATF2"),]
rownames(selectedGDNF)="GDNF"
selectedGDNF=data.frame(t(selectedGDNF))

merged=cbind(selectedGDNF, selectedRET)

#Box plots RET
ggplot(merged, aes(x = response, y = RET, fill=response))+geom_boxplot(outlier.shape = NA) + scale_fill_manual(name= "Response",values=c("grey40", "grey"))+ xlab("Response") + theme_classic()+  ylab("RET") + geom_signif(comparisons = list(c("Sensitive", "Resistant")), map_signif_level=TRUE, step_increase = 0.05, tip_length = 0.01, test="t.test")+ theme(text=element_text(size=15, face="bold"),axis.text = element_text(size=15, face="bold"),axis.text.x = element_text(size=13))+theme(legend.title = element_text(face="bold", size=14),legend.text = element_text(size=14),legend.position = "right")
bact <- compare_means(RET~response, data = merged,method = 't.test')
bact

ggplot(merged, aes(x = responsetype, y = RET, fill=response))+geom_boxplot(outlier.shape = NA) + scale_fill_manual(name= "Response",values=c("grey40", "grey"))+ xlab("Response Type") + theme_classic()+  ylab("RET") + geom_signif(comparisons = list(c("Sensitive", "Adaptive"), c("Sensitive", "Intrinsic"), c("Adaptive", "Intrinsic")), map_signif_level=TRUE, step_increase = 0.05, tip_length = 0.01, test="t.test")+ theme(text=element_text(size=15, face="bold"),axis.text = element_text(size=15, face="bold"),axis.text.x = element_text(size=13))+theme(legend.title = element_text(face="bold", size=14),legend.text = element_text(size=14),legend.position = "right")
bact <- compare_means(RET~responsetype, data = merged,method = 't.test')
bact

ggplot(merged[grepl("Adaptive|PMP", merged$double),], aes(x = double, y = RET, fill=response))+geom_boxplot(outlier.shape = NA) + scale_fill_manual(name= "Response",values=c("grey40", "grey"))+ xlab("Response Type") + theme_classic()+  ylab("RET") + geom_signif(comparisons = list(c("PMP", "Adaptive")), map_signif_level=TRUE, step_increase = 0.05, tip_length = 0.01, test="t.test")+ theme(text=element_text(size=15, face="bold"),axis.text = element_text(size=15, face="bold"),axis.text.x = element_text(size=13))+theme(legend.title = element_text(face="bold", size=14),legend.text = element_text(size=14),legend.position = "right")+scale_x_discrete(labels=c('Adaptive', 'Patient Matched Preresistant'))
bact <- compare_means(RET~double, data = merged[grepl("Adaptive|PMP", merged$double),],method = 't.test')
bact


#Box plots GDNF
ggplot(merged, aes(x = response, y = GDNF, fill=response))+geom_boxplot(outlier.shape = NA) + scale_fill_manual(name= "Response",values=c("grey40", "grey"))+ xlab("Response") + theme_classic()+  ylab("GDNF") + geom_signif(comparisons = list(c("Sensitive", "Resistant")), map_signif_level=TRUE, step_increase = 0.05, tip_length = 0.01, test="t.test")+ theme(text=element_text(size=15, face="bold"),axis.text = element_text(size=15, face="bold"),axis.text.x = element_text(size=13))+theme(legend.title = element_text(face="bold", size=14),legend.text = element_text(size=14),legend.position = "right")
bact <- compare_means(GDNF~response, data = merged,method = 't.test')
bact

ggplot(merged, aes(x = responsetype, y = GDNF, fill=response))+geom_boxplot(outlier.shape = NA) + scale_fill_manual(name= "Response",values=c("grey40", "grey"))+ xlab("Response Type") + theme_classic()+  ylab("GDNF") + geom_signif(comparisons = list(c("Sensitive", "Adaptive"), c("Sensitive", "Intrinsic"), c("Adaptive", "Intrinsic")), map_signif_level=TRUE, step_increase = 0.05, tip_length = 0.01, test="t.test")+ theme(text=element_text(size=15, face="bold"),axis.text = element_text(size=15, face="bold"),axis.text.x = element_text(size=13))+theme(legend.title = element_text(face="bold", size=14),legend.text = element_text(size=14),legend.position = "right")
bact <- compare_means(GDNF~responsetype, data = merged,method = 't.test')
bact

ggplot(merged[grepl("Adaptive|PMP", merged$double),], aes(x = double, y = GDNF, fill=response))+geom_boxplot(outlier.shape = NA) + scale_fill_manual(name= "Response",values=c("grey40", "grey"))+ xlab("Response Type") + theme_classic()+  ylab("GDNF") + geom_signif(comparisons = list(c("PMP", "Adaptive")), map_signif_level=TRUE, step_increase = 0.05, tip_length = 0.01, test="t.test")+ theme(text=element_text(size=15, face="bold"),axis.text = element_text(size=15, face="bold"),axis.text.x = element_text(size=13))+theme(legend.title = element_text(face="bold", size=14),legend.text = element_text(size=14),legend.position = "right")+scale_x_discrete(labels=c('Adaptive', 'Patient Matched Preresistant'))
bact <- compare_means(GDNF~double, data = merged[grepl("Adaptive|PMP", merged$double),],method = 't.test')
bact


#Scatter plots GDNF vs RET
quantile(merged$RET) #determine x-intercept for 75th percentile cutoff
quantile(merged$GDNF)#determine x-intercept for 75th percentile cutoff

ggplot(merged,aes(x=RET,y=GDNF, color=response))+geom_point(size=3)+scale_color_manual(name="Response", labels=c('Resistant', 'Sensitive'), values=c("grey40", "grey"))+ylab("GDNF")+theme(axis.title = element_text(size=15, face="bold"),axis.text = element_text(size=15, face="bold"), legend.title = element_text(face="bold"))+ geom_hline(yintercept = 5.62, color = "grey4", size=0.5)+geom_vline(xintercept = 5.71, color = "grey4", size=0.5)

ggplot(merged,aes(x=RET,y=GDNF, color=double))+geom_point(size=3)+scale_color_manual(name="Response Type", labels=c("Adaptive","Intrinsic","Patient-Matched Preresistant","Sensitive"), values=c("violetred4", "skyblue", "tomato", "grey"))+ylab("GDNF")+theme(axis.title = element_text(size=15, face="bold"),axis.text = element_text(size=15, face="bold"), legend.title = element_text(face="bold"))+ geom_hline(yintercept = 5.62, color = "grey4", size=0.5)+geom_vline(xintercept = 5.71, color = "grey4", size=0.5)

#Bar plot for high RET/GDNF grouping
merged2=merged
merged2$groupboth=ifelse((merged2$RET>5.71530547)&(merged2$GDNF>5.620480),"High RET/GDNF","Non-high RET/GDNF")
merged2$double[merged2$double == 'PMP'] = 'Patient-Matched Preresistant'

ggplot(merged2, aes(x=groupboth, fill=double)) + geom_bar(position="fill", stat="count")  + scale_fill_manual(values=c("violetred4", "skyblue3","tomato","grey"),name="RET  & GDNF Co-Expression Group") + xlab("") + ylab("Percent of Samples") + theme(axis.title = element_text(size=12, face="bold"),axis.text = element_text(size=12, face="bold"), legend.title = element_text(face="bold"))+ scale_y_continuous(labels = scales::percent_format(accuracy = 1))+coord_flip()

#Gene set scoring
#Load gene sets from TCGA analysis
GO_file = "HIGHRETGDNF_top40.gmt"
data2=read.gmt(GO_file)
up=data2[grepl("UPREG",data2$term),]
up=up$gene
down=data2[grepl("DOWNREG",data2$term),]
down=down$gene

#Run score allocation
rankData <- rankGenes(dataz)
scoredf <- simpleScore(rankData, upSet = up, downSet = down)
scoredf$catt=groups$type
scoredf$group=groups$subtype
scoredf$double=groups$double

#Up set scores
ggplot(scoredf, aes(x = catt, y = UpScore, fill=catt))+geom_boxplot(outlier.shape = NA)+ scale_fill_manual(name= "Sample Response",values=c("grey40", "grey"))+ xlab("Response Subtype") + theme_classic()+  ylab("UPREGULATED_HIGHRETGDNF Score") + geom_signif(comparisons = list(c("Sensitive", "Resistant")), map_signif_level=TRUE, step_increase = 0.05, tip_length = 0.01, test="t.test")+ theme(text=element_text(size=15, face="bold"),axis.text = element_text(size=15, face="bold"),axis.text.x = element_text(size=12),axis.title = element_text(size=12),legend.title = element_text(face="bold", size=14),legend.text = element_text(size=14),legend.position = "none")
bact <- compare_means(UpScore~catt, data = scoredf,method = 't.test')
bact

ggplot(scoredf, aes(x = group, y = UpScore, fill=catt))+geom_boxplot(outlier.shape = NA)+ scale_fill_manual(name= "Sample Response",values=c("grey40", "grey"))+ xlab("Response Subtype") + theme_classic()+  ylab("UPREGULATED_HIGHRETGDNF Score") + geom_signif(comparisons = list(c("Adaptive", "Sensitive"),c("Intrinsic","Sensitive"), c("Adaptive","Intrinsic")), map_signif_level=TRUE, step_increase = 0.05, tip_length = 0.01, test="t.test")+ theme(text=element_text(size=15, face="bold"),axis.text = element_text(size=15, face="bold"),axis.text.x = element_text(size=12),axis.title = element_text(size=12),legend.title = element_text(face="bold", size=14),legend.text = element_text(size=14),legend.position = "none")
bact <- compare_means(UpScore~group, data = scoredf,method = 't.test')
bact

ggplot(scoredf[grepl("Adaptive|PMP", scoredf$double),], aes(x = double, y = UpScore, fill=catt))+geom_boxplot(outlier.shape = NA)+ scale_fill_manual(name= "Sample Response",values=c("grey40", "grey"))+ xlab("Response Subtype") + theme_classic()+  ylab("UPREGULATED_HIGHRETGDNF Score") + geom_signif(comparisons = list(c("Adaptive", "PMP")), map_signif_level=TRUE, step_increase = 0.05, tip_length = 0.01, test="t.test")+ theme(text=element_text(size=15, face="bold"),axis.text = element_text(size=15, face="bold"),axis.text.x = element_text(size=12),axis.title = element_text(size=12),legend.title = element_text(face="bold", size=14),legend.text = element_text(size=14),legend.position = "none")+scale_x_discrete(labels=c('Adaptive', 'Patient Matched Preresistant'))
bact <- compare_means(UpScore~double, data = scoredf[grepl("Adaptive|PMP", scoredf$double),],method = 't.test')
bact

#Down set scores
ggplot(scoredf, aes(x = catt, y = DownScore, fill=catt))+geom_boxplot(outlier.shape = NA)+ scale_fill_manual(name= "Sample Response",values=c("grey40", "grey"))+ xlab("Response Subtype") + theme_classic()+  ylab("DOWNREGULATED_HIGHRETGDNF Score") + geom_signif(comparisons = list(c("Sensitive", "Resistant")), map_signif_level=TRUE, step_increase = 0.05, tip_length = 0.01, test="t.test")+ theme(text=element_text(size=15, face="bold"),axis.text = element_text(size=15, face="bold"),axis.text.x = element_text(size=12),axis.title = element_text(size=12),legend.title = element_text(face="bold", size=14),legend.text = element_text(size=14),legend.position = "none")
bact <- compare_means(DownScore~catt, data = scoredf,method = 't.test')
bact

ggplot(scoredf, aes(x = group, y = DownScore, fill=catt))+geom_boxplot(outlier.shape = NA)+ scale_fill_manual(name= "Sample Response",values=c("grey40", "grey"))+ xlab("Response Subtype") + theme_classic()+  ylab("DOWNREGULATED_HIGHRETGDNF Score") + geom_signif(comparisons = list(c("Adaptive", "Sensitive"),c("Intrinsic","Sensitive"), c("Adaptive","Intrinsic")), map_signif_level=TRUE, step_increase = 0.05, tip_length = 0.01, test="t.test")+ theme(text=element_text(size=15, face="bold"),axis.text = element_text(size=15, face="bold"),axis.text.x = element_text(size=12),axis.title = element_text(size=12),legend.title = element_text(face="bold", size=14),legend.text = element_text(size=14),legend.position = "none")
bact <- compare_means(DownScore~group, data = scoredf,method = 't.test')
bact

ggplot(scoredf[grepl("Adaptive|PMP", scoredf$double),], aes(x = double, y = DownScore, fill=catt))+geom_boxplot(outlier.shape = NA)+ scale_fill_manual(name= "Sample Response",values=c("grey40", "grey"))+ xlab("Response Subtype") + theme_classic()+  ylab("DOWNREGULATED_HIGHRETGDNF Score") + geom_signif(comparisons = list(c("Adaptive", "PMP")), map_signif_level=TRUE, step_increase = 0.05, tip_length = 0.01, test="t.test")+ theme(text=element_text(size=15, face="bold"),axis.text = element_text(size=15, face="bold"),axis.text.x = element_text(size=12),axis.title = element_text(size=12),legend.title = element_text(face="bold", size=14),legend.text = element_text(size=14),legend.position = "none")+scale_x_discrete(labels=c('Adaptive', 'Patient Matched Preresistant'))
bact <- compare_means(DownScore~double, data = scoredf[grepl("Adaptive|PMP", scoredf$double),],method = 't.test')
bact

