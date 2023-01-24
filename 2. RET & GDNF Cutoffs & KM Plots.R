#Data preparation for ER+ patients only
library(fabricatr)
library(survminer)
library(survival)


#Load clinical data
survivaldata=read.csv("ClinicalDataCombined.csv",header=TRUE, row.names=1)

#ER+ patient data selection
ERclinical=subset(clinical, clinical$ER =="Positive")
ERpts=rownames(ERclinical)
ER=select(allgenecounts, ERpts)
  
#Select RET data among ER+ samples
RET=subset(ER, rownames(ER) == "RET")
RET=as.data.frame(t(RET))
RET=RET[rownames(ERclinical),,drop=FALSE]
RET=cbind(RET,ERclinical)

##--------CUTOFFS (for each cutoff indicated, obtain the groups data frame and input to the survival plotting section. to proceed to the next cutoff, run the code pertaining to that cutoff then run the survival plotting section again)
df=RET
#Above 75% vs below
splitRET <- split_quantile(x = df$RET, type = 4)
df$group <- splitRET
ord=df[order(df$RET),]
groups <- ord
GRP <- ifelse(groups$group == "4", "Above", "Below")
groups$categ <- GRP
groups <- groups[c(5)]
#Above 50% vs below
splitRET <- split_quantile(x = df$RET, type = 4)
df$group <- splitRET
ord=df[order(df$RET),]
groups <- ord
GRP <- ifelse(groups$group == "4", "Above", ifelse(groups$group == "3","Above","Below"))
groups$categ <- GRP
groups <- groups[c(5)]
#Above 25% vs below
splitRET <- split_quantile(x = df$RET, type = 4)
df$group <- splitRET
ord=df[order(df$RET),]
groups <- ord
GRP <- ifelse(groups$group == "4", "Above", ifelse(groups$group == "3","Above",ifelse(groups$group == "2","Above","Below")))
groups$categ <- GRP
groups <- groups[c(5)]
#Above 75% vs below 25%
splitRET <- split_quantile(x = df$RET, type = 4)
df$group <- splitRET
ord=df[order(df$RET),]
groups <- ord
GRP <- ifelse(groups$group == "4", "Above", ifelse(groups$group == "1","Below","Neither"))
groups$categ <- GRP
groups=subset(groups,groups$categ !="Neither")
groups <- groups[c(5)]
#Above 90% vs below
splitRET <- split_quantile(x = df$RET, type = 10)
df$group <- splitRET
ord=df[order(df$RET),]
groups <- ord
GRP <- ifelse(groups$group == "10", "Above", "Below")
groups$categ <- GRP
groups <- groups[c(5)]
#Above 90% vs below 10%
splitRET <- split_quantile(x = df$RET, type = 10)
df$group <- splitRET
ord=df[order(df$RET),]
groups <- ord
GRP <- ifelse(groups$group == "10", "Above", ifelse(groups$group == "1","Below","Neither"))
groups$categ <- GRP
groups=subset(groups,groups$categ !="Neither")
groups <- groups[c(5)]

###Survival plotting
pts=rownames(groups)
groups$Category <- groups$categ
groups$Category <- factor(groups$Category, levels = c("Below", "Above")) ##modify groups accordingly, 1st group is reference
Category <- groups$Category
survival_group <- survivaldata[pts,]

#Overall survival
survcolumns <- survival_group[c(29,30)]
survival <- as.matrix(survcolumns)
survv <-sub("0:LIVING", "0", survival)
survv2 <-sub("1:DECEASED", "1", survv)
survivaldf <- as.data.frame(survv2)
OS_MONTHS <- survivaldf$OS_MONTHS
OS_STATUS <- survivaldf$OS_STATUS
monthz <- as.numeric(OS_MONTHS)
status <- as.numeric(OS_STATUS)
fit <- survfit(Surv(monthz, status) ~ Category, data = survivaldf)
res.cox <- coxph(Surv(monthz, status) ~ as.factor(Category), data = survivaldf)##hazard ratio determination
summary(res.cox, conf.int=FALSE) 
ggsurvplot(fit, data = survivaldf, pval=TRUE, linetype = 1, break.x.by = 25, xlab = "Time (months)", xlim=c(0,100),conf.int = FALSE, risk.table = FALSE, tables.theme = theme_cleantable(), palette = c("#2E9FDF", "#E7B800"),tables.y.text = FALSE, font.subtitle = c(10, "black"))

#Progression free survival
survcolumns <- survival_group[c(35,36)]
survival <- as.matrix(survcolumns)
survv <-sub("0:CENSORED", "0", survival)
survv2 <-sub("1:PROGRESSION", "1", survv)
survivaldf <- as.data.frame(survv2)
PFS_MONTHS <- survivaldf$PFS_MONTHS
PFS_STATUS <- survivaldf$PFS_STATUS
monthz <- as.numeric(PFS_MONTHS)
status <- as.numeric(PFS_STATUS)
fit <- survfit(Surv(monthz, status) ~ Category, data = survivaldf)
res.cox <- coxph(Surv(monthz, status) ~ as.factor(Category), data = survivaldf)
summary(res.cox, conf.int=FALSE)
nrow(survivaldf)
ggsurvplot(fit, data = survivaldf, pval=TRUE, linetype = 1, break.x.by = 25, xlab = "Time (months)", xlim=c(0,100),conf.int = FALSE, risk.table = FALSE, tables.theme = theme_cleantable(), palette = c("#2E9FDF", "#E7B800"),tables.y.text = FALSE, font.subtitle = c(10, "black"))

#Disease free survival
survcolumns <- survival_group[c(33,34)]
survival <- as.matrix(survcolumns)
survv <-sub("0:DiseaseFree", "0", survival)
survv2 <-sub("1:Recurred/Progressed", "1", survv)
survivaldf <- as.data.frame(survv2)
DFS_MONTHS <- survivaldf$DFS_MONTHS
DFS_STATUS <- survivaldf$DFS_STATUS
monthz <- as.numeric(DFS_MONTHS)
status <- as.numeric(DFS_STATUS)
fit <- survfit(Surv(monthz, status) ~ Category, data = survivaldf)
res.cox <- coxph(Surv(monthz, status) ~ as.factor(Category), data = survivaldf)
summary(res.cox, conf.int=FALSE)
nrow(survivaldf)
ggsurvplot(fit, data = survivaldf, pval=TRUE, linetype = 1, break.x.by = 25, xlab = "Time (months)", xlim=c(0,100),conf.int = FALSE, risk.table = FALSE, tables.theme = theme_cleantable(), palette = c("#2E9FDF", "#E7B800"),tables.y.text = FALSE, font.subtitle = c(10, "black"))

#Disease specific survival
survcolumns <- survival_group[c(31,32)]
survival <- as.matrix(survcolumns)
survv <-sub("0:ALIVE OR DEAD TUMOR FREE", "0", survival)
survv2 <-sub("1:DEAD WITH TUMOR", "1", survv)
survivaldf <- as.data.frame(survv2)
DSS_MONTHS <- survivaldf$DSS_MONTHS
DSS_STATUS <- survivaldf$DSS_STATUS
monthz <- as.numeric(DSS_MONTHS)
status <- as.numeric(DSS_STATUS)
fit <- survfit(Surv(monthz, status) ~ Category, data = survivaldf)
res.cox <- coxph(Surv(monthz, status) ~ as.factor(Category), data = survivaldf)
summary(res.cox, conf.int=FALSE)
nrow(survivaldf)
ggsurvplot(fit, data = survivaldf, pval=TRUE, linetype = 1, break.x.by = 25, xlab = "Time (months)", xlim=c(0,100),conf.int = FALSE, risk.table = FALSE, tables.theme = theme_cleantable(), palette = c("#2E9FDF", "#E7B800"),tables.y.text = FALSE, font.subtitle = c(10, "black"))




###RET & GDNF co-expression KM plots
#Select RET data among ER+ samples
RET=subset(ER, rownames(ER) == "RET")
RET=as.data.frame(t(RET))
RET=RET[rownames(ERclinical),,drop=FALSE]
RET=cbind(RET,ERclinical)

#Select GDNF data among ER+ samples
GDNF=subset(ER, rownames(ER) == "ATF2")
GDNF=as.data.frame(t(GDNF))
names(GDNF)="GDNF"
GDNF=GDNF[rownames(ERclinical),,drop=FALSE]

RETGDNF=cbind(RET,GDNF)

##--------CUTOFFS (similar to above)
df=RETGDNF
#Above 75% vs below
splitRET <- split_quantile(x = df$RET, type = 4)
df$groupRET <- splitRET
splitGDNF <- split_quantile(x = df$GDNF, type = 4)
df$groupGDNF <- splitGDNF
ord=df[order(df$RET),]
groups <- ord
GRP <- ifelse(groups$groupRET == "4" & groups$groupGDNF == "4", "Above", "Below")
groups$categ <- GRP
groups <- groups[c(7)]

###Survival plotting
pts=rownames(groups)
groups$Category <- groups$categ
groups$Category <- factor(groups$Category, levels = c("Below", "Above"))
Category <- groups$Category
survival_group <- survivaldata[pts,]


#Overall survival
survcolumns <- survival_group[c(29,30)]
survival <- as.matrix(survcolumns)
survv <-sub("0:LIVING", "0", survival)
survv2 <-sub("1:DECEASED", "1", survv)
survivaldf <- as.data.frame(survv2)
OS_MONTHS <- survivaldf$OS_MONTHS
OS_STATUS <- survivaldf$OS_STATUS
monthz <- as.numeric(OS_MONTHS)
status <- as.numeric(OS_STATUS)
fit <- survfit(Surv(monthz, status) ~ Category, data = survivaldf)
res.cox <- coxph(Surv(monthz, status) ~ as.factor(Category), data = survivaldf)##hazard ratio determination
summary(res.cox, conf.int=FALSE) 
ggsurvplot(fit, data = survivaldf, pval=TRUE, linetype = 1, break.x.by = 25, xlab = "Time (months)", conf.int = FALSE, xlim=c(0, 100), risk.table = FALSE, palette = c("#00BFC4", "#F8766D"), tables.theme = theme_cleantable(), tables.y.text = FALSE, font.subtitle = c(10, "black"))

#Progression free survival
survcolumns <- survival_group[c(35,36)]
survival <- as.matrix(survcolumns)
survv <-sub("0:CENSORED", "0", survival)
survv2 <-sub("1:PROGRESSION", "1", survv)
survivaldf <- as.data.frame(survv2)
PFS_MONTHS <- survivaldf$PFS_MONTHS
PFS_STATUS <- survivaldf$PFS_STATUS
monthz <- as.numeric(PFS_MONTHS)
status <- as.numeric(PFS_STATUS)
fit <- survfit(Surv(monthz, status) ~ Category, data = survivaldf)
res.cox <- coxph(Surv(monthz, status) ~ as.factor(Category), data = survivaldf)
summary(res.cox, conf.int=FALSE)
nrow(survivaldf)
ggsurvplot(fit, data = survivaldf, pval=TRUE, linetype = 1, break.x.by = 25, xlab = "Time (months)", conf.int = FALSE, xlim=c(0, 100), risk.table = FALSE, palette = c("#00BFC4", "#F8766D"), tables.theme = theme_cleantable(), tables.y.text = FALSE, font.subtitle = c(10, "black"))

#Disease free survival
survcolumns <- survival_group[c(33,34)]
survival <- as.matrix(survcolumns)
survv <-sub("0:DiseaseFree", "0", survival)
survv2 <-sub("1:Recurred/Progressed", "1", survv)
survivaldf <- as.data.frame(survv2)
DFS_MONTHS <- survivaldf$DFS_MONTHS
DFS_STATUS <- survivaldf$DFS_STATUS
monthz <- as.numeric(DFS_MONTHS)
status <- as.numeric(DFS_STATUS)
fit <- survfit(Surv(monthz, status) ~ Category, data = survivaldf)
res.cox <- coxph(Surv(monthz, status) ~ as.factor(Category), data = survivaldf)
summary(res.cox, conf.int=FALSE)
nrow(survivaldf)
ggsurvplot(fit, data = survivaldf, pval=TRUE, linetype = 1, break.x.by = 25, xlab = "Time (months)", conf.int = FALSE, xlim=c(0, 100), risk.table = FALSE, palette = c("#00BFC4", "#F8766D"), tables.theme = theme_cleantable(), tables.y.text = FALSE, font.subtitle = c(10, "black"))

#Disease specific survival
survcolumns <- survival_group[c(31,32)]
survival <- as.matrix(survcolumns)
survv <-sub("0:ALIVE OR DEAD TUMOR FREE", "0", survival)
survv2 <-sub("1:DEAD WITH TUMOR", "1", survv)
survivaldf <- as.data.frame(survv2)
DSS_MONTHS <- survivaldf$DSS_MONTHS
DSS_STATUS <- survivaldf$DSS_STATUS
monthz <- as.numeric(DSS_MONTHS)
status <- as.numeric(DSS_STATUS)
fit <- survfit(Surv(monthz, status) ~ Category, data = survivaldf)
res.cox <- coxph(Surv(monthz, status) ~ as.factor(Category), data = survivaldf)
summary(res.cox, conf.int=FALSE)
nrow(survivaldf)
ggsurvplot(fit, data = survivaldf, pval=TRUE, linetype = 1, break.x.by = 25, xlab = "Time (months)", conf.int = FALSE, xlim=c(0, 100), risk.table = FALSE, palette = c("#00BFC4", "#F8766D"), tables.theme = theme_cleantable(), tables.y.text = FALSE, font.subtitle = c(10, "black"))