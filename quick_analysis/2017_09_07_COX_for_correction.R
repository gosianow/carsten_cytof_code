rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(survival)
library(gridExtra)
library(reshape2)
library(data.table)
library(knitr)
library(ROCR)
library(grid)
library(OptimalCutpoints)
library(lubridate)
library(plotrix)
library(rmeta)
library(broom)
library(survminer)
library(devtools)
library(ggcorrplot)

setwd("/data/Sabrina/CYTOF_revision/COX/")

cytof <- read.delim("/data/Sabrina/CYTOF_revision/2017_09_06_updated_info_CyTOF_patients.txt", comment.char="#", na.strings = c("NA"," ",""))


cytof$batch <- as.factor(cytof$batches)
cytof$start.treatment.date <- dmy(cytof$start.treatment.date)
cytof$date.last.treatment <- dmy(cytof$date.last.treatment)
cytof$date.progression <- dmy(cytof$date.progression...comments.)
cytof$last.seen.in.clinic.date <- dmy(cytof$last.seen.in.clinic.date)
cytof$date.death <- dmy(cytof$date.death)
cytof$RESPONSE <- as.factor(cytof$RESPONSE)
cytof$general.stage <- as.factor(cytof$general.stage)
cytof$pT <- as.factor(cytof$pT)
cytof$pT.ulceration <- as.factor(cytof$pT.ulceration)
cytof$N.stage <- as.factor(cytof$N.stage)
cytof$M.stage <- as.factor(cytof$M.stage)
cytof$last.ipi <- dmy(cytof$last.ipi)
cytof$BRAF.status <- factor(cytof$BRAF.status, levels = c("wt", "mut"))
cytof$NRAS.status <- factor(cytof$NRAS.status, levels = c("wt", "mut"))
cytof$cKIT.status <- factor(cytof$cKIT.status, levels = c("wt", "mut"))


# select CYTOF baseline 
CYTOFlist <- filter(cytof, use == "CyTOF")
CYTOFlist_base <- filter(CYTOFlist, TP == "baseline")

#list factors with batch
list.1 <- colnames(select(cytof, sex, age.sampling, prev.ipi, prev.chemo, prev.radio, prev.target, BRAF.status, NRAS.status, cKIT.status, leukocytes.count, LDH, ANC, ALC, ANC.ALC, S100, hemoglobin, hematocrit, erythro, MCV, MCH, MCHC, RDW, thrombocytes, monocytes, eosino, baso, IG., sodium, potatium, urea, creatinin, eGFR, bilirubin, protein, albumin, AST, ALT, GGT, alk.phosph, CRP, TSH, brain.met, liver.met, lung.met, bone.met, Classical.Monocytes, batch, general.stage, pT, pT.ulceration, N.stage, N.size, M.stage, nb.days.last.ipi)) 


##############################  CyTOF cohort - BASELINE samples only  ##############################
####UNIVARIATE  

COX_uni_CYTOF_base_m <- lapply(list.1, function(x) {                  
  cox <- coxph(Surv(PFS.months, PFS.status) ~ get(x), CYTOFlist_base) 
  cox1 <- tidy(cox)
  cbind(cox1, 
        HRate = round(signif(cox1$estimate,3),3),
        LowerCI = signif(cox1$conf.low,3), 
        UpperCI = signif(cox1$conf.high,3), 
        variables = eval(x), 
        Pvalue = signif(cox1$p.value,3))
})

COX2 <- do.call(rbind, COX_uni_CYTOF_base_m)
COX2 <- select(COX2, variables, term, HRate, LowerCI, UpperCI, Pvalue)
COX2$HRate.exp <- round(signif(exp(COX2$HRate),3),3)
COX2$LowerCI.exp <- signif(exp(COX2$LowerCI),3)
COX2$UpperCI.exp <- signif(exp(COX2$UpperCI),3)
COX2 <- arrange(COX2, Pvalue)

# filter P val < 0.05 and assign to COX_uni_CYTOF_signP
COX_uni_CYTOFm_base_signP <- filter(COX2, Pvalue <= 0.05)

COX2$term <- gsub("get\\(x\\)([0-9a-zA-Z]*)", "\\1", COX2$term)
COX2.1 <- unite(COX2, variables, variables, term, sep=" ")


write.table(COX2.1, "/data/Sabrina/CYTOF_revision/COX/tables/Cox_CYTOF_baseline_univariate.txt", sep="\t", row.names = F)


# Forest plot univariate CyTOF baseline 
COX_uni_CYTOFm_base_forest <- fread("/data/Sabrina/CYTOF_revision/COX/tables/Cox_CYTOF_baseline_univariate.txt", showProgress = T)
COX_uni_CYTOFm_base_forest <- na.omit(COX_uni_CYTOFm_base_forest)
COX_uni_CYTOFm_base_forest <- filter(COX_uni_CYTOFm_base_forest, !grepl("NaN", Pvalue))

var.names <- read.delim("/data/Sabrina/CYTOF_revision/list_variables_names_units_17_07_17.txt", comment.char="#", na.strings = c("NA"," ",""))
var.names <- select(var.names, variables, clean.names)
COX_uni_CYTOFm_base_forest <- left_join(COX_uni_CYTOFm_base_forest, var.names, by="variables")

COX_uni_CYTOFm_base_forest$variables <- as.character(COX_uni_CYTOFm_base_forest$clean.names)

a <- rep(NA, length(COX_uni_CYTOFm_base_forest$HRate))
b <- rep(NA, length(COX_uni_CYTOFm_base_forest$HRate))
c <- rep(NA, length(COX_uni_CYTOFm_base_forest$HRate))
d <- rep(NA, length(COX_uni_CYTOFm_base_forest$HRate))
e <- rep(NA, length(COX_uni_CYTOFm_base_forest$HRate))


COX_uni_CYTOFm_base_forest <- mutate (COX_uni_CYTOFm_base_forest, a, b, c, d, e)
COX_uni_CYTOFm_base_forest <- select(COX_uni_CYTOFm_base_forest, variables, a, HRate, b, LowerCI, c, UpperCI, d, Pvalue, e)

setnames(COX_uni_CYTOFm_base_forest, c("LowerCI", "UpperCI", "a", "b", "c", "d", "e"), 
         c("LCL(95%)", "UCL(95%)", "                                                         ", "   ", "   ", "   ", "   "))

pdf("COXmodel_CyTOF.univariate.pdf", width = 12, height = 18 )
forestplot(rbind(colnames(COX_uni_CYTOFm_base_forest), as.matrix(COX_uni_CYTOFm_base_forest)), 
           c(NA,COX_uni_CYTOFm_base_forest$'HRate'),
           c(NA,COX_uni_CYTOFm_base_forest$'LCL(95%)'),
           c(NA,COX_uni_CYTOFm_base_forest$'UCL(95%)'), 
           zero=0, 
           col=meta.colors(box="red",line="darkblue", summary="royalblue"), 
           is.summary = rep(F, 12), 
           title("UNIVARIATE Cox regression model - DISCOVERY COHORT, baseline samples only"),
           xlab="Hazard Rate", xlog=F, xticks=c(seq(-10,10, by=2)), 
           boxsize = 0.6, clip= c(-10, 10))
dev.off()

######### MULTIVARIATE

signP_list <- COX_uni_CYTOFm_base_signP$variables
signP_list1 <- paste(signP_list, collapse = " + ")

COX_multi_CYTOFm_base <- coxph(Surv(PFS.months, PFS.status) ~ prev.target + thrombocytes + potatium, data = CYTOFlist_base)
COX_multi_CYTOFm_base_t <- tidy(COX_multi_CYTOFm_base)
COX_multi_CYTOFm_base2 <- cbind(COX_multi_CYTOFm_base_t, 
                                    variables = COX_multi_CYTOFm_base_t$term,
                                    HRate = round(signif(COX_multi_CYTOFm_base_t$estimate,3),3),
                                    LowerCI = round(signif(COX_multi_CYTOFm_base_t$conf.low,3),3), 
                                    UpperCI = round(signif(COX_multi_CYTOFm_base_t$conf.high,3),3), 
                                    Pvalue = round(signif(COX_multi_CYTOFm_base_t$p.value,3),3))
COX_multi_CYTOFm_baseline <- select(COX_multi_CYTOFm_base2, variables, HRate, LowerCI, UpperCI, Pvalue)

COX_multi_CYTOFm_baseline$HRate.exp <- round(signif(exp(COX_multi_CYTOFm_baseline$HRate),3),3)
COX_multi_CYTOFm_baseline$LowerCI.exp <- signif(exp(COX_multi_CYTOFm_baseline$LowerCI),3)
COX_multi_CYTOFm_baseline$UpperCI.exp <- signif(exp(COX_multi_CYTOFm_baseline$UpperCI),3)

write.table(COX_multi_CYTOFm_baseline, "/data/Sabrina/CYTOF_revision/COX/tables/Cox_CYTOF_baseline_multivariate.txt", sep="\t", row.names = F)

## Forest plot Multivariate COX regression CYTOF baseline
COX_multi_CYTOFm_base_forest <- fread("/data/Sabrina/CYTOF_revision/COX/tables/Cox_CYTOF_baseline_multivariate.txt", showProgress = T)
COX_multi_CYTOFm_base_forest <- arrange(COX_multi_CYTOFm_base_forest, Pvalue)
COX_multi_CYTOFm_base_forest <- left_join(COX_multi_CYTOFm_base_forest, var.names, by="variables")

COX_multi_CYTOFm_base_forest$variables <- as.character(COX_multi_CYTOFm_base_forest$clean.names)

a <- rep(NA, length(COX_multi_CYTOFm_base_forest$HRate))
b <- rep(NA, length(COX_multi_CYTOFm_base_forest$HRate))
c <- rep(NA, length(COX_multi_CYTOFm_base_forest$HRate))
d <- rep(NA, length(COX_multi_CYTOFm_base_forest$HRate))
e <- rep(NA, length(COX_multi_CYTOFm_base_forest$HRate))


COX_multi_CYTOFm_base_forest <- mutate (COX_multi_CYTOFm_base_forest, a, b, c, d, e)
COX_multi_CYTOFm_base_forest <- select(COX_multi_CYTOFm_base_forest, variables, a, HRate, b, LowerCI, c, UpperCI, d, Pvalue, e)

setnames(COX_multi_CYTOFm_base_forest, c("LowerCI", "UpperCI", "a", "b", "c", "d", "e"), 
         c("LCL(95%)", "UCL(95%)", "                                                      ", "   ", "   ", "   ", "   "))

pdf("COXmodel_CyTOF.multivariate.pdf", width = 12, height = 18 )
forestplot(rbind(colnames(COX_multi_CYTOFm_base_forest), as.matrix(COX_multi_CYTOFm_base_forest)), 
           c(NA,COX_multi_CYTOFm_base_forest$'HRate'),
           c(NA,COX_multi_CYTOFm_base_forest$'LCL(95%)'),
           c(NA,COX_multi_CYTOFm_base_forest$'UCL(95%)'), 
           zero=0, 
           col=meta.colors(box="red",line="darkblue", summary="royalblue"), 
           is.summary = rep(F, 12), 
           title("MULTIVARIATE Cox regression model - DISCOVERY COHORT, baseline samples only"),
           xlab="Hazard Rate", xlog=F, xticks=c(seq(-4,4, by=1)), 
           boxsize = 0.6, clip= c(-5, 5))
dev.off()



############################# FACS cohort ##################################################
# select FACS 
FACSlist <- filter(cytof, use == "FACS")

list.2 <- colnames(select(cytof, sex, age.sampling, prev.ipi, prev.chemo, prev.radio, prev.target, BRAF.status, NRAS.status, cKIT.status, leukocytes.count, LDH, ANC, ALC, ANC.ALC, S100, hemoglobin, hematocrit, erythro, MCV, MCH, MCHC, RDW, thrombocytes, monocytes, eosino, baso, IG., sodium, potatium, urea, creatinin, eGFR, bilirubin, protein, albumin, AST, ALT, GGT, alk.phosph, CRP, TSH, brain.met, liver.met, lung.met, bone.met, Classical.Monocytes, general.stage, pT, pT.ulceration, N.stage, N.size, M.stage, nb.days.last.ipi)) 


####UNIVARIATE  
COX_uni_FACSm <- lapply(list.2, function(x) {                  
  cox <- coxph(Surv(PFS.months, PFS.status) ~ get(x), FACSlist) 
  cox1 <- tidy(cox)
  cbind(cox1, 
        HRate = round(signif(cox1$estimate,3),3),
        LowerCI = signif(cox1$conf.low,3), 
        UpperCI = signif(cox1$conf.high,3), 
        variables = eval(x), 
        Pvalue = signif(cox1$p.value,3))
})

COX2 <- do.call(rbind, COX_uni_FACSm)
COX2 <- select(COX2, variables, term, HRate, LowerCI, UpperCI, Pvalue)
COX2$HRate.exp <- round(signif(exp(COX2$HRate),3),3)
COX2$LowerCI.exp <- signif(exp(COX2$LowerCI),3)
COX2$UpperCI.exp <- signif(exp(COX2$UpperCI),3)
COX2 <- arrange(COX2, Pvalue)

# filter P val < 0.05 and assign to COX_uni_CYTOF_signP
COX_uni_FACSm_signP <- filter(COX2, Pvalue <= 0.05)

COX2$term <- gsub("get\\(x\\)([0-9a-zA-Z]*)", "\\1", COX2$term)
COX2.1 <- unite(COX2, variables, variables, term, sep=" ")


write.table(COX2.1, "/data/Sabrina/CYTOF_revision/COX/tables/Cox_FACS_univariate.txt", sep="\t", row.names = F)


# Forest plot univariate FACS
COX_uni_FACSm_forest <- fread("/data/Sabrina/CYTOF_revision/COX/tables/Cox_FACS_univariate.txt", showProgress = T)
COX_uni_FACSm_forest <- na.omit(COX_uni_FACSm_forest)
COX_uni_FACSm_forest <- filter(COX_uni_FACSm_forest, !grepl("NaN", Pvalue))

var.names <- read.delim("/data/Sabrina/CYTOF_revision/list_variables_names_units_17_07_17.txt", comment.char="#", na.strings = c("NA"," ",""))
var.names <- select(var.names, variables, clean.names)
COX_uni_FACSm_forest <- left_join(COX_uni_FACSm_forest, var.names, by="variables")

COX_uni_FACSm_forest$variables <- as.character(COX_uni_FACSm_forest$clean.names)

a <- rep(NA, length(COX_uni_FACSm_forest$HRate))
b <- rep(NA, length(COX_uni_FACSm_forest$HRate))
c <- rep(NA, length(COX_uni_FACSm_forest$HRate))
d <- rep(NA, length(COX_uni_FACSm_forest$HRate))
e <- rep(NA, length(COX_uni_FACSm_forest$HRate))


COX_uni_FACSm_forest <- mutate (COX_uni_FACSm_forest, a, b, c, d, e)
COX_uni_FACSm_forest <- select(COX_uni_FACSm_forest, variables, a, HRate, b, LowerCI, c, UpperCI, d, Pvalue, e)

setnames(COX_uni_FACSm_forest, c("LowerCI", "UpperCI", "a", "b", "c", "d", "e"), 
         c("LCL(95%)", "UCL(95%)", "                                                         ", "   ", "   ", "   ", "   "))

pdf("COXmodel_FACS.univariate.pdf", width = 12, height = 18 )
forestplot(rbind(colnames(COX_uni_FACSm_forest), as.matrix(COX_uni_FACSm_forest)), 
           c(NA,COX_uni_FACSm_forest$'HRate'),
           c(NA,COX_uni_FACSm_forest$'LCL(95%)'),
           c(NA,COX_uni_FACSm_forest$'UCL(95%)'), 
           zero=0, 
           col=meta.colors(box="red",line="darkblue", summary="royalblue"), 
           is.summary = rep(F, 12), 
           title("UNIVARIATE Cox regression model - VALIDATION COHORT"),
           xlab="Hazard Rate", xlog=F, xticks=c(seq(-10,10, by=2)), 
           boxsize = 0.6, clip= c(-10, 10))
dev.off()

######## MULTIVARIATE

signP_list <- COX_uni_FACSm_signP$variables
signP_list1 <- paste(signP_list, collapse = " + ")

COX_multi_FACSm <- coxph(Surv(PFS.months, PFS.status) ~ Classical.Monocytes + alk.phosph + thrombocytes + LDH + CRP + IG. + prev.chemo + liver.met + M.stage + hemoglobin + albumin + baso, data = FACSlist)
COX_multi_FACSm_t <- tidy(COX_multi_FACSm)
COX_multi_FACSm_2 <- cbind(COX_multi_FACSm_t, 
                                variables = COX_multi_FACSm_t$term,
                                HRate = round(signif(COX_multi_FACSm_t$estimate,3),3),
                                LowerCI = round(signif(COX_multi_FACSm_t$conf.low,3),3), 
                                UpperCI = round(signif(COX_multi_FACSm_t$conf.high,3),3), 
                                Pvalue = round(signif(COX_multi_FACSm_t$p.value,3),3))
COX_multi_FACSm <- select(COX_multi_FACSm_2, variables, HRate, LowerCI, UpperCI, Pvalue)

COX_multi_FACSm$HRate.exp <- round(signif(exp(COX_multi_FACSm$HRate),3),3)
COX_multi_FACSm$LowerCI.exp <- signif(exp(COX_multi_FACSm$LowerCI),3)
COX_multi_FACSm$UpperCI.exp <- signif(exp(COX_multi_FACSm$UpperCI),3)

write.table(COX_multi_FACSm, "/data/Sabrina/CYTOF_revision/COX/tables/Cox_FACS_multivariate.txt", sep="\t", row.names = F)

## Forest plot Multivariate analysis FACS 
COX_multi_FACSm_forest <- fread("/data/Sabrina/CYTOF_revision/COX/tables/Cox_FACS_multivariate.txt", showProgress = T)
COX_multi_FACSm_forest <- arrange(COX_multi_FACSm_forest, Pvalue)
COX_multi_FACSm_forest <- left_join(COX_multi_FACSm_forest, var.names, by="variables")

COX_multi_FACSm_forest$variables <- as.character(COX_multi_FACSm_forest$clean.names)

a <- rep(NA, length(COX_multi_FACSm_forest$HRate))
b <- rep(NA, length(COX_multi_FACSm_forest$HRate))
c <- rep(NA, length(COX_multi_FACSm_forest$HRate))
d <- rep(NA, length(COX_multi_FACSm_forest$HRate))
e <- rep(NA, length(COX_multi_FACSm_forest$HRate))


COX_multi_FACSm_forest <- mutate (COX_multi_FACSm_forest, a, b, c, d, e)
COX_multi_FACSm_forest <- select(COX_multi_FACSm_forest, variables, a, HRate, b, LowerCI, c, UpperCI, d, Pvalue, e)

setnames(COX_multi_FACSm_forest, c("LowerCI", "UpperCI", "a", "b", "c", "d", "e"), 
         c("LCL(95%)", "UCL(95%)", "                                                      ", "   ", "   ", "   ", "   "))

pdf("COXmodel_FACS.multivariate.pdf", width = 12, height = 18 )
forestplot(rbind(colnames(COX_multi_FACSm_forest), as.matrix(COX_multi_FACSm_forest)), 
           c(NA,COX_multi_FACSm_forest$'HRate'),
           c(NA,COX_multi_FACSm_forest$'LCL(95%)'),
           c(NA,COX_multi_FACSm_forest$'UCL(95%)'), 
           zero=0, 
           col=meta.colors(box="red",line="darkblue", summary="royalblue"), 
           is.summary = rep(F, 12), 
           title("MULTIVARIATE Cox regression model - VALIDATION COHORT"),
           xlab="Hazard Rate", xlog=F, xticks=c(seq(-20,20, by=5)), 
           boxsize = 0.6, clip= c(-20, 20))
dev.off()

#####################  COX model combining FACS and CYTOF cohorts, using normalised monocytes values ########################

# select CYTOF+FACS baseline 
list_base <- filter(cytof, TP == "baseline")
list_baseC <- filter(list_base, use == "CyTOF")
list_baseF <- filter(list_base, use == "FACS")

list_baseC <- mutate(list_baseC, normalized.class.monos = (Classical.Monocytes-mean(Classical.Monocytes))/sd(Classical.Monocytes))
list_baseF <- mutate(list_baseF, normalized.class.monos = (Classical.Monocytes-mean(Classical.Monocytes))/sd(Classical.Monocytes))

list_base_norm <-bind_rows(list_baseC, list_baseF)

# test monocytes normalisation

ggplot(data=list_base_norm, aes(x=use, y=Classical.Monocytes))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.3, aes(color=batch))+
  stat_compare_means(method = "t.test",label.x = 0.6, label.y = 20)+
  ggtitle("Classical monocytes frequency in FACS and CyTOF") +
  xlab(" ")+
  ylab("Classical monocytes frequency")
ggsave("class_mono_freq_C+F.pdf")

ggplot(data=list_base_norm, aes(x=use, y=normalized.class.monos))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.3, aes(color=batch))+
  stat_compare_means(method = "t.test",label.x = 0.6, label.y = 2.5)+
  ggtitle("Normalized Classical monocytes frequency in FACS and CyTOF") +
  xlab(" ")+
  ylab("Normalized Classical monocytes frequency")
ggsave("NORMALIZED_class_mono_freq_C+F.pdf")

#without batch and with normalized classical monos
list.3 <- colnames(select(cytof, sex, age.sampling, prev.ipi, prev.chemo, prev.radio, prev.target, BRAF.status, NRAS.status, cKIT.status, leukocytes.count, LDH, ANC, ALC, ANC.ALC, S100, hemoglobin, hematocrit, erythro, MCV, MCH, MCHC, RDW, thrombocytes, monocytes, eosino, baso, IG., sodium, potatium, urea, creatinin, eGFR, bilirubin, protein, albumin, AST, ALT, GGT, alk.phosph, CRP, TSH, brain.met, liver.met, lung.met, bone.met, normalized.class.monos, general.stage, pT, pT.ulceration, N.stage, N.size, M.stage, nb.days.last.ipi)) 


####UNIVARIATE  
COX_uni_base_m <- lapply(list.3, function(x) {                  
  cox <- coxph(Surv(PFS.months, PFS.status) ~ get(x), list_base_norm) 
  cox1 <- tidy(cox)
  cbind(cox1, 
        HRate = round(signif(cox1$estimate,3),3),
        LowerCI = signif(cox1$conf.low,3), 
        UpperCI = signif(cox1$conf.high,3), 
        variables = eval(x), 
        Pvalue = signif(cox1$p.value,3))
})

COX2 <- do.call(rbind, COX_uni_base_m)
COX2 <- select(COX2, variables, term, HRate, LowerCI, UpperCI, Pvalue)
COX2$HRate.exp <- round(signif(exp(COX2$HRate),3),3)
COX2$LowerCI.exp <- signif(exp(COX2$LowerCI),3)
COX2$UpperCI.exp <- signif(exp(COX2$UpperCI),3)
COX2 <- arrange(COX2, Pvalue)

# filter P val < 0.05 and assign to COX_uni_CYTOF_signP
COX_uni_base_signP <- filter(COX2, Pvalue <= 0.05)

COX2$term <- gsub("get\\(x\\)([0-9a-zA-Z]*)", "\\1", COX2$term)
COX2.1 <- unite(COX2, variables, variables, term, sep=" ")


write.table(COX2.1, "/data/Sabrina/CYTOF_revision/COX/tables/Cox_CYTOF_FACS_baseline_univariate.txt", sep="\t", row.names = F)


#Forest plot univariate FACS + CyTOF baseline only
COX_uni_base_forest <- fread("/data/Sabrina/CYTOF_revision/COX/tables/Cox_CYTOF_FACS_baseline_univariate.txt", showProgress = T)
COX_uni_base_forest <- na.omit(COX_uni_base_forest)
COX_uni_base_forest <- filter(COX_uni_base_forest, !grepl("NaN", Pvalue))

var.names <- read.delim("/data/Sabrina/CYTOF_revision/list_variables_names_units_17_07_17.txt", comment.char="#", na.strings = c("NA"," ",""))
var.names <- select(var.names, variables, clean.names)
COX_uni_base_forest <- left_join(COX_uni_base_forest, var.names, by="variables")

COX_uni_base_forest$variables <- as.character(COX_uni_base_forest$clean.names)

a <- rep(NA, length(COX_uni_base_forest$HRate))
b <- rep(NA, length(COX_uni_base_forest$HRate))
c <- rep(NA, length(COX_uni_base_forest$HRate))
d <- rep(NA, length(COX_uni_base_forest$HRate))
e <- rep(NA, length(COX_uni_base_forest$HRate))


COX_uni_base_forest <- mutate (COX_uni_base_forest, a, b, c, d, e)
COX_uni_base_forest <- select(COX_uni_base_forest, variables, a, HRate, b, LowerCI, c, UpperCI, d, Pvalue, e)

setnames(COX_uni_base_forest, c("LowerCI", "UpperCI", "a", "b", "c", "d", "e"), 
         c("LCL(95%)", "UCL(95%)", "                                                         ", "   ", "   ", "   ", "   "))

pdf("COXmodel_CyTOF_FACS.univariate.pdf", width = 12, height = 18 )
forestplot(rbind(colnames(COX_uni_base_forest), as.matrix(COX_uni_base_forest)), 
           c(NA,COX_uni_base_forest$'HRate'),
           c(NA,COX_uni_base_forest$'LCL(95%)'),
           c(NA,COX_uni_base_forest$'UCL(95%)'), 
           zero=0, 
           col=meta.colors(box="red",line="darkblue", summary="royalblue"), 
           is.summary = rep(F, 12), 
           title("UNIVARIATE Cox regression model - BOTH COHORT, baseline samples only"),
           xlab="Hazard Rate", xlog=F, xticks=c(seq(-10,10, by=2)), 
           boxsize = 0.6, clip= c(-10, 10))
dev.off()

############ MULTIVARIATE

signP_list <- COX_uni_base_signP$variables
signP_list1 <- paste(signP_list, collapse = " + ")

COX_multi_base <- coxph(Surv(PFS.months, PFS.status) ~ normalized.class.monos + thrombocytes + LDH + liver.met + alk.phosph + IG. + M.stage + potatium, data = list_base_norm)
COX_multi_base_t <- tidy(COX_multi_base)
COX_multi_base2 <- cbind(COX_multi_base_t, 
                                variables = COX_multi_base_t$term,
                                HRate = round(signif(COX_multi_base_t$estimate,3),3),
                                LowerCI = round(signif(COX_multi_base_t$conf.low,3),3), 
                                UpperCI = round(signif(COX_multi_base_t$conf.high,3),3), 
                                Pvalue = round(signif(COX_multi_base_t$p.value,3),3))
COX_multi_baseline <- select(COX_multi_base2, variables, HRate, LowerCI, UpperCI, Pvalue)

COX_multi_baseline$HRate.exp <- round(signif(exp(COX_multi_baseline$HRate),3),3)
COX_multi_baseline$LowerCI.exp <- signif(exp(COX_multi_baseline$LowerCI),3)
COX_multi_baseline$UpperCI.exp <- signif(exp(COX_multi_baseline$UpperCI),3)

write.table(COX_multi_baseline, "/data/Sabrina/CYTOF_revision/COX/tables/Cox_CYTOF_FACS_baseline_multivariate.txt", sep="\t", row.names = F)

## Forest plot Multivariate analysis FACS + CYTOF baseline only
COX_multi_baseline <- fread("/data/Sabrina/CYTOF_revision/COX/tables/Cox_CYTOF_FACS_baseline_multivariate.txt", showProgress = T)
COX_multi_baseline <- arrange(COX_multi_baseline, Pvalue)
COX_multi_baseline <- left_join(COX_multi_baseline, var.names, by="variables")

COX_multi_baseline$variables <- as.character(COX_multi_baseline$clean.names)

a <- rep(NA, length(COX_multi_baseline$HRate))
b <- rep(NA, length(COX_multi_baseline$HRate))
c <- rep(NA, length(COX_multi_baseline$HRate))
d <- rep(NA, length(COX_multi_baseline$HRate))
e <- rep(NA, length(COX_multi_baseline$HRate))


COX_multi_baseline <- mutate (COX_multi_baseline, a, b, c, d, e)
COX_multi_baseline <- select(COX_multi_baseline, variables, a, HRate, b, LowerCI, c, UpperCI, d, Pvalue, e)

setnames(COX_multi_baseline, c("LowerCI", "UpperCI", "a", "b", "c", "d", "e"), 
         c("LCL(95%)", "UCL(95%)", "                                                      ", "   ", "   ", "   ", "   "))

pdf("COXmodel_CyTOF_FACS.multivariate.pdf", width = 12, height = 18 )
forestplot(rbind(colnames(COX_multi_baseline), as.matrix(COX_multi_baseline)), 
           c(NA,COX_multi_baseline$'HRate'),
           c(NA,COX_multi_baseline$'LCL(95%)'),
           c(NA,COX_multi_baseline$'UCL(95%)'), 
           zero=0, 
           col=meta.colors(box="red",line="darkblue", summary="royalblue"), 
           is.summary = rep(F, 12), 
           title("MULTIVARIATE Cox regression model - BOTH COHORT, baseline samples only"),
           xlab="Hazard Rate", xlog=FALSE, xticks=c(seq(-5,5, by=1)), 
           boxsize = 0.6, clip= c(-5, 5))
dev.off()


######### TEST STEPWISE MODEL

test_C <- select(CYTOFlist_base, PFS.months, PFS.status, Classical.Monocytes,normalized.class.monos, sex, age.sampling, prev.ipi, prev.chemo, prev.radio, leukocytes.count, LDH, ANC, ALC, ANC.ALC, S100,  monocytes )
test_F <- select(FACSlist, PFS.months, PFS.status, Classical.Monocytes,normalized.class.monos, sex, age.sampling, prev.ipi, prev.chemo, prev.radio, leukocytes.count, LDH, ANC, ALC, ANC.ALC, S100,  monocytes )
test_CF <- select(list_base_norm, PFS.months, PFS.status, Classical.Monocytes,normalized.class.monos, sex, age.sampling, prev.ipi, prev.chemo, prev.radio, leukocytes.count, LDH, ANC, ALC, ANC.ALC, S100,  monocytes )


COX_CyTOF_STEP1 <- coxph(Surv(PFS.months,PFS.status) ~ .  ,data= test_C)
COX_CyTOF_STEP2 <- step(COX_CyTOF_STEP1, direction="both")

COX_FACS_STEP1 <- coxph(Surv(PFS.months,PFS.status) ~ .  ,data= test_F)
COX_FACS_STEP2 <- step(COX_FACS_STEP1, direction="both")

COX_CyTOF_FACS_STEP1 <- coxph(Surv(PFS.months,PFS.status) ~ .  ,data= test_CF)
COX_CyTOF_FACS_STEP2 <- step(COX_CyTOF_FACS_STEP1, direction="both")

