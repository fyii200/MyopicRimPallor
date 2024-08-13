##############################################################################################################
##                                              Fabian Yii                                                  ##
##                                          fabian.yii@ed.ac.uk                                             ##
##############################################################################################################
rm(list=ls())
library(sjPlot)
library(ggplot2)
library(lme4)
library(car)
library(dplyr)
library(tidyr)
library(lqmm)
library(RColorBrewer)
source("helperFunctions.R")
options(scipen=999)

## Set up
# # Move up one directory to the and set it as working directory 
# setwd("..") 
# Read full data (137016 eyes of 68508 participants)
d  <- read.csv(paste0("data", .Platform$file.sep, "linked_cleaned_data_long_all.csv")) 
# Blood pressure data
BP <- read.csv(paste0("data", .Platform$file.sep, "UKB_cleaned_data_bloodMeasuresAndBMI_all.csv")) 
# Retinal parameters data
RP <- read.csv(paste0("data", .Platform$file.sep, "fundusDerivedDataReturnToUKbiobank.csv"))
names(RP)[1:2] <- c("id", "fundus")
# Merge datasets
d <- merge(d, BP, all.x=TRUE, by=c("id", "age", "YOB", "sex", "townsend") )
d <- merge(d, RP, all.x=TRUE, by=c("id", "fundus"))
# Remove eyes without fundus photographs (135787 eyes of 68508 participants left)
d <- d[d$fundus != "", ]
# Calculate mean arterial blood pressure: (0.33 × systolic BP) + (0.67 × diastolic BP), as defined
# in https://www.sciencedirect.com/science/article/pii/S0002939404000364?via%3Dihub#section.0010
d$meanBP <- (0.33 * d$SysBP) + (0.67 * d$DiasBP)
# Categorise participants into "White" and "non-White"
d <- d[d$ethnic!="" & d$ethnic!="Do not know" & d$ethnic!="Prefer not to answer", ] # remove participants with missing ethnicity info (978 eyes removed; 134809 eyes of 68012 participants left)
d$ethnicGroup <- ifelse(d$ethnic=="British" | d$ethnic=="Irish" | d$ethnic=="White" | d$ethnic=="Any other white background", "White", "non-White")
d$ethnicGroup <- factor(d$ethnicGroup, levels=c("White", "non-White"))

############################### Part 1: select cohort for subsequent analysis ###############################
## Remove participants with systemic conditions
# Hypertension (35082 eyes removed)
# 99727 eyes of 50299 participants left
d <- subset(d, hypertension==FALSE)
# Diabetes (2506 eyes removed)
# 97221 eyes of 49035 participants left
d <- subset(d, diabetes==FALSE)
# Heart conditions (2159 eyes removed)
# 95062 eyes of 47944 participants left
d$heartConditions <- TRUE
d[which(d$heartFailure==FALSE & d$ischaemicHeartDisease==FALSE & d$multipleValvularHeartDisease==FALSE & d$myocardialInfarction==FALSE & d$cardiomyopathy==FALSE & d$otherHeartDiseases==FALSE), ]$heartConditions <- FALSE
d <- subset(d, heartConditions==FALSE)
## Remove eyes with posterior ocular conditions other than glaucoma
# Chorioretinal diseases (1238 eyes)
# 93824 eyes of 47313 participants left
d <- subset(d, chorioretinalDiseases==FALSE)
# Scleral or globe diseases (314 eyes)
# 93510 eyes of 47153 participants left
d <- subset(d, scleralDiseases==FALSE & globeDiseases==FALSE)
# Optic nerve disorders other than glaucoma (135 eyes)
# 93375 eyes of 47083 participants left
d <- subset(d, opticNerveDiseases==FALSE)
# Glaucoma (1215 eyes)
# 92160 eyes of 46464 participants left
d <- subset(d, glaucoma==FALSE)

## Remove eyes with poor or missing VA (37255 eyes)
# 54905 eyes of 34284 participants left
d <- subset(d, VA<=0 & !is.na(VA))

## Remove eyes without refractive error measurements (638 eyes)
# 54267 eyes of 33953 participants left
d <- subset(d, !is.na(SER))

## Remove eyes with missing or zero IOP (1560 eyes)
# 52707 eyes of 33093 participants left
d <- subset(d, !is.na(IOP) & IOP!=0)

## Remove eyes that fail quality control (18757 eyes)
# 33950 eyes of 24057 participants left
d <- subset(d, palProcError==0 & palRejectThres==0)
excluded <- subset(d, palProcError==1 | palRejectThres==1)

## Compare the characteristics of those passing and failing quality control
# Age
t.test(d[!duplicated(d$id),]$age, excluded[!duplicated(excluded$id),]$age)
# Sex
table(d[!duplicated(d$id),]$sex); table(excluded[!duplicated(excluded$id),]$sex)
observed_table <- matrix(c(13927, 10130, 8218, 6474), nrow = 2, ncol = 2, byrow = T)
rownames(observed_table) <- c('Included', 'Excluded')
colnames(observed_table) <- c('Female', 'Male')
chisq.test(observed_table)
# SER
t.test(d$SER, excluded$SER)
# Ethnic group
table(d[!duplicated(d$id),]$ethnicGroup); table(excluded[!duplicated(excluded$id),]$ethnicGroup)
observed_table <- matrix(c(22048, 2009, 13537, 1155), nrow = 2, ncol = 2, byrow = T)
rownames(observed_table) <- c('Included', 'Excluded')
colnames(observed_table) <- c('White', 'non-White')
chisq.test(observed_table)
# Intraocular pressure
t.test(d$IOP, excluded$IOP)
# Mean blood pressure
t.test(d[!duplicated(d$id),]$meanBP, excluded[!duplicated(excluded$id),]$meanBP)
################################################# Part 1 COMPLETE #################################################


#################################### Part 2: association between SER and pallor ###################################
## NRR pallor vs SER using linear mixed-effects models (LMM)
# More positive pallor values = more pale
OLSglobal <- lmer(palG ~ SER*ethnicGroup + age + sex + IOP + meanBP + (1|id), d); tab_model(OLSglobal, digits=5)
OLS_T     <- lmer(palT ~ SER*ethnicGroup + age + sex + IOP + meanBP + (1|id), d); tab_model(OLS_T, digits=5)
OLS_TI    <- lmer(palTI ~ SER*ethnicGroup + age + sex + IOP + meanBP + (1|id), d); tab_model(OLS_TI, digits=5)
OLS_NI    <- lmer(palNI ~ SER*ethnicGroup + age + sex + IOP + meanBP + (1|id), d); tab_model(OLS_NI, digits=5)
OLS_N     <- lmer(palN ~ SER*ethnicGroup + age + sex + IOP + meanBP + (1|id), d); tab_model(OLS_N, digits=5)
OLS_NS    <- lmer(palNS ~ SER*ethnicGroup + age + sex + IOP + meanBP + (1|id), d); tab_model(OLS_NS, digits=5)
OLS_TS    <- lmer(palTS ~ SER*ethnicGroup + age + sex + IOP + meanBP + (1|id), d); tab_model(OLS_TS, digits=5)
# Check model assumptions (normality and homoscedasticity of residuals)
LMMmodels <- c(OLSglobal, OLS_T, OLS_TS, OLS_TI, OLS_N, OLS_NS, OLS_NI)
for(model in LMMmodels){qqnorm(resid(model)); qqline(resid(model)); plot(fitted(model),residuals(model))}
# Dummy dataframe to save LMM predictions
zones  <- c("palG", "palT", "palTS", "palTI", "palN", "palNS", "palNI")
newDat <- expand.grid(SER=seq(-20, 10, by=1), ethnicGroup=unique(d$ethnicGroup), age=mean(d$age), sex="Female", IOP=mean(d$IOP,na.rm=TRUE), meanBP=mean(d$meanBP,na.rm=TRUE), zone=zones, pallor=NA)
# Save LMM predictions based on the dummy dataframe created above
for(i in 1:length(zones)){
  m <- LMMmodels[i][[1]]
  newDat[which(newDat$zone==zones[i]),]$pallor <- predict(m, subset(newDat, zone==zones[i]), re.form=NA) }

## Secondary analysis: association of OD-fovea distance and vessel concavity with NRR pallor, controlling for SER
tab_model(lmer(palT~adj_OD_fovea_dist + concavity_artery + concavity_vein + SER + age + sex + ethnicGroup + IOP + meanBP + (1|id), d), digits=5)
tab_model(lmer(palTI~adj_OD_fovea_dist + concavity_artery + concavity_vein + SER + age + sex + ethnicGroup + IOP + meanBP + (1|id), d), digits=5)
tab_model(lmer(palNI~adj_OD_fovea_dist + concavity_artery + concavity_vein + SER + age + sex + ethnicGroup + IOP + meanBP + (1|id), d), digits=5)
tab_model(lmer(palN~adj_OD_fovea_dist + concavity_artery + concavity_vein + SER + age + sex + ethnicGroup + IOP + meanBP + (1|id), d), digits=5)
tab_model(lmer(palNS~adj_OD_fovea_dist + concavity_artery + concavity_vein + SER + age + sex + ethnicGroup + IOP + meanBP + (1|id), d), digits=5)
tab_model(lmer(palTS~adj_OD_fovea_dist + concavity_artery + concavity_vein + SER + age + sex + ethnicGroup + IOP + meanBP + (1|id), d), digits=5)
################################################# Part 2 COMPLETE #################################################


################################################### Part 3: plots #################################################
## Set common plot theme
myTheme <- theme_linedraw() + theme(panel.grid.minor=element_blank(),
                                    panel.grid.major=element_blank(),
                                    panel.spacing.x=unit(0.1, "lines"),
                                    strip.background=element_rect(colour=NA, fill="gray90"),
                                    strip.text=element_text(colour="black", face="bold"),
                                    panel.border=element_rect(colour=NA, fill=alpha("gray", 0.1)),
                                    plot.margin=margin(.2,.5,.2,.2, "cm"),
                                    axis.title.x=element_text(margin=margin(t=15)),
                                    axis.title.y=element_text(margin=margin(r=15)),
                                    axis.text=element_text(size=10),
                                    legend.title=element_text(face="bold", size=10)) 

## Global NRR and disc pallor vs SER in eyes without glaucoma
d %>% pivot_longer(cols=c(palG, palDisc), names_to="zone", values_to="pallor") %>% 
  ggplot(aes(x=SER, y=pallor, color=ethnicGroup, size=ethnicGroup)) + myTheme + geom_point(alpha=0.3) +
  geom_line(data=subset(newDat, zone=="palDisc"), linewidth=1) +
  geom_line(data=subset(newDat, zone=="palG"), linewidth=1) +
  scale_color_manual("White", values=c("black", "red"), labels=c("Yes", "No")) +
  scale_size_manual(values=c(0.7, 1.1), guide = "none") +
  xlab("Spherical equivalent refraction (D)") + ylab("Global pallor") +
  facet_wrap(~factor(zone, levels=c("palG", "palDisc")), nrow=2, scales="free_y", labeller=as_labeller(c("palDisc"="Optic disc", "palG"="Neuroretinal rim")))
ggsave(paste0("figures", .Platform$file.sep, "globalPallorVsSER.pdf"), width=8, height=7)

## Sectoral NRR pallor profile plot in eyes without glaucoma
# Put pallor derived from different sectors in separate rows
dLong              <- d %>% pivot_longer(cols=palT:palTS,names_to = "zone", values_to = "pallor")
dLong$zone         <- factor(dLong$zone, levels=c("palT", "palTI", "palNI", "palN", "palNS", "palTS"))
dLong$facetGroup   <- interaction(dLong$ethnicGroup, dLong$zone)
labs               <- c("-14.5", "-10", "-5.5", "-1", "3.5", "8")
dLong              <- dLong %>% mutate(SERgroup=cut(SER, breaks=c(min(d$SER), -6, -3, 0, 3, 6.02, max(d$SER)), include.lowest=TRUE, labels=labs))
dLong$SERgroup     <- as.numeric(levels(dLong$SERgroup))[dLong$SERgroup]
# Plot
dLong %>% ggplot(aes(x=zone)) + 
  geom_line(aes(y=pallor, group=id, alpha=ethnicGroup, col=ethnicGroup)) + 
  labs(x="", y="Sectoral neuroretinal rim pallor") +
  scale_x_discrete(labels=c("T", "TI", "NI", "N", "NS", "TS")) +
  scale_color_manual(name="White", labels=c("Yes", "No"), values=c("gray80", "darkgreen")) +
  scale_alpha_manual(values=c(0.06, 0.11), guide="none") +
  guides(col=guide_legend(override.aes=list(size=5, linewidth=1.3))) +
  myTheme + theme(panel.border=element_rect(colour=NA, fill=NA), legend.position="top")
ggsave(paste0("figures", .Platform$file.sep, "sectoralPallorProfile.pdf"), width=6, height=5.5)

## Sectoral NRR pallor vs SER in eyes without glaucoma
# Specify facet labels
facetTitles=c("White.palT"="Temporal", "White.palTI"="Temporal inferior", "White.palNI"="Nasal inferior",
              "White.palN"="Nasal", "White.palNS"="Nasal superior", "White.palTS"="Temporal superior",
              "non-White.palT"="Temporal", "non-White.palTI"="Temporal inferior", "non-White.palNI"="Nasal inferior",
              "non-White.palN"="Nasal", "non-White.palNS"="Nasal superior", "non-White.palTS"="Temporal superior")
# Define facet group
newDat$facetGroup <- paste0(newDat$ethnicGroup, ".", newDat$zone)
# Plot
dLong %>% ggplot() + labs(x="", y="Sectoral neuroretinal rim pallor") +
  ## Violin plots ##
  geom_violin(aes(x=SERgroup, y=pallor, fill=factor(SERgroup)), color=NA, alpha=0.6) +
  geom_violin_quantiles(aes(x=SERgroup, y=pallor, fill=factor(SERgroup)), quantiles=0.5, linesize=1, col="gray45", show.legend=F) +
  ## Temporal ##
  geom_line(aes(x=SER, y=pallor), subset(newDat, facetGroup=="White.palT"), size=1.2) +
  geom_line(aes(x=SER, y=pallor), subset(newDat, facetGroup=="non-White.palT"), size=1.2) +
  ## Temporal inferior ##
  geom_line(aes(x=SER, y=pallor), subset(newDat, facetGroup=="White.palTI"), size=1.2) +
  geom_line(aes(x=SER, y=pallor), subset(newDat, facetGroup=="non-White.palTI"), size=1.2) +
  ## Temporal superior ##
  geom_line(aes(x=SER, y=pallor), subset(newDat, facetGroup=="White.palTS"), size=1.2) +
  geom_line(aes(x=SER, y=pallor), subset(newDat, facetGroup=="non-White.palTS"), size=1.2) +
  ## Nasal ##
  geom_line(aes(x=SER, y=pallor), subset(newDat, facetGroup=="White.palN"), size=1.2) +
  geom_line(aes(x=SER, y=pallor), subset(newDat, facetGroup=="non-White.palN"), size=1.2) +
  ## Nasal inferior ##
  geom_line(aes(x=SER, y=pallor), subset(newDat, facetGroup=="White.palNI"), size=1.2) +
  geom_line(aes(x=SER, y=pallor), subset(newDat, facetGroup=="non-White.palNI"), size=1.2) +
  ## Nasal superior ##
  geom_line(aes(x=SER, y=pallor), subset(newDat, facetGroup=="White.palNS"), size=1.2) +
  geom_line(aes(x=SER, y=pallor), subset(newDat, facetGroup=="non-White.palNS"), size=1.2) +
  ## Legends ##
  scale_fill_manual(name="Spherical equivalent refraction (D)", labels=c("< -6", "-6 to -3", "-3 to 0", "0 to +3", "+3 to +6", "> +6"), values=alpha(brewer.pal(11,"Spectral")[c(1,3,5,7,9,10)], 0.7)) +
  ## Subplots ##
  facet_wrap(~factor(facetGroup, levels=c("White.palT", "White.palTI", "White.palNI", "White.palN", "White.palNS", "White.palTS", "non-White.palT", "non-White.palTI", "non-White.palNI", "non-White.palN", "non-White.palNS", "non-White.palTS")), labeller=as_labeller(facetTitles), scales="fixed", nrow=2) + myTheme +
  guides(fill=guide_legend(nrow=1), color=guide_legend(override.aes=list(size=2))) + theme(legend.position="bottom", legend.margin=margin(-25, 5, 0, 5), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_y_continuous(sec.axis=sec_axis(~ ., name = "WHITE                                                                         NON-WHITE", breaks=NULL, labels=NULL))
ggsave(paste0("figures", .Platform$file.sep, "sectoralPallorVsSER.pdf"), width=8.2, height=7.8)  
################################################# Part 3 COMPLETE ################################################
  











